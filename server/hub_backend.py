"""Literature search backend for the snp2cell hub viewer (PubMed direct).

Calls NCBI Entrez E-utilities and Europe PMC REST directly — NO LLM web search.
Citations are authoritative (real PMIDs, real titles, real DOIs). LLM is only
invoked as an OPTIONAL second step to synthesize the retrieved abstracts into
a short summary; it never invents citations.

Endpoints:
    GET  /literature/<type>/<key>     → return cached JSON if exists, else 404
    POST /literature/<type>/<key>     → run a fresh PubMed search; return result
    GET  /literature/list             → list cached entries
    GET  /healthz                     → liveness check

`type` ∈ {"snp", "tf", "gene"}.  `key` = rsid or symbol.

Body for POST (JSON):
    {
      "context": {
        "disease": "AF",
        "target_cs": "MyocardialSleeveCells",
        "rsid": "rs6702619",          # for type=snp
        "chr": "1", "pos": 99580690,   # for type=snp (extra search terms)
        "symbol": "PALMD",             # for type=tf/gene
        "extra_terms": "atrial myocyte" # optional extra free-text
      },
      "force": false,                   # set true to bypass cache
      "synthesize": false,              # set true to ask claude/codex for a
                                        # short synthesis on top of retrieved
                                        # abstracts (no web search; pure
                                        # summarization of provided text)
      "max_results": 20,
      "min_year": 1990
    }

Run alongside the http.server hosting the hub:
    python scripts/snp2cell_hub_backend.py --port 8766
    # optional --ncbi-api-key=<KEY> to raise rate-limit from 3/s to 10/s
    # optional --email=<contact> for NCBI compliance
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import subprocess
import sys
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, urlencode
import urllib.request
import urllib.error

log = logging.getLogger(__name__)

ALLOWED_TYPES = {"snp", "tf", "gene"}
DEFAULT_PORT = 8766

# ---- AG1 (AlphaGenome) integration --------------------------------------
AG1_PYTHON = os.environ.get(
    "AG1_PYTHON",
    "/Users/takahiro/miniforge3/envs/scanpy/bin/python",
)
AG1_SCRIPT = Path(os.environ.get(
    "AG1_SCRIPT",
    "/Users/takahiro/Desktop/project/hypersampling/ag1/pipeline_final/integration/ag1.py",
))
AG1_VARIANT_RE = re.compile(r"^(?:chr)?([0-9XYM]+)_(\d+)_([ACGT]+)_([ACGT]+)$")
# CORS allow-list. The deployed Cloudflare Pages URL is read from env
# (CF_PAGES_ORIGIN) so we don't bake a fixed domain into the image.
ALLOWED_ORIGINS = tuple(filter(None, [
    "http://localhost:8765", "http://127.0.0.1:8765",
    "http://localhost:8000", "http://127.0.0.1:8000",
    os.environ.get("ALLOWED_ORIGIN"),                               # primary deploy URL (e.g. https://178-105-162-190.nip.io)
    os.environ.get("CF_PAGES_ORIGIN"),                              # legacy: Cloudflare Pages URL
    os.environ.get("CUSTOM_ORIGIN"),                                # additional custom origin
]))

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# ---- In-process rate limiter ------------------------------------------------
# Bucket per (client_ip, endpoint) → list of recent request timestamps.
# Cleared lazily on each check.
import threading as _rl_threading
_RATE_WINDOW_SEC = 60
_RATE_LIMITS = {
    "/chat": int(os.environ.get("CHAT_RATE_PER_MIN", "20")),
    # Add more endpoints here if needed; absent endpoints are unlimited.
}
_rl_buckets: dict[tuple[str, str], list[float]] = {}
_rl_lock = _rl_threading.Lock()


def rate_limit_check(client_ip: str, endpoint: str) -> tuple[bool, int]:
    """Return (allowed, retry_after_sec). retry_after_sec is 0 when allowed."""
    cap = _RATE_LIMITS.get(endpoint)
    if cap is None:
        return True, 0
    now = time.time()
    cutoff = now - _RATE_WINDOW_SEC
    key = (client_ip, endpoint)
    with _rl_lock:
        bucket = [t for t in _rl_buckets.get(key, []) if t > cutoff]
        if len(bucket) >= cap:
            oldest = bucket[0]
            return False, max(1, int(_RATE_WINDOW_SEC - (now - oldest)))
        bucket.append(now)
        _rl_buckets[key] = bucket
    return True, 0

# Disease label → preferred PubMed query expansion
DISEASE_TERMS = {
    "AF": "(atrial fibrillation[MeSH Terms] OR atrial fibrillation OR atrial flutter)",
    "CAD": "(coronary artery disease[MeSH Terms] OR coronary heart disease OR myocardial infarction[MeSH Terms])",
    "AVS": "(aortic valve stenosis[MeSH Terms] OR calcific aortic stenosis OR aortic valve calcification)",
}


def safe_key(key: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_:.\-]", "_", key)[:64]
    return cleaned or hashlib.sha1(key.encode()).hexdigest()[:12]


def cache_path(hub_dir: Path, type_: str, key: str, context_extra: str = "") -> Path:
    safe = safe_key(key)
    if context_extra:
        h = hashlib.sha1(context_extra.encode()).hexdigest()[:8]
        safe = f"{safe}__{h}"
    return hub_dir / "literature" / type_ / f"{safe}.json"


# ---- query builders ---------------------------------------------------------

def build_queries(type_: str, ctx: dict) -> list[tuple[str, str]]:
    """Return a list of (label, query) candidates from strict to loose.
    The backend tries each in order until at least 1 PubMed hit is returned.
    Cell_state is NEVER added as an AND clause — it is not a MeSH/indexed term
    and would zero-out the result set. It's used only for the synthesis step.
    """
    disease_id = ctx.get("disease", "")
    disease_term = DISEASE_TERMS.get(disease_id, disease_id)
    extra = ctx.get("extra_terms", "")

    candidates: list[tuple[str, str]] = []

    if type_ == "snp":
        rsid = ctx.get("rsid", "")
        chrom = ctx.get("chr", "")
        pos = ctx.get("pos", "")
        if rsid and rsid.startswith("rs"):
            if disease_term:
                candidates.append(("rsid+disease", f"{rsid}[All Fields] AND {disease_term}"))
            # Loose: just rsid (catches all literature mentioning the variant)
            candidates.append(("rsid", f"{rsid}[All Fields]"))
        elif chrom and pos:
            # No rs# — use coordinate as soft hint
            if disease_term:
                candidates.append(("region+disease",
                                   f"({chrom}p[All Fields] OR {chrom}q[All Fields]) AND {disease_term}"))

    elif type_ in ("tf", "gene"):
        sym = ctx.get("symbol", "")
        if sym:
            sym_clause = f"({sym}[Title/Abstract] OR {sym}[Gene Symbol])"
            if disease_term:
                candidates.append((f"{type_}+disease", f"{sym_clause} AND {disease_term}"))
            candidates.append((type_, sym_clause))
    else:
        raise ValueError(f"unknown type {type_!r}")

    if extra:
        # Tighten the most specific query with the extra terms
        if candidates:
            label, q = candidates[0]
            candidates.insert(0, (f"{label}+extra", f"{q} AND ({extra})"))
    return candidates


# ---- HTTP helper ------------------------------------------------------------

def _http_get(url: str, params: dict, timeout: int = 20,
              api_key: str | None = None, email: str | None = None) -> bytes:
    if api_key:
        params = {**params, "api_key": api_key}
    if email:
        params = {**params, "email": email, "tool": "snp2cell-hub"}
    qs = urlencode(params)
    req = urllib.request.Request(url + "?" + qs,
                                 headers={"User-Agent": "snp2cell-hub-backend/2.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


# ---- PubMed -----------------------------------------------------------------

def pubmed_search(query: str, max_results: int = 20, *,
                  min_year: int | None = None,
                  api_key: str | None = None,
                  email: str | None = None) -> dict:
    """E-utilities esearch → esummary. Returns {citations: [...], query: str}."""
    if min_year:
        q = f"{query} AND {min_year}:3000[dp]"
    else:
        q = query
    # esearch
    raw = _http_get(NCBI_ESEARCH, {
        "db": "pubmed", "term": q, "retmax": str(max_results),
        "retmode": "json", "sort": "relevance",
    }, api_key=api_key, email=email)
    es = json.loads(raw).get("esearchresult", {})
    ids = es.get("idlist", [])
    if not ids:
        return {"citations": [], "query": q, "count": 0}
    # esummary
    raw = _http_get(NCBI_ESUMMARY, {
        "db": "pubmed", "id": ",".join(ids), "retmode": "json",
    }, api_key=api_key, email=email)
    summ = json.loads(raw).get("result", {})
    out = []
    for pid in ids:
        e = summ.get(pid, {})
        if not e:
            continue
        # Extract DOI from articleids
        doi = None
        for ai in e.get("articleids") or []:
            if ai.get("idtype") == "doi":
                doi = ai.get("value")
                break
        authors = ", ".join(a.get("name", "") for a in (e.get("authors") or [])[:3])
        if len(e.get("authors", []) or []) > 3:
            authors += " et al."
        out.append({
            "pmid": pid,
            "doi": doi,
            "title": e.get("title", "").strip("."),
            "authors": authors,
            "journal": e.get("source"),
            "year": int(e.get("pubdate", "0").split()[0]) if e.get("pubdate") else None,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
        })
    return {"citations": out, "query": q, "count": int(es.get("count") or len(ids))}


def europe_pmc_search(query: str, max_results: int = 10) -> dict:
    """Europe PMC search — broader coverage, includes preprints. Optional."""
    raw = _http_get(EUROPE_PMC, {
        "query": query,
        "format": "json",
        "pageSize": str(max_results),
        "resultType": "lite",
    })
    res = json.loads(raw).get("resultList", {}).get("result", [])
    out = []
    for r in res:
        pmid = r.get("pmid")
        doi = r.get("doi")
        out.append({
            "pmid": pmid,
            "doi": doi,
            "title": r.get("title", "").strip("."),
            "authors": r.get("authorString", ""),
            "journal": r.get("journalTitle"),
            "year": int(r["pubYear"]) if r.get("pubYear") else None,
            "url": (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid
                    else f"https://europepmc.org/article/{r.get('source')}/{r.get('id')}"),
            "source": r.get("source"),
        })
    return {"citations": out, "query": query, "count": len(out)}


# ---- LLM synthesis (optional) ----------------------------------------------

class AgentError(Exception):
    pass


# BYOK keys can appear in provider SDK exception strings (Anthropic, Google).
# Redact them before surfacing any AgentError text to clients or logs.
_API_KEY_RE = re.compile(r"(sk-ant-[A-Za-z0-9_\-]{8,}|AIza[A-Za-z0-9_\-]{20,})")


def _scrub_key(s) -> str:
    """Replace anything that looks like an Anthropic / Google API key with
    '[redacted-key]'. Apply to every exception we wrap into AgentError so the
    user's BYOK secret never reaches the SSE/JSON response or journalctl."""
    try:
        return _API_KEY_RE.sub("[redacted-key]", str(s))
    except Exception:
        return "[unprintable]"


def _atomic_write_text(path: Path, content: str) -> None:
    """Write content to `path` via a temp file + os.replace so a daemon
    thread killed mid-flush cannot leave a half-written file that crashes
    a live read. The replace is atomic on POSIX."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}.{os.urandom(4).hex()}")
    try:
        tmp.write_text(content)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try: tmp.unlink()
            except OSError: pass


def _load_json_safe(path: Path) -> dict | None:
    """Read a JSON file and return its parsed dict, OR None if the file is
    missing, unreadable, or truncated/corrupt. Used for any cache file that
    a concurrent writer might be touching."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        log.warning("cache read failed for %s: %s", path, e)
        return None


def _call_anthropic_sdk(prompt: str, *, api_key: str,
                         allowed_tools: list[str] | None = None,
                         timeout: int = 300,
                         model: str = "claude-sonnet-4-6",
                         max_tokens: int = 4096) -> str:
    """Call Anthropic API directly via SDK; return concatenated text content.

    api_key: user-supplied BYOK Anthropic key (sk-ant-...). NEVER stored or
    logged on this server — passed through to Anthropic for this request only.
    allowed_tools: if it contains "WebSearch", enable Anthropic's built-in
    web_search server-side tool (handles its own retrieval + tool-use loop).
    """
    try:
        from anthropic import Anthropic, APIError, AuthenticationError, RateLimitError
    except ImportError:
        raise AgentError("anthropic SDK not installed (pip install anthropic)")
    if not api_key or not api_key.startswith("sk-ant-"):
        raise AgentError("invalid Anthropic API key (must start with 'sk-ant-')")

    client = Anthropic(api_key=api_key, timeout=timeout)
    tools = []
    if allowed_tools and any(t in ("WebSearch", "web_search") for t in allowed_tools):
        # Tool-use round-trips re-send the full prompt + accumulated
        # tool_use/tool_result blocks each turn, which compounds against
        # the org's per-minute input-token cap (30k/min Anthropic Tier 1,
        # 1M/min Gemini 2.0 free) AND lengthens the user's wait. Keep
        # budgets minimal: 2 searches + 3 fetches = enough for 2 strong
        # citations with verification. Tightened from 3/5 → 2/3 after
        # free-tier users were exhausting TPM in a handful of chats.
        tools.append({"type": "web_search_20250305",
                      "name": "web_search",
                      "max_uses": 2})
        tools.append({"type": "web_fetch_20250910",
                      "name": "web_fetch",
                      "max_uses": 3})

    try:
        resp = client.messages.create(
            model=model,
            max_tokens=max_tokens,
            tools=tools or [],
            messages=[{"role": "user", "content": prompt}],
        )
    except AuthenticationError as e:
        raise AgentError(f"anthropic auth failed (check your API key): {_scrub_key(e)}")
    except RateLimitError as e:
        raise AgentError(f"anthropic rate-limited (your key has hit Anthropic's RPM/TPM cap): {_scrub_key(e)}")
    except APIError as e:
        raise AgentError(f"anthropic API error: {_scrub_key(e)}")

    text_parts = [b.text for b in resp.content if getattr(b, "type", None) == "text"]
    return "".join(text_parts).strip()


def _call_gemini_sdk(prompt: str, *, api_key: str,
                      allowed_tools: list[str] | None = None,
                      timeout: int = 300,
                      model: str | None = None,
                      max_tokens: int = 4096) -> str:
    """Call Google Gemini API; return text. Mirrors _call_anthropic_sdk.

    api_key: BYOK Google API key (AIzaSy...) from aistudio.google.com.
    allowed_tools: ['WebSearch'] enables Gemini's google_search grounding
    (Google's equivalent of Anthropic's web_search built-in tool).
    """
    try:
        from google import genai
        from google.genai import types
    except ImportError:
        raise AgentError("google-genai SDK not installed (pip install google-genai)")
    if not api_key:
        raise AgentError("invalid Google API key (empty)")

    # Default: gemini-2.0-flash. Free tier gets 15 RPM / 1M TPM vs 2.5-flash's
    # 10 RPM / 250k TPM, and we see far fewer 503s from 2.0. User can override
    # via GEMINI_MODEL env var if they want 2.5-flash quality on a paid plan.
    model = model or os.environ.get("GEMINI_MODEL", "gemini-2.0-flash")
    client = genai.Client(api_key=api_key)
    tools_arg = []
    if allowed_tools and any(t in ("WebSearch", "web_search") for t in allowed_tools):
        tools_arg.append(types.Tool(google_search=types.GoogleSearch()))
        # url_context lets Gemini fetch a specific URL (e.g. the PubMed page
        # for a candidate PMID) to verify reachability + title before citing.
        try:
            tools_arg.append(types.Tool(url_context=types.UrlContext()))
        except AttributeError:
            pass  # older google-genai versions without UrlContext type

    try:
        resp = client.models.generate_content(
            model=model,
            contents=prompt,
            config=types.GenerateContentConfig(
                max_output_tokens=max_tokens,
                tools=tools_arg or None,
            ),
        )
    except Exception as e:
        emsg = str(e).lower()
        if "api key" in emsg or "permission" in emsg or "auth" in emsg or "invalid argument" in emsg:
            raise AgentError(f"gemini auth failed (check your Google API key): {_scrub_key(e)}")
        if "quota" in emsg or "rate" in emsg or "exceeded" in emsg or "resource_exhausted" in emsg:
            raise AgentError(f"gemini quota/rate exceeded: {_scrub_key(e)}")
        raise AgentError(f"gemini API error: {_scrub_key(e)}")
    return (resp.text or "").strip()


def _call_llm_sdk(prompt: str, *, api_key: str, **kwargs) -> str:
    """Dispatch to the provider implied by the api_key prefix.
       'sk-ant-...' -> Anthropic ; otherwise -> Google Gemini.
    """
    if api_key and api_key.startswith("sk-ant-"):
        return _call_anthropic_sdk(prompt, api_key=api_key, **kwargs)
    return _call_gemini_sdk(prompt, api_key=api_key, **kwargs)


def _stream_llm_sdk(prompt: str, *, api_key: str,
                     allowed_tools: list[str] | None = None,
                     timeout: int = 300,
                     model: str | None = None,
                     max_tokens: int = 4096):
    """Yield (kind, payload) tuples from the LLM streaming endpoint.

      kind == 'text' : payload is a string of text deltas
      kind == 'tool' : payload is {'name': 'web_search'|'web_fetch'|...,
                                   'input': {<input dict>}}

    Dispatch by api_key prefix (sk-ant-* → Anthropic, else Gemini).
    """
    if api_key and api_key.startswith("sk-ant-"):
        try:
            from anthropic import Anthropic, RateLimitError, AuthenticationError, APIError
        except ImportError:
            raise AgentError("anthropic SDK not installed")
        client = Anthropic(api_key=api_key, timeout=timeout)
        tools = []
        if allowed_tools and any(t in ("WebSearch", "web_search") for t in allowed_tools):
            # Same budget as non-streaming path; multi-turn tool use
            # compounds against per-minute input-token caps. See
            # _call_anthropic_sdk for the reasoning.
            tools.append({"type": "web_search_20250305", "name": "web_search", "max_uses": 2})
            tools.append({"type": "web_fetch_20250910", "name": "web_fetch", "max_uses": 3})
        try:
            with client.messages.stream(
                model=model or "claude-sonnet-4-6",
                max_tokens=max_tokens,
                tools=tools or [],
                messages=[{"role": "user", "content": prompt}],
            ) as stream:
                # Iterate raw events so we can surface tool_use blocks to
                # the UI ("Searching PubMed for ..." etc.) in addition to
                # the text deltas. Server-side tools (web_search,
                # web_fetch) appear as content_block_start events with
                # block.type == 'server_tool_use'.
                for event in stream:
                    et = getattr(event, "type", None)
                    if et == "content_block_start":
                        block = getattr(event, "content_block", None)
                        bt = getattr(block, "type", None)
                        if bt in ("server_tool_use", "tool_use"):
                            yield ("tool", {
                                "name": getattr(block, "name", "tool"),
                                "input": getattr(block, "input", None) or {},
                            })
                    elif et == "content_block_delta":
                        delta = getattr(event, "delta", None)
                        dt = getattr(delta, "type", None)
                        if dt == "text_delta":
                            text = getattr(delta, "text", "") or ""
                            if text:
                                yield ("text", text)
        except RateLimitError as e:
            # Map the per-minute input-token cap (30k/min on Tier 1) to a
            # user-friendly hint that suggests waiting or switching key.
            raise AgentError(
                "anthropic rate-limited (your Anthropic key hit its per-minute "
                "input-token cap, typically 30,000/min on Tier 1). "
                "Wait ~60s and retry, upgrade the Anthropic plan, or switch "
                f"to a Gemini key in the BYOK panel. Raw: {_scrub_key(e)}"
            )
        except AuthenticationError as e:
            raise AgentError(f"anthropic auth failed (check your API key): {_scrub_key(e)}")
        except APIError as e:
            raise AgentError(f"anthropic API error: {_scrub_key(e)}")
    else:
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            raise AgentError("google-genai SDK not installed")
        client = genai.Client(api_key=api_key)
        tools_arg = []
        if allowed_tools and any(t in ("WebSearch", "web_search") for t in allowed_tools):
            tools_arg.append(types.Tool(google_search=types.GoogleSearch()))
            try:
                tools_arg.append(types.Tool(url_context=types.UrlContext()))
            except AttributeError:
                pass
        # Try requested model; on 503/UNAVAILABLE fall back to gemini-2.0-flash
        # (more stable backend) before giving up. Each attempt has its own
        # token budget; if the first emits any text we don't fall back.
        # Default is gemini-2.0-flash (free tier 15 RPM / 1M TPM, fewer 503s).
        # On overload, try 2.5-flash as the slower-but-higher-quality fallback.
        candidates = [model or os.environ.get("GEMINI_MODEL", "gemini-2.0-flash"),
                      "gemini-2.5-flash"]
        seen = set()
        last_err: Exception | None = None
        for mdl in candidates:
            if mdl in seen:
                continue
            seen.add(mdl)
            emitted_any = False
            grounding_announced = False
            try:
                for chunk in client.models.generate_content_stream(
                    model=mdl,
                    contents=prompt,
                    config=types.GenerateContentConfig(
                        max_output_tokens=max_tokens,
                        tools=tools_arg or None,
                    ),
                ):
                    # Announce grounding tool use the first time we see
                    # web_search queries in the chunk metadata, so the UI
                    # can show "Searching PubMed for ...".
                    if not grounding_announced:
                        try:
                            gm = (chunk.candidates[0].grounding_metadata
                                  if getattr(chunk, "candidates", None) else None)
                            queries = getattr(gm, "web_search_queries", None) if gm else None
                            if queries:
                                grounding_announced = True
                                yield ("tool", {
                                    "name": "google_search",
                                    "input": {"query": " / ".join(queries[:3])},
                                })
                        except (AttributeError, IndexError, TypeError):
                            pass
                    if chunk.text:
                        emitted_any = True
                        yield ("text", chunk.text)
                return                      # ok, stream finished
            except Exception as e:
                last_err = e
                emsg = str(e).lower()
                if emitted_any:
                    # Mid-stream failure: re-raise (can't switch mid-flight).
                    if "api key" in emsg or "permission" in emsg or "auth" in emsg or "invalid argument" in emsg:
                        raise AgentError(f"gemini auth failed: {_scrub_key(e)}")
                    if "quota" in emsg or "rate" in emsg or "exceeded" in emsg or "resource_exhausted" in emsg:
                        raise AgentError(f"gemini quota/rate exceeded: {_scrub_key(e)}")
                    raise AgentError(f"gemini stream interrupted: {_scrub_key(e)}")
                # Pre-stream failure: decide whether to fall back to next model.
                if "503" in emsg or "unavailable" in emsg or "overloaded" in emsg:
                    log.warning("gemini %s 503 unavailable, trying next model", mdl)
                    continue
                if "api key" in emsg or "permission" in emsg or "auth" in emsg or "invalid argument" in emsg:
                    raise AgentError(f"gemini auth failed (check your Google API key): {_scrub_key(e)}")
                if "quota" in emsg or "rate" in emsg or "exceeded" in emsg or "resource_exhausted" in emsg:
                    raise AgentError(f"gemini quota/rate exceeded: {_scrub_key(e)}")
                # Unknown error — don't loop, surface it.
                raise AgentError(f"gemini API error: {_scrub_key(e)}")
        # All candidates exhausted with 503.
        raise AgentError(f"gemini service unavailable (tried {candidates}). "
                          f"Either retry in a minute or switch to an Anthropic key. "
                          f"Last error: {last_err}")


# Backwards-compat alias for legacy call sites (will be removed after audit).
def _run_claude(prompt: str, *, api_key: str | None = None,
                allowed_tools: list[str] | None = None,
                timeout: int = 480, effort: str = "low") -> str:
    if not api_key:
        raise AgentError("BYOK Anthropic API key required (none supplied)")
    return _call_anthropic_sdk(prompt, api_key=api_key,
                                allowed_tools=allowed_tools, timeout=timeout)


def _strip_json(text: str) -> dict | list | None:
    """Extract the first balanced JSON object or array from arbitrary text.

    Scans for `{` or `[`, then bracket-counts (string-aware, escape-aware)
    until the matching close. Stops at the first parseable JSON value,
    ignoring any trailing markdown / source links the agent may append.
    """
    if not text:
        return None
    # Strip markdown code fences if the whole response is wrapped in them.
    # Gemini in particular frequently wraps JSON in ```json ... ``` and the
    # previous non-greedy regex broke on JSON containing nested ] or }.
    fence_strip = re.sub(r"^\s*```(?:json|JSON)?\s*\n?", "", text)
    fence_strip = re.sub(r"\n?\s*```\s*$", "", fence_strip)
    if fence_strip != text:
        # Fenced response: try parsing the stripped body directly first.
        try:
            return json.loads(fence_strip.strip())
        except json.JSONDecodeError:
            text = fence_strip   # fall through to bracket counter

    for start in range(len(text)):
        open_ch = text[start]
        if open_ch not in "{[":
            continue
        close_ch = "}" if open_ch == "{" else "]"
        depth = 0
        in_str = False
        escape = False
        i = start
        while i < len(text):
            c = text[i]
            if escape:
                escape = False
            elif in_str:
                if c == "\\":
                    escape = True
                elif c == '"':
                    in_str = False
            else:
                if c == '"':
                    in_str = True
                elif c == open_ch:
                    depth += 1
                elif c == close_ch:
                    depth -= 1
                    if depth == 0:
                        cand = text[start:i + 1]
                        try:
                            return json.loads(cand)
                        except json.JSONDecodeError:
                            break
            i += 1
    return None


def synthesize_via_claude(citations: list, ctx: dict, *,
                          anthropic_api_key: str,
                          timeout: int = 90) -> str:
    """LLM synthesis WITHOUT web search — only summarises retrieved citations.
    Requires BYOK Anthropic key from caller."""
    if not citations:
        return ""
    payload = json.dumps([{"title": c.get("title"), "year": c.get("year"),
                           "journal": c.get("journal"), "authors": c.get("authors"),
                           "pmid": c.get("pmid")} for c in citations[:10]], indent=2)
    prompt = (
        "You are summarizing literature for a cardiac genetics researcher. "
        f"Disease: {ctx.get('disease','')}, cell state: {ctx.get('target_cs','')}. "
        "Given the following PubMed citations (titles + metadata only — do NOT "
        "search the web, only synthesize what is reasonable from the titles), "
        "write a 2-4 sentence synthesis. Be cautious — if a title is ambiguous "
        "or off-topic, say so. Avoid speculation. Return plain text only.\n\n"
        f"Citations JSON:\n{payload}"
    )
    return _call_anthropic_sdk(prompt, api_key=anthropic_api_key,
                                allowed_tools=None, timeout=timeout)


def agent_search(type_: str, ctx: dict, *,
                 anthropic_api_key: str,
                 timeout: int = 300, max_results: int = 6) -> dict:
    """Run Claude with WebSearch tool to gather literature.
    Requires BYOK Anthropic key.

    The agent is instructed to return JSON with **real** PMIDs/DOIs only;
    we cross-check returned PMIDs against PubMed via the local fetcher to
    catch any hallucinated identifiers. Items that fail verification are
    marked unverified.
    """
    if type_ == "snp":
        target = ctx.get("rsid") or f"chr{ctx.get('chr')}:{ctx.get('pos')}"
        kind_phrase = f"GWAS variant {target}"
    elif type_ in ("tf", "gene"):
        target = ctx.get("symbol", "")
        kind_phrase = f"{type_.upper()} {target}"
    else:
        raise ValueError(f"unknown type {type_}")

    disease = ctx.get("disease", "")
    cs = ctx.get("target_cs", "")

    prompt = (
        f"You are a cardiac genetics literature curator. Use **at most 3 "
        f"WebSearch calls** to find papers about **{kind_phrase}** in the "
        f"context of **{disease}** (cardiac cell state: {cs}). "
        f"Prefer the PubMed website (pubmed.ncbi.nlm.nih.gov) for "
        f"deterministic PMIDs. After 3 searches, stop and synthesize. "
        f"\n\nReturn ONLY a JSON object (no markdown fence, no prose). "
        f"For every citation, provide the REAL PMID or DOI exactly as "
        f"it appears on PubMed. Do NOT fabricate identifiers. "
        f"If you cannot find a PMID, leave it null but include the DOI.\n"
        f"Schema:\n"
        f"{{\n"
        f'  "citations": [\n'
        f'    {{"pmid": str|null, "doi": str|null, "title": str, '
        f'"authors": str, "year": int, "journal": str, '
        f'"key_finding": str, "relevance": str}}\n'
        f'  ],\n'
        f'  "summary": "2-4 sentence synthesis tailored to {disease}/{cs}",\n'
        f'  "broader_context": "1-2 sentences on preprints/news/reviews '
        f'not in PubMed (optional, can be empty string)"\n'
        f"}}\n\n"
        f"Max {max_results} citations. Skip papers that do not specifically "
        f"mention {target}. Be terse — do not write essays in the fields."
    )
    raw = _call_anthropic_sdk(prompt, api_key=anthropic_api_key,
                                allowed_tools=["WebSearch"],
                                timeout=timeout)
    obj = _strip_json(raw)
    if not isinstance(obj, dict):
        return {"_raw": raw, "_parse_error": True, "citations": []}
    return obj


# Mapping from cell_state -> broader cell_type (used to expand user phrases
# like "ValveFB" / "Valve Fibroblasts" into the full set of cell_states the
# viewer carries).
STATE_TO_TYPE: dict[str, str] = {
    "VentricularCardiomyocytesLeft": "VentricularCardiomyocytes",
    "VentricularCardiomyocytesLeftStressed": "VentricularCardiomyocytes",
    "VentricularCardiomyocytesRight": "VentricularCardiomyocytes",
    "VentricularCardiomyocytesSeptal": "VentricularCardiomyocytes",
    "VentricularCardiomyocytesRightStressed": "VentricularCardiomyocytes",
    "ValveFibroblastsImmune": "ValveFibroblasts",
    "ValveFibroblastsCardiacSkeleton": "ValveFibroblasts",
    "ValveFibroblastsBaseline": "ValveFibroblasts",
    "SmoothMuscleCellsPericytesIntermediate": "SmoothMuscleCells",
    "SmoothMuscleCellsGreatArtery": "SmoothMuscleCells",
    "SmoothMuscleCellsCoronaryArtery": "SmoothMuscleCells",
    "SmoothMuscleCellsArterial": "SmoothMuscleCells",
    "SmoothMuscleCellsAtrial": "SmoothMuscleCells",
    "PericytesVentricular": "Pericytes",
    "PericytesAtrial": "Pericytes",
    "NeuralCells": "NeuralCells",
    "NeuralCellsNerveAssociated": "NeuralCells",
    "SchwannCells": "NeuralCells",
    "NerveFibroblastsEndoneurial": "NerveFibroblasts",
    "NerveFibroblastsPerineurial": "NerveFibroblasts",
    "MacrophagesLYVE1pos": "MyeloidCells",
    "Mono_Macrophages": "MyeloidCells",
    "MacrophagesCycling": "MyeloidCells",
    "MacrophagesCXCL8pos": "MyeloidCells",
    "DendriticCells": "MyeloidCells",
    "MacrophagesLYVE1posAGBL4pos": "MyeloidCells",
    "Monocytes": "MyeloidCells",
    "MacrophagesLipidAssociated": "MyeloidCells",
    "MastCells": "MastCells",
    "TCellsCD8pos": "LymphoidCellsNonB",
    "TCellsCD4pos": "LymphoidCellsNonB",
    "NaturalKillerCellsCD56hi": "LymphoidCellsNonB",
    "NaturalKillerCellsCD16hi": "LymphoidCellsNonB",
    "TCellsCD4posRegulatory": "LymphoidCellsNonB",
    "BPlasmaCells": "LymphoidCellsB",
    "BCells": "LymphoidCellsB",
    "LymphaticEndothelialCells": "LymphaticEndothelialCells",
    "FibroblastsAtrial": "Fibroblasts",
    "FibroblastsPCOLCE2": "Fibroblasts",
    "FibroblastsCD44": "Fibroblasts",
    "FibroblastsAPOD": "Fibroblasts",
    "FibroblastsVascular": "Fibroblasts",
    "FibroblastsVentricular": "Fibroblasts",
    "FibroblastsActivated": "Fibroblasts",
    "FibroblastsCXCL8": "Fibroblasts",
    "EpicardialCells": "EpicardialCells",
    "EndothelialCellsCapillary": "EndothelialCells",
    "EndothelialCellsVenous": "EndothelialCells",
    "EndothelialCellsArterial": "EndothelialCells",
    "EndothelialCellsArterialLarge": "EndothelialCells",
    "EndothelialCellsNOVA1Capillary": "EndothelialCells",
    "EndothelialCellsNOVA1Arterial": "EndothelialCells",
    "EndothelialCellsNOVA1Venous": "EndothelialCells",
    "EndocardialCells": "EndocardialCells",
    "VentricularConductionSystemDistal": "CardiacConductionSystem",
    "VentricularConductionSystemProximal": "CardiacConductionSystem",
    "AtrialConductionSystem": "CardiacConductionSystem",
    "PacemakerCells": "CardiacConductionSystem",
    "MyocardialSleeveCells": "CardiacConductionSystem",
    "AtrialCardiomyocytesRight": "AtrialCardiomyocytes",
    "AtrialCardiomyocytesLeft": "AtrialCardiomyocytes",
    "AtrialCardiomyocytesLeftStressed": "AtrialCardiomyocytes",
    "AtrialCardiomyocytesRightStressed": "AtrialCardiomyocytes",
    "Adipocytes": "Adipocytes",
    "AdipocytesAKR1C1": "Adipocytes",
}


def _type_to_states() -> dict[str, list[str]]:
    """Reverse mapping (cell_type -> sorted list of cell_states)."""
    out: dict[str, list[str]] = {}
    for state, ctype in STATE_TO_TYPE.items():
        out.setdefault(ctype, []).append(state)
    for ctype in out:
        out[ctype].sort()
    return out


CHAT_SYSTEM_PROMPT = """You are an interactive assistant embedded in a Cytoscape eGRN viewer
for cardiac GWAS analysis. The user is exploring a concentric network: outer
ring = top-N seed genes, middle ring = peaks (regulatory regions), inner
ring = TFs. Each peak→gene and TF→peak edge carries a set of source
cell_types (which cardiac cell types the eGRN connection was inferred in).

CITATION POLICY — **INLINE, BUT MODEST IN COUNT**.
You will be given a **VERIFIED_CITATIONS** list (real PubMed entries
retrieved deterministically before you ran). The new rule for `message`:

  * **Cite every factual / mechanistic claim inline** with one PMID.
    A "claim" = any sentence that asserts biology, drug effect, GWAS
    finding, mechanism, statistic, association, etc. that a reader
    could plausibly want to verify. Example shape:
      "TBX5 inactivation causes spontaneous atrial fibrillation
       (PMID 27009272) and disrupts the atrial enhancer landscape
       (PMID 30385751)."
    Each clause that makes an independent claim gets its own PMID.
    Two claims in one sentence → two PMIDs, separated by ;.
    Pure descriptions of the current viewer state (e.g. "TBX5 is a
    seed gene with GWAS-z 4304") are NOT claims and need no citation.
  * Cite ONLY entries from VERIFIED_CITATIONS by their PMID
    (e.g. "PMID 27009272"). The post-processor turns these into
    coloured cards with the title + journal + year.
  * Do NOT fabricate PMIDs, DOIs, or titles. Do NOT pull citations
    from your training data.
  * **Repeat citations are fine** — if the same PMID supports several
    claims, cite it each time. The user wants to see exactly which
    paper backs which statement, not a bibliography at the end.
  * **Target 2-5 unique citations per response.** Inline citations on
    the 2-5 most critical claims; for secondary supporting points you
    may either re-use one of those PMIDs or leave the sentence
    uncited. Do NOT try to cite every clause — that bloats the reply
    and wastes the tool-call budget. The post-processor dedups for
    the side panel.
  * In the JSON `citations` field, include one entry per UNIQUE PMID
    you cited, each with `pmid`, `key_finding` (one-line takeaway).
  * **If VERIFIED_CITATIONS is empty** OR no entry covers a specific
    claim, **YOU MUST actively use WebSearch / google_search (up to
    15 calls)** to find real PubMed PMIDs for the user's question.
    Put those PMIDs inline in the `message` AND in the JSON
    `citations` field — the backend will re-verify each one via NCBI
    esummary. Verified-but-off-grounding hits get an orange tag in
    the UI (positive signal: you found something the deterministic
    NCBI prefetch missed). Hallucinated/unreachable PMIDs get
    flagged red.

  * **MANDATORY SELF-VERIFICATION (after WebSearch, before citing)**:
    For EVERY PMID you intend to cite that did NOT come from
    VERIFIED_CITATIONS, you MUST re-access the PubMed page yourself
    to confirm it actually resolves:
      1. Call `web_fetch` (Anthropic) / `url_context` (Gemini) on
         `https://pubmed.ncbi.nlm.nih.gov/<pmid>/` for each candidate
         PMID.
      2. Confirm: (a) the page returns content (not 404 / "not
         found"), and (b) the title on the page is on-topic for the
         claim you are making.
      3. If either check fails → DROP that PMID. Do NOT cite it.
         Try another candidate or rephrase as a hypothesis.
      4. Only PMIDs that survived this self-fetch step go into the
         inline citation AND `citations[]`.
    Budget your tool calls: roughly half for search, half for fetch.
    The user has explicitly asked for this re-access step — skipping
    it is a policy violation.
  * **FORBIDDEN FALLBACK**: Do NOT tell the user "search PubMed
    yourself"; do NOT paste a `pubmed.ncbi.nlm.nih.gov/?term=...`
    URL; do NOT say "VERIFIED_CITATIONS is empty so no papers can be
    cited". Those are punts. The user has explicitly said the agent
    should do the PubMed search and return verified PMIDs. Only if
    WebSearch ALSO returns nothing relevant may you phrase a claim
    as a hypothesis — and you must say so explicitly ("WebSearch
    returned no directly supporting paper; hypothesis only").
  * Do NOT cite from training data without WebSearch confirmation.

The graph_state JSON contains the COMPLETE current subgraph:
  • seed_genes — outer ring, full list with chr / gwas_z / target_de_z / specificity
  • peaks — middle ring, chr / pos / gwas_z / snp_overlap / snp_rsids
  • tfs — inner ring, gwas_z / specificity
  • snps — every SNP overlapping a visible peak, with OT credible set L2G + GTEx
  • edges — every visible edge as {s, t, k, cts}
        s = source id, t = target id, k = 1 for TF2R / 0 for R2G,
        cts = array of indices into cell_type_table (e.g. [0,3] means edge
        is supported by cell_type_table[0] and cell_type_table[3]).
  • cell_type_table — array of cell_type strings indexed by edges[].cts.
You can therefore answer questions like "which TFs are cell-type-specific to
EndothelialCells?" by scanning edges + cell_type_table.

You can:
  1. EXPLAIN concepts, paper findings, gene/TF biology to the user.
  2. RETURN ACTIONS that modify the viewer to make the user's question
     visible: highlight nodes, highlight paths, change top-N, focus on a
     seed, trigger PubMed search, etc.

CELL TYPE vs CELL STATE — IMPORTANT.
The viewer carries 65+ cell STATES (e.g. ValveFibroblastsImmune,
ValveFibroblastsBaseline, ValveFibroblastsCardiacSkeleton), which roll
up to broader cell TYPES (e.g. all three above → ValveFibroblasts).
A CELL_TYPE_TO_STATES mapping is provided in the prompt below.
When the user mentions a cell TYPE (or a common abbreviation like
"ValveFB", "SMC", "VCM", "ACM", "EC"), do NOT silently pick one state.

THE MOST COMMON QUERY PATTERN IS: "gene/TF in cell_type" (e.g.
"TBX5 in ValveFB", "PALMD in EC", "PDGFD in SMC"). For this pattern,
your workflow MUST be:
  Step 1. Expand the cell_type the user mentioned into ALL related
          cell_STATES via CELL_TYPE_TO_STATES.
  Step 2. Search ALL those cell_states (using CROSS_PAYLOAD_GENE_HITS
          when present) for the gene/TF, ordered by disease priority:
            (a) Search the CURRENT graph_state.disease FIRST.
            (b) Then the other diseases (CAD / AF / AVS) in turn.
  Step 3. Report the consolidated finding to the user — list which
          (disease × cell_state) combinations show signal, including
          inner_min / GWAS-z / DE-z values from CROSS_PAYLOAD_GENE_HITS.
  Step 4. Pick the SINGLE BEST hit (highest inner_min) and emit
          set_disease + set_target_cs to navigate there. If the current
          disease already has signal, prefer to stay in the current
          disease even if another disease has a slightly stronger hit.
  Step 5. After switching, emit highlight_nodes for the gene/TF so it
          is visible on the new network.

Common short aliases users use:
   "ValveFB" / "valve fibroblast"   → ValveFibroblasts*
   "SMC" / "smooth muscle"          → SmoothMuscleCells*
   "VCM" / "ventricular CM"         → VentricularCardiomyocytes*
   "ACM" / "atrial CM"              → AtrialCardiomyocytes*
   "EC" / "endothelial"             → EndothelialCells* (incl NOVA1 subset)
   "Myeloid" / "macrophage"         → MyeloidCells* (macrophages + DC + mono)
   "CCS" / "conduction system"      → CardiacConductionSystem*

PROACTIVE ROUTING — DO NOT ASK FOR PERMISSION.
When the user's question is grounded in a disease or cell_state that does
NOT match graph_state.disease / graph_state.target_cs, you MUST:
  (a) Emit a set_disease action FIRST if the disease differs (CAD/AF/AVS)
  (b) Emit a set_target_cs action to a sensible match from
      graph_state.manifest_target_cs_list — pick the closest cell_state for
      the user's intent, e.g.:
         "SMC" / "vascular smooth muscle" / "coronary artery SMC"
              → SmoothMuscleCellsCoronaryArtery
         "valve fibroblast" → ValveFibroblastsCardiacSkeleton
         "endothelial" → EndothelialCellsArterial (or specific subtype)
         "sleeve cells" / "pulmonary vein sleeve"
              → MyocardialSleeveCells
  (c) THEN write your message. The viewer applies switch-actions before
      rendering highlights, so node IDs you mention will refer to the NEW
      subgraph.
Never tell the user "the current viewer is loaded for X" — just switch to
the relevant context first and answer in that context. Mention briefly at
the end which disease/cell_state you routed to (e.g. "(switched to CAD /
SmoothMuscleCellsCoronaryArtery)") so the user knows.

CRITICAL — ACTIONS ARE EXECUTED, PROSE IS NOT.
Anything you write in `message` is just text shown to the user. The viewer
only changes state when you emit the corresponding entry in `actions`.
If your message claims you "switched", "updated", "navigated to",
"focused on", "highlighted", "set top-N", etc., you MUST include the
matching action object in `actions`. Examples:
  * Message says "Switched to AF / AtrialCardiomyocytes" → REQUIRES
      `{"type":"set_disease","args":{"disease":"AF"}}` AND
      `{"type":"set_target_cs","args":{"cs":"AtrialCardiomyocytes"}}` in actions.
  * Message says "Highlighted TBX5 and PLN" → REQUIRES
      `{"type":"highlight_nodes","args":{"ids":["TBX5","PLN"]}}` in actions.
A message that promises a state change without the matching action is a
bug — the viewer will not update and the user sees an inconsistent UI.

Use WebSearch + WebFetch tools when the user asks for literature context.
ALWAYS reference REAL PMIDs / paper titles — never fabricate. When you cite a
paper, include a clickable URL like https://pubmed.ncbi.nlm.nih.gov/<PMID>/
directly in the `message` field so the user can open it.

Respond with ONE JSON object (no markdown fence, no prose before/after):
{
  "message": "<your text response. Refer to papers by PMID like 'PMID 27009272' — the backend will replace these with clickable verified citation cards.>",
  "actions": [<list of action objects, empty list OK>],
  "citations": [
    {"pmid": "27009272", "key_finding": "<1-line takeaway>"}
  ]
}

Action schema:
  {"type": "highlight_nodes", "args": {"ids": ["TBX5", "PLN", "chr4:..."]}}
      → adds magenta outline + un-fades these nodes
  {"type": "highlight_path", "args": {"from": "TBX5", "to": "PLN"}}
      → traces shortest path in the eGRN, highlights nodes + edges along it
  {"type": "highlight_cell_type_edges", "args": {"cell_type": "AtrialCardiomyocytes"}}
      → highlights every edge whose cell_types contains this label
  {"type": "set_view", "args": {"top_n": 50, "peaks_per_seed": 6,
                                  "tfs_per_peak": 4,
                                  "edge_highlight_cell_types": ["CardiacConductionSystem"]}}
      → updates viewer parameters and rebuilds the network (parameters are
        clamped: top_n ≤ 200, peaks_per_seed ≤ 10, tfs_per_peak ≤ 6)
  {"type": "set_disease", "args": {"disease": "CAD"}}
      → switches the entire payload to a different EFO/disease, reloads
        the DE-z matrix and the target-cs payload. Valid values: see
        graph_state.manifest_disease_list (typically AF / CAD / AVS).
        Common phrasing hints:
          AF  ← "atrial fibrillation", "AFib", "AF"
          CAD ← "coronary artery disease", "coronary disease", "CAD",
                "myocardial infarction"
          AVS ← "aortic valve stenosis", "calcific aortic stenosis", "AVS"
        Emit this BEFORE set_target_cs when both change.
  {"type": "set_target_cs", "args": {"cs": "VentricularCardiomyocytesLeft"}}
      → switches the network's target cell_state (reloads the payload).
        Use this whenever the user asks about a *different* cell state
        than graph_state.target_cs. cs MUST be one of the values in
        graph_state.manifest_target_cs_list — choose the closest match if
        the user uses an informal name (e.g. "sleeve cells" →
        MyocardialSleeveCells).
  {"type": "set_ref_cs", "args": {"cs_list": ["PacemakerCells",
                                                "AtrialConductionSystem"]}}
      → switches specificity contrast to "vs selected" with these cs as the
        reference set. Requires ≥2 entries (the target_cs itself is
        excluded automatically).
  {"type": "focus_on_seed", "args": {"id": "PLN"}}
      → hides all seeds except this one and its peak/TF neighborhood
  {"type": "search_pubmed", "args": {"type": "tf|gene|snp",
                                       "key": "TBX5", "mode": "pubmed|agent|hybrid"}}
      → triggers a literature search for this entity; result will be shown
        in the detail panel when the user clicks
  {"type": "reset_highlights", "args": {}}
      → clears all chat-driven highlights

Use IDs EXACTLY as they appear in graph_state. For peaks the ID is the
chrN:start-end string. Be conservative: 3-10 nodes in highlight_nodes is
plenty; avoid spamming.

IMPORTANT: Whenever your `message` text mentions any specific gene symbol,
TF, peak region (chrN:start-end), or rsid that EXISTS in graph_state,
ALSO include that ID in a `highlight_nodes` action so the user sees it
on the network. The viewer will additionally auto-extract obvious matches
from your prose, but explicit actions are preferred.
"""


def build_chat_prompt(message: str, graph_state: dict,
                      history: list, max_chars: int = 120000,
                      grounding_citations: list | None = None,
                      snp_cross_hits: list | None = None,
                      gene_cross_hits: list | None = None) -> str:
    """Assemble a single-shot prompt for `claude --print`.

    The viewer ships the full visible subgraph as graph_state. We forward
    it verbatim — Claude has a 200K-token context window so even 50-100 KB
    of JSON fits comfortably. We only truncate if the payload exceeds a
    sane safety cap.
    """
    state_json = json.dumps(graph_state, separators=(",", ":"))
    if len(state_json) > max_chars:
        # Truncate gracefully (drop large arrays from the end) by emitting
        # a header summary and the heaviest arrays first if size is excessive.
        truncated = {**graph_state}
        for key in ("edges", "peaks", "tfs", "snps", "seed_genes"):
            if len(json.dumps(truncated, separators=(",", ":"))) <= max_chars:
                break
            if isinstance(truncated.get(key), list):
                truncated[key] = truncated[key][: max(50, len(truncated[key]) // 2)]
        truncated["_truncated"] = True
        state_json = json.dumps(truncated, separators=(",", ":"))

    parts = [CHAT_SYSTEM_PROMPT,
             "\n## CELL_TYPE_TO_STATES (use to expand user phrases)\n"
             f"```json\n{json.dumps(_type_to_states(), indent=1)}\n```\n",
             f"\n## Current viewer state\n```json\n{state_json}\n```\n"]
    if grounding_citations:
        cit_compact = json.dumps([{
            "pmid": c.get("pmid"), "title": c.get("title"),
            "authors": c.get("authors"), "year": c.get("year"),
            "journal": c.get("journal"), "search_term": c.get("search_term"),
        } for c in grounding_citations], indent=1)
        parts.append("\n## VERIFIED_CITATIONS (cite ONLY these PMIDs)\n"
                     f"```json\n{cit_compact}\n```\n")
    else:
        parts.append("\n## VERIFIED_CITATIONS\n(empty — no grounded "
                     "papers retrieved; do not cite specific PMIDs)\n")
    # Cross-payload gene/TF hits — user asked about a gene that exists in
    # OTHER (disease, cs) payloads. Surface them so agent can answer or route.
    if gene_cross_hits:
        out_of_view = [h for h in gene_cross_hits if not h["in_current_view"]]
        if out_of_view:
            payload = json.dumps(gene_cross_hits, indent=1)
            parts.append(
                "\n## CROSS_PAYLOAD_GENE_HITS — IMPORTANT\n"
                f"The user mentioned gene/TF symbol(s) that exist in OTHER "
                f"(disease, cs) payloads (not the current view "
                f"{graph_state.get('disease')}/{graph_state.get('target_cs')}). "
                f"Each hit shows where that gene is a seed (high inner_min of "
                f"GWAS-z × target DE-z) or a TF. Use this to ANSWER the user "
                f"directly — you can already report which cs has significant "
                f"expression without seeing it in the current view. If they ask "
                f"to see it, emit set_disease + set_target_cs to the top hit:\n"
                f"```json\n{payload}\n```\n"
            )

    # Cross-payload SNP hits — user mentioned a SNP that lives outside the
    # currently-loaded subgraph. Agent should switch to the best (disease, cs)
    # combo before answering.
    if snp_cross_hits:
        out_of_view = [h for h in snp_cross_hits if not h["in_current_view"]]
        in_view = [h for h in snp_cross_hits if h["in_current_view"]]
        if out_of_view:
            payload = json.dumps(snp_cross_hits, indent=1)
            parts.append(
                "\n## CROSS_PAYLOAD_SNP_HITS — IMPORTANT\n"
                f"The user mentioned SNP(s) that exist in OTHER (disease, cs) "
                f"payloads. Current view: disease={graph_state.get('disease')} "
                f"target_cs={graph_state.get('target_cs')}. Lookup result:\n"
                f"```json\n{payload}\n```\n"
                "ACTION REQUIRED: pick the most biologically relevant hit "
                "(preference order: matching disease keyword in user message > "
                "most CAD-relevant cell_state for CAD SNPs / etc.). "
                "Emit set_disease + set_target_cs FIRST, then answer in the "
                "new context. If the SNP appears in multiple cs of the same "
                "disease, prefer the one with the strongest biological "
                "rationale for the user's question.\n"
            )
        elif in_view:
            parts.append(
                "\n## SNP IN CURRENT VIEW\nThe user's SNP is already loaded "
                "in the current payload — no payload switch needed.\n"
            )
    if history:
        parts.append("\n## Recent chat history\n")
        for m in history[-6:]:
            role = m.get("role", "user")
            content = (m.get("content", "") or "")[:600]
            parts.append(f"**{role}**: {content}\n")
    parts.append(f"\n## User asks now\n{message}\n\n"
                 "Return your JSON response. Remember: cite only PMIDs "
                 "from VERIFIED_CITATIONS, in the `citations` array, with "
                 "a 1-line `key_finding`. Inline mentions in the message "
                 "should reference them as 'PMID 12345' (no markdown link "
                 "— the viewer renders the citation card).")
    return "\n".join(parts)


def _extract_search_terms(message: str, graph_state: dict, max_terms: int = 6) -> list[str]:
    """Pull obvious search anchors from the user's question."""
    terms: list[str] = []
    seen = set()
    def _add(t):
        if t and t not in seen: seen.add(t); terms.append(t)

    # rsids
    for m in re.finditer(r"\brs\d+\b", message):
        _add(m.group(0))
    # cell_state names (camel-cased so exact match safe)
    for cs in graph_state.get("all_cs_in_de_z") or []:
        if cs in message:
            _add(cs)
    # gene/TF symbols in the message that exist as nodes in graph_state
    node_ids: set[str] = set()
    for arr in (graph_state.get("seed_genes", []), graph_state.get("tfs", [])):
        for n in arr or []:
            nid = n.get("id")
            if nid: node_ids.add(nid)
    # Case-insensitive symbol match (uppercase before checking the index).
    for m in re.finditer(r"\b[A-Za-z][A-Za-z0-9-]{1,}\b", message):
        sym = m.group(0).upper()
        if sym in node_ids:
            _add(sym)
    return terms[:max_terms]


def fetch_grounding_citations(terms: list[str], disease_id: str, *,
                              api_key: str | None = None,
                              email: str | None = None,
                              max_per_term: int = 3) -> list[dict]:
    """Pre-fetch real PubMed citations to ground the agent — STRICT mode."""
    disease_term = DISEASE_TERMS.get(disease_id, disease_id)
    out: list[dict] = []
    seen_pmids: set[str] = set()
    for t in terms:
        # Skip very common one-letter / disease tokens that would over-match
        if len(t) < 3 and not t.startswith("rs"):
            continue
        q = f"({t}) AND {disease_term}" if disease_term else t
        try:
            res = pubmed_search(q, max_results=max_per_term,
                                api_key=api_key, email=email)
        except Exception as e:
            log.warning("grounding search failed for %s: %s", t, e)
            continue
        for c in res.get("citations", []):
            pmid = c.get("pmid")
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                out.append({**c, "search_term": t})
        time.sleep(0.34)
    return out


def verify_citations_in_response(obj: dict, grounding: list[dict], *,
                                  api_key: str | None = None,
                                  email: str | None = None) -> dict:
    """Cross-check every PMID in obj.citations against grounding + esummary.
    Annotate each entry with `verified=True/False` and `off_grounding=True`
    if real but not in our pre-fetched set."""
    grounding_pmids = {c["pmid"]: c for c in grounding if c.get("pmid")}
    cits = obj.get("citations") or []
    if not isinstance(cits, list):
        cits = []
    # Verify any PMIDs not in grounding via esummary
    unknown = [c["pmid"] for c in cits
               if c.get("pmid") and str(c["pmid"]) not in grounding_pmids]
    real_titles: dict[str, dict] = {}
    if unknown:
        try:
            raw = _http_get(NCBI_ESUMMARY, {
                "db": "pubmed", "id": ",".join(str(p) for p in unknown),
                "retmode": "json",
            }, api_key=api_key, email=email)
            summ = json.loads(raw).get("result", {})
            for pid in unknown:
                if str(pid) in summ:
                    real_titles[str(pid)] = summ[str(pid)]
        except Exception as e:
            log.warning("verification esummary failed: %s", e)
    # Pre-compute confidence inputs (user votes + agent recurrence) so every
    # output row can carry a confidence band.
    cred_agg = _credible_aggregate()
    out = []
    for c in cits:
        pmid = str(c.get("pmid") or "")
        if not pmid:
            out.append({**c, "verified": False, "hallucinated": True}); continue
        if pmid in grounding_pmids:
            g = grounding_pmids[pmid]
            row = {**g, "verified": True,
                    "key_finding": c.get("key_finding") or g.get("title")}
        elif pmid in real_titles:
            r = real_titles[pmid]
            row = {
                "pmid": pmid,
                "title": r.get("title", "").strip("."),
                "authors": ", ".join(a.get("name", "") for a in (r.get("authors") or [])[:3]) +
                            (" et al." if (r.get("authors") or []) and len(r["authors"]) > 3 else ""),
                "journal": r.get("source"),
                "year": int(r.get("pubdate", "0").split()[0]) if r.get("pubdate") else None,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "verified": True,
                "off_grounding": True,
                "key_finding": c.get("key_finding"),
            }
        else:
            row = {
                "pmid": pmid,
                "title": c.get("title"),
                "verified": False,
                "hallucinated": True,
                "key_finding": c.get("key_finding"),
            }
        # Attach user-vote counts + recurrence-derived confidence band.
        agg = cred_agg.get(pmid, {})
        up = agg.get("up", 0); down = agg.get("down", 0)
        row["credible_up"] = up
        row["credible_down"] = down
        row["confidence"] = _credible_confidence(up, down, c.get("_recurrence", 0))
        # `from_cache` / `curated` are set on the grounding entry when we
        # inject from the agent-references cache. Propagate so the UI can
        # show the 🧠 curated badge + inline summary.
        g_match = grounding_pmids.get(pmid, {})
        if isinstance(c, dict) and c.get("from_cache"):
            row["from_cache"] = True
        elif g_match.get("from_cache"):
            row["from_cache"] = True
        cur = _load_curated(pmid)
        if g_match.get("curated") or cur:
            row["curated"] = True
        if cur:
            row["curated_summary"] = cur.get("summary")
            row["curated_key_findings"] = cur.get("key_findings", [])
        out.append(row)
    return out


def extend_grounding_with_cache(grounding: list[dict], *,
                                 disease: str, cs: str,
                                 genes: list[str] | None = None,
                                 max_add: int = 8) -> list[dict]:
    """Append high-confidence cached agent references (matching the given
    context) to the deterministic NCBI grounding list. Dedups by pmid;
    NCBI prefetch entries always win. Returns the extended list."""
    cached = _agent_refs_by_context(disease, cs, genes)
    if not cached:
        return grounding
    cred_agg = _credible_aggregate()
    seen = {c.get("pmid") for c in grounding if c.get("pmid")}
    added = 0
    for pmid, r in sorted(cached.items(),
                          key=lambda kv: kv[1]["recurrence"], reverse=True):
        if pmid in seen:
            continue
        agg = cred_agg.get(pmid, {"up": 0, "down": 0})
        conf = _credible_confidence(agg["up"], agg["down"], r["recurrence"])
        if conf != "high":
            continue   # only auto-inject high-confidence cached refs
        # If we have a curated AI summary for this PMID (vote-curated
        # prior knowledge), prefer it over the raw key_finding so the
        # model sees the most informative snippet.
        cur = _load_curated(pmid)
        if cur:
            entry_kf = cur.get("summary") or r["key_finding"]
            kfs = cur.get("key_findings") or []
            if kfs:
                entry_kf = (entry_kf + " · " +
                              " ; ".join(str(k) for k in kfs[:3]))[:700]
        else:
            entry_kf = r["key_finding"]
        grounding.append({
            "pmid": pmid,
            "title": r["title"],
            "journal": r["journal"],
            "year": r["year"],
            "key_finding": entry_kf,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "from_cache": True,
            "curated": bool(cur),
            "_recurrence": r["recurrence"],
        })
        seen.add(pmid); added += 1
        if added >= max_add:
            break
    if added:
        log.info("grounding extended with %d high-conf cached refs (disease=%s, cs=%s)",
                 added, disease, cs)
    return grounding


def chat_with_agent(message: str, graph_state: dict, history: list,
                    *, anthropic_api_key: str,
                    use_websearch: bool = True, timeout: int = 300,
                    ncbi_api_key: str | None = None,
                    ncbi_email: str | None = None,
                    hub_dir: Path | None = None) -> dict:
    # 1. RAG grounding — pre-fetch real PubMed citations
    terms = _extract_search_terms(message, graph_state)
    disease_id = graph_state.get("disease", "")
    cs_id = graph_state.get("target_cs", "")
    grounding = fetch_grounding_citations(terms, disease_id,
                                           api_key=ncbi_api_key, email=ncbi_email)
    # 1a. Extend grounding with high-confidence cached agent refs (same
    #     disease + cs context, or sharing a search term).
    grounding = extend_grounding_with_cache(grounding,
                                              disease=disease_id, cs=cs_id,
                                              genes=terms)
    log.info("chat grounding: %d citations (incl cache) from %d terms (%s)",
             len(grounding), len(terms), terms)

    # 1b. Cross-payload SNP + gene lookup — find any SNP/gene/TF mentioned
    #     in the user message that lives in OTHER (disease, cs) combos.
    snp_hits, gene_hits = [], []
    if hub_dir is not None:
        snp_hits = lookup_snps_in_message(message, hub_dir,
                                           current_disease=graph_state.get("disease", ""),
                                           current_cs=graph_state.get("target_cs", ""))
        gene_hits = lookup_genes_in_message(message, hub_dir,
                                              current_disease=graph_state.get("disease", ""),
                                              current_cs=graph_state.get("target_cs", ""))
        if snp_hits or gene_hits:
            log.info("chat cross-lookup: %d SNP hits, %d gene rows (%d genes)",
                     len(snp_hits), len(gene_hits), len({h["symbol"] for h in gene_hits}))

    # 2. Build prompt with verified citation list + cross-payload hits attached
    prompt = build_chat_prompt(message, graph_state, history,
                                grounding_citations=grounding,
                                snp_cross_hits=snp_hits,
                                gene_cross_hits=gene_hits)
    tools = ["WebSearch"] if use_websearch else None
    raw = _call_llm_sdk(prompt, api_key=anthropic_api_key,
                          allowed_tools=tools, timeout=timeout)
    obj = _strip_json(raw)
    if not isinstance(obj, dict):
        return {"message": "(Agent did not return parseable JSON.)",
                "actions": [], "citations": [], "grounding_count": len(grounding),
                "_raw": raw}

    # 3. Verify every PMID the agent returned
    obj["citations"] = verify_citations_in_response(obj, grounding,
                                                     api_key=ncbi_api_key, email=ncbi_email)
    obj["grounding_count"] = len(grounding)
    obj["grounding_terms"] = terms
    # 4. Persist every verified PMID with the current context so the next
    #    chat in this context can reuse them as grounding.
    ctx = {"disease": disease_id, "cs": graph_state.get("target_cs", ""),
           "genes": terms}
    _agent_refs_append(obj["citations"], ctx)
    # Curation is now USER-TRIGGERED via POST /citations/curate (cf. the
    # 📝 Summarize button in the citation card). We do not auto-fire it
    # here because each summary costs ~1.5k BYOK tokens; on Tier 1 plans
    # that can push the user over the per-minute cap during active chat.
    return obj


def chat_with_agent_stream(message: str, graph_state: dict, history: list,
                            *, anthropic_api_key: str,
                            use_websearch: bool = True, timeout: int = 300,
                            ncbi_api_key: str | None = None,
                            ncbi_email: str | None = None,
                            hub_dir: Path | None = None):
    """Generator: yield (event_name, data_dict) tuples.

    Event sequence:
      ('chunk', {'text': '...'})    repeated, as LLM tokens arrive
      ('done',  {full response object})  once at the end
    """
    terms = _extract_search_terms(message, graph_state)
    disease_id = graph_state.get("disease", "")
    cs_id = graph_state.get("target_cs", "")
    grounding = fetch_grounding_citations(terms, disease_id,
                                           api_key=ncbi_api_key, email=ncbi_email)
    grounding = extend_grounding_with_cache(grounding,
                                              disease=disease_id, cs=cs_id,
                                              genes=terms)
    log.info("chat-stream grounding: %d citations (incl cache) from %d terms (%s)",
             len(grounding), len(terms), terms)
    snp_hits, gene_hits = [], []
    if hub_dir is not None:
        snp_hits = lookup_snps_in_message(message, hub_dir,
                                           current_disease=graph_state.get("disease", ""),
                                           current_cs=graph_state.get("target_cs", ""))
        gene_hits = lookup_genes_in_message(message, hub_dir,
                                              current_disease=graph_state.get("disease", ""),
                                              current_cs=graph_state.get("target_cs", ""))
    prompt = build_chat_prompt(message, graph_state, history,
                                grounding_citations=grounding,
                                snp_cross_hits=snp_hits,
                                gene_cross_hits=gene_hits)
    tools = ["WebSearch"] if use_websearch else None
    raw_parts = []
    for kind, payload in _stream_llm_sdk(prompt, api_key=anthropic_api_key,
                                          allowed_tools=tools, timeout=timeout):
        if kind == "text":
            raw_parts.append(payload)
            yield ("chunk", {"text": payload})
        elif kind == "tool":
            # Surface server-tool use ("Searching PubMed for ...",
            # "Fetching pubmed.ncbi.nlm.nih.gov/<pmid>/") to the UI as
            # an SSE 'tool' event so the chat bubble can show an
            # informative pill while the model is mid-tool-call.
            yield ("tool", payload)
    raw = "".join(raw_parts)
    obj = _strip_json(raw)
    if not isinstance(obj, dict):
        yield ("done", {"message": "(Agent did not return parseable JSON.)",
                         "actions": [], "citations": [],
                         "grounding_count": len(grounding), "_raw": raw})
        return
    obj["citations"] = verify_citations_in_response(obj, grounding,
                                                     api_key=ncbi_api_key, email=ncbi_email)
    obj["grounding_count"] = len(grounding)
    obj["grounding_terms"] = terms
    ctx = {"disease": disease_id, "cs": cs_id, "genes": terms}
    _agent_refs_append(obj["citations"], ctx)
    yield ("done", obj)
    # Note: curation (PubMed abstract → AI summary → /var/lib/heartgrn/
    # curated_knowledge/<pmid>.json) is USER-TRIGGERED via the 📝
    # Summarize button → POST /citations/curate. We intentionally do
    # NOT auto-summarize here to keep the user's BYOK token spend
    # explicit.


def verify_pmids(citations: list, api_key: str | None = None,
                  email: str | None = None) -> list:
    """Hit esummary for each claimed PMID; mark items as unverified if
    the PMID does not exist or returns no title."""
    out = []
    pmids = [c.get("pmid") for c in citations if c.get("pmid")]
    if not pmids:
        return [{**c, "verified": False} for c in citations]
    try:
        raw = _http_get(NCBI_ESUMMARY, {
            "db": "pubmed", "id": ",".join(str(p) for p in pmids), "retmode": "json",
        }, api_key=api_key, email=email)
        summ = json.loads(raw).get("result", {})
    except Exception:
        summ = {}
    for c in citations:
        pmid = c.get("pmid")
        if pmid and str(pmid) in summ:
            real = summ[str(pmid)]
            out.append({**c,
                        "verified": True,
                        "title": real.get("title", c.get("title")).strip("."),
                        "journal": real.get("source", c.get("journal")),
                        "year": int(real.get("pubdate", "0").split()[0]) if real.get("pubdate") else c.get("year"),
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        else:
            out.append({**c, "verified": False})
    return out


# ---- Cross-payload SNP index ------------------------------------------------
import threading as _threading_for_index
_snp_index_cache: dict | None = None
_snp_index_lock = _threading_for_index.Lock()


def build_snp_index(hub_dir: Path) -> dict:
    """One-shot scan of every `hub/payloads/*.json` building lookup maps:
       by_variant_id  "8_30423317_G_A"  → [{disease, cs, peak, rsid, pip, trait}]
       by_rsid        "rs2820315"       → [same as above]
       by_gene        "MYO9B"           → [{disease, cs, role, gwas_z, target_de_z, cat}]
    """
    idx: dict = {"by_variant_id": {}, "by_rsid": {}, "by_gene": {},
                  "n_payloads": 0, "n_snps": 0, "n_gene_rows": 0,
                  "ts": __import__("datetime").datetime.utcnow().isoformat() + "Z"}
    pd_ = hub_dir / "payloads"
    if not pd_.exists():
        return idx
    for pf in sorted(pd_.glob("*.json")):
        try:
            p = json.loads(pf.read_text())
        except Exception as e:
            log.warning("snp_index: failed to parse %s: %s", pf.name, e); continue
        disease = p.get("disease"); cs = p.get("target_cs")
        idx["n_payloads"] += 1
        # SNPs
        for s in p.get("snps", []) or []:
            chr_ = s.get("chr"); pos = s.get("pos")
            ref = s.get("ref"); alt = s.get("alt")
            # Skip rows with any missing locus field — without this guard the
            # f-string below produced a literal "None_None_None_None" key that
            # accumulated every malformed row and polluted cross-payload search.
            if not (chr_ and pos and ref and alt):
                continue
            vid = f"{chr_}_{pos}_{ref}_{alt}"
            rsid = (s.get("rsid") or "").strip()
            entry = {"disease": disease, "cs": cs,
                      "peak": s.get("peak"), "rsid": rsid,
                      "pip": s.get("pip"), "trait": s.get("trait")}
            idx["by_variant_id"].setdefault(vid, []).append(entry)
            if rsid.startswith("rs"):
                idx["by_rsid"].setdefault(rsid, []).append(entry)
            idx["n_snps"] += 1
        # Seed genes — these appear because they passed the top-N inner_min filter
        for seed in p.get("seed_genes", []) or []:
            sid = seed.get("id")
            if not sid: continue
            idx["by_gene"].setdefault(sid, []).append({
                "disease": disease, "cs": cs, "role": "seed",
                "gwas_z": seed.get("gwas_z"),
                "target_de_z": seed.get("target_de_z"),
                "inner_min": seed.get("inner_min"),
            })
            idx["n_gene_rows"] += 1
        # TFs
        for tf in p.get("tfs", []) or []:
            tid = tf.get("id")
            if not tid: continue
            idx["by_gene"].setdefault(tid, []).append({
                "disease": disease, "cs": cs, "role": "tf",
                "gwas_z": tf.get("gwas_z"),
                "target_de_z": tf.get("target_de_z"),
            })
            idx["n_gene_rows"] += 1
    log.info("snp_index built: %d payloads, %d snp rows, %d gene rows, "
             "%d uniq variants, %d uniq gene symbols",
              idx["n_payloads"], idx["n_snps"], idx["n_gene_rows"],
              len(idx["by_variant_id"]), len(idx["by_gene"]))
    return idx


def get_snp_index(hub_dir: Path) -> dict:
    global _snp_index_cache
    with _snp_index_lock:
        if _snp_index_cache is None:
            _snp_index_cache = build_snp_index(hub_dir)
        return _snp_index_cache


def lookup_genes_in_message(message: str, hub_dir: Path,
                             current_disease: str = "",
                             current_cs: str = "",
                             max_per_gene: int = 10) -> list[dict]:
    """Extract uppercase gene symbols from message → cross-payload hits.
    For each gene, return up to `max_per_gene` rows sorted by target_de_z
    descending so the agent can route to the cs with strongest expression.
    """
    idx = get_snp_index(hub_dir)
    by_gene = idx.get("by_gene") or {}
    # Common false-positive tokens (acronyms, citation tags, etc.)
    SKIP = {"GWAS","SNP","DNA","RNA","AF","CAD","AVS","PMID","DOI","URL",
            "JSON","API","MHz","TBD","TFs","JACC","PMC","PNAS","NIH","OK",
            "ALL","MeSH","TF","FDR","HGNC","GTEx","VEP","UTR","CDS","EUR",
            "EAS","AFR","SD","CI","HR","BP","UK","LD","PIP"}
    # Case-insensitive match so "tbx5" / "Tbx5" / "TBX5" all hit the index.
    sym_re = re.compile(r"\b([A-Za-z][A-Za-z0-9-]{1,})\b")
    seen: set[str] = set()
    hits: list[dict] = []
    for m in sym_re.finditer(message):
        sym = m.group(1).upper()
        if sym in seen or sym in SKIP or len(sym) < 3:
            continue
        seen.add(sym)
        rows = by_gene.get(sym) or []
        # Sort: max target_de_z first (so the cs with strongest expression bubbles up)
        rows_sorted = sorted(
            rows,
            key=lambda r: (r.get("target_de_z") if r.get("target_de_z") is not None else -1e18),
            reverse=True,
        )[:max_per_gene]
        for h in rows_sorted:
            hits.append({"symbol": sym, **h,
                         "in_current_view": (h["disease"] == current_disease
                                              and h["cs"] == current_cs)})
    return hits


def lookup_snps_in_message(message: str, hub_dir: Path,
                            current_disease: str = "",
                            current_cs: str = "") -> list[dict]:
    """Extract variant_id / rsid tokens from message → cross-payload hits."""
    idx = get_snp_index(hub_dir)
    seen: set[str] = set()
    hits: list[dict] = []
    # variant_id pattern: 4_110637255_G_C or chr4_110637255_G_C
    vid_re = re.compile(r"\b(?:chr)?([0-9XYM]+)_(\d+)_([ACGT]+)_([ACGT]+)\b")
    for m in vid_re.finditer(message):
        vid = f"{m.group(1)}_{m.group(2)}_{m.group(3)}_{m.group(4)}"
        if vid in seen: continue
        seen.add(vid)
        for h in idx["by_variant_id"].get(vid, []):
            hits.append({"matched_token": vid, **h})
    # rsid
    rs_re = re.compile(r"\brs\d+\b")
    for m in rs_re.finditer(message):
        rs = m.group(0)
        if rs in seen: continue
        seen.add(rs)
        for h in idx["by_rsid"].get(rs, []):
            hits.append({"matched_token": rs, **h})
    # Categorise: in-current vs elsewhere
    for h in hits:
        h["in_current_view"] = (h["disease"] == current_disease
                                  and h["cs"] == current_cs)
    return hits


# ---- AG1 subprocess runner --------------------------------------------------

import threading
_ag1_locks: dict[str, threading.Lock] = {}
_ag1_locks_guard = threading.Lock()


def _ag1_dir_for(hub_dir: Path, variant_id: str) -> Path:
    return hub_dir / "ag1_cache" / variant_id


def _read_json_safe(p: Path) -> dict | None:
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def _ag1_status_path(hub_dir: Path, variant_id: str) -> Path:
    return _ag1_dir_for(hub_dir, variant_id) / "status.json"


def _ag1_write_status(hub_dir: Path, variant_id: str, **kw) -> dict:
    d = _ag1_dir_for(hub_dir, variant_id); d.mkdir(parents=True, exist_ok=True)
    existing = _read_json_safe(d / "status.json") or {}
    existing.update(kw)
    existing["ts"] = __import__("datetime").datetime.utcnow().isoformat() + "Z"
    (d / "status.json").write_text(json.dumps(existing, indent=2))
    return existing


def _ag1_summarise_csv(csv_path: Path, top_n: int = 200) -> dict:
    """Return rows sorted by abs_score (descending) from atac/rna CSV.
    Default keeps top 200 — viewer wraps in scrollable container."""
    if not csv_path.exists():
        return {"top": []}
    import csv as _csv
    rows = []
    with csv_path.open() as f:
        for row in _csv.DictReader(f):
            try:
                rows.append({
                    "cell_state": row.get("cell_state"),
                    "gene_name":  row.get("gene_name"),
                    "raw_score":  float(row.get("raw_score") or 0),
                    "abs_score":  float(row.get("abs_score") or 0),
                    "track_name": row.get("track_name"),
                })
            except Exception:
                pass
    rows.sort(key=lambda r: r["abs_score"], reverse=True)
    return {"top": rows[:top_n], "n_rows": len(rows)}


def _ag1_score_subprocess(hub_dir: Path, variant_id: str) -> None:
    """Run ag1.py score in background; write status updates + summary."""
    out_dir = _ag1_dir_for(hub_dir, variant_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("[ag1] scoring %s → %s", variant_id, out_dir)
    cmd = [AG1_PYTHON, str(AG1_SCRIPT), "score",
           "--variant", variant_id, "--out", str(out_dir)]
    log_path = out_dir / "score.log"
    try:
        with log_path.open("w") as f:
            res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                                  timeout=1500)
        if res.returncode != 0:
            _ag1_write_status(hub_dir, variant_id, status="error",
                              error=f"exit {res.returncode}",
                              log=str(log_path.relative_to(hub_dir)))
            return
        atac_csv = out_dir / f"{variant_id}_atac.csv"
        rna_csv  = out_dir / f"{variant_id}_rna.csv"
        summary = {
            "atac": _ag1_summarise_csv(atac_csv),
            "rna":  _ag1_summarise_csv(rna_csv),
            "atac_csv": str(atac_csv.relative_to(hub_dir)),
            "rna_csv":  str(rna_csv.relative_to(hub_dir)),
            "log":      str(log_path.relative_to(hub_dir)),
        }
        _ag1_write_status(hub_dir, variant_id, status="done", **summary)
        log.info("[ag1] done %s (%d atac, %d rna)", variant_id,
                 summary["atac"].get("n_rows", 0), summary["rna"].get("n_rows", 0))
    except subprocess.TimeoutExpired:
        _ag1_write_status(hub_dir, variant_id, status="error", error="timeout")
    except Exception as e:
        _ag1_write_status(hub_dir, variant_id, status="error", error=str(e))


def _ag1_render_subprocess(hub_dir: Path, variant_id: str,
                            cell_states: list[str], window: int) -> dict:
    """Run ag1.py render (synchronous — fast once score cache exists)."""
    out_dir = _ag1_dir_for(hub_dir, variant_id); out_dir.mkdir(parents=True, exist_ok=True)
    cs_args = []
    for cs in cell_states:
        cs_args += ["--cell-state", cs]
    cmd = [AG1_PYTHON, str(AG1_SCRIPT), "render",
           "--variant", variant_id, "--window", str(window),
           "--out", str(out_dir)] + cs_args
    log_path = out_dir / "render.log"
    try:
        with log_path.open("w") as f:
            res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                                  timeout=1200)
        if res.returncode != 0:
            return {"error": f"exit {res.returncode}",
                    "log": str(log_path.relative_to(hub_dir))}
        # ag1 render writes track_<variant_id>.{pdf,png} into out_dir
        pngs = sorted(out_dir.glob("*.png"),
                      key=lambda p: p.stat().st_mtime, reverse=True)
        if not pngs:
            return {"error": "no PNG produced",
                    "log": str(log_path.relative_to(hub_dir))}
        return {"png": str(pngs[0].relative_to(hub_dir)),
                "cell_states": cell_states, "window": window}
    except subprocess.TimeoutExpired:
        return {"error": "render timeout"}
    except Exception as e:
        return {"error": str(e)}


# ---- credible citation set (user-curated trusted PMIDs) ---------------------
CREDIBLE_CITATIONS_FILE = Path(os.environ.get(
    "CREDIBLE_CITATIONS_FILE",
    "/var/lib/heartgrn/credible_citations.jsonl",
))
CREDIBLE_IP_SALT = os.environ.get("CREDIBLE_IP_SALT", "credible-default-salt")


def _credible_path() -> Path:
    """Return write-target; fall back to /tmp if /var/lib path is read-only."""
    p = CREDIBLE_CITATIONS_FILE
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    except (PermissionError, OSError):
        return Path("/tmp/credible_citations.jsonl")


def _hash_ip(ip: str) -> str:
    """One-way hash so we can dedup without storing raw client IPs."""
    return hashlib.sha256(f"{ip}|{CREDIBLE_IP_SALT}".encode()).hexdigest()[:16]


# ---- agent reference cache (tagged PMID store) -----------------------------
AGENT_REFS_FILE = Path(os.environ.get(
    "AGENT_REFS_FILE",
    "/var/lib/heartgrn/agent_references.jsonl",
))


def _agent_refs_path() -> Path:
    p = AGENT_REFS_FILE
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    except (PermissionError, OSError):
        return Path("/tmp/agent_references.jsonl")


def _agent_refs_append(citations: list, context: dict) -> None:
    """Append every cited PMID with the current chat context as tags.
    Tags are how we filter cached refs for re-use on the next chat in the
    same context (disease + cs + overlapping genes)."""
    if not citations or not isinstance(citations, list):
        return
    p = _agent_refs_path()
    try:
        with p.open("a") as fh:
            for c in citations:
                pmid = str(c.get("pmid") or "").strip()
                if not pmid or not pmid.isdigit():
                    continue
                if c.get("hallucinated") or not c.get("verified"):
                    continue   # only persist real, verified PMIDs
                entry = {
                    "pmid": pmid,
                    "title": (c.get("title") or "")[:250],
                    "journal": (c.get("journal") or "")[:80],
                    "year": c.get("year"),
                    "key_finding": (c.get("key_finding") or "")[:300],
                    "tags": {
                        "disease": context.get("disease", ""),
                        "cs": context.get("cs", ""),
                        "genes": (context.get("genes") or [])[:20],
                    },
                    "ts": __import__("datetime").datetime.utcnow().isoformat() + "Z",
                }
                fh.write(json.dumps(entry, separators=(",", ":")) + "\n")
    except OSError as e:
        log.warning("agent_refs append failed: %s", e)


def _agent_refs_by_context(disease: str, cs: str,
                            genes: list[str] | None = None) -> dict:
    """Aggregate cached PMIDs matching the given context.
    Returns {pmid: {pmid, title, journal, year, key_finding, recurrence,
                     tags_seen, last_ts}}.
    Match rule: same disease AND (same cs OR overlapping gene)."""
    p = _agent_refs_path()
    if not p.exists():
        return {}
    gene_set = set((g or "").upper() for g in (genes or []) if g)
    out: dict = {}
    try:
        with p.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    e = json.loads(line)
                except json.JSONDecodeError:
                    continue
                pmid = str(e.get("pmid") or "")
                if not pmid:
                    continue
                tags = e.get("tags") or {}
                t_disease = tags.get("disease", "")
                t_cs = tags.get("cs", "")
                t_genes = set((g or "").upper() for g in (tags.get("genes") or []))
                if t_disease != disease:
                    continue
                if cs and t_cs != cs and not (gene_set and t_genes & gene_set):
                    continue
                r = out.setdefault(pmid, {
                    "pmid": pmid,
                    "title": e.get("title", ""),
                    "journal": e.get("journal", ""),
                    "year": e.get("year"),
                    "key_finding": e.get("key_finding", ""),
                    "recurrence": 0,
                    "last_ts": "",
                    "tags_seen": [],
                })
                r["recurrence"] += 1
                ts = e.get("ts", "")
                if ts > r["last_ts"]:
                    r["last_ts"] = ts
                r["tags_seen"].append(tags)
    except OSError:
        return {}
    return out


# ---- curated knowledge (PMID → vote-curated AI summary) --------------------
CURATED_KNOWLEDGE_DIR = Path(os.environ.get(
    "CURATED_KNOWLEDGE_DIR",
    "/var/lib/heartgrn/curated_knowledge",
))


def _curated_knowledge_dir() -> Path:
    try:
        CURATED_KNOWLEDGE_DIR.mkdir(parents=True, exist_ok=True)
        return CURATED_KNOWLEDGE_DIR
    except (PermissionError, OSError):
        return Path("/tmp/curated_knowledge")


def _curated_path(pmid: str) -> Path:
    return _curated_knowledge_dir() / f"{pmid}.json"


def _load_curated(pmid: str) -> dict | None:
    return _load_json_safe(_curated_path(pmid))


def _fetch_pubmed_abstract(pmid: str, *, api_key: str | None = None,
                            email: str | None = None) -> str:
    """One-shot efetch for a PubMed PMID. Returns the plain-text abstract
    (with title), or empty string on failure."""
    try:
        raw = _http_get(NCBI_EFETCH, {
            "db": "pubmed", "id": pmid,
            "rettype": "abstract", "retmode": "text",
        }, api_key=api_key, email=email)
        return (raw or "").strip()
    except Exception as e:
        log.warning("efetch abstract failed for PMID %s: %s", pmid, e)
        return ""


def _summarize_and_persist(pmid: str, citation_meta: dict, context: dict,
                            *, api_key: str,
                            ncbi_api_key: str | None = None,
                            ncbi_email: str | None = None) -> None:
    """Background worker. Skip if a summary already exists; otherwise fetch
    the PubMed abstract, ask the LLM for a 3-4 sentence digest + key findings,
    and persist as /var/lib/heartgrn/curated_knowledge/<pmid>.json. Failures
    are logged but never raised — this runs in a daemon thread and must not
    crash the parent."""
    if not pmid or _load_curated(pmid):
        return
    abstract = _fetch_pubmed_abstract(pmid, api_key=ncbi_api_key, email=ncbi_email)
    if not abstract or len(abstract) < 60:
        log.warning("curated skipped (no abstract) PMID=%s", pmid)
        return
    prompt = (
        "Summarize the following PubMed abstract for a cardiac gene-regulatory-"
        "network research context in 3-4 sentences. Focus on: the key biological "
        "finding, the gene(s) / variant(s) involved, the cell type or tissue, and "
        "the clinical relevance to cardiovascular disease. Then list 3-5 short "
        "bullet 'key findings' (≤ 14 words each).\n\n"
        "OUTPUT FORMAT — strict: emit a RAW JSON object, nothing else. "
        "Do NOT wrap it in ```json … ``` or any other markdown fences. "
        "Do NOT add a preamble (`Here is the summary:`) or trailing prose. "
        "First character of your response must be `{` and last must be `}`.\n"
        '{"summary": "...", "key_findings": ["...", "..."]}\n\n'
        f"PMID: {pmid}\n"
        f"Title: {citation_meta.get('title','')}\n"
        f"Journal: {citation_meta.get('journal','')}  Year: {citation_meta.get('year','')}\n\n"
        f"Abstract:\n{abstract[:7000]}"
    )
    try:
        # 1500 tokens gives ample headroom; Gemini in particular emits
        # verbose JSON with markdown fences that ate budget at 700.
        raw = _call_llm_sdk(prompt, api_key=api_key, max_tokens=1500, timeout=90)
    except Exception as e:
        log.warning("curate summarize LLM call failed PMID=%s: %s", pmid, e)
        return
    parsed = _strip_json(raw)
    if not isinstance(parsed, dict):
        log.warning("curate summary unparseable for PMID=%s; raw_len=%d head=%r tail=%r",
                     pmid, len(raw or ""), (raw or "")[:400], (raw or "")[-200:])
        return
    summary = (parsed.get("summary") or "").strip()
    key_findings = parsed.get("key_findings") or []
    if not isinstance(key_findings, list):
        key_findings = []
    if not summary:
        log.warning("curate summary empty for PMID=%s", pmid)
        return
    entry = {
        "pmid": pmid,
        "title": citation_meta.get("title"),
        "journal": citation_meta.get("journal"),
        "year": citation_meta.get("year"),
        "abstract": abstract[:10000],
        "summary": summary,
        "key_findings": [str(x)[:200] for x in key_findings][:8],
        "first_high_ts": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "context": context,
    }
    try:
        _atomic_write_text(_curated_path(pmid), json.dumps(entry, indent=2))
        log.info("curated knowledge saved PMID=%s len(summary)=%d kf=%d",
                  pmid, len(summary), len(entry["key_findings"]))
    except OSError as e:
        log.warning("curate persist failed PMID=%s: %s", pmid, e)


def _credible_confidence(up: int, down: int, recurrence: int = 0) -> str:
    """Confidence band. HUMAN VOTES are the gatekeeper for 'high' — agent
    recurrence alone cannot promote a PMID to high.

       flagged → users net-down-voted (down > up)
       high    → human net up-votes >= 2  (recurrence is NOT sufficient
                                            because a deterministic LLM can
                                            cycle a PMID across context-
                                            matched chats without any human
                                            signal; Codex flagged this as a
                                            gameable path that silently
                                            leaks into next chat's grounding)
       medium  → net == 1, or recurrence >= 2
       low     → otherwise (no votes, single agent cite)
    """
    up = int(up or 0); down = int(down or 0); recurrence = int(recurrence or 0)
    if down > up:
        return "flagged"
    net = up - down
    if net >= 2:
        return "high"
    if net == 1 or recurrence >= 2:
        return "medium"
    return "low"


def _credible_aggregate() -> dict:
    """Aggregate the append-only jsonl into per-PMID up/down votes.
    Latest vote per (pmid, ip_hash) wins, so users can change their mind."""
    p = _credible_path()
    if not p.exists():
        return {}
    # Read all entries, keep the latest per (pmid, ip_hash) pair.
    latest: dict = {}  # (pmid, ip_hash) -> entry
    try:
        with p.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:
                    continue
                pmid = entry.get("pmid")
                ip_hash = entry.get("ip_hash", "")
                if not pmid:
                    continue
                key = (pmid, ip_hash)
                prev = latest.get(key)
                if prev is None or entry.get("ts", "") >= prev.get("ts", ""):
                    latest[key] = entry
    except OSError:
        return {}
    # Aggregate per pmid.
    out: dict = {}
    for (pmid, _ip), entry in latest.items():
        vote = entry.get("vote", "up")
        e = out.setdefault(pmid, {"up": 0, "down": 0, "last_ts": "",
                                   "contexts": [], "title": ""})
        if vote == "up":
            e["up"] += 1
        elif vote == "down":
            e["down"] += 1
        e["last_ts"] = max(e["last_ts"], entry.get("ts", ""))
        if entry.get("context") and entry["context"] not in e["contexts"]:
            e["contexts"].append(entry["context"])
        if entry.get("title") and not e["title"]:
            e["title"] = entry["title"]
    return out


# ---- HTTP handler -----------------------------------------------------------

class LitHandler(BaseHTTPRequestHandler):
    server_version = "snp2cell-hub-backend/2.0"

    @property
    def hub_dir(self) -> Path:
        return self.server.hub_dir  # type: ignore[attr-defined]

    def _set_cors(self):
        origin = self.headers.get("Origin", "")
        if origin in ALLOWED_ORIGINS:
            self.send_header("Access-Control-Allow-Origin", origin)
            self.send_header("Vary", "Origin")
        elif not ALLOWED_ORIGINS or os.environ.get("ALLOW_ANY_ORIGIN") == "1":
            self.send_header("Access-Control-Allow-Origin", "*")
        # else: no ACAO header → browser blocks (intentional)
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, X-API-Key")
        self.send_header("Access-Control-Max-Age", "86400")

    def _json(self, status: int, obj: object):
        body = json.dumps(obj, default=str).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self._set_cors()
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_response(204); self._set_cors(); self.end_headers()

    def log_message(self, fmt, *args):
        log.info("%s - %s", self.address_string(), fmt % args)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/healthz":
            self._json(200, {"status": "ok", "hub_dir": str(self.hub_dir),
                             "backend": "pubmed-direct/2.0"})
            return
        if path == "/" or path == "":
            # Friendly landing page so a stray browser visitor sees something
            # useful instead of bare 404. Backend has no UI of its own.
            self._json(200, {
                "service": "HeartGRN Viewer — backend API",
                "status": "ok",
                "info": "This is the API endpoint. The interactive viewer lives elsewhere.",
                "endpoints": {
                    "GET /healthz": "liveness check",
                    "GET /literature/{type}/{key}": "cached PubMed/Europe PMC results (type ∈ snp|tf|gene)",
                    "POST /literature/{type}/{key}": "run PubMed lookup; mode=agent/hybrid + synthesize=true require X-API-Key (sk-ant-...)",
                    "POST /chat": "interactive agent (requires X-API-Key BYOK Anthropic key)",
                    "GET /ag1/status/{variant}": "AlphaGenome score job status",
                    "POST /ag1/score/{variant}": "kick off AG1 score subprocess",
                    "POST /ag1/render/{variant}": "render AG1 figure",
                },
                "frontend": "https://heartgrn-viewer-online.pages.dev (Cloudflare Pages, password-gated staging)",
                "repo": "https://github.com/blackjtaka/heartgrn-viewer-online",
            })
            return
        m = re.fullmatch(r"/ag1/status/(.+)", path)
        if m:
            variant_id = m.group(1)
            status = _read_json_safe(_ag1_status_path(self.hub_dir, variant_id))
            if not status:
                self._json(404, {"status": "idle"}); return
            # Lazy re-summarise: older caches stored only 8 rows; new default is
            # 200. If CSVs exist locally and the existing summary is small,
            # re-build from CSV (no AG1 re-inference) and persist.
            if status.get("status") == "done":
                need_resum = (len((status.get("atac") or {}).get("top") or []) < 50
                              or len((status.get("rna") or {}).get("top") or []) < 50)
                if need_resum:
                    out_dir = _ag1_dir_for(self.hub_dir, variant_id)
                    atac_csv = out_dir / f"{variant_id}_atac.csv"
                    rna_csv  = out_dir / f"{variant_id}_rna.csv"
                    if atac_csv.exists() and rna_csv.exists():
                        status["atac"] = _ag1_summarise_csv(atac_csv)
                        status["rna"]  = _ag1_summarise_csv(rna_csv)
                        _ag1_write_status(self.hub_dir, variant_id, **status)
                        log.info("ag1 re-summarised %s (top→%d/%d)", variant_id,
                                 status["atac"].get("n_rows", 0),
                                 status["rna"].get("n_rows", 0))
                # Attach the list of pre-rendered REF/ALT/Δ tracks (PNGs) so
                # the frontend can offer them as plot buttons in beta mode.
                # Filename convention: {variant_id}_{cell_state}_w{window}.png
                out_dir = _ag1_dir_for(self.hub_dir, variant_id)
                fig_re = re.compile(rf"^{re.escape(variant_id)}_(.+)_w(\d+)\.png$")
                figs = []
                if out_dir.exists():
                    for png in sorted(out_dir.glob(f"{variant_id}_*.png")):
                        mm = fig_re.match(png.name)
                        if not mm:
                            continue
                        figs.append({
                            "cell_state": mm.group(1),
                            "window": int(mm.group(2)),
                            "png_url": f"ag1_cache/{variant_id}/{png.name}",
                        })
                status["figures"] = figs
            self._json(200, status); return
        # GET /citations/credible — list user-curated trusted PMIDs (aggregated)
        if path == "/citations/credible":
            agg = _credible_aggregate()
            self._json(200, {"credible": agg, "total_pmids": len(agg)}); return
        m = re.fullmatch(r"/citations/credible/(\d{1,10})", path)
        if m:
            pmid = m.group(1)
            agg = _credible_aggregate()
            self._json(200, {"pmid": pmid, **(agg.get(pmid, {"count": 0}))}); return

        # /snp/lookup/<rsid_or_variant_id> — debug helper
        m = re.fullmatch(r"/snp/lookup/(.+)", path)
        if m:
            tok = m.group(1)
            idx = get_snp_index(self.hub_dir)
            hits = (idx["by_rsid"].get(tok) or
                    idx["by_variant_id"].get(tok) or [])
            self._json(200, {"token": tok, "hits": hits,
                              "n_payloads_indexed": idx.get("n_payloads"),
                              "n_uniq_variants": len(idx.get("by_variant_id") or {})});
            return
        # /gene/lookup/<symbol> — debug helper for cross-payload gene index
        m = re.fullmatch(r"/gene/lookup/(.+)", path)
        if m:
            sym = m.group(1)
            idx = get_snp_index(self.hub_dir)
            rows = (idx.get("by_gene") or {}).get(sym) or []
            rows_sorted = sorted(rows,
                key=lambda r: (r.get("target_de_z") if r.get("target_de_z") is not None else -1e18),
                reverse=True)
            self._json(200, {"symbol": sym, "n_hits": len(rows_sorted),
                              "top": rows_sorted[:20],
                              "n_uniq_gene_symbols": len(idx.get("by_gene") or {})});
            return
        # /gene/search/<q> — case-insensitive substring search across all
        # gene symbols in the cross-payload index. Returns flat hits, each
        # tagged with its symbol so the frontend can render a single list.
        m = re.fullmatch(r"/gene/search/(.+)", path)
        if m:
            q = m.group(1).upper()
            idx = get_snp_index(self.hub_dir)
            by_gene = idx.get("by_gene") or {}
            # Match shorter symbols first (more likely the exact intent),
            # cap symbol count to avoid giant payloads.
            syms = [s for s in by_gene if q in s.upper()]
            syms.sort(key=lambda s: (0 if s.upper() == q else 1, len(s), s))
            flat = []
            for sym in syms[:30]:
                for r in by_gene[sym][:8]:
                    flat.append({**r, "symbol": sym})
            # Rank: exact symbol match first, then by target_de_z desc.
            flat.sort(key=lambda r: (0 if r["symbol"].upper() == q else 1,
                                       -(r.get("target_de_z") or -1e18)))
            self._json(200, {"q": q,
                              "n_symbols": len(syms),
                              "matches": flat[:40]})
            return
        if path == "/literature/list":
            entries = []
            for t in ALLOWED_TYPES:
                tdir = self.hub_dir / "literature" / t
                if not tdir.exists():
                    continue
                for f in tdir.glob("*.json"):
                    entries.append({"type": t, "key": f.stem,
                                    "mtime": f.stat().st_mtime})
            self._json(200, {"entries": entries})
            return
        m = re.fullmatch(r"/literature/(\w+)/(.+)", path)
        if m:
            type_, key = m.group(1), m.group(2)
            if type_ not in ALLOWED_TYPES:
                self._json(400, {"error": f"unknown type {type_!r}"}); return
            p = cache_path(self.hub_dir, type_, key)
            if not p.exists():
                self._json(404, {"error": "not cached", "key": key}); return
            self._json(200, json.loads(p.read_text())); return
        self._json(404, {"error": "not found"})

    def do_POST(self):
        path = urlparse(self.path).path

        # /ag1/score|render/<variant_id> — disabled in beta deployment.
        # On-demand AG1 inference requires the AG1 Python script (only present
        # on the developer's local workstation), so neither score nor render
        # can run server-side. Cached AG1 results stay queryable via GET
        # /ag1/status/<variant_id>.
        m = re.fullmatch(r"/ag1/(score|render)/(.+)", path)
        if m:
            action, variant_id = m.group(1), m.group(2)
            self._json(503, {
                "error": "ag1_disabled_in_beta",
                "action": action,
                "variant_id": variant_id,
                "message": ("AG1 on-demand inference is disabled in this beta. "
                            "Only the pre-cached representative variants under "
                            "hub/ag1_cache/ are available, queryable via "
                            "GET /ag1/status/<variant_id>."),
            }); return

        # /chat: interactive Ask Agent (BYOK — user supplies X-API-Key header)
        if path == "/chat":
            # CSRF defense: Origin header must match ALLOWED_ORIGINS.
            # Browsers always send Origin on cross-origin POST; same-origin
            # also fine. Curl / scripts without Origin are blocked unless
            # ALLOW_ANY_ORIGIN=1 is set (e.g. for testing).
            req_origin = self.headers.get("Origin", "").strip()
            if ALLOWED_ORIGINS and os.environ.get("ALLOW_ANY_ORIGIN") != "1":
                if not req_origin or req_origin not in ALLOWED_ORIGINS:
                    log.warning("chat blocked: origin=%r not in ALLOWED_ORIGINS=%r",
                                req_origin, ALLOWED_ORIGINS)
                    self._json(403, {"error": "forbidden_origin",
                                     "message": "Origin header missing or not allowed"})
                    return
            # Per-IP rate limit. Cloudflare forwards the real client IP via CF-Connecting-IP;
            # respect X-Forwarded-For as a fallback, otherwise use the socket peer.
            client_ip = (self.headers.get("CF-Connecting-IP")
                         or (self.headers.get("X-Forwarded-For", "").split(",")[0].strip() or None)
                         or self.client_address[0])
            allowed, retry_after = rate_limit_check(client_ip, "/chat")
            if not allowed:
                self.send_response(429)
                self.send_header("Retry-After", str(retry_after))
                self.send_header("Content-Type", "application/json")
                self._set_cors()
                body = json.dumps({"error": "rate_limited",
                                   "message": f"Too many requests; retry after {retry_after}s",
                                   "retry_after": retry_after}).encode()
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return

            user_key = self.headers.get("X-API-Key", "").strip()
            if not user_key:
                self._json(401, {"error": "byok_required",
                                 "message": "Provide an API key in X-API-Key header (Anthropic sk-ant-... or Google AIzaSy...)"}); return
            if not (user_key.startswith("sk-ant-") or user_key.startswith("AIza")):
                self._json(401, {"error": "invalid_key_format",
                                 "message": "API key must start with 'sk-ant-' (Anthropic) or 'AIza' (Google Gemini)"}); return
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length) if length else b""
            try:
                body = json.loads(raw) if raw else {}
            except json.JSONDecodeError as e:
                self._json(400, {"error": f"bad json body: {e}"}); return
            msg = (body.get("message") or "").strip()
            if not msg:
                self._json(400, {"error": "empty message"}); return
            graph_state = body.get("graph_state", {})
            history = body.get("history", []) or []
            use_websearch = bool(body.get("use_websearch", True))
            try:
                resp = chat_with_agent(msg, graph_state, history,
                                       anthropic_api_key=user_key,
                                       use_websearch=use_websearch,
                                       ncbi_api_key=self.server.ncbi_api_key,  # type: ignore[attr-defined]
                                       ncbi_email=self.server.ncbi_email,       # type: ignore[attr-defined]
                                       hub_dir=self.hub_dir)
            except AgentError as e:
                # Surface auth/quota errors to user clearly so they can fix their key
                err_lower = str(e).lower()
                if "auth" in err_lower or "invalid" in err_lower:
                    self._json(401, {"error": "anthropic_auth_failed", "detail": str(e)}); return
                if "rate" in err_lower or "limit" in err_lower:
                    self._json(429, {"error": "anthropic_rate_limited", "detail": str(e)}); return
                self._json(502, {"error": "agent failed", "detail": str(e)}); return
            # NEVER log the user_key; only log non-secret metadata
            log.info("chat → %d actions (%s, %s)",
                     len(resp.get("actions", [])),
                     graph_state.get("disease", "?"),
                     graph_state.get("target_cs", "?"))
            self._json(200, resp)
            return

        # POST /chat/stream — same as /chat but returns Server-Sent Events
        # so the frontend can render text as it arrives. Final 'done' event
        # carries the post-processed response (verified citations, actions).
        if path == "/chat/stream":
            req_origin = self.headers.get("Origin", "").strip()
            if ALLOWED_ORIGINS and os.environ.get("ALLOW_ANY_ORIGIN") != "1":
                if not req_origin or req_origin not in ALLOWED_ORIGINS:
                    log.warning("chat-stream blocked: origin=%r not in ALLOWED_ORIGINS=%r",
                                req_origin, ALLOWED_ORIGINS)
                    self._json(403, {"error": "forbidden_origin",
                                      "got_origin": req_origin,
                                      "allowed": list(ALLOWED_ORIGINS)})
                    return
            client_ip = (self.headers.get("CF-Connecting-IP")
                         or (self.headers.get("X-Forwarded-For", "").split(",")[0].strip() or None)
                         or self.client_address[0])
            allowed, retry_after = rate_limit_check(client_ip, "/chat")
            if not allowed:
                self._json(429, {"error": "rate_limited",
                                  "message": f"Too many requests; retry after {retry_after}s",
                                  "retry_after": retry_after}); return
            user_key = self.headers.get("X-API-Key", "").strip()
            if not user_key or not (user_key.startswith("sk-ant-") or user_key.startswith("AIza")):
                self._json(401, {"error": "byok_required",
                                  "message": "API key must start with 'sk-ant-' or 'AIza'"}); return
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length) if length else b""
            try:
                body = json.loads(raw) if raw else {}
            except json.JSONDecodeError as e:
                self._json(400, {"error": f"bad json body: {e}"}); return
            msg = (body.get("message") or "").strip()
            if not msg:
                self._json(400, {"error": "empty message"}); return
            graph_state = body.get("graph_state", {})
            history = body.get("history", []) or []
            use_websearch = bool(body.get("use_websearch", True))
            # SSE headers. Use Connection: close so the client reader sees
            # a definitive end-of-stream when we return from this handler,
            # instead of waiting for more data on a keep-alive socket and
            # tripping the 120s AbortController on the frontend.
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Connection", "close")
            self.send_header("X-Accel-Buffering", "no")
            self._set_cors()
            self.end_headers()
            # Also signal to BaseHTTPRequestHandler that this connection
            # must NOT be reused after this response.
            self.close_connection = True
            try:
                for event_name, payload in chat_with_agent_stream(
                        msg, graph_state, history,
                        anthropic_api_key=user_key,
                        use_websearch=use_websearch,
                        ncbi_api_key=self.server.ncbi_api_key,   # type: ignore[attr-defined]
                        ncbi_email=self.server.ncbi_email,       # type: ignore[attr-defined]
                        hub_dir=self.hub_dir):
                    sse = f"event: {event_name}\ndata: {json.dumps(payload, separators=(',', ':'))}\n\n"
                    self.wfile.write(sse.encode("utf-8"))
                    self.wfile.flush()
            except AgentError as e:
                err_lower = str(e).lower()
                # Distinguish provider 503/overload from real rate-limit so
                # the UI can suggest "Google is busy, retry" vs "you hit your
                # quota, wait or upgrade". 503 / unavailable / overloaded are
                # provider-side capacity issues; 429 / rate / quota are user-
                # side or per-key caps.
                if "auth" in err_lower or "invalid" in err_lower:
                    code = "auth_failed"
                elif "unavailable" in err_lower or "overloaded" in err_lower or "503" in err_lower:
                    code = "service_unavailable"
                elif "rate" in err_lower or "limit" in err_lower or "quota" in err_lower or "429" in err_lower:
                    code = "rate_limited"
                else:
                    code = "agent_failed"
                sse = f"event: error\ndata: {json.dumps({'error': code, 'detail': str(e)})}\n\n"
                try:
                    self.wfile.write(sse.encode("utf-8")); self.wfile.flush()
                except Exception:
                    pass
            except Exception as e:
                # Any non-AgentError exception (network blip, NCBI timeout,
                # serialization error, etc.) — surface it as an SSE error
                # event so the client doesn't just see an empty stream.
                import traceback
                log.error("chat-stream unexpected error: %s\n%s",
                          e, traceback.format_exc())
                sse = f"event: error\ndata: {json.dumps({'error': 'server_error', 'detail': str(e)[:300]})}\n\n"
                try:
                    self.wfile.write(sse.encode("utf-8")); self.wfile.flush()
                except Exception:
                    pass
            log.info("chat-stream → done")
            return

        # POST /citations/credible — user votes a PMID as credible (up) or
        # not credible (down). Append-only jsonl; same (pmid, IP) can flip
        # vote later — latest wins on aggregate.
        if path == "/citations/credible":
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length) if length else b""
            try:
                body = json.loads(raw) if raw else {}
            except json.JSONDecodeError:
                self._json(400, {"error": "bad json"}); return
            pmid = str(body.get("pmid") or "").strip()
            if not pmid or not re.fullmatch(r"\d{1,10}", pmid):
                self._json(400, {"error": "valid pmid required"}); return
            vote = str(body.get("vote") or "up").strip().lower()
            if vote not in ("up", "down"):
                self._json(400, {"error": "vote must be 'up' or 'down'"}); return
            context = str(body.get("context") or "").strip()[:300]
            title = str(body.get("title") or "").strip()[:300]
            client_ip = (self.headers.get("CF-Connecting-IP")
                         or (self.headers.get("X-Forwarded-For", "").split(",")[0].strip() or None)
                         or self.client_address[0])
            entry = {
                "pmid": pmid,
                "ip_hash": _hash_ip(client_ip),
                "vote": vote,
                "ts": __import__("datetime").datetime.utcnow().isoformat() + "Z",
                "context": context,
                "title": title,
            }
            p = _credible_path()
            try:
                with p.open("a") as fh:
                    fh.write(json.dumps(entry, separators=(",", ":")) + "\n")
            except OSError as e:
                self._json(500, {"error": f"write failed: {e}"}); return
            agg = _credible_aggregate().get(pmid, {"up": 0, "down": 0})
            log.info("credible vote pmid=%s vote=%s up=%d down=%d",
                     pmid, vote, agg.get("up", 0), agg.get("down", 0))
            self._json(200, {"pmid": pmid, "vote": vote,
                             "up": agg.get("up", 0), "down": agg.get("down", 0)})
            return

        # POST /citations/curate — user-triggered summarization. Fetches
        # the PubMed abstract and asks the BYOK LLM for a 3-4 sentence
        # digest + key findings, persists to /var/lib/heartgrn/
        # curated_knowledge/<pmid>.json. Synchronous (~10-60s); the
        # caller is expected to show a spinner. Body: {"pmid": "...",
        # "title?": "...", "context?": {...}}.
        if path == "/citations/curate":
            client_ip = (self.headers.get("CF-Connecting-IP")
                         or (self.headers.get("X-Forwarded-For", "").split(",")[0].strip() or None)
                         or self.client_address[0])
            allowed, retry_after = rate_limit_check(client_ip, "/chat")
            if not allowed:
                self._json(429, {"error": "rate_limited",
                                  "message": f"Too many requests; retry after {retry_after}s",
                                  "retry_after": retry_after}); return
            user_key = self.headers.get("X-API-Key", "").strip()
            if not user_key or not (user_key.startswith("sk-ant-") or user_key.startswith("AIza")):
                self._json(401, {"error": "byok_required",
                                  "message": "API key must start with 'sk-ant-' or 'AIza'"}); return
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length) if length else b""
            try:
                body = json.loads(raw) if raw else {}
            except json.JSONDecodeError:
                self._json(400, {"error": "bad json"}); return
            pmid = str(body.get("pmid") or "").strip()
            if not pmid or not re.fullmatch(r"\d{1,10}", pmid):
                self._json(400, {"error": "valid pmid required"}); return
            existing = _load_curated(pmid)
            force = bool(body.get("force"))
            if existing and not force:
                self._json(200, {"status": "exists", "pmid": pmid,
                                  "summary": existing.get("summary"),
                                  "key_findings": existing.get("key_findings", []),
                                  "title": existing.get("title")})
                return
            citation_meta = {
                "pmid": pmid,
                "title": str(body.get("title") or "")[:300],
                "journal": str(body.get("journal") or "")[:80],
                "year": body.get("year"),
            }
            context = body.get("context") or {}
            if not isinstance(context, dict):
                context = {}
            try:
                _summarize_and_persist(
                    pmid, citation_meta, context,
                    api_key=user_key,
                    ncbi_api_key=self.server.ncbi_api_key,   # type: ignore[attr-defined]
                    ncbi_email=self.server.ncbi_email)        # type: ignore[attr-defined]
            except Exception as e:
                self._json(500, {"error": f"summarize failed: {e}"}); return
            saved = _load_curated(pmid)
            if not saved:
                self._json(500, {"error": "summary not saved (LLM may have returned unparseable JSON; check server logs)"}); return
            log.info("curate POST → saved pmid=%s by ip=%s", pmid, client_ip)
            self._json(200, {"status": "saved", "pmid": pmid,
                              "summary": saved.get("summary"),
                              "key_findings": saved.get("key_findings", []),
                              "title": saved.get("title")})
            return

        m = re.fullmatch(r"/literature/(\w+)/(.+)", path)
        if not m:
            self._json(404, {"error": "not found"}); return
        type_, key = m.group(1), m.group(2)
        if type_ not in ALLOWED_TYPES:
            self._json(400, {"error": f"unknown type {type_!r}"}); return
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length) if length else b""
        try:
            body = json.loads(raw) if raw else {}
        except json.JSONDecodeError as e:
            self._json(400, {"error": f"bad json body: {e}"}); return
        ctx = body.get("context", {})
        ctx.setdefault("rsid" if type_ == "snp" else "symbol", key)
        force = bool(body.get("force", False))
        synthesize = bool(body.get("synthesize", False))
        mode = body.get("mode", "pubmed")           # "pubmed" | "agent" | "hybrid"
        max_results = int(body.get("max_results", 15))
        min_year = body.get("min_year")

        if mode not in ("pubmed", "agent", "hybrid"):
            self._json(400, {"error": f"invalid mode {mode!r}"}); return

        context_extra = f"{ctx.get('disease','')}__{ctx.get('target_cs','')}__{mode}"
        cp = cache_path(self.hub_dir, type_, key, context_extra)
        if cp.exists() and not force:
            cached_blob = _load_json_safe(cp)
            if cached_blob is not None:
                self._json(200, {"cached": True, **cached_blob}); return
            # corrupt/truncated cache → fall through to a fresh build

        if mode == "agent":
            # BYOK required for LLM-driven literature search
            user_key = self.headers.get("X-API-Key", "").strip()
            if not user_key or not user_key.startswith("sk-ant-"):
                self._json(401, {"error": "byok_required_for_agent_mode",
                                 "message": "agent/hybrid modes require X-API-Key header (Anthropic key starting sk-ant-)"}); return
            try:
                agent_out = agent_search(type_, ctx,
                                          anthropic_api_key=user_key,
                                          max_results=max_results)
            except AgentError as e:
                self._json(502, {"error": "agent failed", "detail": str(e)}); return
            cits = agent_out.get("citations", []) if isinstance(agent_out, dict) else []
            verified = verify_pmids(cits,
                                    api_key=self.server.ncbi_api_key,   # type: ignore[attr-defined]
                                    email=self.server.ncbi_email)        # type: ignore[attr-defined]
            entry = {
                "type": type_, "key": key, "context": ctx, "mode": mode,
                "agent_citations": verified,
                "agent_summary": agent_out.get("summary") if isinstance(agent_out, dict) else None,
                "agent_broader": agent_out.get("broader_context") if isinstance(agent_out, dict) else None,
                "agent_raw": agent_out.get("_raw") if isinstance(agent_out, dict) and agent_out.get("_parse_error") else None,
                "ts": __import__("datetime").datetime.utcnow().isoformat() + "Z",
            }
            _atomic_write_text(cp, json.dumps(entry, indent=2))
            n_verified = sum(1 for c in verified if c.get("verified"))
            log.info("agent → %s (%d cit, %d verified)",
                     cp.relative_to(self.hub_dir), len(verified), n_verified)
            self._json(200, {"cached": False, **entry})
            return

        try:
            candidates = build_queries(type_, ctx)
            if not candidates:
                self._json(400, {"error": "could not build query — no rsid/symbol in context"}); return
            api_key = self.server.ncbi_api_key  # type: ignore[attr-defined]
            email = self.server.ncbi_email      # type: ignore[attr-defined]
            pm = None
            tried: list[dict] = []
            chosen_label = None
            chosen_query = None
            for label, q in candidates:
                log.info("PubMed try [%s/%s] %s: %s", type_, key, label, q)
                res = pubmed_search(q, max_results=max_results,
                                    min_year=min_year, api_key=api_key,
                                    email=email)
                tried.append({"label": label, "query": q, "count": res["count"]})
                time.sleep(0.34 if not api_key else 0.11)
                if res["citations"]:
                    pm = res
                    chosen_label = label
                    chosen_query = q
                    break
            if pm is None:
                pm = {"citations": [], "count": 0, "query": candidates[-1][1]}
                chosen_query = candidates[-1][1]
                chosen_label = candidates[-1][0]
            # Europe PMC search with the loosest query (most likely to find hits)
            epmc_q = candidates[-1][1]
            try:
                epmc = europe_pmc_search(epmc_q, max_results=8)
            except Exception as e:
                log.warning("Europe PMC failed: %s", e)
                epmc = {"citations": [], "error": str(e)}
        except urllib.error.HTTPError as e:
            self._json(502, {"error": "pubmed http error",
                              "detail": f"{e.code} {e.reason}"}); return
        except Exception as e:
            self._json(500, {"error": "search failed", "detail": str(e)}); return

        # Dedup PubMed vs EuropePMC by PMID
        pmid_seen = {c["pmid"] for c in pm["citations"] if c.get("pmid")}
        extra_epmc = [c for c in epmc["citations"]
                      if c.get("pmid") not in pmid_seen]

        synthesis = ""
        synth_err = None
        if synthesize and pm["citations"]:
            user_key = self.headers.get("X-API-Key", "").strip()
            if user_key and user_key.startswith("sk-ant-"):
                try:
                    synthesis = synthesize_via_claude(pm["citations"], ctx,
                                                       anthropic_api_key=user_key)
                except AgentError as e:
                    synth_err = str(e)
            else:
                synth_err = "byok_required: synthesis needs X-API-Key (sk-ant-...)"

        entry = {
            "type": type_, "key": key, "context": ctx, "mode": mode,
            "query": chosen_query,
            "query_label": chosen_label,
            "tried_queries": tried,
            "pubmed": {"citations": pm["citations"], "count_total": pm["count"]},
            "europe_pmc": {"citations": extra_epmc[:6]},
            "synthesis": synthesis,
            "synthesis_error": synth_err,
            "ts": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        }
        # Hybrid: run agent on top of PubMed results to gather preprints/news (BYOK)
        if mode == "hybrid":
            user_key = self.headers.get("X-API-Key", "").strip()
            if not user_key or not user_key.startswith("sk-ant-"):
                entry["agent_error"] = "byok_required: hybrid mode needs X-API-Key (sk-ant-...)"
            else:
                try:
                    agent_out = agent_search(type_, ctx,
                                              anthropic_api_key=user_key,
                                              max_results=max_results,
                                              timeout=180)
                    cits = agent_out.get("citations", []) if isinstance(agent_out, dict) else []
                    verified = verify_pmids(cits,
                                            api_key=self.server.ncbi_api_key,   # type: ignore[attr-defined]
                                            email=self.server.ncbi_email)        # type: ignore[attr-defined]
                    entry["agent_citations"] = verified
                    entry["agent_summary"] = agent_out.get("summary") if isinstance(agent_out, dict) else None
                    entry["agent_broader"] = agent_out.get("broader_context") if isinstance(agent_out, dict) else None
                except AgentError as e:
                    entry["agent_error"] = _scrub_key(e)
        _atomic_write_text(cp, json.dumps(entry, indent=2))
        log.info("→ %s (PubMed %d, Europe PMC extra %d)",
                 cp.relative_to(self.hub_dir),
                 len(pm["citations"]), len(extra_epmc))
        self._json(200, {"cached": False, **entry})


# ---- entry point ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--hub-dir", type=Path,
                    default=Path(__file__).resolve().parent.parent / "hub")
    ap.add_argument("--bind", default="127.0.0.1")
    ap.add_argument("--ncbi-api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"),
                    help="email for NCBI E-utilities compliance")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    hub_dir = args.hub_dir.resolve()
    if not hub_dir.exists():
        sys.exit(f"hub_dir does not exist: {hub_dir}")
    (hub_dir / "literature").mkdir(exist_ok=True)
    for t in ALLOWED_TYPES:
        (hub_dir / "literature" / t).mkdir(exist_ok=True)

    server = ThreadingHTTPServer((args.bind, args.port), LitHandler)
    server.hub_dir = hub_dir                          # type: ignore[attr-defined]
    server.ncbi_api_key = args.ncbi_api_key           # type: ignore[attr-defined]
    server.ncbi_email = args.email                    # type: ignore[attr-defined]

    log.info("Literature backend (PubMed direct) listening on http://%s:%d "
             "(hub=%s, ncbi_api_key=%s)",
             args.bind, args.port, hub_dir, "set" if args.ncbi_api_key else "no")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        log.info("shutdown"); server.server_close()


if __name__ == "__main__":
    main()
