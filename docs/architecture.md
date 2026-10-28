# CardioNav (HeartGRN Viewer) — online deployment architecture

## High-level diagram

```
┌──────────────────────────────────────────────────────────────────┐
│ User browser                                                     │
│  ├─ Cytoscape graph viewer (no key)                              │
│  └─ Agent chat panel (X-API-Key from localStorage)               │
└──────────────────────────────────────────────────────────────────┘
        │ HTTPS  https://www.heartcellatlas.org/cardionav/...
        ▼
┌──────────────────────────────────────────────────────────────────┐
│ Sanger load balancer — strips the /cardionav/ prefix             │
└──────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────┐
│ Single host                                                      │
│                                                                  │
│  nginx  (deploy/nginx.conf.sample)                               │
│   ├─ TLS (Let's Encrypt / certbot)                               │
│   ├─ HTTP Basic Auth on everything except /healthz               │
│   ├─ /            → static files from /opt/heartgrn/hub          │
│   └─ /chat, /chat/stream, /ag1/, /literature/, /gene/, /snp/,    │
│      /citations/, /healthz → proxy to 127.0.0.1:8766             │
│                                                                  │
│  hub_backend.py  (systemd: heartgrn-backend.service, :8766)      │
│   ├─ /chat, /chat/stream  ── LLM tool-use loop (BYOK)            │
│   ├─ /literature/{snp,tf,gene}/*  ── NCBI E-utilities / EPMC     │
│   ├─ /ag1/*  ── AlphaGenome scores (cache, optional on-demand)   │
│   ├─ /gene/*, /snp/*  ── cross-payload lookup                    │
│   └─ /citations/*  ── user-curated credible-PMID votes           │
└──────────────────────────────────────────────────────────────────┘
        │ HTTPS (with the user's key)
        ▼
┌──────────────────────────────────────────────────────────────────┐
│ Anthropic API (default claude-sonnet-4-6)                        │
│   or Google Gemini API (default gemini-2.5-flash)                │
└──────────────────────────────────────────────────────────────────┘
```

The frontend calls the backend on the same origin. `hub/assets/config.js`
sets `API_BASE` to `"/cardionav"` when the page was reached under that path,
and to `""` otherwise, so the same files work behind the Sanger load balancer
and when the host is reached directly.

## Key principles

1. **BYOK only** — the user provides their own Anthropic (`sk-ant-…`) or Gemini (`AIza…`) key for `/chat`; the provider is chosen from the key prefix
2. **Never persist user keys** — the backend forwards the key for one request and discards it; key-shaped strings are redacted from error messages
3. **Graph is open** — Cytoscape exploration requires no key (but the whole site sits behind Basic Auth during the preview)
4. **Citation grounding** — `/literature` queries NCBI E-utilities and Europe PMC server-side, so cited PMIDs are real records rather than LLM output
5. **No database** — state lives in files: literature cache, AG1 cache, and the credible-citation log (`CREDIBLE_CITATIONS_FILE`, default `/var/lib/heartgrn/credible_citations.jsonl`)

## Security

- `X-API-Key` header (not cookie) for the user's LLM key
- HTTP Basic Auth (nginx, `/etc/nginx/heartgrn.htpasswd`) on every path except `/healthz`
- CORS / Origin check: only origins listed in `ALLOWED_ORIGIN` / `CUSTOM_ORIGIN` (plus localhost) are accepted on POST
- Rate limits on `/chat`: nginx `limit_req` zone `chat` (recommended 30 req/min per IP, burst 10) and the backend's own per-IP limit (`CHAT_RATE_PER_MIN`, default 20)
- Response headers set by nginx: `X-Content-Type-Options`, `X-Frame-Options: DENY`, `Referrer-Policy: no-referrer`, HSTS
- No `Content-Security-Policy` header is set at present
- HTTP → HTTPS redirect (certbot-managed)
