#!/usr/bin/env python3
"""
Enrich each payload SNP with the credibleSet study_id from
hub/snp_annotations.json. No external API hits — pure local join.

For each SNP s in each payload:
  variant_id = f"{s.chr}_{s.pos}_{s.ref}_{s.alt}"
  cs_list    = annotations[variant_id].ot_variant.credibleSets
  match cs_list rows to s.trait by substring / word overlap.
  - exactly 1 best match  → s.study_id = cs.study_id
  - >1 tied best matches  → s.study_id_candidates = [...]
                            s.study_id = None
  - 0 matches             → s.study_id = None
                            s.study_id_note = "no_trait_match"

Atomic per-file write (temp file + os.replace). Backups not kept since
git history is the rollback.
"""
import argparse, json, os, re, sys, time
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
PAYLOAD_DIR = ROOT / "hub/payloads"
ANN_FILE    = ROOT / "hub/snp_annotations.json"


def _norm(s: str) -> str:
    """Lowercase, collapse whitespace, drop parenthetical suffix."""
    s = (s or "").lower().strip()
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _match_score(payload_trait_norm: str, cs_trait_norm: str) -> int:
    """0=no match, 1=word overlap, 2=substring containment (preferred)."""
    if not payload_trait_norm or not cs_trait_norm:
        return 0
    if payload_trait_norm in cs_trait_norm or cs_trait_norm in payload_trait_norm:
        return 2
    pw = {w for w in re.findall(r"[a-z]+", payload_trait_norm) if len(w) > 3}
    cw = {w for w in re.findall(r"[a-z]+", cs_trait_norm) if len(w) > 3}
    return 1 if (pw & cw) else 0


def _atomic_write(path: Path, content: str) -> None:
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}.{os.urandom(4).hex()}")
    try:
        tmp.write_text(content)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try: tmp.unlink()
            except OSError: pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="don't write files")
    args = ap.parse_args()

    if not ANN_FILE.exists():
        sys.exit(f"FATAL: {ANN_FILE} not found")
    if not PAYLOAD_DIR.exists():
        sys.exit(f"FATAL: {PAYLOAD_DIR} not found")

    print(f"Loading annotations from {ANN_FILE.name}...")
    ann_data = json.loads(ANN_FILE.read_text())
    ann = ann_data.get("annotations") or {}
    print(f"  {len(ann)} variant annotations loaded")
    print()

    files = sorted(PAYLOAD_DIR.glob("*.json"))
    print(f"Scanning {len(files)} payloads in {PAYLOAD_DIR}...")
    print()

    total_snps = 0
    total_with_study = 0
    total_ambiguous = 0
    total_no_match = 0
    total_no_ann = 0
    per_outcome = Counter()
    t0 = time.time()

    for pf in files:
        try:
            p = json.loads(pf.read_text())
        except Exception as e:
            print(f"  ! parse error {pf.name}: {e}")
            continue

        snps = p.get("snps") or []
        if not snps:
            continue
        changed = False

        for s in snps:
            total_snps += 1
            chr_ = s.get("chr"); pos = s.get("pos")
            ref = s.get("ref"); alt = s.get("alt")
            if not (chr_ and pos and ref and alt):
                continue   # build_snp_index now skips these too
            vid = f"{chr_}_{pos}_{ref}_{alt}"
            entry = ann.get(vid)
            if not entry:
                # No annotation row at all — leave a marker so we don't keep
                # re-trying on every run.
                if s.get("study_id") is None and "study_id" not in s:
                    s["study_id"] = None
                    s["study_id_note"] = "no_annotation"
                    changed = True
                total_no_ann += 1
                per_outcome["no_annotation"] += 1
                continue
            cs_list = ((entry.get("ot_variant") or {}).get("credibleSets") or [])
            if not cs_list:
                s["study_id"] = None
                s["study_id_note"] = "no_credible_sets"
                changed = True
                total_no_match += 1
                per_outcome["no_credible_sets"] += 1
                continue
            p_trait = _norm(s.get("trait"))
            best_score = 0
            best_rows: list[dict] = []
            for cs in cs_list:
                sc = _match_score(p_trait, _norm(cs.get("trait")))
                if sc > best_score:
                    best_score = sc
                    best_rows = [cs]
                elif sc == best_score and sc > 0:
                    best_rows.append(cs)

            if best_score == 0:
                # No trait match at all — pick the only credibleSet if there's
                # exactly one, otherwise leave ambiguous
                if len(cs_list) == 1:
                    sid = cs_list[0].get("study_id")
                    s["study_id"] = sid
                    s["study_id_note"] = "single_cs_no_trait_match"
                    per_outcome["single_cs_only"] += 1
                else:
                    s["study_id"] = None
                    s["study_id_candidates"] = [c.get("study_id") for c in cs_list]
                    s["study_id_note"] = "no_trait_match"
                    total_no_match += 1
                    per_outcome["no_trait_match"] += 1
                changed = True
                continue

            if len(best_rows) == 1:
                sid = best_rows[0].get("study_id")
                s["study_id"] = sid
                # Re-write the trait we matched against (so the agent can
                # see the canonical OT trait string if it differs from the
                # original GWAS sumstat label).
                cs_trait = best_rows[0].get("trait")
                if cs_trait and cs_trait != s.get("trait"):
                    s["study_trait_ot"] = cs_trait
                total_with_study += 1
                per_outcome[f"unique_score{best_score}"] += 1
                changed = True
            else:
                # Multiple cs tied at the same score — record all candidates
                s["study_id"] = None
                s["study_id_candidates"] = [c.get("study_id") for c in best_rows]
                s["study_id_note"] = f"ambiguous_n{len(best_rows)}"
                total_ambiguous += 1
                per_outcome[f"ambiguous_score{best_score}"] += 1
                changed = True

        if changed and not args.dry_run:
            _atomic_write(pf, json.dumps(p, separators=(",", ":")))

    print("=" * 60)
    print(f"Done in {time.time() - t0:.1f}s")
    print(f"  Total SNP rows:        {total_snps}")
    print(f"  Unique study_id set:   {total_with_study}")
    print(f"  Ambiguous:             {total_ambiguous}")
    print(f"  No trait match:        {total_no_match}")
    print(f"  No annotation row:     {total_no_ann}")
    print()
    print("Per-outcome counts:")
    for k in sorted(per_outcome.keys()):
        print(f"  {k:30} {per_outcome[k]}")
    if args.dry_run:
        print("\n(dry-run — no files modified)")


if __name__ == "__main__":
    main()
