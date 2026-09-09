#!/usr/bin/env python3
"""
Turn the stage-1 inventory into a worklist for the stage-2 array job.

One line per virus: taxid<TAB>slug<TAB>count. Slugs are filesystem-safe and
deduplicated, since NCBI species names collide once punctuation is stripped.

Viruses are ordered by sequence count descending, so if the array is cut short
the viruses with the most data are already done.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def slugify(name: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    return re.sub(r"_+", "_", s)[:48] or "virus"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", type=Path,
                    default=Path("data/panviral/virus_inventory.json"))
    ap.add_argument("--counts", type=Path,
                    default=Path("data/panviral/cache/counts.jsonl"),
                    help="fallback when the inventory has not been written yet")
    ap.add_argument("--out", type=Path, default=Path("data/panviral/worklist.tsv"))
    ap.add_argument("--min-count", type=int, default=150)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--skip-existing", action="store_true",
                    help="drop viruses that already have a manifest.json")
    args = ap.parse_args()

    # Merge both sources rather than preferring one. The inventory JSON is
    # rewritten only when stage 1 finishes, so a smoke-test run can leave a
    # stale file behind -- preferring it silently produced an empty worklist,
    # which would have made the array a no-op. Taking the union keyed on tax_id
    # is correct whether stage 1 has finished, is still running, or was
    # interrupted partway.
    merged: dict[int, dict] = {}
    sources = []

    if args.inventory.exists():
        try:
            for r in json.loads(args.inventory.read_text()).get("viruses", []):
                merged[int(r["tax_id"])] = r
            sources.append(f"{args.inventory} ({len(merged)})")
        except Exception as e:  # noqa: BLE001 - malformed/partial file
            sources.append(f"{args.inventory} (unreadable: {e})")

    if args.counts.exists():
        n = 0
        for line in args.counts.open():
            try:
                r = json.loads(line)
            except Exception:  # noqa: BLE001 - partial final line
                continue
            tid = int(r["tax_id"])
            if r.get("count", 0) > merged.get(tid, {}).get("count", -1):
                merged[tid] = r
                n += 1
        sources.append(f"{args.counts} ({n} new/updated)")

    if not merged:
        raise SystemExit("no inventory and no counts cache; run stage 1 first")

    src = " + ".join(sources)
    rows = [r for r in merged.values() if r.get("count", 0) >= args.min_count]
    rows.sort(key=lambda r: -r["count"])

    seen: set[str] = set()
    out_lines = []
    for r in rows:
        slug = slugify(r["name"])
        if slug in seen:
            slug = f"{slug}_{r['tax_id']}"
        seen.add(slug)
        if args.skip_existing and (args.out.parent / slug / "manifest.json").exists():
            continue
        out_lines.append(f"{r['tax_id']}\t{slug}\t{r['count']}")
        if args.limit and len(out_lines) >= args.limit:
            break

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"source: {src}")
    print(f"wrote {args.out} with {len(out_lines):,} viruses "
          f"(>= {args.min_count} sequences)")
    for line in out_lines[:15]:
        tid, slug, n = line.split("\t")
        print(f"  {int(n):>10,}  {slug}")


if __name__ == "__main__":
    main()
