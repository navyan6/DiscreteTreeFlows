#!/usr/bin/env python3
"""
Group the freshly pulled VOC-origin sequences into trees whose roots the model
has never seen (see results/voc_threat_panel/LEAKAGE.md).

Differs from prepare_covid_geo.py in one important way. That script chunks a
country's sequences sequentially by date, so a chunk is a narrow time slice; here
each tree must span from *before* the variant emerged to *after*, otherwise
augur's inferred root is already a VOC and there is nothing to forecast. So each
tree is a temporally stratified sample over the whole emergence window, with the
target lineage capped at `--target-frac` of leaves.

Building several disjoint trees per population (`--n-trees`) gives several
independent roots per VOC, which is what "run new diverse roots" needs.

Input : data/covid_voc_roots/raw/{pop}_covid_seqs.fasta  (+ _metadata.csv)
Output: data/covid_voc_roots/{pop}/vocroots{pop}_group_NNN.fasta (+ .csv)

Note the Spike CDS is not extracted here -- nextclade only exists on Betty, so
that step runs there (scripts/slurm_covid_voc_roots_pipeline.sh) before this
script's grouping is consumed by run_all_groups.py.

Example:
  python scripts/prepare_covid_voc_roots.py --n-trees 3 --group-size 300
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MIN_GROUP = 60


def load_meta(meta_path: Path) -> dict[str, dict]:
    with meta_path.open() as fh:
        return {r["acc_base"]: r for r in csv.DictReader(fh)}


def load_seqs(fasta: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    acc = None
    buf: list[str] = []
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split()[0].split(".")[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def stratified_trees(
    recs: list[dict], n_trees: int, size: int, target_frac: float, rng: random.Random
) -> list[list[dict]]:
    """Split records into n_trees disjoint, temporally stratified samples."""
    target = [r for r in recs if r["role"] == "target"]
    bg = [r for r in recs if r["role"] != "target"]
    rng.shuffle(target)
    rng.shuffle(bg)

    n_t = min(len(target), int(size * target_frac))
    n_b = size - n_t
    trees: list[list[dict]] = []
    ti = bi = 0
    for _ in range(n_trees):
        if len(target) - ti < n_t // 2 or len(bg) - bi < n_b // 2:
            break
        chunk = target[ti : ti + n_t] + bg[bi : bi + n_b]
        ti += n_t
        bi += n_b
        if len(chunk) < MIN_GROUP:
            break
        chunk.sort(key=lambda r: r["date"])
        trees.append(chunk)

    # Sparse populations (India, Peru) can't fill even one stratified tree; take
    # everything as a single tree rather than dropping the VOC from the panel.
    if not trees and len(recs) >= MIN_GROUP:
        chunk = sorted(recs, key=lambda r: r["date"])
        trees.append(chunk)
    return trees


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default="data/covid_voc_roots")
    ap.add_argument("--n-trees", type=int, default=3, help="independent roots per population")
    ap.add_argument("--group-size", type=int, default=300)
    ap.add_argument("--target-frac", type=float, default=0.4, help="max leaf frac of target lineage")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--require-spike", action="store_true",
                    help="skip populations lacking {pop}_spike.fasta (use on Betty)")
    ap.add_argument("--seen-accessions", type=Path,
                    default=ROOT / "results" / "voc_threat_panel" / "seen_accessions.txt")
    args = ap.parse_args()

    base = ROOT / args.base
    raw = base / "raw"
    seen = {
        line.strip().split(".")[0]
        for line in args.seen_accessions.read_text().splitlines()
        if line.strip()
    } if args.seen_accessions.exists() else set()

    rng = random.Random(args.seed)
    manifest = []
    for fasta in sorted(raw.glob("*_covid_seqs.fasta")):
        pop = fasta.name.replace("_covid_seqs.fasta", "")
        meta = load_meta(raw / f"{pop}_metadata.csv")

        # The pipeline trains on the Spike CDS, so prefer covid_extract_spike.py's
        # output. That needs nextclade (Betty only); locally we fall back to whole
        # genomes just so the grouping logic can be exercised.
        spike = raw / f"{pop}_spike.fasta"
        if spike.exists():
            seqs = load_seqs(spike)
        elif args.require_spike:
            print(f"!! {pop}: missing {spike.name}; run covid_extract_spike.py first")
            continue
        else:
            print(f"   {pop}: no {spike.name}, using whole genomes (LOCAL DRY RUN ONLY)")
            seqs = load_seqs(fasta)

        recs = []
        collisions = 0
        for acc, seq in seqs.items():
            if acc in seen:
                collisions += 1
                continue
            m = meta.get(acc)
            if not m or len(m["date"]) < 10:
                continue
            recs.append({**m, "seq": seq})
        recs.sort(key=lambda r: r["date"])
        if collisions:
            print(f"!! {pop}: dropped {collisions} accessions that ARE in a training tree")

        trees = stratified_trees(recs, args.n_trees, args.group_size, args.target_frac, rng)
        out_dir = base / pop
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"vocroots{pop}"
        for gi, chunk in enumerate(trees, start=1):
            with (out_dir / f"{prefix}_group_{gi:03d}.fasta").open("w") as fh:
                for r in chunk:
                    fh.write(f">{r['acc_base']},{r['date']}\n{r['seq']}\n")
            with (out_dir / f"{prefix}_group_{gi:03d}.csv").open("w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["name", "date"])
                for r in chunk:
                    w.writerow([r["acc_base"], r["date"]])
            n_t = sum(1 for r in chunk if r["role"] == "target")
            manifest.append(
                {
                    "population": pop,
                    "prefix": prefix,
                    "gid": gi,
                    "data_dir": str(out_dir.relative_to(ROOT)),
                    "n_leaves": len(chunk),
                    "n_target_lineage": n_t,
                    "date_min": chunk[0]["date"],
                    "date_max": chunk[-1]["date"],
                    "pango_counts": dict(
                        sorted(
                            ((k, v) for k, v in _counts(chunk).items()),
                            key=lambda kv: -kv[1],
                        )
                    ),
                }
            )
            print(
                f"  {pop} g{gi:03d}: {len(chunk)} leaves "
                f"({n_t} target) {chunk[0]['date']}..{chunk[-1]['date']}"
            )
        if not trees:
            print(f"  {pop}: too few records for a tree ({len(recs)})")

    out = base / "voc_roots_manifest.json"
    out.write_text(json.dumps({"n_trees": len(manifest), "trees": manifest}, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)} ({len(manifest)} trees)")
    print("Next on Betty: sbatch scripts/slurm_covid_voc_roots_pipeline.sh")


def _counts(chunk: list[dict]) -> dict[str, int]:
    c: dict[str, int] = defaultdict(int)
    for r in chunk:
        c[r["pango"]] += 1
    return dict(c)


if __name__ == "__main__":
    main()
