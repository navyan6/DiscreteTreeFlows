#!/usr/bin/env python3
"""
Sequence-level leakage check for the new VOC-origin roots.

Accession-level disjointness (already verified: 0 overlap) is necessary but not
sufficient. SARS-CoV-2 Spike sequences repeat constantly -- early B.1 genomes
differ from Wuhan-Hu-1 by little more than D614G -- so a "new" root can carry a
Spike string byte-identical to something the model trained on. That is not
automatically disqualifying under the tree-level standard in LEAKAGE.md, but it
has to be measured and reported rather than discovered by a reviewer.

Reports, per new tree:
  root_seq_in_train   - is the inferred root's Spike identical to a train Spike?
  leaf_identity_frac  - fraction of leaves whose Spike is identical to a train Spike
  novel_leaf_frac     - fraction the model has demonstrably never seen

This also becomes the gate for any GISAID top-up, where accession matching breaks
outright: GISAID uses EPI_ISL ids and many records are deposited in both
databases, so only sequence identity can detect the overlap.

Run on Betty after slurm_covid_voc_roots_pipeline.sh completes:
  python scripts/check_voc_roots_overlap.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from voc_threat_lib import seq_at  # noqa: E402


def iter_fasta(path: Path):
    name = None
    buf: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def digest(seq: str) -> str:
    return hashlib.sha1(seq.replace("-", "").upper().encode()).hexdigest()


def is_internal(name: str) -> bool:
    return name.startswith("NODE_") or name in ("root", "root|root")


def parse_mut(m: str) -> tuple[int, str]:
    i = 1
    while i < len(m) and m[i].isdigit():
        i += 1
    return int(m[1:i]), m[i:]


def load_roles(raw_dir: Path) -> dict[str, str]:
    import csv

    roles: dict[str, str] = {}
    for p in raw_dir.glob("*_metadata.csv"):
        with p.open() as fh:
            for r in csv.DictReader(fh):
                roles[r["acc_base"]] = r.get("role", "")
    return roles


def load_pop_vocs(summary: Path) -> dict[str, str]:
    if not summary.exists():
        return {}
    return {
        s["population"]: s.get("voc_id")
        for s in json.loads(summary.read_text())
        if s.get("voc_id")
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-dir", type=Path, default=ROOT / "data" / "covid" / "train")
    ap.add_argument("--roots-base", type=Path, default=ROOT / "data" / "covid_voc_roots")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "voc_threat_panel" / "voc_roots_seq_overlap.json")
    args = ap.parse_args()

    train_hashes: set[str] = set()
    n_train = 0
    for p in sorted(args.train_dir.glob("group_*_anc_aa.fasta")):
        for name, seq in iter_fasta(p):
            train_hashes.add(digest(seq))
            n_train += 1
    print(f"indexed {n_train} training sequences -> {len(train_hashes)} distinct Spike strings")
    if not train_hashes:
        print("no training sequences found; check --train-dir")
        return

    # Aggregate identity is dominated by the pre-emergence background, where a
    # shared Spike is unavoidable and harmless. The number that decides whether a
    # tree can carry a recovery claim is whether the *VOC-carrying* leaves are novel.
    roles = load_roles(args.roots_base / "raw")
    vocs = load_pop_vocs(args.roots_base / "raw" / "pull_summary.json")
    panel = json.loads((ROOT / "benchmarks" / "voc_threat_panel.json").read_text())
    bundles = {
        v["id"]: ([parse_mut(m) for m in v["defining_bundle"]],
                  int(v.get("bundle_min_hits", len(v["defining_bundle"]))))
        for v in panel["vocs"]
    }

    results = []
    print(f"\n{'population':<20}{'tree':<10}{'leaves':>7}{'all':>7}{'bundle':>8}"
          f"{'n_bndl':>8}  root_seen")
    for fasta in sorted(args.roots_base.glob("*/group_*_anc_aa.fasta")):
        pop = fasta.parent.name
        recs = list(iter_fasta(fasta))
        leaves = {n: s for n, s in recs if not is_internal(n)}
        internals = {n: s for n, s in recs if is_internal(n)}
        if not leaves:
            continue
        hit = sum(1 for s in leaves.values() if digest(s) in train_hashes)

        voc = vocs.get(pop)
        bundle, min_hits = bundles.get(voc, ([], 0))
        # seq_at maps reference positions through an alignment, so deletion-
        # carrying lineages (Alpha, Beta, Delta, BA.1) are read at the right
        # residue rather than a shifted one. See SPIKE_FRAME_BUG.md.
        bundle_leaves = {
            n: s for n, s in leaves.items()
            if bundle and sum(1 for pos, aa in bundle if seq_at(s, pos) == aa) >= min_hits
        }
        b_hit = sum(1 for s in bundle_leaves.values() if digest(s) in train_hashes)
        b_frac = round(b_hit / len(bundle_leaves), 4) if bundle_leaves else None

        root_seq = internals.get("NODE_0000000")  # augur's root label
        results.append(
            {
                "population": pop,
                "voc": voc,
                "tree": fasta.name,
                "n_leaves": len(leaves),
                "leaf_identity_frac": round(hit / len(leaves), 4),
                "n_bundle_leaves": len(bundle_leaves),
                "bundle_leaf_identity_frac": b_frac,
                "novel_bundle_leaves": len(bundle_leaves) - b_hit,
                "root_seq_in_train": (digest(root_seq) in train_hashes) if root_seq else None,
                "usable_for_recovery_claim": bool(
                    bundle_leaves and (len(bundle_leaves) - b_hit) > 0
                ),
            }
        )
        gid = fasta.name.split("_")[1]
        bs = f"{b_frac:.3f}" if b_frac is not None else "   -  "
        print(
            f"  {pop:<18}{gid:<10}{len(leaves):>7}{hit / len(leaves):>7.3f}{bs:>8}"
            f"{len(bundle_leaves):>8}  {results[-1]['root_seq_in_train']}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {"n_train_seqs": n_train, "n_distinct_train_spikes": len(train_hashes),
             "trees": results},
            indent=2,
        )
        + "\n"
    )
    print(f"\nwrote {args.out.relative_to(ROOT)}")
    if not results:
        print("no trees yet -- run slurm_covid_voc_roots_pipeline.sh first")


if __name__ == "__main__":
    main()
