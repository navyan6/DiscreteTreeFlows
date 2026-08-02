#!/usr/bin/env python3
"""
Post-process PhylaFlow sampled trees into a TreeSBM topology-prior pool.

PhylaFlow is posterior-basin transport in BHV space for a fixed alignment — it is
NOT root-conditioned forward generation. This adapter never reimplements the
flow; it only consumes trees produced by PhylaFlow's own sampling scripts and
writes anonymized bare newicks for TopologyPriorMethod (row: phylaflow_adapted).

Typical upstream (run inside PhylaFlow's env + data layout; see EXTERNAL.md):
    python scripts/evaluate_per_dataset_sample_kl.py \\
        --config configs/...yaml --checkpoint /path/to.ckpt \\
        --sample-config ... --output-dir samples/treesbm_N16 \\
        --num-samples 50 --num-datasets 8 --dump-trees

Then (treesbm or PhylaFlow env; needs ete3):
    python phylaflow_sample.py \\
        --input-dir samples/treesbm_N16/tree_dumps \\
        --ntips 16 --n-samples 300 \\
        --out /path/to/DiscreteTreeFlows/benchmarks/external_pools/sampled/phylaflow_N16.nwk

Also accepts a single .nwk (one topology per line) or a .json dump with a
`sampled_trees` list.

COPY into $LABHOME/baselines/PhylaFlow/ if convenient; paths are absolute-safe.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

try:
    from ete3 import Tree as EteTree
except ImportError as e:
    raise SystemExit(
        "phylaflow_sample.py needs ete3 (PhylaFlow env or: pip install ete3)\n"
        f"ImportError: {e}"
    ) from e


def _anonymize(newick: str, ntips: int) -> str | None:
    """Bare topology with leaves relabeled 0..ntips-1; None if leaf count ≠ ntips."""
    tree = EteTree(str(newick), format=1)
    leaves = tree.get_leaves()
    if len(leaves) != ntips:
        return None
    for i, leaf in enumerate(leaves):
        leaf.name = str(i)
    # format=9: topology + leaf names only (no branch lengths / internal labels)
    return tree.write(format=9)


def _collect_from_json(path: Path) -> list[str]:
    payload = json.loads(path.read_text())
    if isinstance(payload, dict):
        trees = payload.get("sampled_trees") or payload.get("trees") or []
        if isinstance(trees, str):
            trees = [trees]
        return [str(t) for t in trees]
    if isinstance(payload, list):
        return [str(t) for t in payload]
    return []


def _collect_from_nwk(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def collect_trees(inputs: list[Path]) -> list[str]:
    out: list[str] = []
    for path in inputs:
        if path.is_dir():
            for p in sorted(path.rglob("*")):
                if p.suffix.lower() == ".json":
                    out.extend(_collect_from_json(p))
                elif p.suffix.lower() in {".nwk", ".newick", ".tree"}:
                    out.extend(_collect_from_nwk(p))
            continue
        if path.suffix.lower() == ".json":
            out.extend(_collect_from_json(path))
        else:
            out.extend(_collect_from_nwk(path))
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Convert PhylaFlow samples → anonymized phylaflow_N{N}.nwk pool")
    ap.add_argument("--input", nargs="+", required=True,
                    help="PhylaFlow tree_dumps/ dir, .json dumps, and/or .nwk files")
    ap.add_argument("--input-dir", default=None,
                    help="Alias for a single dump directory (same as --input DIR)")
    ap.add_argument("--ntips", type=int, required=True)
    ap.add_argument("--n-samples", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True,
                    help="Output path, e.g. .../external_pools/sampled/phylaflow_N16.nwk")
    args = ap.parse_args()

    paths = [Path(p) for p in args.input]
    if args.input_dir:
        paths.append(Path(args.input_dir))
    for p in paths:
        if not p.exists():
            raise SystemExit(f"Input not found: {p}")

    raw = collect_trees(paths)
    if not raw:
        raise SystemExit(
            "No trees found. Produce dumps with PhylaFlow's "
            "evaluate_per_dataset_sample_kl.py --dump-trees (see EXTERNAL.md)."
        )

    accepted: list[str] = []
    skipped = 0
    for nw in raw:
        anon = _anonymize(nw, args.ntips)
        if anon is None:
            skipped += 1
            continue
        accepted.append(anon)

    if not accepted:
        raise SystemExit(
            f"Parsed {len(raw)} trees but none had exactly N={args.ntips} leaves "
            f"(skipped={skipped}). Check PhylaFlow case leaf counts."
        )

    rng = random.Random(args.seed)
    if len(accepted) > args.n_samples:
        accepted = rng.sample(accepted, args.n_samples)
    # de-dup while preserving order (identical topologies collapse)
    uniq, seen = [], set()
    for t in accepted:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    # if too few unique after dedup, allow repeats from accepted to hit n_samples
    while len(uniq) < min(args.n_samples, len(accepted)):
        pick = rng.choice(accepted)
        uniq.append(pick)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(uniq[: args.n_samples]) + "\n")
    print(
        f"wrote {min(len(uniq), args.n_samples)} topologies "
        f"(from {len(raw)} raw, skipped_wrong_N={skipped}) -> {out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
