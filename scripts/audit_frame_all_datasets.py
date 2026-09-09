#!/usr/bin/env python3
"""
How far does the reference-frame problem reach across the paper's datasets?

Two numbers decide it per dataset:

  root_len_modal_frac  - fraction of trees whose root is the modal (reference)
                         length. TreeSBM only substitutes, so generated leaves
                         inherit the root's frame; a root off-reference puts an
                         entire generated tree in the wrong frame.
  leaf_frame_eq_root   - fraction of observed leaves in the same frame as their
                         root. Column-wise generated-vs-observed comparison is
                         only meaningful for these.

A dataset with no indels shows ~1.0 for both and needs no correction.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

SHIFT_HAMMING = 50

LAB = "/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data"

DATASETS = [
    ("SARS-CoV-2 Spike (test/Brazil)", "data/covid/test"),
    ("SARS-CoV-2 Spike (train)", "data/covid/train"),
    ("Influenza H3N2 HA", "data/h3n2/test"),
    ("Influenza H1N1 HA", "data/h1n1/test"),
    ("HIV Env (geo test)", f"{LAB}/hiv_geo/test"),
    ("HIV Env (test)", f"{LAB}/hiv/test"),
    ("Antibody clones", f"{LAB}/ab_clones/test"),
]


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


def is_leaf(name: str) -> bool:
    return not name.startswith("NODE_") and name not in ("root", "root|root")


def audit(data_dir: Path) -> dict | None:
    paths = sorted(data_dir.glob("group_*_anc_aa.fasta"))
    if not paths:
        return None
    root_lens = Counter()
    leaf_lens = Counter()
    n_leaves = same_frame = 0
    per_tree = []
    for path in paths:
        recs = list(iter_fasta(path))
        internals = {k: v for k, v in recs if not is_leaf(k)}
        leaves = [v for k, v in recs if is_leaf(k)]
        if not leaves or not internals:
            continue
        root = internals.get("NODE_0000000") or next(iter(internals.values()))
        root_lens[len(root)] += 1
        ok = 0
        for s in leaves:
            leaf_lens[len(s)] += 1
            if len(s) == len(root) and sum(1 for a, b in zip(s, root) if a != b) <= SHIFT_HAMMING:
                ok += 1
        n_leaves += len(leaves)
        same_frame += ok
        per_tree.append({"tree": path.name, "root_len": len(root),
                         "frac_same_frame": round(ok / len(leaves), 4)})

    modal_len = leaf_lens.most_common(1)[0][0] if leaf_lens else None
    n_trees = sum(root_lens.values())
    return {
        "data_dir": str(data_dir),
        "n_trees": n_trees,
        "n_leaves": n_leaves,
        "modal_seq_len": modal_len,
        "root_len_modal_frac": round(root_lens.get(modal_len, 0) / max(n_trees, 1), 4),
        "leaf_len_modal_frac": round(leaf_lens.get(modal_len, 0) / max(n_leaves, 1), 4),
        "leaf_frame_eq_root": round(same_frame / max(n_leaves, 1), 4),
        "root_len_hist": dict(root_lens.most_common(6)),
        "per_tree": per_tree,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "voc_threat_panel" / "frame_audit_all_datasets.json")
    args = ap.parse_args()

    out = {}
    print(f"{'dataset':<34}{'trees':>7}{'modal_len':>11}{'roots@modal':>13}"
          f"{'leaves@modal':>14}{'leaf_frame=root':>17}")
    for label, rel in DATASETS:
        d = Path(rel) if rel.startswith('/') else ROOT / rel
        res = audit(d)
        if res is None:
            print(f"{label:<34}  (no group_*_anc_aa.fasta under {rel})")
            continue
        out[label] = res
        print(f"{label:<34}{res['n_trees']:>7}{str(res['modal_seq_len']):>11}"
              f"{res['root_len_modal_frac']:>13.3f}{res['leaf_len_modal_frac']:>14.3f}"
              f"{res['leaf_frame_eq_root']:>17.3f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
