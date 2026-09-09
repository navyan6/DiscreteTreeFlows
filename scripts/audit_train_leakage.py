#!/usr/bin/env python3
"""
Count VOC-bundle carriers in the COVID training split, using alignment-aware
mutation calling (voc_threat_lib.seq_at maps reference positions through indels).

Supersedes the first pass, which indexed sequences directly by reference
position and therefore undercounted every deletion-carrying lineage -- Alpha,
Beta, Delta and Omicron BA.1. See SPIKE_FRAME_BUG.md.

  python scripts/audit_train_leakage.py --split train
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from voc_threat_lib import load_panel, muts_present  # noqa: E402

TARGETS = ["Gamma", "Beta", "Delta", "Omicron_BA1", "Alpha", "Lambda", "Mu", "Zeta"]


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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data" / "covid" / "train")
    ap.add_argument("--split", default="train")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "voc_threat_panel" / "train_leakage_audit_fixed.json")
    args = ap.parse_args()

    panel = load_panel()
    vocs = {v["id"]: v for v in panel["vocs"]}
    stats = {
        t: {"groups_with_bundle": 0, "leaves_with_bundle": 0, "groups": {}} for t in TARGETS
    }
    total_leaves = 0
    groups = sorted(args.data_dir.glob("group_*_anc_aa.fasta"))
    print(f"scanning {len(groups)} {args.split} trees …")

    for n, path in enumerate(groups, 1):
        gid = int(path.name.split("_")[1])
        leaves = [s for name, s in iter_fasta(path) if is_leaf(name)]
        total_leaves += len(leaves)
        for t in TARGETS:
            voc = vocs[t]
            bundle = voc["defining_bundle"]
            min_hits = int(voc.get("bundle_min_hits", len(bundle)))
            hits = sum(1 for s in leaves if len(muts_present(s, bundle)) >= min_hits)
            if hits:
                stats[t]["groups_with_bundle"] += 1
                stats[t]["leaves_with_bundle"] += hits
                stats[t]["groups"][gid] = [hits, len(leaves)]
        if n % 50 == 0:
            print(f"  {n}/{len(groups)}")

    print(f"\n{args.split}: {len(groups)} trees, {total_leaves} leaves\n")
    print(f"{'VOC':<14}{'trees':>7}{'carrier leaves':>16}{'frac':>9}")
    for t in TARGETS:
        s = stats[t]
        print(
            f"{t:<14}{s['groups_with_bundle']:>7}{s['leaves_with_bundle']:>16}"
            f"{s['leaves_with_bundle'] / max(total_leaves, 1):>9.4f}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {"split": args.split, "n_trees": len(groups), "train_leaves": total_leaves,
             "stats": stats},
            indent=2,
        )
        + "\n"
    )
    print(f"\nwrote {args.out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
