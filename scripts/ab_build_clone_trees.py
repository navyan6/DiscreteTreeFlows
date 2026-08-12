#!/usr/bin/env python3
"""
Table 5 — MAFFT + FastTree (AA) + midpoint root + leaf-root ASR stub.

For each clone FASTA under clones/by_clone/*.fasta:
  1. MAFFT --auto
  2. FastTree -lg (protein)
  3. Midpoint root (dendropy)
  4. Root sequence = majority-vote AA at aligned columns among leaves
     (germline ASR placeholder until IgPhyML / true germline reconstruct)

Writes data layout compatible with TreeDataset-ish groups:
  out/{train,val,test}/group_XXX.{fasta,nwk}  (after split script moves them)
  or out/all/group_XXX.* before split.

Usage:
  python scripts/ab_build_clone_trees.py \\
    --clones-dir data/ab_t5_500k/clones \\
    --out-dir data/ab_t5_500k/trees_all \\
    --max-clones 400
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(cmd: list[str], capture: bool = False) -> str:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed: {' '.join(cmd)}\n{r.stderr[-1500:]}")
    return r.stdout if capture else ""


def majority_root(aligned_fasta: Path) -> str:
    seqs = []
    header = None
    parts: list[str] = []
    with aligned_fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    seqs.append("".join(parts))
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            seqs.append("".join(parts))
    if not seqs:
        return ""
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        counts: dict[str, int] = {}
        for s in seqs:
            if i >= len(s):
                continue
            aa = s[i]
            if aa in "-.":
                continue
            counts[aa] = counts.get(aa, 0) + 1
        if not counts:
            out.append("-")
        else:
            out.append(max(counts, key=counts.get))
    return "".join(out).replace("-", "")


def midpoint_root_newick(nwk_text: str) -> str:
    try:
        import dendropy
    except ImportError as e:
        raise SystemExit("dendropy required for midpoint root") from e
    tree = dendropy.Tree.get(data=nwk_text, schema="newick")
    tree.reroot_at_midpoint(update_bipartitions=False)
    return tree.as_string(schema="newick").strip()


def process_clone(
    fasta: Path,
    out_dir: Path,
    group_idx: int,
    mafft: str,
    fasttree: str,
) -> dict | None:
    g = f"group_{group_idx:03d}"
    work = out_dir / g
    work.mkdir(parents=True, exist_ok=True)
    aligned = work / f"{g}_aligned.fasta"
    nwk_raw = work / f"{g}_unrooted.nwk"
    nwk = work / f"{g}.nwk"
    leaf_fa = work / f"{g}.fasta"
    root_fa = work / f"{g}_root.fasta"

    # copy leaves
    shutil.copy(fasta, leaf_fa)

    if not shutil.which(mafft):
        raise SystemExit(f"mafft not found ({mafft})")
    stdout = run([mafft, "--auto", "--quiet", str(leaf_fa)], capture=True)
    aligned.write_text(stdout)

    ft = shutil.which(fasttree) or shutil.which("FastTree") or shutil.which("fasttree")
    if not ft:
        raise SystemExit("FastTree not found on PATH")
    # protein LG
    r = subprocess.run([ft, "-lg", str(aligned)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-1000:])
    nwk_raw.write_text(r.stdout.strip() + "\n")
    rooted = midpoint_root_newick(r.stdout)
    nwk.write_text(rooted + "\n")

    root_seq = majority_root(aligned)
    root_fa.write_text(f">ROOT\n{root_seq}\n")
    # also write anc_aa placeholder = root
    (work / f"{g}_anc_aa.fasta").write_text(f">NODE_ROOT\n{root_seq}\n")

    n_leaves = sum(1 for line in leaf_fa.open() if line.startswith(">"))
    return {
        "group": g,
        "clone_fasta": str(fasta),
        "n_leaves": n_leaves,
        "root_len": len(root_seq),
        "nwk": str(nwk),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clones-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--max-clones", type=int, default=400)
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--fasttree", default="FastTree")
    args = ap.parse_args()

    by_clone = args.clones_dir / "by_clone"
    fastas = sorted(by_clone.glob("*.fasta"))
    if not fastas:
        raise SystemExit(f"No clone FASTAs under {by_clone}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for i, fa in enumerate(fastas[: args.max_clones]):
        print(f"[{i+1}/{min(len(fastas), args.max_clones)}] {fa.name}", flush=True)
        try:
            rec = process_clone(fa, args.out_dir, i, args.mafft, args.fasttree)
            if rec:
                results.append(rec)
        except Exception as e:
            print(f"  FAIL {fa.name}: {e}", flush=True)

    summary = {
        "n_clone_fastas": len(fastas),
        "n_built": len(results),
        "out_dir": str(args.out_dir),
        "rooting": "midpoint",
        "asr": "majority_leaf_vote_placeholder",
        "note": "Germline ASR / IgPhyML not wired — replace root before FINAL Table 5.",
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
