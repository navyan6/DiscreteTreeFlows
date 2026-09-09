#!/usr/bin/env python3
"""
Build an H3N2 (or H1N1) named-clade trajectory figure parallel to COVID Gamma.

Picks a held-out test root whose observed descendants hit a named clade label
(from nextclade / GISAID clade fields when present; else uses antigenic-site
signature hits), finds the closest TreeSBM generated leaf, and writes:
  - mutation path table (root → gen leaf)
  - antigenic-site hit summary
  - markdown stats (mirrors results/covid_tree_viz/group_003_gamma_figure_stats.md)

Usage:
  python scripts/build_flu_variant_path_figure.py \
    --virus h3n2 \
    --ckpt checkpoints/h3n2_v3_lit_hotspot/best.pt \
    --data data/h3n2/test \
    --gen-dir results/h3n2_tree_viz \
    --out-dir results/h3n2_variant_path \
    --clade-name "3C.2a1b.2a.2"
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def hamming(a: str, b: str) -> int:
    L = min(len(a), len(b))
    return sum(1 for i in range(L) if a[i] != b[i]) + abs(len(a) - len(b))


def muts(parent: str, child: str) -> list[str]:
    L = min(len(parent), len(child))
    out = []
    for i in range(L):
        if parent[i] != child[i] and parent[i] != "-" and child[i] != "-":
            out.append(f"{parent[i]}{i+1}{child[i]}")
    return out


# H3 antigenic sites (Koel / WHO HA1 numbering ≈ our MSA columns for aligned HA)
H3_ANTIGENIC = {
    145, 155, 156, 158, 159, 189, 193, 225, 226, 227, 228
}


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    cur = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(cur).upper()
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
        if name is not None:
            seqs[name] = "".join(cur).upper()
    return seqs


def pick_named_clade_tip(obs: dict[str, str], clade_name: str) -> tuple[str, str] | None:
    """Return (tip_id, seq) if header/id contains clade_name; else longest tip."""
    clade_l = clade_name.lower()
    for k, v in obs.items():
        if clade_l in k.lower() or clade_l in k.replace("_", ".").lower():
            return k, v
    # Fallback: use tip with median Hamming to root later — return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--virus", default="h3n2", choices=["h3n2", "h1n1"])
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--gen-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--clade-name", default="3C.2a1b.2a.2")
    ap.add_argument("--group", default=None, help="e.g. group_046; auto-pick if unset")
    ap.add_argument("--ckpt", default=None)
    args = ap.parse_args()

    virus = args.virus
    data = args.data or Path(f"data/{virus}/test")
    gen_dir = args.gen_dir or Path(f"results/{virus}_tree_viz")
    out_dir = args.out_dir or Path(f"results/{virus}_variant_path")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover groups with both observed anc_aa and generated leaves
    groups = []
    for fa in sorted(data.glob("group_*_anc_aa.fasta")):
        g = fa.name.replace("_anc_aa.fasta", "")
        groups.append(g)
    if args.group:
        groups = [args.group] if args.group in groups or True else groups

    report = {
        "virus": virus,
        "clade_name": args.clade_name,
        "candidates": [],
        "selected": None,
    }

    best = None
    for g in groups[:40]:
        anc = data / f"{g}_anc_aa.fasta"
        if not anc.is_file():
            continue
        seqs = load_fasta(anc)
        # root often NODE_0000000 or first internal
        root_id = next((k for k in seqs if "NODE_0000000" in k or k.endswith("_root")), None)
        if root_id is None:
            # TreeTime: often the first NODE_* with most children — use shortest tip-distance proxy
            tips = {k: v for k, v in seqs.items() if not k.startswith("NODE_")}
            if not tips:
                continue
            root_id = min(tips, key=lambda k: len(tips[k]))  # weak fallback
            # Prefer an internal NODE if present
            nodes = [k for k in seqs if k.startswith("NODE_")]
            if nodes:
                root_id = sorted(nodes)[0]
        root_seq = seqs[root_id]
        tips = {k: v for k, v in seqs.items() if not str(k).startswith("NODE_")}
        named = pick_named_clade_tip(tips, args.clade_name)
        if named is None and tips:
            # pick tip farthest from root (most diverged) as "variant-like"
            tip_id = max(tips, key=lambda k: hamming(root_seq, tips[k]))
            tip_seq = tips[tip_id]
            named_source = "max_div_tip"
        elif named is not None:
            tip_id, tip_seq = named
            named_source = "clade_match"
        else:
            continue

        # Look for generated leaves fasta under gen_dir
        gen_fa = None
        for cand in [
            gen_dir / f"{g}_gen_leaves.fasta",
            gen_dir / f"{g}_matched_leaves.fasta",
            gen_dir / g / "gen_leaves.fasta",
        ]:
            if cand.is_file():
                gen_fa = cand
                break
        gen_leaves = load_fasta(gen_fa) if gen_fa else {}
        if not gen_leaves:
            report["candidates"].append(
                {
                    "group": g,
                    "named_source": named_source,
                    "tip": tip_id,
                    "obs_hamming": hamming(root_seq, tip_seq),
                    "gen": None,
                    "note": "no generated leaves yet — run slurm_h3n2_tree_viz.sh first",
                }
            )
            continue

        closest_id, closest_d = None, 10**9
        for gid, gs in gen_leaves.items():
            d = hamming(gs, tip_seq)
            if d < closest_d:
                closest_d, closest_id = d, gid
        closest_seq = gen_leaves[closest_id]
        path_muts = muts(root_seq, closest_seq)
        anti_hits = [m for m in path_muts if int(re.sub(r"[A-Z]", "", m) or 0) in H3_ANTIGENIC]
        cand = {
            "group": g,
            "named_source": named_source,
            "tip": tip_id,
            "obs_hamming": hamming(root_seq, tip_seq),
            "gen_id": closest_id,
            "gen_hamming_to_tip": closest_d,
            "n_path_muts": len(path_muts),
            "antigenic_hits": anti_hits,
            "path_muts_head": path_muts[:40],
        }
        report["candidates"].append(cand)
        score = closest_d  # lower better
        if best is None or score < best[0]:
            best = (score, cand, root_seq, tip_seq, closest_seq, path_muts)

    if best is None:
        # Still write a scaffold markdown pointing at tree_viz
        md = out_dir / f"{virus}_variant_path_PENDING.md"
        md.write_text(
            f"# {virus.upper()} named-clade path — PENDING generations\n\n"
            f"Target clade: `{args.clade_name}`\n\n"
            f"No generated leaves found under `{gen_dir}`.\n"
            f"Run `sbatch scripts/slurm_h3n2_tree_viz.sh` (or H1 analogue), then re-run this script.\n"
        )
        (out_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"Wrote pending scaffold {md}")
        return

    _, cand, root_seq, tip_seq, gen_seq, path_muts = best
    report["selected"] = cand
    (out_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    g = cand["group"]
    md = out_dir / f"{g}_{virus}_variant_path.md"
    lines = [
        f"# {virus.upper()} {g} — named-clade trajectory (parallel to COVID Gamma)",
        "",
        f"Target clade / VOC-analogue: **{args.clade_name}** (source={cand['named_source']})",
        "",
        "## Observed tip",
        f"- tip: `{cand['tip']}`",
        f"- Hamming(root → tip): **{cand['obs_hamming']}**",
        "",
        "## Closest TreeSBM leaf",
        f"- gen id: `{cand['gen_id']}`",
        f"- Hamming(gen → tip): **{cand['gen_hamming_to_tip']}**",
        f"- path mutations (root→gen): **{cand['n_path_muts']}**",
        f"- antigenic-site hits: {', '.join(cand['antigenic_hits']) or '(none)'}",
        "",
        "### Path mutations (first 40)",
        "",
        "```",
        ", ".join(path_muts[:40]),
        "```",
        "",
        f"Ckpt hint: `{args.ckpt or f'checkpoints/{virus}_*_lit_hotspot/best.pt'}`",
        f"Full JSON: `{out_dir / 'report.json'}`",
        "",
    ]
    md.write_text("\n".join(lines) + "\n")
    # Write sequences for plotting pipelines
    with (out_dir / f"{g}_seqs.fasta").open("w") as f:
        f.write(f">root\n{root_seq}\n")
        f.write(f">obs_tip\n{tip_seq}\n")
        f.write(f">gen_closest\n{gen_seq}\n")
    print(f"Wrote {md}")
    print(json.dumps(cand, indent=2))


if __name__ == "__main__":
    main()
