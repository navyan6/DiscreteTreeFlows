#!/usr/bin/env python3
"""
Screen COVID tree groups for VOC signature content (Track A).

Scores each group against benchmarks/voc_threat_panel.json:
  - root vs defining bundle (prospective if root lacks full bundle)
  - leaf signature / acquired recall
  - best matching VOC

Examples:
  # Brazil geo-test (held-out)
  python scripts/screen_voc_threat_trees.py \\
    --data-dir data/covid/test --split test \\
    --out results/voc_threat_panel/brazil_screen.json

  # Origin-country train groups (SA / India / UK / Peru / Colombia)
  python scripts/screen_voc_threat_trees.py \\
    --data-dir data/covid/train --split train \\
    --groups 300,301,302,303,190,191,192,323,289,34 \\
    --out results/voc_threat_panel/origin_screen.json

  # Primary panel only
  python scripts/screen_voc_threat_trees.py --data-dir data/covid/test \\
    --vocs Gamma,Beta,Delta,Omicron_BA1,Alpha --out-md results/voc_threat_panel/brazil_screen.md
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from voc_threat_lib import (  # noqa: E402
    date_range_from_csv,
    find_root_seq,
    group_paths,
    leaf_seqs,
    load_fasta_seqs,
    load_panel,
    score_group_against_voc,
    voc_by_id,
)


def discover_groups(data_dir: Path) -> list[int]:
    ids = []
    for p in sorted(data_dir.glob("group_*_anc_aa.fasta")):
        try:
            ids.append(int(p.name.split("_")[1]))
        except (IndexError, ValueError):
            continue
    return ids


def csv_for_split(paths: dict[str, Path], split: str) -> Path | None:
    key = {"test": "csv", "train": "csv_train", "val": "csv_val"}.get(split, "csv")
    p = paths.get(key)
    if p and p.exists():
        return p
    for k in ("csv", "csv_train", "csv_val"):
        if paths[k].exists():
            return paths[k]
    return None


def screen_one(
    data_dir: Path,
    gid: int,
    vocs: list[dict],
    split: str,
) -> dict:
    paths = group_paths(data_dir, gid)
    if not paths["anc"].exists():
        return {"gid": gid, "error": f"missing {paths['anc']}"}

    seqs = load_fasta_seqs(paths["anc"])
    root_id, root_seq = find_root_seq(seqs, paths["nwk"] if paths["nwk"].exists() else None)
    leaves = leaf_seqs(seqs)
    # drop root from leaves if present
    leaves = {k: v for k, v in leaves.items() if k != root_id}

    csv_path = csv_for_split(paths, split)
    dmin, dmax = date_range_from_csv(csv_path) if csv_path else (None, None)

    per_voc = {}
    for voc in vocs:
        per_voc[voc["id"]] = score_group_against_voc(root_seq, leaves, voc)

    # best VOC: prefer prospective + high leaf exact recall
    ranked = sorted(
        per_voc.values(),
        key=lambda s: (
            1 if s["prospective"] else 0,
            s["voc_exact_mut_recall"],
            s["voc_bundle_frac"],
            s["max_sig_on_leaf"],
        ),
        reverse=True,
    )
    best = ranked[0] if ranked else None

    return {
        "gid": gid,
        "split": split,
        "data_dir": str(data_dir),
        "root_id": root_id,
        "root_len": len(root_seq),
        "n_leaves": len(leaves),
        "date_min": dmin,
        "date_max": dmax,
        "best_voc": best["voc_id"] if best else None,
        "best_prospective": best["prospective"] if best else None,
        "best_exact_recall": best["voc_exact_mut_recall"] if best else None,
        "voc": per_voc,
    }


def write_md(rows: list[dict], path: Path, panel_ids: list[str]) -> None:
    lines = [
        "# VOC threat tree screen",
        "",
        "| gid | dates | n_leaves | best_voc | prospective | exact_recall | bundle_frac | root_bundle |",
        "|-----|-------|----------|----------|-------------|--------------|-------------|-------------|",
    ]
    for r in rows:
        if r.get("error"):
            lines.append(f"| {r['gid']} | — | — | ERROR | — | — | — | {r['error']} |")
            continue
        best_id = r["best_voc"]
        v = r["voc"].get(best_id, {}) if best_id else {}
        lines.append(
            f"| {r['gid']:03d} | {r.get('date_min')}…{r.get('date_max')} | {r['n_leaves']} | "
            f"{best_id} | {r.get('best_prospective')} | "
            f"{r.get('best_exact_recall', 0):.2f} | {v.get('voc_bundle_frac', 0):.3f} | "
            f"{v.get('root_bundle_hits', [])} |"
        )
    lines.extend(["", "## Per-VOC highlights (prospective + recall≥0.3)", ""])
    for vid in panel_ids:
        lines.append(f"### {vid}")
        hits = []
        for r in rows:
            if r.get("error"):
                continue
            s = r["voc"].get(vid)
            if not s:
                continue
            if s["prospective"] and s["voc_exact_mut_recall"] >= 0.3:
                hits.append(
                    f"- group {r['gid']:03d} ({r.get('date_min')}…{r.get('date_max')}): "
                    f"exact_recall={s['voc_exact_mut_recall']:.2f} "
                    f"bundle_frac={s['voc_bundle_frac']:.3f} "
                    f"acquired_recall={s['voc_acquired_mut_recall']:.2f}"
                )
        lines.extend(hits or ["- (none)"])
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--split", choices=["train", "val", "test"], default="test")
    ap.add_argument("--groups", type=str, default="", help="Comma list; default=all with anc_aa")
    ap.add_argument("--vocs", type=str, default="", help="Comma VOC ids; default=primary_panel")
    ap.add_argument("--out", type=Path, default=Path("results/voc_threat_panel/screen.json"))
    ap.add_argument("--out-md", type=Path, default=None)
    ap.add_argument(
        "--candidates-out",
        type=Path,
        default=None,
        help="Write ranked case-study candidates JSON for Betty submit",
    )
    args = ap.parse_args()

    panel = load_panel(args.panel)
    if args.vocs:
        voc_ids = [x.strip() for x in args.vocs.split(",") if x.strip()]
    else:
        voc_ids = list(panel.get("primary_panel") or [v["id"] for v in panel["vocs"]])
    vocs = [voc_by_id(panel, v) for v in voc_ids]

    if args.groups:
        gids = [int(x) for x in args.groups.split(",") if x.strip()]
    else:
        gids = discover_groups(args.data_dir)

    rows = []
    for gid in gids:
        print(f"screening group {gid:03d} …", flush=True)
        rows.append(screen_one(args.data_dir, gid, vocs, args.split))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "data_dir": str(args.data_dir),
        "split": args.split,
        "voc_ids": voc_ids,
        "n_groups": len(rows),
        "groups": rows,
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.out}")

    md_path = args.out_md or args.out.with_suffix(".md")
    write_md(rows, md_path, voc_ids)
    print(f"wrote {md_path}")

    # Ranked candidates: prospective, high recall, prefer held-out / primary
    cands = []
    for r in rows:
        if r.get("error"):
            continue
        for vid, s in r["voc"].items():
            if not s["prospective"]:
                continue
            if s["voc_exact_mut_recall"] < 0.25 and s["voc_bundle_frac"] < 0.05:
                continue
            cands.append(
                {
                    "gid": r["gid"],
                    "split": r["split"],
                    "data_dir": r["data_dir"],
                    "voc_id": vid,
                    "date_min": r.get("date_min"),
                    "date_max": r.get("date_max"),
                    "n_leaves": r["n_leaves"],
                    "exact_recall": s["voc_exact_mut_recall"],
                    "acquired_recall": s["voc_acquired_mut_recall"],
                    "bundle_frac": s["voc_bundle_frac"],
                    "root_bundle_hits": s["root_bundle_hits"],
                    "score": (
                        s["voc_exact_mut_recall"]
                        + 0.5 * s["voc_acquired_mut_recall"]
                        + 0.25 * s["voc_bundle_frac"]
                        + (0.1 if r["split"] == "test" else 0.0)
                    ),
                }
            )
    cands.sort(key=lambda x: x["score"], reverse=True)
    cout = args.candidates_out or args.out.parent / f"candidates_{args.split}.json"
    cout.write_text(json.dumps({"candidates": cands[:50]}, indent=2) + "\n")
    print(f"wrote {cout} ({min(50, len(cands))} candidates)")

    # also CSV summary
    csv_path = args.out.with_suffix(".csv")
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "gid", "split", "date_min", "date_max", "n_leaves",
                "voc_id", "prospective", "exact_recall", "acquired_recall",
                "bundle_frac", "root_n_sig", "root_bundle",
            ],
        )
        w.writeheader()
        for r in rows:
            if r.get("error"):
                continue
            for vid, s in r["voc"].items():
                w.writerow(
                    {
                        "gid": r["gid"],
                        "split": r["split"],
                        "date_min": r.get("date_min"),
                        "date_max": r.get("date_max"),
                        "n_leaves": r["n_leaves"],
                        "voc_id": vid,
                        "prospective": s["prospective"],
                        "exact_recall": f"{s['voc_exact_mut_recall']:.4f}",
                        "acquired_recall": f"{s['voc_acquired_mut_recall']:.4f}",
                        "bundle_frac": f"{s['voc_bundle_frac']:.4f}",
                        "root_n_sig": s["root_n_sig"],
                        "root_bundle": ";".join(s["root_bundle_hits"]),
                    }
                )
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
