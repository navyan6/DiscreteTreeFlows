#!/usr/bin/env python3
"""
Evaluate 2026 BDBV outbreak variant recovery on held-out test trees.

Targets (full L 1-based aa, FJ217161 / EBOV Makona for D759G):
  L:K1738N, L:Q1770R  (primary BDBV 2026)
  L:D759G             (EBOV sanity on pan / EBOV-only runs)

Also reports mut/site/aa|hit stratified by literature hotspot mask vs background.

Usage:
  python scripts/eval_bdbv_outbreak_variants.py \\
    --checkpoint checkpoints/bdbv_v1_mutrec/best.pt \\
    --data data/bdbv_temporal/test \\
    --max-seq-len 900 \\
    --out checkpoints/eval_bdbv_2026_variants_bdbv_v1_mutrec.json
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import torch
from Bio import AlignIO
from transformers import AutoTokenizer, EsmForMaskedLM

from scripts.bdbv_common import REF_EBOV_ACCESSION, REF_FJ217161, load_window_config
from scripts.eval_single_tree import (
    AA_VOCAB,
    generate_tree,
    get_leaves,
    load_models,
    positional_recovery,
    seq_identity,
)
from src.dataset import TreeDataset
from src.treeencoder.plm_embeddings import ESM2Embedder


def _mean(xs: list[float]) -> float:
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else float("nan")


def _ref_pos_to_msa_col(msa, ref_row_idx: int, aa_1based: int) -> int | None:
    seq = str(msa[ref_row_idx].seq)
    pos = 0
    for i, c in enumerate(seq):
        if c == "-":
            continue
        pos += 1
        if pos == aa_1based:
            return i
    return None


def _msa_col_to_window_idx(col: int, cfg: dict) -> int | None:
    aa_start = cfg["aa_start"]
    aa_end = cfg["aa_end"]
    if col < aa_start or col >= aa_end:
        return None
    return col - aa_start


def resolve_variant_indices(
    variants: list[dict],
    cfg: dict,
    msa_path: Path,
) -> list[dict]:
    msa = AlignIO.read(msa_path, "fasta")
    ref_idx = 0
    for i, rec in enumerate(msa):
        rid = rec.id
        if REF_FJ217161.split(".")[0] in rid or REF_EBOV_ACCESSION.split(".")[0] in rid:
            ref_idx = i
            break
    out = []
    for v in variants:
        aa = int(v["aa_ref"])
        ref_acc = v.get("ref", REF_FJ217161)
        if REF_EBOV_ACCESSION.split(".")[0] in ref_acc:
            for i, rec in enumerate(msa):
                if REF_EBOV_ACCESSION.split(".")[0] in rec.id or rec.id.startswith("KM034562"):
                    ref_idx = i
                    break
        col = _ref_pos_to_msa_col(msa, ref_idx, aa)
        wi = _msa_col_to_window_idx(col, cfg) if col is not None else None
        out.append({**v, "msa_col": col, "window_idx": wi})
    return out


def site_metrics(root_seq: str, gt_seq: str, gen_seq: str, wi: int) -> dict:
    if wi is None or wi >= min(len(root_seq), len(gt_seq), len(gen_seq)):
        return {
            "gt_mutates": False,
            "site_recall": float("nan"),
            "aa_acc_given_hit": float("nan"),
            "root_aa": "",
            "gt_aa": "",
            "gen_aa": "",
        }
    r, g, s = root_seq[wi], gt_seq[wi], gen_seq[wi]
    gt_mut = r != g
    site_hit = s != r
    aa_hit = s == g if gt_mut else float("nan")
    return {
        "gt_mutates": gt_mut,
        "site_recall": float(site_hit) if gt_mut else float("nan"),
        "aa_acc_given_hit": float(aa_hit) if gt_mut and site_hit else (
            float(aa_hit) if gt_mut else float("nan")
        ),
        "root_aa": r,
        "gt_aa": g,
        "gen_aa": s,
    }


def any_leaf_site_metrics(root_seq: str, gt_seq: str, gen_seqs: list[str], wi: int) -> dict:
    best = site_metrics(root_seq, gt_seq, root_seq, wi)
    if wi is None:
        return best
    gt_mut = root_seq[wi] != gt_seq[wi] if wi < len(root_seq) and wi < len(gt_seq) else False
    if not gt_mut:
        best["gt_mutates"] = False
        return best
    any_site = False
    any_aa = False
    for gs in gen_seqs:
        if wi >= len(gs):
            continue
        if gs[wi] != root_seq[wi]:
            any_site = True
        if gs[wi] == gt_seq[wi]:
            any_aa = True
    best["site_recall_any_descendant"] = float(any_site)
    best["aa_acc_any_descendant"] = float(any_aa)
    return best


def stratified_recovery(
    root_seq: str,
    gt_seq: str,
    gen_seq: str,
    lit_mask: torch.Tensor | None,
) -> dict:
    rec = positional_recovery(root_seq, gt_seq, gen_seq)
    if lit_mask is None or len(lit_mask) != len(root_seq):
        return {
            "lit_mut_recovery": float("nan"),
            "bg_mut_recovery": float("nan"),
            "lit_site_recall": float("nan"),
            "bg_site_recall": float("nan"),
        }
    L = len(root_seq)
    lit_mut = lit_hit = lit_aa = lit_n = 0
    bg_mut = bg_hit = bg_aa = bg_n = 0
    for i in range(L):
        if root_seq[i] == gt_seq[i]:
            continue
        if i < len(lit_mask) and lit_mask[i]:
            lit_n += 1
            if gen_seq[i] != root_seq[i]:
                lit_hit += 1
                if gen_seq[i] == gt_seq[i]:
                    lit_aa += 1
        else:
            bg_n += 1
            if gen_seq[i] != root_seq[i]:
                bg_hit += 1
                if gen_seq[i] == gt_seq[i]:
                    bg_aa += 1
    return {
        "lit_mut_recovery": (lit_hit / lit_n * (lit_aa / lit_hit)) if lit_n and lit_hit else float("nan"),
        "bg_mut_recovery": (bg_hit / bg_n * (bg_aa / bg_hit)) if bg_n and bg_hit else float("nan"),
        "lit_site_recall": lit_hit / lit_n if lit_n else float("nan"),
        "bg_site_recall": bg_hit / bg_n if bg_n else float("nan"),
        "lit_aa_acc_given_hit": lit_aa / lit_hit if lit_hit else float("nan"),
        "bg_aa_acc_given_hit": bg_aa / bg_hit if bg_hit else float("nan"),
        "lit_mut_sites": lit_n,
        "bg_mut_sites": bg_n,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--data", default="data/bdbv_temporal/test")
    ap.add_argument("--max-seq-len", type=int, default=900)
    ap.add_argument("--mutation-rate-scale", type=float, default=0.5)
    ap.add_argument("--branch-rate-scale", type=float, default=6.0)
    ap.add_argument("--n-steps", type=int, default=100)
    ap.add_argument("--max-leaves", type=int, default=200)
    ap.add_argument("--max-trees", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--variants-json",
        default="results/bdbv_l_mask/mut_hotspot_mask_eval_variants.json",
    )
    ap.add_argument(
        "--lit-mask",
        default="results/bdbv_l_mask/mut_hotspot_mask_lit.pt",
    )
    ap.add_argument(
        "--msa",
        default="results/bdbv_l_conservation/l_msa.fasta",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_window_config()
    if args.max_seq_len != cfg.get("max_seq_len"):
        args.max_seq_len = int(cfg["max_seq_len"])

    var_path = ROOT / args.variants_json
    variants_raw = json.loads(var_path.read_text())["variants"]
    resolved = resolve_variant_indices(variants_raw, cfg, ROOT / args.msa)

    lit_mask = None
    lit_path = ROOT / args.lit_mask
    if lit_path.is_file():
        blob = torch.load(lit_path, map_location="cpu", weights_only=False)
        lit_mask = blob.get("mut_hotspot_mask") or blob.get("mask")
        if isinstance(lit_mask, torch.Tensor):
            lit_mask = lit_mask.bool()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    node_enc, tree_enc, rate_heads, col_entropy = load_models(
        args.checkpoint, device, args.max_seq_len
    )
    embedder = ESM2Embedder(device=device)
    model_id = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    esm_model = EsmForMaskedLM.from_pretrained(model_id).to(device)
    esm_model.eval()
    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB], dtype=torch.long
    )

    ds = TreeDataset(str(ROOT / args.data), max_seq_len=args.max_seq_len)
    n = min(len(ds.groups), args.max_trees)
    rng = random.Random(args.seed)

    variant_stats: dict[str, dict] = {
        v["name"]: {"site_recall": [], "aa_acc_given_hit": [], "aa_any_desc": []}
        for v in resolved
    }
    lit_strat = {"lit_mut_recovery": [], "bg_mut_recovery": [], "lit_site_recall": [], "bg_site_recall": []}
    per_tree = []

    for i in range(n):
        random.seed(args.seed + i)
        torch.manual_seed(args.seed + i)
        batch = ds[i]
        root_id = batch["node_ids"][batch["root_index"]]
        root_seq = batch["seqs"][root_id]
        gt_leaves = [n for n in batch["node_ids"] if batch["is_leaf"][n]]
        if not gt_leaves:
            continue
        try:
            gen = generate_tree(
                root_seq,
                args.n_steps,
                args.max_seq_len,
                args.branch_rate_scale,
                args.max_leaves,
                args.mutation_rate_scale,
                node_enc,
                tree_enc,
                rate_heads,
                embedder,
                tokenizer,
                esm_model,
                aa_token_ids,
                device,
                col_entropy=col_entropy,
            )
        except Exception as e:
            print(f"[{i+1}/{n}] ERROR: {e}")
            continue

        gen_seqs = [gen.node_seqs[g] for g in get_leaves(gen)]
        gt_sample = rng.sample(gt_leaves, min(5, len(gt_leaves)))
        tree_row = {"tree": i, "variants": {}}

        for v in resolved:
            wi = v.get("window_idx")
            if wi is None:
                tree_row["variants"][v["name"]] = {"window_idx": None, "note": "outside window"}
                continue
            recalls, aa_hits, any_aa = [], [], []
            for gl in gt_sample:
                gt_seq = batch["seqs"][gl]
                best = None
                best_id = -1.0
                for gs in gen_seqs:
                    idv = seq_identity(gt_seq, gs)
                    if idv > best_id:
                        best_id, best = idv, gs
                if best is None:
                    continue
                m = site_metrics(root_seq, gt_seq, best, wi)
                ad = any_leaf_site_metrics(root_seq, gt_seq, gen_seqs, wi)
                if m["gt_mutates"]:
                    if m["site_recall"] == m["site_recall"]:
                        recalls.append(m["site_recall"])
                    if m["aa_acc_given_hit"] == m["aa_acc_given_hit"]:
                        aa_hits.append(m["aa_acc_given_hit"])
                    if ad.get("aa_acc_any_descendant") == ad.get("aa_acc_any_descendant"):
                        any_aa.append(ad["aa_acc_any_descendant"])
            tree_row["variants"][v["name"]] = {
                "window_idx": wi,
                "site_recall_mean": _mean(recalls),
                "aa_acc_given_hit_mean": _mean(aa_hits),
                "aa_any_descendant_mean": _mean(any_aa),
                "n_gt_mut": len(recalls),
            }
            if recalls:
                variant_stats[v["name"]]["site_recall"].extend(recalls)
            if aa_hits:
                variant_stats[v["name"]]["aa_acc_given_hit"].extend(aa_hits)
            if any_aa:
                variant_stats[v["name"]]["aa_any_desc"].extend(any_aa)

        for gl in gt_sample[:1]:
            gt_seq = batch["seqs"][gl]
            best = max(gen_seqs, key=lambda gs: seq_identity(gt_seq, gs), default=None)
            if best:
                st = stratified_recovery(root_seq, gt_seq, best, lit_mask)
                for k in lit_strat:
                    if st[k] == st[k]:
                        lit_strat[k].append(st[k])

        per_tree.append(tree_row)

    summary = {
        "checkpoint": args.checkpoint,
        "data": args.data,
        "max_seq_len": args.max_seq_len,
        "n_trees": n,
        "variants": {},
        "lit_stratification": {k: _mean(v) for k, v in lit_strat.items()},
        "resolved_sites": resolved,
    }
    for name, stats in variant_stats.items():
        summary["variants"][name] = {
            "site_recall": _mean(stats["site_recall"]),
            "aa_acc_given_hit": _mean(stats["aa_acc_given_hit"]),
            "aa_any_descendant": _mean(stats["aa_any_desc"]),
            "n_samples": len(stats["site_recall"]),
        }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {"summary": summary, "per_tree": per_tree}
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
