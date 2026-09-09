#!/usr/bin/env python3
"""Fixed-topology sequence paint (viral dest bake-off).

Observed Newick + BL are held fixed. Only the substitution process changes:
  neutral     — flat AA CTMC Q (same generator as NeutralBD, no free BD)
  codon_gy94  — GY94 collapsed to 20 AA (no learned residual)
  treesbm     — paper/new ckpt RateHeads residual on that R0
Optional --nt-fire-scale boosts fire at stall motifs (eval only).

Does not write checkpoints. Does not overwrite paper ckpt dirs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import torch
from Bio import SeqIO

from benchmarks.metrics.sequences import positional_recovery
from benchmarks.metrics.trees import children_map
from src.bridge.mutation_sample import mutate_sequence_independent
from src.bridge.nt_stall_fire import boost_fire_at_stalls
from src.dataset import fill_missing_node_seqs, parse_newick
from src.r0_backends import SubstitutionMatrixR0Backend, build_r0_backend
from src.tree_state import TreeState
from src.treeencoder.edges import build_edges
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.structural_features import compute_structural_features


def _load_tree(data_dir: Path, g: int) -> TreeState:
    root_id, node_ids, edges, bls = parse_newick(str(data_dir / f"group_{g:03d}_rooted.nwk"))
    seqs = {
        rec.id: str(rec.seq)
        for rec in SeqIO.parse(data_dir / f"group_{g:03d}_anc_aa.fasta", "fasta")
    }
    seqs = fill_missing_node_seqs(root_id, edges, seqs)
    seqs = {nid: seqs[nid] for nid in node_ids}
    has_ch = {p for p, _ in edges}
    leaves = [nid for nid in node_ids if nid not in has_ch]
    return TreeState(
        node_ids=node_ids,
        root_id=root_id,
        edges=edges,
        branch_lengths=bls,
        node_seqs=seqs,
        active_leaves=leaves,
    )


def _mean(xs):
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else float("nan")


def paint_r0(tree: TreeState, backend, mrs: float, nt_scale: float, max_len: int) -> dict[str, str]:
    cm = children_map(tree)
    seqs = {tree.root_id: tree.node_seqs[tree.root_id]}
    order, stack = [], [tree.root_id]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(cm.get(n, []))
    for p in order:
        for c in cm.get(p, []):
            bl = float(tree.branch_lengths.get((p, c), 0.0))
            parent = seqs[p]
            log_R = backend.log_mutation_rates([parent], max_len)[0]
            if nt_scale > 0:
                log_R = boost_fire_at_stalls(log_R, parent, nt_scale)
            seqs[c] = mutate_sequence_independent(
                log_R, parent, min(len(parent), max_len), bl, mrs,
            )
    return seqs


def paint_treesbm(
    tree: TreeState,
    *,
    node_enc,
    tree_enc,
    rate_heads,
    embedder,
    get_lm_logits,
    tokenizer,
    esm_model,
    aa_token_ids,
    device,
    col_entropy,
    mrs: float,
    nt_scale: float,
    max_len: int,
) -> dict[str, str]:
    from src.bridge.fitness_tilt import tilt_log_R0_by_fitness
    from src.bridge.losses import _build_seq_indices

    cm = children_map(tree)
    seqs = {tree.root_id: tree.node_seqs[tree.root_id]}
    order, stack = [], [tree.root_id]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(cm.get(n, []))
    for p in order:
        for c in cm.get(p, []):
            bl = float(tree.branch_lengths.get((p, c), 0.0))
            parent = seqs[p]
            mini = TreeState.root_only(parent)
            nids = mini.node_ids
            n2i = {nid: i for i, nid in enumerate(nids)}
            act = list(mini.active_leaves)
            aidx = [n2i[v] for v in act]
            struct = compute_structural_features(mini, n2i).to(device)
            lap = compute_laplacian_pe(mini, n2i, 8, device=device)
            eidx, _, eattr = build_edges(mini, n2i)
            eidx = eidx.to(device)
            bl_t = eattr.squeeze(-1).to(device)
            plm = embedder.embed_sequences([mini.node_seqs[nid] for nid in nids]).to(device)
            log_R0 = get_lm_logits(
                tokenizer, esm_model, aa_token_ids, [parent], max_len, device,
            )
            log_R0 = tilt_log_R0_by_fitness(
                log_R0,
                beta=float(getattr(rate_heads, "_fitness_beta", 0.0) or 0.0),
                score=getattr(rate_heads, "_fitness_score", "log_R0"),
                mode=getattr(rate_heads, "_fitness_tilt_mode", "site_local"),
            )
            aa_idx = None
            if getattr(rate_heads, "needs_aa_indices", False):
                aa_idx = _build_seq_indices([parent], max_len, device)
            with torch.no_grad():
                h = node_enc(plm, struct, lap)
                H, _ = tree_enc(
                    h, nids, {nid: 0.0 for nid in nids}, eidx, bl_t, t_scalar=0.0,
                )
                out = rate_heads(
                    H, aidx, log_R0,
                    site_entropy=col_entropy,
                    aa_indices=aa_idx,
                    log_pssm=getattr(rate_heads, "_train_log_pssm", None),
                )
            log_R = out["log_R_theta_mut"][0]
            if nt_scale > 0:
                log_R = boost_fire_at_stalls(log_R.cpu(), parent, nt_scale).to(log_R.device)
            seqs[c] = mutate_sequence_independent(
                log_R, parent, min(len(parent), max_len), bl, mrs,
            )
    return seqs


def score_leaves(tree: TreeState, gen: dict[str, str]) -> dict:
    root = tree.node_seqs[tree.root_id]
    recs = []
    for lid in tree.active_leaves:
        gt = tree.node_seqs.get(lid, "")
        g = gen.get(lid, "")
        if not gt or not g:
            continue
        recs.append(positional_recovery(root, gt, g))
    return {
        "n_leaves": len(recs),
        "mut_recovery": _mean([r["mut_recovery"] for r in recs]),
        "site_recall": _mean([r["site_recall"] for r in recs]),
        "aa_acc_given_hit": _mean([r["aa_acc_given_hit"] for r in recs]),
        "cons_retention": _mean([r["cons_retention"] for r in recs]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-seq-len", type=int, required=True)
    ap.add_argument("--checkpoint", default=None, help="TreeSBM ckpt (read-only)")
    ap.add_argument("--methods", default="neutral,codon_gy94,treesbm")
    ap.add_argument("--mutation-rate-scale", type=float, default=1.0)
    ap.add_argument("--nt-fire-scale", type=float, default=0.0)
    ap.add_argument("--max-trees", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data_dir = Path(args.data)
    groups = sorted(
        int(p.stem.split("_")[1])
        for p in data_dir.glob("group_*_rooted.nwk")
        if (data_dir / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists()
    )[: args.max_trees]

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    mrs = args.mutation_rate_scale
    want_tsbm = any(m.startswith("treesbm") for m in methods)
    tsbm = None
    if want_tsbm:
        if not args.checkpoint:
            raise SystemExit("treesbm method needs --checkpoint")
        from scripts.eval_single_tree import get_lm_logits, load_models
        from src.treeencoder.plm_embeddings import ESM2Embedder
        from transformers import AutoTokenizer, EsmForMaskedLM

        device = "cuda" if torch.cuda.is_available() else "cpu"
        node_enc, tree_enc, rate_heads, col_ent = load_models(
            args.checkpoint, device, args.max_seq_len
        )
        tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        esm = EsmForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device).eval()
        aa_ids = torch.tensor(
            [tok.convert_tokens_to_ids(a) for a in "ACDEFGHIKLMNPQRSTVWY"],
            dtype=torch.long,
        )
        tsbm = dict(
            node_enc=node_enc, tree_enc=tree_enc, rate_heads=rate_heads,
            embedder=ESM2Embedder(device=device), get_lm_logits=get_lm_logits,
            tokenizer=tok, esm_model=esm, aa_token_ids=aa_ids, device=device,
            col_entropy=col_ent, max_len=args.max_seq_len, mrs=mrs,
        )

    backends = {
        "neutral": SubstitutionMatrixR0Backend(model="neutral"),
        "codon_gy94": build_r0_backend("codon_gy94"),
    }

    per_method = {m: [] for m in methods}
    for g in groups:
        tree = _load_tree(data_dir, g)
        for m in methods:
            nt = args.nt_fire_scale if m.endswith("_ntfire") else 0.0
            base = m.replace("_ntfire", "")
            if base in backends:
                gen = paint_r0(tree, backends[base], mrs, nt, args.max_seq_len)
            elif base == "treesbm":
                gen = paint_treesbm(tree, nt_scale=nt, **tsbm)
            else:
                raise SystemExit(f"unknown method {m}")
            rec = score_leaves(tree, gen)
            rec["group"] = g
            per_method[m].append(rec)
        print(f"group {g:03d} done", flush=True)

    summary = {}
    for m, rows in per_method.items():
        summary[m] = {
            k: _mean([r[k] for r in rows])
            for k in ("mut_recovery", "site_recall", "aa_acc_given_hit", "cons_retention")
        }
        summary[m]["n_trees"] = len(rows)
        print(m, summary[m], flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "data": str(data_dir),
        "checkpoint": args.checkpoint,
        "mrs": mrs,
        "nt_fire_scale": args.nt_fire_scale,
        "summary": summary,
        "per_tree": per_method,
    }, indent=2) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
