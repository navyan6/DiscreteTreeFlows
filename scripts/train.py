#!/usr/bin/env python3
"""
Bridge matching training loop for TreeSBM (Algorithm 1).

Each step:
  1. Load T1 from dataset, compute PLM embeddings once.
  2. Sample t ~ U(0, t_max_clip).
  3. Construct T_t via SampleBridgeState (Algorithm 2).
  4. Rebuild TreeState / structural features / Laplacian PE for T_t.
  5. Run NodeEncoder and TreeEncoder (with time conditioning) to get H_t.
  6. Run RateHeads for active leaves of T_t.
  7. Compute bridge matching loss vs. T1 targets.

Usage:
    python scripts/train.py --data data/train --epochs 100 --lr 1e-4
"""

import argparse
import math
import random
from pathlib import Path

import torch
import torch.nn as nn
import json
from collections import defaultdict
from torch.utils.data import DataLoader, Subset

ROOT = Path(__file__).parent.parent
import sys; sys.path.insert(0, str(ROOT))

from src.dataset import TreeDataset
from src.tree_state import TreeState
from src.treeencoder.plm_embeddings import ESM2Embedder
from src.treeencoder.node_encoder import NodeEncoder
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.edges import build_edges
from src.networks import TreeEncoder, RateHeads
from src.bridge.sample_bridge_state import sample_bridge_state
from src.bridge.losses import bridge_losses


def compute_site_entropy_from_log_probs(log_R0_mut: torch.Tensor) -> torch.Tensor:
    """Per-position Shannon entropy in nats from [n_active, L, 20] log-probs/logits."""
    log_probs = torch.log_softmax(log_R0_mut, dim=-1)
    probs = log_probs.exp()
    return -(probs * log_probs).sum(dim=-1)


AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"


def compute_empirical_column_entropy(dataset, max_seq_len: int) -> torch.Tensor:
    """
    Per-column Shannon entropy over ALL training-tree sequences, normalized to
    [0,1] by log(20). This is the "which sites historically vary" hotspot signal
    (antigenic sites -- e.g. spike RBD, D614G -- score high), computed ONCE from
    the training alignment and then used identically at train and generation so
    the two never diverge. Columns with no observed residue get 0.

    Uses only `dataset` (the TRAIN split) -- never val/test -- so no leakage.
    """
    import numpy as np

    tbl = np.full(256, 255, dtype=np.int16)
    for i, aa in enumerate(AA_VOCAB):
        tbl[ord(aa)] = i  # non-AA bytes (X, -, *, ...) stay 255 -> excluded

    counts = np.zeros((max_seq_len, 20), dtype=np.int64)
    for i in range(len(dataset)):
        for seq in dataset[i]["seqs"].values():
            if not seq:
                continue
            # latin-1 is a 1:1 byte map (never drops/shifts a column); protein
            # chars are all <128 so this equals ascii for them.
            arr = tbl[np.frombuffer(seq[:max_seq_len].encode("latin-1"), dtype=np.uint8)]
            valid = arr < 20
            if valid.any():
                pos = np.nonzero(valid)[0]
                np.add.at(counts, (pos, arr[valid]), 1)

    total = counts.sum(axis=1, keepdims=True)
    probs = counts / np.maximum(total, 1)
    safe_log = np.log(np.maximum(probs, 1e-12))   # avoid log(0); masked by probs=0 anyway
    ent = -(probs * safe_log).sum(axis=1)          # nats
    ent = np.clip(ent / math.log(20.0), 0.0, 1.0)
    ent[total[:, 0] == 0] = 0.0
    return torch.tensor(ent, dtype=torch.float32)


def forward_bridge_step(
    batch: dict,
    t: float,
    node_enc: NodeEncoder,
    tree_enc: TreeEncoder,
    rate_heads: RateHeads,
    device: str,
    lap_dim: int = 8,
    max_seq_len: int = 566,
    lambda_top: float = 0.1,
    lambda_br: float = 0.1,
    lambda_stop: float = 0.1,
    lambda_pll: float = 0.01,
    lambda_mut: float = 5.0,
    bridge_c: float = 1.0,
    use_site_entropy: bool = False,
    use_entropy_loss_weighting: bool = False,
    entropy_weight_alpha: float = 1.0,
    entropy_weight_floor: float = 1.0,
    entropy_is_normalized: bool = False,
    col_entropy: torch.Tensor | None = None,
    embedder: ESM2Embedder | None = None,
) -> tuple[dict | None, int]:
    """
    One forward pass of Algorithm 1.
    or
    Returns (losses_dict, n_active_leaves).
    """
    node_ids= batch["node_ids"]              # list[str], N
    node_times_t    = batch["node_times"]            # [N] tensor
    seqs            = batch["seqs"]                  # dict[str, str]
    edges           = batch["edges"]                 # list[(parent, child)]
    branch_lengths  = batch["branch_lengths"]        # dict[(str,str), float]
    root_id         = node_ids[batch["root_index"]]

    node_times_dict = {nid: node_times_t[i].item() for i, nid in enumerate(node_ids)}

    if batch.get("plm_embeddings") is not None:
        plm_T1 = batch["plm_embeddings"].to(device)              # [N, 320] cached
    elif embedder is not None:
        sequences = [seqs[nid] for nid in node_ids]
        with torch.no_grad():
            plm_T1 = embedder.embed_sequences(sequences).to(device)
    else:
        raise RuntimeError("No PLM embeddings: run scripts/precompute_plm.py first")
    plm_map = {nid: i for i, nid in enumerate(node_ids)}

    T_t = sample_bridge_state(
        t=t,
        node_ids=node_ids,
        node_times_dict=node_times_dict,
        edges=edges,
        branch_lengths=branch_lengths,
        seqs=seqs,
        root_id=root_id,
    )

    node_ids_t      = T_t["node_ids_t"]
    edges_t         = T_t["edges_t"]
    branch_lengths_t = T_t["branch_lengths_t"]
    seqs_t          = T_t["seqs_t"]
    active_leaves_t = T_t["active_leaves_t"]

    if len(node_ids_t) == 0 or len(active_leaves_t) == 0:
        return None, 0

    tree_t = TreeState(
        node_ids=node_ids_t,
        root_id=root_id,
        edges=edges_t,
        branch_lengths=branch_lengths_t,
        node_seqs=seqs_t,
        active_leaves=active_leaves_t,
    )
    node_to_idx_t = {nid: i for i, nid in enumerate(node_ids_t)}

    #compute necessary features

    struct_t = compute_structural_features(tree_t, node_to_idx_t).to(device)
    lap_t = compute_laplacian_pe(tree_t, node_to_idx_t, lap_dim, device=device)

    plm_t = torch.stack([plm_T1[plm_map[nid]] for nid in node_ids_t]).to(device)  # [N_t, 320]

    
    h_t = node_enc(plm_t.to(device), struct_t.to(device), lap_t.to(device))  # [N_t, 128]

    # edge tensors
    edge_index_t, _, edge_attr_t = build_edges(tree_t, node_to_idx_t)
    edge_index_t = edge_index_t.to(device)
    branch_lens_t = edge_attr_t.squeeze(-1).to(device) 

    # treeencoder + time
    H_t, _ = tree_enc(
        h_t, node_ids_t, node_times_dict, edge_index_t, branch_lens_t, t_scalar=t
    ) 

    # R0 must be computed before rate_heads so the mutation head can condition on it
    active_idx_t = [node_to_idx_t[nid] for nid in active_leaves_t]
    if batch.get("log_ref_mut_rates") is not None:
        log_R0_mut = torch.stack([
            batch["log_ref_mut_rates"][plm_map[nid]] for nid in active_leaves_t
        ]).to(device)                                      # [n_active, 566, 20]
    else:
        log_R0_mut = torch.zeros(len(active_leaves_t), max_seq_len, 20, device=device)

    site_entropy = None
    ent_is_norm = entropy_is_normalized
    if use_site_entropy or use_entropy_loss_weighting:
        if col_entropy is not None:
            # Empirical column entropy from the training alignment, already
            # normalized to [0,1]. A [L] vector; RateHeads / bridge_losses
            # broadcast it over active leaves. Used identically here and at
            # generation (loaded from the checkpoint) so the two never diverge.
            site_entropy = col_entropy.to(device=log_R0_mut.device, dtype=log_R0_mut.dtype)
            ent_is_norm = True
        else:
            # Fallback: ESM self-entropy H(softmax(log_R0)) per position.
            site_entropy = compute_site_entropy_from_log_probs(log_R0_mut)
            if entropy_is_normalized:
                site_entropy = (site_entropy / math.log(20.0)).clamp(min=0.0, max=1.0)

    out = rate_heads(H_t, active_idx_t, log_R0_mut, site_entropy=site_entropy)
    # out["log_R_theta_mut"] = log_R0 + c_θ, computed inside RateHeads

    losses = bridge_losses(
        log_R_theta_mut=out["log_R_theta_mut"],
        log_R_theta_branch=out["branching_rate"],
        branch_length_pred=out["branch_length"],
        stop_prob=out["stop_prob"],
        log_R0_mut=log_R0_mut,
        seqs_t=[seqs_t[nid] for nid in active_leaves_t],
        active_leaves=active_leaves_t,
        T1_mut_targets=T_t["T1_mut_targets"],
        T1_child_counts=T_t["T1_child_counts"],
        T1_child_bls=T_t["T1_child_bls"],
        t=t,
        max_seq_len=max_seq_len,
        lambda_top=lambda_top,
        lambda_br=lambda_br,
        lambda_stop=lambda_stop,
        lambda_pll=lambda_pll,
        lambda_mut=lambda_mut,
        bridge_c=bridge_c,
        device=device,
        site_entropy=site_entropy,
        use_entropy_loss_weighting=use_entropy_loss_weighting,
        entropy_weight_alpha=entropy_weight_alpha,
        entropy_weight_floor=entropy_weight_floor,
        entropy_is_normalized=ent_is_norm,
    )

    return losses, len(active_leaves_t)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",        default="data/train")
    parser.add_argument("--val-data",    default=None,
                        help="Pre-split val dir, e.g. data/h3n2/val (temporal) or "
                             "data/covid/val (geographic). If set with --test-data, "
                             "bypasses the random subtype split.")
    parser.add_argument("--test-data",   default=None,
                        help="Pre-split test dir, e.g. data/h3n2/test or data/covid/test.")
    parser.add_argument("--epochs",      type=int,   default=100)
    parser.add_argument("--lr",          type=float, default=1e-4)
    parser.add_argument("--val-frac",    type=float, default=0.1)
    parser.add_argument("--test-frac",   type=float, default=0.1)
    parser.add_argument("--ckpt-dir",    default="checkpoints")
    parser.add_argument("--seed",        type=int,   default=42)
    parser.add_argument("--t-max",       type=float, default=0.95,
                        help="Max t to sample (avoids 1/(1-t) blow-up near t=1)")
    parser.add_argument("--lambda-top",  type=float, default=0.1)
    parser.add_argument("--lambda-br",   type=float, default=0.1)
    parser.add_argument("--lambda-stop", type=float, default=0.1)
    parser.add_argument("--lambda-pll",  type=float, default=0.01)
    parser.add_argument("--lambda-mut",  type=float, default=5.0,
                        help="Upweight loss at mutating positions (T_t_aa != T1_aa) within L_rate")
    parser.add_argument("--bridge-c",    type=float, default=1.0,
                        help="Reference resampling rate c in the conditional bridge target "
                             "(kappa = exp(-c(1-t))); larger = sharper terminal pull earlier")
    parser.add_argument("--per-site-pos-emb", action="store_true",
                        help="Add a learned positional embedding to the mutation head so "
                             "c_theta can act per-site (attacks the recovery ceiling). "
                             "Changes the architecture -> needs a fresh checkpoint.")
    parser.add_argument("--use-site-entropy", action="store_true",
                        help="Inject per-position Shannon entropy into the mutation head.")
    parser.add_argument("--use-entropy-loss-weighting", action="store_true",
                        help="Weight mutating-position bridge loss by per-position entropy.")
    parser.add_argument("--entropy-weight-alpha", type=float, default=1.0,
                        help="Slope for entropy-based mutation-loss weights.")
    parser.add_argument("--entropy-weight-floor", type=float, default=1.0,
                        help="Positive baseline added to entropy-based mutation-loss weights.")
    parser.add_argument("--entropy-is-normalized", action="store_true",
                        help="Treat provided entropy values as already normalized to [0,1].")
    parser.add_argument("--entropy-source", choices=["esm_self", "empirical"], default="esm_self",
                        help="Which per-position entropy feeds --use-site-entropy / "
                             "--use-entropy-loss-weighting. 'esm_self': H(softmax(log_R0)) "
                             "per position (free, self-consistent, = ESM uncertainty). "
                             "'empirical': Shannon entropy of each TRAIN-alignment column "
                             "(antigenic hotspots; computed once, saved in the checkpoint, "
                             "and reused at generation so train/inference never diverge).")
    parser.add_argument("--max-seq-len", type=int,   default=566)
    parser.add_argument("--patience",    type=int,   default=30,
                        help="Early stopping: stop if val loss doesn't improve for this many epochs")
    parser.add_argument("--n-t-samples", type=int,  default=4,
                        help="Number of t values sampled per tree per epoch")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from {ckpt-dir}/best.pt (weights, optimizer, scheduler, "
                             "epoch, best_val, patience_counter) instead of starting fresh. "
                             "Needed for long runs that may hit the SLURM time limit.")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    # ── data: pre-split dirs (temporal or geographic) OR random subtype-binned split
    if args.val_data and args.test_data:
        print("Pre-split data: loading train/val/test from separate dirs")
        train_ds = TreeDataset(args.data,      max_seq_len=args.max_seq_len)
        val_ds   = TreeDataset(args.val_data,  max_seq_len=args.max_seq_len)
        test_ds  = TreeDataset(args.test_data, max_seq_len=args.max_seq_len)
        dataset  = train_ds  # used for the PLM-cache probe / export below
        Path(args.ckpt_dir).mkdir(exist_ok=True)
        print(f"Total — Train: {len(train_ds)}  Val: {len(val_ds)}  Test: {len(test_ds)}")
    else:
        dataset = TreeDataset(args.data, max_seq_len=args.max_seq_len)

        def subtype_of(group: int) -> str:
            if   1   <= group <= 48:  return "h3n2"
            elif 49  <= group <= 55:  return "swine"
            elif group == 56:         return "avian"
            elif 57  <= group <= 106: return "h1n1_ha"
            elif 107 <= group <= 156: return "h1n1_na"
            elif 157 <= group <= 206: return "h1n1_2015"
            elif 207 <= group <= 238: return "fluB_yam"
            elif 239 <= group <= 284: return "fluB_vic"
            return "unknown"

        bins: dict[str, list[int]] = defaultdict(list)
        for i in range(len(dataset)):
            bins[subtype_of(dataset.groups[i])].append(i)

        rng = random.Random(args.seed)
        train_idx, val_idx, test_idx = [], [], []
        for subtype in sorted(bins):
            indices = bins[subtype]
            rng.shuffle(indices)
            n = len(indices)
            if n < 3:
                train_idx += indices
                print(f"  {subtype}: {n} tree(s) → all train (too small to split)")
                continue
            n_t  = max(1, int(n * args.test_frac))
            n_v  = max(1, int(n * args.val_frac))
            n_tr = n - n_t - n_v
            train_idx += indices[:n_tr]
            val_idx   += indices[n_tr:n_tr + n_v]
            test_idx  += indices[n_tr + n_v:]
            print(f"  {subtype}: {n} trees → {n_tr} train / {n_v} val / {n_t} test")

        train_ds = Subset(dataset, train_idx)
        val_ds   = Subset(dataset, val_idx)
        test_ds  = Subset(dataset, test_idx)
        print(f"Total — Train: {len(train_idx)}  Val: {len(val_idx)}  Test: {len(test_idx)}")

        # save split indices so the held-out test set is always recoverable
        split_path = Path(args.ckpt_dir) / "split_indices.json"
        split_path.parent.mkdir(exist_ok=True)
        with open(split_path, "w") as f:
            json.dump({"train": train_ds.indices, "val": val_ds.indices, "test": test_ds.indices}, f)

    # batch_size=1 (trees vary in node count — no collation)
    train_loader = DataLoader(train_ds, batch_size=1, shuffle=True,  collate_fn=lambda x: x[0])
    val_loader   = DataLoader(val_ds,   batch_size=1, shuffle=False, collate_fn=lambda x: x[0])
    test_loader  = DataLoader(test_ds,  batch_size=1, shuffle=False, collate_fn=lambda x: x[0])

    # ── models 
    # Only instantiate ESM2 if PLM caches are missing (fallback)
    first_batch = dataset[0]
    if first_batch.get("plm_embeddings") is None:
        print("WARNING: No PLM cache found — run scripts/precompute_plm.py for faster training")
        embedder = ESM2Embedder(device=device)
    else:
        embedder = None
        print("PLM embeddings cached — skipping ESM2 at training time")

    node_enc   = NodeEncoder(d_plm=320, d_struct=3, d_laplacian=8, d_node=128).to(device)
    tree_enc   = TreeEncoder(d_model=128, n_layers=4, n_heads=8, dropout=0.1).to(device)
    rate_heads = RateHeads(d_model=128, max_seq_len=args.max_seq_len,
                           use_pos_emb=args.per_site_pos_emb,
                           use_site_entropy=args.use_site_entropy).to(device)

    # Empirical column entropy (computed once from the TRAIN split) if requested;
    # otherwise None -> forward_bridge_step falls back to ESM self-entropy.
    col_entropy = None
    if args.entropy_source == "empirical" and (args.use_site_entropy or args.use_entropy_loss_weighting):
        print("Computing empirical column entropy from the training alignment...")
        col_entropy = compute_empirical_column_entropy(train_ds, args.max_seq_len).to(device)
        print(f"  col_entropy: [{col_entropy.numel()}]  mean={col_entropy.mean():.3f}  "
              f"max={col_entropy.max():.3f}  nonzero={(col_entropy > 0).sum().item()}")

    params = (
        list(node_enc.parameters()) +
        list(tree_enc.parameters()) +
        list(rate_heads.parameters())
    )
    optimizer = torch.optim.AdamW(params, lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs, eta_min=1e-6
    )

    ckpt_dir = Path(args.ckpt_dir)
    ckpt_dir.mkdir(exist_ok=True)
    best_val = float("inf")
    patience_counter = 0
    start_epoch = 1

    resume_path = ckpt_dir / "best.pt"
    if args.resume and resume_path.exists():
        ckpt = torch.load(resume_path, map_location=device, weights_only=False)
        node_enc.load_state_dict(ckpt["node_enc"])
        tree_enc.load_state_dict(ckpt["tree_enc"])
        rate_heads.load_state_dict(ckpt["rate_heads"])
        optimizer.load_state_dict(ckpt["optimizer"])
        if "scheduler" in ckpt:
            scheduler.load_state_dict(ckpt["scheduler"])
        best_val = ckpt["val_loss"]
        patience_counter = ckpt.get("patience_counter", 0)
        start_epoch = ckpt["epoch"] + 1
        print(f"Resumed from {resume_path}: epoch {ckpt['epoch']}, "
              f"best_val={best_val:.4f}, patience_counter={patience_counter}")

    # ── training loop
    for epoch in range(start_epoch, args.epochs + 1):
        node_enc.train(); tree_enc.train(); rate_heads.train()
        train_loss = 0.0
        n_steps = 0
        loss_breakdown = {"L_rate": 0.0, "L_mut": 0.0, "L_cons": 0.0, "L_top": 0.0, "L_br": 0.0, "L_stop": 0.0, "L_pll": 0.0,
                          "mean_mut_entropy": 0.0, "mean_mut_weight": 0.0, "max_mut_weight": 0.0,
                          "L_br_pred_std": 0.0, "L_br_target_std": 0.0}

        for batch in train_loader:
            for _ in range(args.n_t_samples):
                t = random.uniform(0.0, args.t_max)
                optimizer.zero_grad()

                losses, n_active = forward_bridge_step(
                    batch, t, node_enc, tree_enc, rate_heads, device,
                    max_seq_len=args.max_seq_len,
                    lambda_top=args.lambda_top,
                    lambda_br=args.lambda_br,
                    lambda_stop=args.lambda_stop,
                    lambda_pll=args.lambda_pll,
                    lambda_mut=args.lambda_mut,
                    bridge_c=args.bridge_c,
                    use_site_entropy=args.use_site_entropy,
                    use_entropy_loss_weighting=args.use_entropy_loss_weighting,
                    entropy_weight_alpha=args.entropy_weight_alpha,
                    entropy_weight_floor=args.entropy_weight_floor,
                    entropy_is_normalized=args.entropy_is_normalized,
                    col_entropy=col_entropy,
                    embedder=embedder,
                )
                if losses is None or n_active == 0:
                    continue
                if torch.isnan(losses["total"]):
                    continue

                losses["total"].backward()
                for p in params:
                    if p.grad is not None:
                        p.grad.nan_to_num_(0.0, 0.0, 0.0)
                nn.utils.clip_grad_norm_(params, max_norm=1.0)
                optimizer.step()

                train_loss += losses["total"].item()
                for k in loss_breakdown:
                    loss_breakdown[k] += losses[k].item()
                n_steps += 1

        if n_steps > 0:
            train_loss /= n_steps
            for k in loss_breakdown:
                loss_breakdown[k] /= n_steps

        # ── validation (fixed t=0.5 for reproducibility) 
        node_enc.eval(); tree_enc.eval(); rate_heads.eval()
        val_loss = 0.0
        n_val_steps = 0
        with torch.no_grad():
            for batch in val_loader:
                losses, n_active = forward_bridge_step(
                    batch, t=0.5,
                    node_enc=node_enc, tree_enc=tree_enc, rate_heads=rate_heads,
                    device=device, max_seq_len=args.max_seq_len,
                    lambda_top=args.lambda_top, lambda_br=args.lambda_br,
                    lambda_stop=args.lambda_stop, lambda_pll=args.lambda_pll,
                    lambda_mut=args.lambda_mut, bridge_c=args.bridge_c,
                    use_site_entropy=args.use_site_entropy,
                    use_entropy_loss_weighting=args.use_entropy_loss_weighting,
                    entropy_weight_alpha=args.entropy_weight_alpha,
                    entropy_weight_floor=args.entropy_weight_floor,
                    entropy_is_normalized=args.entropy_is_normalized,
                    col_entropy=col_entropy,
                    embedder=embedder,
                )
                if losses is None or n_active == 0:
                    continue
                val_loss += losses["total"].item()
                n_val_steps += 1

        if n_val_steps > 0:
            val_loss /= n_val_steps

        scheduler.step()
        lr = scheduler.get_last_lr()[0]
        print(
            f"Epoch {epoch:03d}  "
            f"train={train_loss:.4f} "
            f"(rate={loss_breakdown['L_rate']:.3f} "
            f"mut={loss_breakdown['L_mut']:.3f} "
            f"cons={loss_breakdown['L_cons']:.3f} "
            f"top={loss_breakdown['L_top']:.3f} "
            f"br={loss_breakdown['L_br']:.6f} "
            f"stop={loss_breakdown['L_stop']:.3f} "
            f"pll={loss_breakdown['L_pll']:.3f})  "
            f"val={val_loss:.4f}  lr={lr:.2e}  "
            f"[mut_entropy={loss_breakdown['mean_mut_entropy']:.3f} "
            f"mut_w={loss_breakdown['mean_mut_weight']:.3f} "
            f"mut_w_max={loss_breakdown['max_mut_weight']:.3f}]  "
            f"[br_pred_std={loss_breakdown['L_br_pred_std']:.6f} "
            f"br_target_std={loss_breakdown['L_br_target_std']:.6f}]"
        )

        if val_loss < best_val:
            best_val = val_loss
            patience_counter = 0
            torch.save({
                "epoch": epoch,
                "node_enc": node_enc.state_dict(),
                "tree_enc": tree_enc.state_dict(),
                "rate_heads": rate_heads.state_dict(),
                "optimizer": optimizer.state_dict(),
                "val_loss": val_loss,
                "config": {
                    "use_pos_emb": args.per_site_pos_emb,
                    "use_site_entropy": args.use_site_entropy,
                    "use_entropy_loss_weighting": args.use_entropy_loss_weighting,
                    "entropy_weight_alpha": args.entropy_weight_alpha,
                    "entropy_weight_floor": args.entropy_weight_floor,
                    "entropy_is_normalized": args.entropy_is_normalized,
                    "entropy_source": args.entropy_source,
                },
                # empirical column-entropy vector [L] (None for esm_self), so
                # generation reuses the exact same signal training saw.
                "col_entropy": col_entropy.cpu() if col_entropy is not None else None,
            }, ckpt_dir / "best.pt")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"\nEarly stopping at epoch {epoch} (no improvement for {args.patience} epochs)")
                break

    print(f"\nBest val loss: {best_val:.4f}  -> {ckpt_dir}/best.pt")

    # ── test evaluation on best checkpoint ──
    ckpt = torch.load(ckpt_dir / "best.pt", map_location=device, weights_only=False)
    node_enc.load_state_dict(ckpt["node_enc"])
    tree_enc.load_state_dict(ckpt["tree_enc"])
    rate_heads.load_state_dict(ckpt["rate_heads"])
    node_enc.eval(); tree_enc.eval(); rate_heads.eval()
    test_loss = 0.0
    n_test_steps = 0
    with torch.no_grad():
        for batch in test_loader:
            losses, n_active = forward_bridge_step(
                batch, t=0.5,
                node_enc=node_enc, tree_enc=tree_enc, rate_heads=rate_heads,
                device=device, max_seq_len=args.max_seq_len,
                lambda_top=args.lambda_top, lambda_br=args.lambda_br,
                lambda_stop=args.lambda_stop, lambda_pll=args.lambda_pll,
                lambda_mut=args.lambda_mut, bridge_c=args.bridge_c,
                use_site_entropy=args.use_site_entropy,
                use_entropy_loss_weighting=args.use_entropy_loss_weighting,
                entropy_weight_alpha=args.entropy_weight_alpha,
                entropy_weight_floor=args.entropy_weight_floor,
                entropy_is_normalized=args.entropy_is_normalized,
                col_entropy=col_entropy,
                embedder=embedder,
            )
            if losses is None or n_active == 0:
                continue
            test_loss += losses["total"].item()
            n_test_steps += 1
    if n_test_steps > 0:
        test_loss /= n_test_steps
    print(f"Test  loss: {test_loss:.4f}  ({n_test_steps} trees)")

    # export embeddings
    print("\nExporting embeddings with trained weights")
    export_embeddings(dataset, node_enc, tree_enc, device, Path(args.data))
    print("Done.")


def export_embeddings(
    dataset: TreeDataset,
    node_enc: NodeEncoder,
    tree_enc: TreeEncoder,
    device: str,
    data_dir: Path,
):

    node_enc.eval()
    tree_enc.eval()

    with torch.no_grad():
        for i in range(len(dataset)):
            batch = dataset[i]
            g = batch["group"]
            out_path = data_dir / f"group_{g:03d}_trained_emb.pt"

            node_ids      = batch["node_ids"]
            node_times_t  = batch["node_times"]
            edges         = batch["edges"]
            branch_lengths = batch["branch_lengths"]
            plm_T1        = batch["plm_embeddings"].to(device)  # [N, 320]

            node_times_dict = {nid: node_times_t[i].item() for i, nid in enumerate(node_ids)}
            root_id = node_ids[batch["root_index"]]

            has_children = {p for p, c in edges}
            tree_T1 = TreeState(
                node_ids=node_ids, root_id=root_id,
                edges=edges, branch_lengths=branch_lengths,
                node_seqs=batch["seqs"],
                active_leaves=[nid for nid in node_ids if nid not in has_children],
            )
            n2i = {nid: j for j, nid in enumerate(node_ids)}

            struct = compute_structural_features(tree_T1, n2i).to(device)
            lap    = compute_laplacian_pe(tree_T1, n2i, 8, device=device)

            # NodeEncoder: fuse PLM + structural + Laplacian → [N, 128]
            node_emb = node_enc(plm_T1, struct, lap)

            # TreeEncoder at t=1.0 (full tree, no bridge sampling)
            edge_index, _, edge_attr = build_edges(tree_T1, n2i)
            edge_index = edge_index.to(device)
            branch_lens = edge_attr.squeeze(-1).to(device)

            ctx_emb, _ = tree_enc(
                node_emb, node_ids, node_times_dict,
                edge_index, branch_lens, t_scalar=1.0,
            )

            torch.save({
                "node_ids":  node_ids,
                "plm":       plm_T1.cpu(),       # [N, 320] raw ESM2
                "node_emb":  node_emb.cpu(),     # [N, 128] NodeEncoder (PLM+struct+lap)
                "ctx_emb":   ctx_emb.cpu(),      # [N, 128] TreeEncoder contextual
            }, out_path)
            print(f"  group_{g:03d}: plm={tuple(plm_T1.shape)}  "
                  f"node_emb={tuple(node_emb.shape)}  ctx_emb={tuple(ctx_emb.shape)}")


if __name__ == "__main__":
    main()
