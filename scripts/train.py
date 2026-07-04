#!/usr/bin/env python3
"""
Bridge matching training loop for TreeSBM (Algorithm 1).

Each step:
  1. Load T1 from dataset, compute PLM embeddings once.
  2. Sample t ~ U(0, t_max_clip).
  3. Construct T_t via SampleBridgeState (Algorithm 2).
  4. Rebuild TreeState / structural features / Laplacian PE for T_t.
  5. Run NodeEncoder → TreeEncoder (with time conditioning) → H_t.
  6. Run RateHeads for active leaves of T_t.
  7. Compute bridge matching loss vs. T1 targets.

Usage:
    python scripts/train.py --data data/train --epochs 100 --lr 1e-4
"""

import argparse
import random
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split

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
    embedder: ESM2Embedder | None = None,
) -> tuple[dict | None, int]:
    """
    One forward pass of Algorithm 1.

    Returns (losses_dict, n_active_leaves).
    Returns (None, 0) if T_t has no usable structure.
    """
    node_ids        = batch["node_ids"]              # list[str], N
    node_times_t    = batch["node_times"]            # [N] tensor
    seqs            = batch["seqs"]                  # dict[str, str]
    edges           = batch["edges"]                 # list[(parent, child)]
    branch_lengths  = batch["branch_lengths"]        # dict[(str,str), float]
    root_id         = node_ids[batch["root_index"]]

    node_times_dict = {nid: node_times_t[i].item() for i, nid in enumerate(node_ids)}

    # ── 1. PLM embeddings for all T1 nodes ───────────────────────────────────
    if batch.get("plm_embeddings") is not None:
        plm_T1 = batch["plm_embeddings"].to(device)              # [N, 320] cached
    elif embedder is not None:
        sequences = [seqs[nid] for nid in node_ids]
        with torch.no_grad():
            plm_T1 = embedder.embed_sequences(sequences).to(device)
    else:
        raise RuntimeError("No PLM embeddings: run scripts/precompute_plm.py first")
    plm_map = {nid: i for i, nid in enumerate(node_ids)}

    # ── 2. Sample intermediate tree T_t ───────────────────────────────────────
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

    # ── 3. Build TreeState for T_t ────────────────────────────────────────────
    tree_t = TreeState(
        node_ids=node_ids_t,
        root_id=root_id,
        edges=edges_t,
        branch_lengths=branch_lengths_t,
        node_seqs=seqs_t,
        active_leaves=active_leaves_t,
    )
    node_to_idx_t = {nid: i for i, nid in enumerate(node_ids_t)}

    # ── 4. Structural features + Laplacian PE for T_t (recomputed on device) ─
    struct_t = compute_structural_features(tree_t, node_to_idx_t).to(device)
    lap_t    = compute_laplacian_pe(tree_t, node_to_idx_t, lap_dim, device=device)

    # ── 5. PLM embeddings for T_t nodes (use T1 embeddings, no re-embedding) ─
    plm_t = torch.stack([plm_T1[plm_map[nid]] for nid in node_ids_t]).to(device)  # [N_t, 320]

    # ── 6. NodeEncoder ────────────────────────────────────────────────────────
    h_t = node_enc(plm_t.to(device), struct_t.to(device), lap_t.to(device))  # [N_t, 128]

    # ── 7. Edge tensors for T_t ───────────────────────────────────────────────
    edge_index_t, _, edge_attr_t = build_edges(tree_t, node_to_idx_t)
    edge_index_t = edge_index_t.to(device)
    branch_lens_t = edge_attr_t.squeeze(-1).to(device)  # [2E_t]

    # ── 8. TreeEncoder with time conditioning ────────────────────────────────
    H_t, _ = tree_enc(
        h_t, node_ids_t, node_times_dict, edge_index_t, branch_lens_t, t_scalar=t
    )  # [N_t, 128]

    # ── 9. RateHeads for active leaves of T_t ─────────────────────────────────
    active_idx_t = [node_to_idx_t[nid] for nid in active_leaves_t]
    out = rate_heads(H_t, active_idx_t)

    # ── 10. Bridge matching loss ───────────────────────────────────────────────
    losses = bridge_losses(
        out=out,
        active_leaves=active_leaves_t,
        T1_seqs=seqs,
        T1_child_counts=T_t["T1_child_counts"],
        T1_child_bls=T_t["T1_child_bls"],
        t=t,
        max_seq_len=max_seq_len,
        lambda_top=lambda_top,
        lambda_br=lambda_br,
        device=device,
    )

    return losses, len(active_leaves_t)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",        default="data/train")
    parser.add_argument("--epochs",      type=int,   default=100)
    parser.add_argument("--lr",          type=float, default=1e-4)
    parser.add_argument("--val-frac",    type=float, default=0.1)
    parser.add_argument("--ckpt-dir",    default="checkpoints")
    parser.add_argument("--seed",        type=int,   default=42)
    parser.add_argument("--t-max",       type=float, default=0.95,
                        help="Max t to sample (avoids 1/(1-t) blow-up near t=1)")
    parser.add_argument("--lambda-top",  type=float, default=0.1)
    parser.add_argument("--lambda-br",   type=float, default=0.1)
    parser.add_argument("--max-seq-len", type=int,   default=566)
    parser.add_argument("--patience",    type=int,   default=30,
                        help="Early stopping: stop if val loss doesn't improve for this many epochs")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    # ── data ──────────────────────────────────────────────────────────────────
    dataset = TreeDataset(args.data, max_seq_len=args.max_seq_len)
    n_val   = max(1, int(len(dataset) * args.val_frac))
    n_train = len(dataset) - n_val
    train_ds, val_ds = random_split(dataset, [n_train, n_val])
    print(f"Train: {n_train}  Val: {n_val}")

    # batch_size=1 (trees vary in node count — no collation)
    train_loader = DataLoader(train_ds, batch_size=1, shuffle=True,  collate_fn=lambda x: x[0])
    val_loader   = DataLoader(val_ds,   batch_size=1, shuffle=False, collate_fn=lambda x: x[0])

    # ── models ────────────────────────────────────────────────────────────────
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
    rate_heads = RateHeads(d_model=128, max_seq_len=args.max_seq_len).to(device)

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

    # ── training loop ─────────────────────────────────────────────────────────
    for epoch in range(1, args.epochs + 1):
        node_enc.train(); tree_enc.train(); rate_heads.train()
        train_loss = 0.0
        n_steps = 0
        loss_breakdown = {"L_seq": 0.0, "L_top": 0.0, "L_br": 0.0}

        for batch in train_loader:
            t = random.uniform(0.0, args.t_max)
            optimizer.zero_grad()

            losses, n_active = forward_bridge_step(
                batch, t, node_enc, tree_enc, rate_heads, device,
                max_seq_len=args.max_seq_len,
                lambda_top=args.lambda_top,
                lambda_br=args.lambda_br,
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

        # ── validation (fixed t=0.5 for reproducibility) ─────────────────────
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
            f"(seq={loss_breakdown['L_seq']:.3f} "
            f"top={loss_breakdown['L_top']:.3f} "
            f"br={loss_breakdown['L_br']:.3f})  "
            f"val={val_loss:.4f}  lr={lr:.2e}"
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
            }, ckpt_dir / "best.pt")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"\nEarly stopping at epoch {epoch} (no improvement for {args.patience} epochs)")
                break

    print(f"\nBest val loss: {best_val:.4f}  →  {ckpt_dir}/best.pt")

    # ── export embeddings for all trees using trained weights ─────────────────
    print("\nExporting embeddings with trained weights...")
    export_embeddings(dataset, node_enc, tree_enc, device, Path(args.data))
    print("Done.")


def export_embeddings(
    dataset: TreeDataset,
    node_enc: NodeEncoder,
    tree_enc: TreeEncoder,
    device: str,
    data_dir: Path,
):
    """
    For each tree in the dataset, save:
      - plm:       [N, 320]  raw ESM2 embeddings (already cached, just copied)
      - node_emb:  [N, 128]  NodeEncoder output (PLM + structural + Laplacian PE fused)
      - ctx_emb:   [N, 128]  TreeEncoder output at t=1.0 (fully contextualized)

    Saved to data_dir/group_NNN_trained_emb.pt
    """
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

            # Build full T1 TreeState for structural features / Laplacian
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
