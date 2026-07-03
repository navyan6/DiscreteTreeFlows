#!/usr/bin/env python3
"""
Load ASR outputs → build TreeState → run graph transformer → save node/edge embeddings.

Usage:
    python scripts/embed_tree.py \
        --tree       data/train/group_001_rooted.nwk \
        --aa-fasta   data/train/ancestral_001_aa.fasta \
        --node-data  data/train/group_001_branch_lengths.json \
        --out        data/train/group_001_embeddings.pt
"""

import argparse
import json
from pathlib import Path

import torch
from Bio import Phylo, SeqIO
from io import StringIO

from src.tree_state import TreeState
from src.treeencoder.plm_embeddings import ESM2Embedder
from src.treeencoder.tree_adapter import tree_state_to_encoder_input
from src.treeencoder.node_encoder import NodeEncoder
from src.treeencoder.attention_mask import build_temporal_attention_mask


# ── tree parsing ──────────────────────────────────────────────────────────────

def parse_newick(nwk_path: str) -> tuple[str, list[tuple[str, str]], dict]:
    """
    Parse a Newick file into (root_id, edges, branch_lengths).
    Internal nodes without names are assigned NODE_XXXXXXX labels to match
    augur's naming convention (augur writes names into the Newick).
    """
    tree = Phylo.read(nwk_path, "newick")

    edges = []
    branch_lengths = {}
    _counter = [0]

    def node_name(clade):
        if clade.name:
            return clade.name
        name = f"NODE_{_counter[0]:07d}"
        clade.name = name
        _counter[0] += 1
        return name

    def walk(parent_clade):
        parent_name = node_name(parent_clade)
        for child in parent_clade.clades:
            child_name = node_name(child)
            bl = child.branch_length if child.branch_length is not None else 0.0
            edges.append((parent_name, child_name))
            branch_lengths[(parent_name, child_name)] = bl
            walk(child)

    root_clade = tree.root
    walk(root_clade)
    root_id = node_name(root_clade)

    # Collect all node IDs in BFS order
    all_ids = [root_id]
    visited = {root_id}
    queue = [root_id]
    children_map = {}
    for p, c in edges:
        children_map.setdefault(p, []).append(c)
    while queue:
        curr = queue.pop(0)
        for ch in children_map.get(curr, []):
            if ch not in visited:
                visited.add(ch)
                all_ids.append(ch)
                queue.append(ch)

    return root_id, all_ids, edges, branch_lengths


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree",      required=True, help="Rooted Newick (.nwk)")
    parser.add_argument("--aa-fasta",  required=True, help="AA FASTA (tips + internal nodes)")
    parser.add_argument("--node-data", required=True, help="augur branch_lengths JSON (contains numdate)")
    parser.add_argument("--out",       required=True, help="Output .pt file")
    parser.add_argument("--lap-dim",   type=int, default=8,  help="Laplacian PE dimension")
    parser.add_argument("--d-node",    type=int, default=128, help="NodeEncoder output dim")
    parser.add_argument("--device",    default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    device = args.device
    print(f"Device: {device}")

    # 1. Parse tree
    print("Parsing tree...")
    root_id, node_ids, edges, branch_lengths = parse_newick(args.tree)
    print(f"  {len(node_ids)} nodes, {len(edges)} edges, root={root_id}")

    # 2. Load AA sequences
    print("Loading AA sequences...")
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.aa_fasta, "fasta")}
    missing = [nid for nid in node_ids if nid not in seqs]
    if missing:
        print(f"  WARNING: {len(missing)} nodes have no sequence — filling with gaps")
        ref_len = len(next(iter(seqs.values())))
        for nid in missing:
            seqs[nid] = "-" * ref_len
    print(f"  {len(seqs)} sequences loaded")

    # 3. Load node times (numdate) from augur branch_lengths JSON
    print("Loading node times...")
    with open(args.node_data) as f:
        node_data = json.load(f)["nodes"]
    node_times = {}
    for nid in node_ids:
        entry = node_data.get(nid, {})
        node_times[nid] = entry.get("numdate", 0.0)
    print(f"  Time range: {min(node_times.values()):.2f} – {max(node_times.values()):.2f}")

    # 4. Build TreeState
    print("Building TreeState...")
    # Identify active leaves (nodes with no children)
    has_children = {p for p, _ in edges}
    active_leaves = [nid for nid in node_ids if nid not in has_children]

    tree_state = TreeState(
        node_ids=node_ids,
        root_id=root_id,
        edges=edges,
        branch_lengths=branch_lengths,
        node_seqs=seqs,
        active_leaves=active_leaves,
    )
    print(f"  {tree_state.n_nodes()} nodes, {tree_state.n_leaves()} leaves")

    # 5. ESM2 embeddings (320d per node)
    print("Computing ESM2 embeddings (this may take a few minutes)...")
    embedder = ESM2Embedder(device=device)
    sequences = [seqs[nid] for nid in node_ids]
    plm_embeddings = embedder.embed_sequences(sequences)  # [N, 320]
    print(f"  PLM embeddings: {plm_embeddings.shape}")

    # 6. Build TreeEncoderInput (structural features + Laplacian PE + edges)
    print("Building encoder input...")
    encoder_input = tree_state_to_encoder_input(
        tree=tree_state,
        node_embeddings=plm_embeddings,
        laplacian_dim=args.lap_dim,
    )

    # 7. Fuse with NodeEncoder → 128d
    print("Running NodeEncoder...")
    node_enc = NodeEncoder(
        d_plm=320,
        d_struct=3,
        d_laplacian=args.lap_dim,
        d_node=args.d_node,
    ).to(device)
    node_enc.eval()

    with torch.no_grad():
        h = node_enc(
            encoder_input.x.to(device),
            encoder_input.structural_features.to(device),
            encoder_input.lap_pe.to(device),
        )  # [N, d_node]

    print(f"  Node embeddings: {h.shape}")

    # 8. Temporal attention mask
    print("Building temporal attention mask...")
    attn_mask = build_temporal_attention_mask(node_ids, node_times).to(device)  # [N, N]

    # 9. Save everything
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "node_ids": node_ids,
        "node_embeddings": h.cpu(),          # [N, d_node]
        "edge_index": encoder_input.edge_index,   # [2, 2E]
        "edge_attr": encoder_input.edge_attr,     # [2E, 1] branch lengths
        "edge_type": encoder_input.edge_type,     # [2E] 0=parent→child, 1=child→parent
        "attn_mask": attn_mask.cpu(),             # [N, N] bool
        "node_times": torch.tensor([node_times[nid] for nid in node_ids]),  # [N]
        "root_index": encoder_input.root_index,
    }
    torch.save(payload, out_path)
    print(f"\nSaved to {out_path}")
    print(f"  node_embeddings: {payload['node_embeddings'].shape}")
    print(f"  edge_index:      {payload['edge_index'].shape}")
    print(f"  attn_mask:       {payload['attn_mask'].shape}")


if __name__ == "__main__":
    main()
