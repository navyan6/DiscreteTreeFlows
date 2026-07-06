"""
PyTorch Dataset over processed phylogenetic tree groups.

Each item is one tree: loads rooted NWK + ancestral AA FASTA + branch_lengths JSON
and returns everything needed for one forward pass.
"""

import json
from pathlib import Path

import torch
from torch.utils.data import Dataset
from Bio import Phylo, SeqIO

from src.tree_state import TreeState
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.edges import build_edges

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}
PAD_IDX = len(AA_VOCAB)  # 20


def aa_seq_to_tensor(seq: str, length: int) -> torch.Tensor:
    """Convert AA string to [length] int tensor. Unknown AAs → PAD_IDX."""
    t = torch.full((length,), PAD_IDX, dtype=torch.long)
    for i, aa in enumerate(seq[:length]):
        t[i] = AA_TO_IDX.get(aa, PAD_IDX)
    return t


def parse_newick(nwk_path: str):
    tree = Phylo.read(nwk_path, "newick")
    edges, branch_lengths = [], {}
    _c = [0]

    def name(clade):
        if clade.name:
            return clade.name
        n = f"NODE_{_c[0]:07d}"
        clade.name = n
        _c[0] += 1
        return n

    def walk(parent):
        pn = name(parent)
        for child in parent.clades:
            cn = name(child)
            bl = child.branch_length or 0.0
            edges.append((pn, cn))
            branch_lengths[(pn, cn)] = bl
            walk(child)

    walk(tree.root)
    root_id = name(tree.root)

    # BFS order
    node_ids = [root_id]
    visited = {root_id}
    children = {}
    for p, c in edges:
        children.setdefault(p, []).append(c)
    queue = [root_id]
    while queue:
        curr = queue.pop(0)
        for ch in children.get(curr, []):
            if ch not in visited:
                visited.add(ch)
                node_ids.append(ch)
                queue.append(ch)

    return root_id, node_ids, edges, branch_lengths


class TreeDataset(Dataset):
    def __init__(self, data_dir: str, laplacian_dim: int = 8, max_seq_len: int = 566):
        self.data_dir = Path(data_dir)
        self.laplacian_dim = laplacian_dim
        self.max_seq_len = max_seq_len

        # Find all complete groups (need all 3 files)
        self.groups = sorted([
            int(p.stem.split("_")[1])
            for p in self.data_dir.glob("group_*_rooted.nwk")
            if (self.data_dir / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists()
            and (self.data_dir / p.name.replace("_rooted.nwk", "_bl.json")).exists()
        ])
        print(f"TreeDataset: {len(self.groups)} trees found")

    def __len__(self):
        return len(self.groups)

    def __getitem__(self, idx: int) -> dict:
        g = self.groups[idx]
        d = self.data_dir

        root_id, node_ids, edges, branch_lengths = parse_newick(
            str(d / f"group_{g:03d}_rooted.nwk")
        )

        seqs = {
            rec.id: str(rec.seq)
            for rec in SeqIO.parse(d / f"group_{g:03d}_anc_aa.fasta", "fasta")
        }
        ref_len = len(next(iter(seqs.values())))
        for nid in node_ids:
            if nid not in seqs:
                seqs[nid] = "-" * ref_len

        with open(d / f"group_{g:03d}_bl.json") as f:
            node_data = json.load(f)["nodes"]
        node_times = {
            nid: node_data.get(nid, {}).get("numdate", 0.0) for nid in node_ids
        }

        has_children = {p for p, _ in edges}
        active_leaves = [nid for nid in node_ids if nid not in has_children]

        tree_state = TreeState(
            node_ids=node_ids, root_id=root_id, edges=edges,
            branch_lengths=branch_lengths, node_seqs=seqs,
            active_leaves=active_leaves,
        )

        node_to_idx = {nid: i for i, nid in enumerate(node_ids)}

        structural = compute_structural_features(tree_state, node_to_idx)
        lap_pe = compute_laplacian_pe(tree_state, node_to_idx, self.laplacian_dim)
        edge_index, edge_type, edge_attr = build_edges(tree_state, node_to_idx)

        # Target: AA sequences as integer tensors [N, max_seq_len]
        targets = torch.stack([
            aa_seq_to_tensor(seqs[nid], self.max_seq_len) for nid in node_ids
        ])

        # Active leaf indices
        leaf_indices = [node_to_idx[nid] for nid in active_leaves]

        # Load cached PLM embeddings if available
        plm_path = d / f"group_{g:03d}_plm.pt"
        if plm_path.exists():
            cached = torch.load(plm_path, weights_only=True)
            plm_embeddings = cached["plm"]  # [N, 320], BFS-ordered to match node_ids
        else:
            plm_embeddings = None

        # Load cached reference mutation log-rates if available
        ref_path = d / f"group_{g:03d}_ref_rates.pt"
        if ref_path.exists():
            log_ref_mut_rates = torch.load(ref_path, weights_only=True)["log_mut_rates"]
        else:
            log_ref_mut_rates = None  # [N, 566, 20] or None

        return {
            "group": g,
            "node_ids": node_ids,
            "node_times": torch.tensor([node_times[nid] for nid in node_ids], dtype=torch.float32),
            "structural_features": structural,        # [N, 3]
            "lap_pe": lap_pe,                         # [N, lap_dim]
            "edge_index": edge_index,                 # [2, 2E]
            "edge_attr": edge_attr,                   # [2E, 1]
            "targets": targets,                       # [N, max_seq_len]
            "leaf_indices": leaf_indices,
            "root_index": node_to_idx[root_id],
            "seqs": seqs,                             # for ESM2 embedding
            # raw graph topology (needed by SampleBridgeState)
            "edges": edges,
            "branch_lengths": branch_lengths,
            # precomputed ESM2 [N, 320] or None if not yet cached
            "plm_embeddings": plm_embeddings,
            # precomputed log R0 mutation rates [N, 566, 20] or None
            "log_ref_mut_rates": log_ref_mut_rates,
        }
