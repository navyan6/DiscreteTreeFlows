#!/usr/bin/env python3
"""
Comprehensive evaluation of a single generated tree.

  A. Sequence quality     — root divergence, leaf diversity, ESM PLL, post-prune survival
  B. Phylogenetic coherence — sibling/parent vs random-pair sequence identity
  C. Tree structure       — bifurcating check, depth, branch lengths, Sackin index, cherries
  D. GT comparison        — depth/branch distributions, best-match seq identity

Usage:
    python scripts/eval_single_tree.py \
        --checkpoint checkpoints/best.pt \
        --data data/train \
        --group 1 \
        --n-steps 50 \
        --max-leaves 200 \
        --branch-rate-scale 6.0
"""

import argparse
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM

from src.dataset import TreeDataset
from src.tree_state import TreeState
from src.treeencoder.node_encoder import NodeEncoder
from src.treeencoder.plm_embeddings import ESM2Embedder
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.edges import build_edges
from src.networks import TreeEncoder, RateHeads

AA_VOCAB  = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}


# ── utilities ─────────────────────────────────────────────────────────────────

def seq_identity(a: str, b: str) -> float:
    L = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:L], b[:L])) / L if L else 0.0


def get_lm_logits(tokenizer, esm_model, aa_token_ids, sequences, max_seq_len, device):
    log_rates = torch.zeros(len(sequences), max_seq_len, 20, dtype=torch.float32, device=device)
    with torch.no_grad():
        tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                           truncation=False).to(device)
        logits = esm_model(**tokens).logits
    seq_lens = tokens["attention_mask"].sum(dim=1)
    for i in range(len(sequences)):
        actual_L = int(seq_lens[i].item()) - 2
        aa_logits = logits[i, 1:actual_L + 1, :][:, aa_token_ids]
        log_probs = F.log_softmax(aa_logits, dim=-1)
        clip = min(actual_L, max_seq_len)
        log_rates[i, :clip, :] = log_probs[:clip]
    return log_rates


def esm_pll_seq(log_R0_i: torch.Tensor, seq: str, max_seq_len: int) -> float:
    vals = [log_R0_i[pos, AA_TO_IDX[aa]].item()
            for pos, aa in enumerate(seq[:max_seq_len]) if aa in AA_TO_IDX]
    return sum(vals) / len(vals) if vals else float("-inf")


def children_map(tree: TreeState) -> dict:
    cm = defaultdict(list)
    for p, c in tree.edges:
        cm[p].append(c)
    return dict(cm)


def get_leaves(tree: TreeState) -> list[str]:
    cm = children_map(tree)
    return [n for n in tree.node_ids if n not in cm]


def bfs_depths(tree: TreeState) -> dict[str, int]:
    depths = {tree.root_id: 0}
    queue = [tree.root_id]
    cm = children_map(tree)
    while queue:
        node = queue.pop(0)
        for child in cm.get(node, []):
            depths[child] = depths[node] + 1
            queue.append(child)
    return depths


def sackin_index(tree: TreeState, leaves: list[str]) -> int:
    depths = bfs_depths(tree)
    return sum(depths.get(l, 0) for l in leaves)


def cherry_count(tree: TreeState) -> int:
    cm = children_map(tree)
    leaf_set = set(get_leaves(tree))
    return sum(1 for ch in cm.values()
               if len(ch) == 2 and all(c in leaf_set for c in ch))


def section(title: str):
    print(f"\n{'='*60}\n  {title}\n{'='*60}")


# ── model loading ─────────────────────────────────────────────────────────────

def load_models(checkpoint, device, max_seq_len):
    ckpt = torch.load(checkpoint, map_location=device, weights_only=False)
    print(f"Checkpoint: epoch {ckpt.get('epoch','?')}  "
          f"val_loss={ckpt.get('val_loss', 0):.4f}")
    node_enc = NodeEncoder(d_plm=320, d_struct=3, d_laplacian=8, d_node=128).to(device)
    tree_enc = TreeEncoder(d_model=128, n_layers=4, n_heads=8, dropout=0.1).to(device)
    r_heads  = RateHeads(d_model=128, max_seq_len=max_seq_len).to(device)
    node_enc.load_state_dict(ckpt["node_enc"])
    tree_enc.load_state_dict(ckpt["tree_enc"])
    r_heads.load_state_dict(ckpt["rate_heads"])
    for m in [node_enc, tree_enc, r_heads]:
        m.eval()
        for p in m.parameters():
            p.requires_grad = False
    return node_enc, tree_enc, r_heads


# ── generation ────────────────────────────────────────────────────────────────

def generate_tree(root_seq, n_steps, max_seq_len, branch_rate_scale, max_leaves, mutation_rate_scale,
                  node_enc, tree_enc, rate_heads, embedder,
                  tokenizer, esm_model, aa_token_ids, device):
    tree = TreeState.root_only(root_seq)
    node_birth_step = {tree.root_id: 0}
    dt = 1.0 / n_steps

    for step in range(n_steps):
        t = step / n_steps
        if not tree.active_leaves:
            break

        node_ids_t  = tree.node_ids
        node_to_idx = {nid: i for i, nid in enumerate(node_ids_t)}
        active_leaves = list(tree.active_leaves)
        active_idx    = [node_to_idx[v] for v in active_leaves]
        node_times_dict = {nid: node_birth_step.get(nid, 0) / n_steps
                           for nid in node_ids_t}

        struct_t = compute_structural_features(tree, node_to_idx).to(device)
        lap_t    = compute_laplacian_pe(tree, node_to_idx, 8, device=device)
        edge_index_t, _, edge_attr_t = build_edges(tree, node_to_idx)
        edge_index_t  = edge_index_t.to(device)
        branch_lens_t = edge_attr_t.squeeze(-1).to(device)
        plm_t = embedder.embed_sequences(
            [tree.node_seqs[nid] for nid in node_ids_t]).to(device)

        active_seqs = [tree.node_seqs[v] for v in active_leaves]
        log_R0_mut  = get_lm_logits(tokenizer, esm_model, aa_token_ids,
                                     active_seqs, max_seq_len, device)
        with torch.no_grad():
            h_t     = node_enc(plm_t, struct_t, lap_t)
            H_t, _  = tree_enc(h_t, node_ids_t, node_times_dict,
                                edge_index_t, branch_lens_t, t_scalar=t)
            out     = rate_heads(H_t, active_idx, log_R0_mut)

        new_node_seqs = dict(tree.node_seqs)

        for i, leaf_id in enumerate(active_leaves):
            seq     = tree.node_seqs[leaf_id]
            seq_len = min(len(seq), max_seq_len)
            new_seq = list(seq)

            for pos in range(seq_len):
                curr_idx = AA_TO_IDX.get(seq[pos], -1)
                if curr_idx < 0:
                    continue
                probs = out["log_R_theta_mut"][i, pos].softmax(-1)
                probs_mut = probs.clone()
                probs_mut[curr_idx] = 0.0
                total = probs_mut.sum().item()
                if total > 0 and torch.rand(1).item() < total * dt * mutation_rate_scale:
                    new_seq[pos] = AA_VOCAB[torch.multinomial(probs_mut / total, 1).item()]
            new_node_seqs[leaf_id] = "".join(new_seq)

            at_cap   = len(tree.active_leaves) >= max_leaves
            lam      = out["branching_rate"][i].item() * branch_rate_scale
            p_branch = 1.0 - math.exp(-max(0.0, lam) * dt)
            n_ch     = 0 if at_cap else (2 if torch.rand(1).item() < p_branch else 0)

            if n_ch > 0:
                child_seqs = [new_node_seqs[leaf_id]] * n_ch
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges, branch_lengths=tree.branch_lengths,
                    node_seqs=new_node_seqs, active_leaves=list(tree.active_leaves))
                tree = tree.branch_node(leaf_id, child_seqs)
                bl_pred = out["branch_length"][i].item()
                new_children = tree.get_children(leaf_id)
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges,
                    branch_lengths={**tree.branch_lengths,
                                    **{(leaf_id, c): bl_pred for c in new_children}},
                    node_seqs=tree.node_seqs,
                    active_leaves=list(tree.active_leaves))
                new_node_seqs = dict(tree.node_seqs)
                for child_id in new_children:
                    node_birth_step.setdefault(child_id, step + 1)

        tree = TreeState(
            node_ids=tree.node_ids, root_id=tree.root_id,
            edges=tree.edges, branch_lengths=tree.branch_lengths,
            node_seqs=new_node_seqs, active_leaves=list(tree.active_leaves))

    return tree


# ── evaluation sections ───────────────────────────────────────────────────────

def eval_sequence_quality(gen_tree, gen_leaves, root_seq, max_seq_len,
                          tokenizer, esm_model, aa_token_ids, device,
                          pll_prune_threshold):
    section("A. SEQUENCE QUALITY")

    # Root-to-leaf divergence
    divs = [1.0 - seq_identity(root_seq, gen_tree.node_seqs[l]) for l in gen_leaves]
    print(f"Root-to-leaf divergence (fraction of positions mutated from root):")
    print(f"  mean={sum(divs)/len(divs):.4f}  "
          f"min={min(divs):.4f}  max={max(divs):.4f}")

    # Leaf-to-leaf pairwise diversity
    sample = random.sample(gen_leaves, min(60, len(gen_leaves)))
    pairs  = [(a, b) for i, a in enumerate(sample) for b in sample[i+1:]]
    pairs  = random.sample(pairs, min(300, len(pairs)))
    leaf_div = [1.0 - seq_identity(gen_tree.node_seqs[a], gen_tree.node_seqs[b])
                for a, b in pairs]
    print(f"Leaf-to-leaf pairwise diversity ({len(pairs)} pairs):")
    print(f"  mean={sum(leaf_div)/len(leaf_div):.4f}  "
          f"min={min(leaf_div):.4f}  max={max(leaf_div):.4f}")

    # ESM PLL
    leaf_seqs = [gen_tree.node_seqs[l] for l in gen_leaves]
    print(f"Computing ESM-2 PLL for {len(gen_leaves)} generated leaves...")
    log_R0 = get_lm_logits(tokenizer, esm_model, aa_token_ids,
                            leaf_seqs, max_seq_len, device)
    plls = [esm_pll_seq(log_R0[i], leaf_seqs[i], max_seq_len)
            for i in range(len(gen_leaves))]
    print(f"ESM PLL (nats/position):")
    print(f"  mean={sum(plls)/len(plls):.4f}  "
          f"min={min(plls):.4f}  max={max(plls):.4f}")

    # Pruning
    surviving = [l for l, p in zip(gen_leaves, plls) if p >= pll_prune_threshold]
    print(f"Post-prune survival (PLL >= {pll_prune_threshold} nats/pos):")
    print(f"  {len(surviving)}/{len(gen_leaves)} leaves survive  "
          f"({100*len(surviving)/len(gen_leaves):.0f}%)")

    return plls, surviving


def eval_phylogenetic_coherence(gen_tree, gen_leaves):
    section("B. PHYLOGENETIC COHERENCE")
    print("Tests whether nearby nodes in the tree are more similar than random pairs.")

    cm      = children_map(gen_tree)
    leaf_set = set(gen_leaves)

    # Sibling leaf pairs (both children of same parent are leaves)
    sibling_ids = []
    for node, children in cm.items():
        if len(children) == 2:
            a, b = children
            if a in leaf_set and b in leaf_set:
                sibling_ids.append(seq_identity(
                    gen_tree.node_seqs[a], gen_tree.node_seqs[b]))

    # Parent-child pairs where child is a leaf
    parent_child_ids = [
        seq_identity(gen_tree.node_seqs[p], gen_tree.node_seqs[c])
        for p, c in gen_tree.edges if c in leaf_set
    ]

    # Random leaf pairs (baseline)
    sample  = random.sample(gen_leaves, min(60, len(gen_leaves)))
    pairs   = [(a, b) for i, a in enumerate(sample) for b in sample[i+1:]]
    pairs   = random.sample(pairs, min(300, len(pairs)))
    random_ids = [seq_identity(gen_tree.node_seqs[a], gen_tree.node_seqs[b])
                  for a, b in pairs]

    def fmt(vals, label):
        if not vals:
            return f"  {label}: N/A"
        m = sum(vals) / len(vals)
        return f"  {label} (n={len(vals):3d}): mean={m:.4f}  min={min(vals):.4f}  max={max(vals):.4f}"

    print(fmt(sibling_ids,      "Sibling leaf pairs    "))
    print(fmt(parent_child_ids, "Parent-child pairs    "))
    print(fmt(random_ids,       "Random leaf pairs     "))

    if sibling_ids and random_ids:
        delta = sum(sibling_ids)/len(sibling_ids) - sum(random_ids)/len(random_ids)
        result = "PASS ✓" if delta > 0 else "FAIL ✗"
        print(f"\n  Sibling vs random delta: {delta:+.4f}  [{result}]")


def eval_tree_structure(gen_tree, gen_leaves):
    section("C. TREE STRUCTURE")

    cm       = children_map(gen_tree)
    internal = [n for n in gen_tree.node_ids if n in cm]

    # Bifurcating
    counts = [len(cm[n]) for n in internal]
    all_bif = all(c == 2 for c in counts)
    print(f"Strictly bifurcating: {all_bif}  (child counts seen: {sorted(set(counts))})")

    # Depth
    depths      = bfs_depths(gen_tree)
    leaf_depths = [depths.get(l, 0) for l in gen_leaves]
    print(f"Leaf depth from root:")
    print(f"  mean={sum(leaf_depths)/len(leaf_depths):.1f}  "
          f"min={min(leaf_depths)}  max={max(leaf_depths)}")

    # Branch lengths
    bls = list(gen_tree.branch_lengths.values())
    if bls:
        print(f"Branch lengths ({len(bls)} edges, all positive: {all(b>0 for b in bls)}):")
        print(f"  mean={sum(bls)/len(bls):.6f}  "
              f"min={min(bls):.6f}  max={max(bls):.6f}")

    # Sackin + cherries
    sak     = sackin_index(gen_tree, gen_leaves)
    cherries = cherry_count(gen_tree)
    print(f"Sackin index:  {sak}  (lower = more balanced)")
    print(f"Cherry count:  {cherries}  (sibling leaf pairs)")

    return leaf_depths, bls


def eval_gt_comparison(gen_tree, gen_leaves, gt_batch):
    section("D. GROUND TRUTH COMPARISON")

    gt_node_ids = gt_batch["node_ids"]
    gt_seqs     = gt_batch["seqs"]
    gt_edges    = gt_batch["edges"]
    gt_bls_dict = gt_batch["branch_lengths"]
    gt_root_id  = gt_node_ids[gt_batch["root_index"]]

    gt_cm      = defaultdict(list)
    for p, c in gt_edges:
        gt_cm[p].append(c)
    gt_leaf_set = set(n for n in gt_node_ids if n not in gt_cm)
    gt_leaves   = list(gt_leaf_set)
    # Filter seqs to node_ids only (FASTA may have extra entries)
    gt_seqs_filtered = {nid: gt_seqs.get(nid, "") for nid in gt_node_ids}

    print(f"GT:  {len(gt_node_ids)} nodes, {len(gt_leaves)} leaves")
    print(f"Gen: {len(gen_tree.node_ids)} nodes, {len(gen_leaves)} leaves")

    # Best-match seq identity: each GT leaf → closest generated leaf
    sample_gt = random.sample(gt_leaves, min(100, len(gt_leaves)))
    best_ids  = [max(seq_identity(gt_seqs[gl], gen_tree.node_seqs[gl2])
                     for gl2 in gen_leaves)
                 for gl in sample_gt]
    print(f"\nBest-match identity (GT leaf → nearest gen leaf, n={len(sample_gt)}):")
    print(f"  mean={sum(best_ids)/len(best_ids):.4f}  "
          f"min={min(best_ids):.4f}  max={max(best_ids):.4f}")

    # Branch length distributions
    gt_bls_vals  = list(gt_bls_dict.values())
    gen_bls_vals = list(gen_tree.branch_lengths.values())
    if gt_bls_vals and gen_bls_vals:
        print(f"\nBranch length distribution:")
        print(f"  GT:  mean={sum(gt_bls_vals)/len(gt_bls_vals):.6f}  "
              f"max={max(gt_bls_vals):.6f}")
        print(f"  Gen: mean={sum(gen_bls_vals)/len(gen_bls_vals):.6f}  "
              f"max={max(gen_bls_vals):.6f}")

    # Depth distribution
    gt_tree_obj = TreeState(
        node_ids=gt_node_ids, root_id=gt_root_id,
        edges=gt_edges, branch_lengths=gt_bls_dict,
        node_seqs=gt_seqs_filtered, active_leaves=list(gt_leaf_set))
    gt_depths     = bfs_depths(gt_tree_obj)
    gt_leaf_depths = [gt_depths.get(l, 0) for l in gt_leaves]
    gen_depths_map = bfs_depths(gen_tree)
    gen_leaf_depths = [gen_depths_map.get(l, 0) for l in gen_leaves]

    print(f"\nLeaf depth from root:")
    print(f"  GT:  mean={sum(gt_leaf_depths)/len(gt_leaf_depths):.1f}  "
          f"min={min(gt_leaf_depths)}  max={max(gt_leaf_depths)}")
    print(f"  Gen: mean={sum(gen_leaf_depths)/len(gen_leaf_depths):.1f}  "
          f"min={min(gen_leaf_depths)}  max={max(gen_leaf_depths)}")

    # Sackin index + cherries
    gt_sak  = sackin_index(gt_tree_obj, gt_leaves)
    gen_sak = sackin_index(gen_tree, gen_leaves)
    gt_ch   = cherry_count(gt_tree_obj)
    gen_ch  = cherry_count(gen_tree)
    print(f"\nSackin index:  GT={gt_sak}  Gen={gen_sak}  "
          f"(ratio={gen_sak/gt_sak:.2f})" if gt_sak > 0 else "")
    print(f"Cherry count:  GT={gt_ch}    Gen={gen_ch}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",          default="checkpoints/best.pt")
    parser.add_argument("--data",                default="data/train")
    parser.add_argument("--group",               type=int, required=True)
    parser.add_argument("--n-steps",             type=int,   default=50)
    parser.add_argument("--max-seq-len",         type=int,   default=566)
    parser.add_argument("--branch-rate-scale",   type=float, default=6.0)
    parser.add_argument("--max-leaves",          type=int,   default=200)
    parser.add_argument("--pll-prune-threshold",  type=float, default=-3.0)
    parser.add_argument("--mutation-rate-scale",  type=float, default=1.0,
                        help="Multiply CTMC mutation rate by this; >1 forces more mutations")
    parser.add_argument("--seed",                type=int,   default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    node_enc, tree_enc, rate_heads = load_models(args.checkpoint, device, args.max_seq_len)
    embedder = ESM2Embedder(device=device)

    model_id = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    esm_model = EsmForMaskedLM.from_pretrained(model_id).to(device)
    esm_model.eval()
    for p in esm_model.parameters():
        p.requires_grad = False
    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB], dtype=torch.long)

    dataset  = TreeDataset(args.data, max_seq_len=args.max_seq_len)
    group_idx = next((i for i in range(len(dataset))
                      if dataset.groups[i] == args.group), None)
    if group_idx is None:
        print(f"Group {args.group} not found"); sys.exit(1)

    gt_batch = dataset[group_idx]
    root_id  = gt_batch["node_ids"][gt_batch["root_index"]]
    root_seq = gt_batch["seqs"][root_id]
    print(f"\nGroup {args.group}: root_len={len(root_seq)}  "
          f"gt_nodes={len(gt_batch['node_ids'])}")
    print(f"Generating ({args.n_steps} steps, max_leaves={args.max_leaves}, "
          f"scale={args.branch_rate_scale})...")

    gen_tree  = generate_tree(
        root_seq, args.n_steps, args.max_seq_len,
        args.branch_rate_scale, args.max_leaves, args.mutation_rate_scale,
        node_enc, tree_enc, rate_heads, embedder,
        tokenizer, esm_model, aa_token_ids, device)
    gen_leaves = get_leaves(gen_tree)
    print(f"Generated: {len(gen_tree.node_ids)} nodes, {len(gen_leaves)} leaves")

    eval_sequence_quality(gen_tree, gen_leaves, root_seq, args.max_seq_len,
                          tokenizer, esm_model, aa_token_ids, device,
                          args.pll_prune_threshold)
    eval_phylogenetic_coherence(gen_tree, gen_leaves)
    eval_tree_structure(gen_tree, gen_leaves)
    eval_gt_comparison(gen_tree, gen_leaves, gt_batch)

    # Save sequences
    out_path = Path("checkpoints") / f"gen_group{args.group}.fasta"
    cm = children_map(gen_tree)
    with open(out_path, "w") as f:
        for nid in gen_tree.node_ids:
            tag = ("root" if nid == gen_tree.root_id
                   else ("leaf" if nid not in cm else "internal"))
            f.write(f">{nid}|{tag}\n{gen_tree.node_seqs[nid]}\n")
    print(f"\nSequences saved to {out_path}")


if __name__ == "__main__":
    main()
