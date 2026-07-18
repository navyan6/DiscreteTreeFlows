#!/usr/bin/env python3
"""
Orchestrator: run every method on the held-out-root generation task, validate
every sample, score both tracks, and write the long-format results.csv that
make_table.py aggregates.

Empirical track  -> per generated sample: RF / quartet / branch-W / terminal-edit
                    vs the single real target subtree (averaged over K, over roots).
Simulated track  -> Tree-JS / Split-JS between the K generated trees and the M
                    reference trees (per regime), plus mean gen-vs-ref distances.

Needs dendropy + pyvolve (baselines), torch/ESM (pLM-prior, TreeSBM); runs on the
cluster. Native methods only unless --checkpoint / ESM available.

    python benchmarks/run_table.py --test-data data/h3n2/test \
        --params benchmarks/results/params.json --checkpoint checkpoints/h3n2_temporal/best.pt \
        --N 16 --K 100 --M 50 --max-roots 100 --out benchmarks/results/results.csv
"""

import argparse
import csv
import json
import random
import time
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from src.tree_state import TreeState
from benchmarks.heldout.build_examples import build_examples
from benchmarks.methods.bd_methods import NeutralBD, EmpiricalBD
from benchmarks.methods.plm_prior import PLMPrior
from benchmarks import validity as V
from benchmarks.metrics.matched import sequence_matched_rf, quartet_distance, terminal_edit_distance
from benchmarks.metrics.branch_lengths import branch_length_wasserstein
from benchmarks.metrics.distributions import tree_js, split_js
from benchmarks.sim.reference import simulate_reference, REGIMES

FIELDS = ["method", "track", "sim_regime", "root_id", "N", "H", "sample_seed",
          "valid", "tree_kl", "split_kl", "rf", "quartet", "branch_w_all",
          "terminal_edit", "runtime"]


def rebuild_target(ex: dict) -> TreeState:
    edges, bls = [], {}
    for k, v in ex["target_branch_lengths"].items():
        p, c = k.split("|"); edges.append((p, c)); bls[(p, c)] = v
    node_seqs = ex["target_node_seqs"]
    node_ids = list(dict.fromkeys([ex["root_id"]] + [n for e in edges for n in e] + list(node_seqs)))
    leaves = [n for n in node_ids if n not in {p for p, _ in edges}]
    return TreeState(node_ids=node_ids, root_id=ex["root_id"], edges=edges,
                     branch_lengths=bls, node_seqs=node_seqs, active_leaves=leaves)


def _qd(a, b):
    try:
        return quartet_distance(a, b)
    except ImportError:
        return float("nan")


def score_empirical(gen: TreeState, target: TreeState) -> dict:
    return {"tree_kl": float("nan"), "split_kl": float("nan"),
            "rf": sequence_matched_rf(gen, target), "quartet": _qd(gen, target),
            "branch_w_all": branch_length_wasserstein(gen, target)["all"],
            "terminal_edit": terminal_edit_distance(gen, target)["mean"]}


def score_simulated(gens: list[TreeState], refs: list[TreeState], seed: int) -> dict:
    rng = random.Random(seed)
    rf, q, bw, te = [], [], [], []
    for g in gens:
        r = rng.choice(refs)
        rf.append(sequence_matched_rf(g, r)); bw.append(branch_length_wasserstein(g, r)["all"])
        te.append(terminal_edit_distance(g, r)["mean"])
        qv = _qd(g, r)
        if qv == qv:
            q.append(qv)
    return {"tree_kl": tree_js(gens, refs), "split_kl": split_js(gens, refs),
            "rf": mean(rf) if rf else float("nan"),
            "quartet": mean(q) if q else float("nan"),
            "branch_w_all": mean(bw) if bw else float("nan"),
            "terminal_edit": mean(te) if te else float("nan")}


def build_methods(args, params, esm):
    methods = [NeutralBD(params["birth"], params["death"]),
               EmpiricalBD(params["birth"], params["death"], model=args.empirical_model)]
    if esm is not None:
        methods.append(PLMPrior(esm.lm_logits, params["birth"], params["death"],
                                params.get("subst_scale", 1.0)))
    if args.checkpoint:
        from benchmarks.methods.treesbm import TreeSBMMethod
        methods.append(TreeSBMMethod(args.checkpoint))
    return methods


class ESM:
    def __init__(self, device, max_len=566):
        import torch
        from transformers import AutoTokenizer, EsmForMaskedLM
        from scripts.eval_single_tree import get_lm_logits, AA_VOCAB
        self._gl, self._t = get_lm_logits, torch
        mid = "facebook/esm2_t6_8M_UR50D"
        self.tok = AutoTokenizer.from_pretrained(mid)
        self.model = EsmForMaskedLM.from_pretrained(mid).to(device).eval()
        self.aa = torch.tensor([self.tok.convert_tokens_to_ids(a) for a in AA_VOCAB])
        self.device, self.max_len = device, max_len

    def lm_logits(self, seq):
        return self._gl(self.tok, self.model, self.aa, [seq], self.max_len, self.device)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-data", default="data/h3n2/test")
    ap.add_argument("--params", default="benchmarks/results/params.json")
    ap.add_argument("--checkpoint", default=None)
    ap.add_argument("--empirical-model", default="JTT")
    ap.add_argument("--N", type=int, nargs="+", default=[16])
    ap.add_argument("--K", type=int, default=100)
    ap.add_argument("--M", type=int, default=50)
    ap.add_argument("--max-roots", type=int, default=100)
    ap.add_argument("--regimes", nargs="+", default=REGIMES)
    ap.add_argument("--no-esm", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="benchmarks/results/results.csv")
    args = ap.parse_args()

    params = json.loads((ROOT / args.params).read_text())
    esm = None if args.no_esm else ESM("cuda" if _cuda() else "cpu")
    methods = build_methods(args, params, esm)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for N in args.N:
            examples = build_examples(ROOT / args.test_data, N, seed=args.seed)[: args.max_roots]
            print(f"N={N}: {len(examples)} roots x {len(methods)} methods")
            for ex in examples:
                target = rebuild_target(ex)
                root_seq, H = ex["root_seq"], ex["H"]
                for method in methods:
                    gens, valids = [], []
                    for k in range(args.K):
                        g = method.generate(root_seq, N, H, seed=args.seed * 10000 + k)
                        vr = V.validate(g.tree, root_seq, N, H)
                        gens.append(g); valids.append(vr)
                    valid_trees = [g.tree for g, vr in zip(gens, valids) if vr["valid"]]
                    base = dict(method=method.name, root_id=ex["root_id"], N=N, H=H,
                                valid=len(valid_trees), sample_seed=args.seed,
                                runtime=sum(g.meta.get("runtime", 0.0) for g in gens))
                    # empirical track (mean over valid samples vs the one true subtree)
                    if valid_trees:
                        es = [score_empirical(t, target) for t in valid_trees]
                        row = {**base, "track": "empirical", "sim_regime": ""}
                        for m in ("rf", "quartet", "branch_w_all", "terminal_edit",
                                  "tree_kl", "split_kl"):
                            vals = [e[m] for e in es if e[m] == e[m]]
                            row[m] = mean(vals) if vals else float("nan")
                        w.writerow(row)
                    # simulated track (per regime)
                    for regime in args.regimes:
                        refs = simulate_reference(root_seq, N, H, regime, args.M,
                                                  params["birth"], params["death"],
                                                  seed=args.seed + hash(ex["root_id"]) % 9973)
                        if valid_trees:
                            ss = score_simulated(valid_trees, refs, seed=args.seed)
                            w.writerow({**base, "track": f"sim_{regime}",
                                        "sim_regime": regime, **ss})
                    f.flush()
    print(f"wrote {out}")


def _cuda():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


if __name__ == "__main__":
    main()
