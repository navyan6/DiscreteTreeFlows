#!/usr/bin/env python3
"""
Generate NeutralBD (and optional pLM-prior) leaf FASTAs for a VOC case root,
then score with eval_voc_threat_recovery.py.

Used as baseline comparison for TreeSBM VOC threat recovery.

Examples:
  python scripts/gen_voc_baseline_leaves.py \\
    --case-dir results/voc_threat_panel/cases/Gamma_g003_test \\
    --voc Gamma --methods neutral_bd --n-leaves 250 --out-dir \\
    results/voc_threat_panel/cases/Gamma_g003_test/baselines
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from voc_threat_lib import find_root_seq, load_fasta_seqs  # noqa: E402


def _write_leaves_fasta(path: Path, root_seq: str, leaf_seqs: dict[str, str]) -> None:
    recs = [SeqRecord(Seq(root_seq), id="root|root", description="baseline_root")]
    for lid, seq in leaf_seqs.items():
        recs.append(SeqRecord(Seq(seq), id=lid, description="baseline_leaf"))
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(recs, path, "fasta")


def _leaves_from_treestate(ts) -> dict[str, str]:
    leaves = {}
    for nid in ts.node_ids:
        if nid == ts.root_id:
            continue
        if not ts.is_leaf(nid):
            continue
        seq = ts.node_seqs.get(nid, "")
        if seq:
            leaves[nid] = seq
    if not leaves:
        raise RuntimeError("no leaf sequences on TreeState")
    return leaves


def gen_neutral(root_seq: str, n_leaves: int, birth: float, death: float, H: float, seed: int):
    from benchmarks.methods.bd_methods import NeutralBD

    m = NeutralBD(birth, death)
    gt = m.generate(root_seq, N=n_leaves, H=H, seed=seed)
    return _leaves_from_treestate(gt.tree)


def gen_plm(root_seq: str, n_leaves: int, birth: float, death: float, H: float, seed: int, device: str):
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM
    from benchmarks.methods.plm_prior import PLMPrior
    from voc_threat_lib import AA_VOCAB

    tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    model = EsmForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device)
    model.eval()

    def lm_logits(seq: str):
        enc = tok(seq, return_tensors="pt", add_special_tokens=True)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            out = model(**enc).logits[0, 1:-1, :]
        idx = [tok.convert_tokens_to_ids(a) for a in AA_VOCAB]
        return out[:, idx].float().cpu()

    m = PLMPrior(lm_logits, birth, death)
    gt = m.generate(root_seq, N=n_leaves, H=H, seed=seed)
    return _leaves_from_treestate(gt.tree)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--voc", required=True)
    ap.add_argument("--methods", nargs="+", default=["neutral_bd"])
    ap.add_argument("--n-leaves", type=int, default=250)
    ap.add_argument("--birth", type=float, default=50.0)
    ap.add_argument("--death", type=float, default=0.0)
    ap.add_argument("--H", type=float, default=1.0, help="BD horizon")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--skip-eval", action="store_true")
    args = ap.parse_args()

    case = args.case_dir
    obs = case / "observed_anc_aa.fasta"
    if not obs.exists():
        raise SystemExit(f"missing {obs}")
    seqs = load_fasta_seqs(obs)
    nwk = case / "observed.nwk"
    rid, root = find_root_seq(seqs, nwk if nwk.exists() else None)
    out_dir = args.out_dir or (case / "baselines")
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {"root_id": rid, "root_len": len(root), "methods": {}}

    for method in args.methods:
        print(f"generating {method} n={args.n_leaves} …", flush=True)
        if method == "neutral_bd":
            leaves = gen_neutral(root, args.n_leaves, args.birth, args.death, args.H, args.seed)
        elif method == "plm_prior":
            leaves = gen_plm(
                root, args.n_leaves, args.birth, args.death, args.H, args.seed, args.device
            )
        else:
            raise SystemExit(f"unknown method {method}")

        fa = out_dir / f"{method}.fasta"
        _write_leaves_fasta(fa, root, leaves)
        print(f"  wrote {fa} ({len(leaves)} leaves)")

        eval_out = out_dir / f"voc_eval_{args.voc}_{method}.json"
        if not args.skip_eval:
            import subprocess

            cmd = [
                sys.executable,
                str(ROOT / "scripts" / "eval_voc_threat_recovery.py"),
                "--voc",
                args.voc,
                "--gen-fasta",
                str(fa),
                "--obs-fasta",
                str(obs),
                "--topk",
                "10",
                "--evescape",
                str(ROOT / "data" / "covid" / "evescape_spike_rbd.pt"),
                "--out",
                str(eval_out),
            ]
            if nwk.exists():
                cmd.extend(["--obs-nwk", str(nwk)])
            subprocess.run(cmd, check=True)
            blob = json.loads(eval_out.read_text())
            g = blob["gen_vs_score_root"]
            summary["methods"][method] = {
                "exact_recall": g["voc_exact_mut_recall"],
                "acquired_recall": g["voc_acquired_mut_recall"],
                "bundle_any": g["voc_bundle_any"],
                "topk_union": blob["topk"]["union_sig"],
                "eval": str(eval_out),
            }

    (out_dir / "baseline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
