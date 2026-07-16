#!/usr/bin/env python3
"""
Precompute raw ESM2 embeddings [N, 320] for all groups in data/train/.
this saves compute time later 
"""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import argparse

import torch
from Bio import SeqIO

from src.dataset import parse_newick
from src.treeencoder.plm_embeddings import ESM2Embedder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/train",
                        help="Dir of group_NNN_rooted.nwk / _anc_aa.fasta (e.g. data/h3n2/train)")
    args = parser.parse_args()
    DATA = ROOT / args.data

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    embedder = ESM2Embedder(device=device)

    groups = sorted([
        int(p.stem.split("_")[1])
        for p in DATA.glob("group_*_rooted.nwk")
        if (DATA / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists()
    ])
    print(f"Found {len(groups)} complete groups\n")

    for g in groups:
        out_path = DATA / f"group_{g:03d}_plm.pt"
        if out_path.exists():
            print(f"[{g:03d}] already cached, skipping")
            continue

        nwk   = DATA / f"group_{g:03d}_rooted.nwk"
        fasta = DATA / f"group_{g:03d}_anc_aa.fasta"

        root_id, node_ids, _, _ = parse_newick(str(nwk))

        seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}
        ref_len = len(next(iter(seqs.values())))
        for nid in node_ids:
            if nid not in seqs:
                seqs[nid] = "-" * ref_len

        sequences = [seqs[nid] for nid in node_ids]

        print(f"[{g:03d}] embedding {len(sequences)} sequences ...", end=" ", flush=True)
        plm = embedder.embed_sequences(sequences, batch_size=32)  # [N, 320]
        print(f"done  shape={tuple(plm.shape)}")

        torch.save({"node_ids": node_ids, "plm": plm.cpu()}, out_path)
        print(f"[{g:03d}] saved → {out_path.name}")

    print("\nAll groups done.")


if __name__ == "__main__":
    main()
