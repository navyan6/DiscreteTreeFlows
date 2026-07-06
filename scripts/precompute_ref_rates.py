#!/usr/bin/env python3
"""
Precompute ESM-2 reference mutation log-rates [N, 566, 20] for all groups.
use in forward pass
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
import torch.nn.functional as F
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForMaskedLM

from src.dataset import parse_newick

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY" 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/train")
    parser.add_argument("--max-seq-len", type=int, default=566)
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args()

    data_dir = ROOT / args.data
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    model_id = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = EsmForMaskedLM.from_pretrained(model_id).to(device)
    model.eval()
    for p in model.parameters():
        p.requires_grad = False

    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB], dtype=torch.long
    ) 

    groups = sorted([
        int(p.stem.split("_")[1])
        for p in data_dir.glob("group_*_rooted.nwk")
        if (data_dir / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists()
    ])
    print(f"Found {len(groups)} complete groups\n")

    for g in groups:
        out_path = data_dir / f"group_{g:03d}_ref_rates.pt"
        if out_path.exists():
            print(f"[{g:03d}] already cached, skipping")
            continue

        root_id, node_ids, _, _ = parse_newick(
            str(data_dir / f"group_{g:03d}_rooted.nwk")
        )
        seqs = {
            rec.id: str(rec.seq)
            for rec in SeqIO.parse(data_dir / f"group_{g:03d}_anc_aa.fasta", "fasta")
        }
        ref_len = len(next(iter(seqs.values())))
        for nid in node_ids:
            if nid not in seqs:
                seqs[nid] = "-" * ref_len
        sequences = [seqs[nid] for nid in node_ids]
        N, L = len(sequences), args.max_seq_len

        print(f"[{g:03d}] {N} sequences ...", end=" ", flush=True)
        log_mut_rates = torch.zeros(N, L, 20, dtype=torch.float32)

        for start in range(0, N, args.batch_size):
            batch_seqs = sequences[start : start + args.batch_size]
            with torch.no_grad():
                tokens = tokenizer(
                    batch_seqs, return_tensors="pt", padding=True, truncation=False
                ).to(device)
                logits = model(**tokens).logits  
            seq_lens = tokens["attention_mask"].sum(dim=1)  

            for i in range(len(batch_seqs)):
                actual_L = int(seq_lens[i].item()) - 2 
                aa_logits = logits[i, 1 : actual_L + 1, :][:, aa_token_ids] 
                log_probs = F.log_softmax(aa_logits, dim=-1)
                clip = min(actual_L, L)
                log_mut_rates[start + i, :clip, :] = log_probs[:clip].cpu()

        torch.save({"node_ids": node_ids, "log_mut_rates": log_mut_rates}, out_path)
        print(f"done  shape={tuple(log_mut_rates.shape)}")

    print("\nAll groups done.")


if __name__ == "__main__":
    main()
