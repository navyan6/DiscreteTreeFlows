#!/usr/bin/env python3
"""
Sample K new topologies from a trained PhyloVAE model.

Not autoregressive like ARTreeFormer: LVMTree.forward() only *evaluates*
existing tree vectors against a given z, it never exposes an unconditional
"sample a new tree" method (main.py only has train/test/rep modes -- `rep`
computes 2D embeddings of *existing* training-set topologies, not new
samples). The actual generative process, traced from LVMTree itself:
    z ~ N(0, I)                                  (the VAE prior)
    latent_logits = decoder(z)                   [K, (ntips-3)*(ntips-1)]
    cond_probs_mat = cond_prob_mat(latent_logits) [K, ntips-3, 2*ntips-4]
then each of the (ntips-3) attachment decisions is drawn *independently*
from its row of cond_probs_mat (all conditioned on the same z, not on each
other -- the continuous latent is what's supposed to capture the
correlation between decisions), and decoded to a tree via vec2tree(), the
exact inverse of the tree2vec() encoding process_data() used to build the
training vectors in the first place. Verified the edge_mask valid-position
arithmetic matches vec2tree's construction step-for-step (row i has
2i+3 valid positions matching an (i+3)-taxon tree's edge count; total
valid positions ((ntips-3)(ntips-1)/2 summed) exactly equals the decoder's
flat output dim) -- not executed against a real checkpoint, so verify the
first sampled pool looks sane (right leaf count, no crashes) before trusting
it at scale.

COPY THIS FILE into the PhyloVAE repo root before running (needs
`from src.latent_tree_model import VAETree`, same as main.py). Run in the
phylovae conda env.

Usage (from PhyloVAE repo root; base.mode/data.* args just need to match
whatever was used for training, only to load the same architecture):
    python phylovae_sample.py \
        base.mode=train data.dataset=DATASET data.rep_id=1 \
        decoder.num_layers=4 decoder.latent_dim=2 objective.batch_size=10 objective.n_particles=32 \
        --checkpoint results/tde/DATASET/rep_1/.../final.pt \
        --ntips 16 --n-samples 300 --out ../pools/phylovae_N16.nwk
"""

import argparse

import numpy as np
import torch
from omegaconf import OmegaConf

from src.latent_tree_model import VAETree
from src.vector_representation import vec2tree


class _FakeDataset:
    # VAETree.__init__ reads self.dataloader.dataset.wts for an (unused-at-
    # sampling-time) entropy bookkeeping constant; sampling itself never
    # touches the dataloader again, so a dummy avoids needing real data here.
    wts = np.array([1.0])


class _FakeLoader:
    dataset = _FakeDataset()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--ntips", type=int, required=True)
    ap.add_argument("--n-samples", type=int, default=300)
    ap.add_argument("--out", required=True)
    ap.add_argument("--key", default="ema", choices=["ema", "model"])
    args, unknown = ap.parse_known_args()

    cfg_file = OmegaConf.load("config.yaml")
    cfg_cmd = OmegaConf.from_cli(unknown)
    cfg = OmegaConf.merge(cfg_file, cfg_cmd)
    cfg.base.device = "cuda" if torch.cuda.is_available() else "cpu"

    taxa = [str(i) for i in range(args.ntips)]
    model = VAETree(taxa, _FakeLoader(), None, cfg=cfg).to(cfg.base.device)
    ckpt = torch.load(args.checkpoint, map_location=cfg.base.device)
    model.load_state_dict(ckpt[args.key])
    model.eval()

    ltm = model.latent_tree_model
    with torch.no_grad():
        samp_z, _ = ltm.sample_z(args.n_samples)              # [K, latent_dim]
        latent_logits = ltm.decoder(samp_z)                   # [K, (ntips-3)*(ntips-1)]
        cond_probs_mat = ltm.cond_prob_mat(latent_logits)      # [K, ntips-3, 2*ntips-4]

    trees = []
    for k in range(args.n_samples):
        vec = [torch.multinomial(cond_probs_mat[k, i], 1).item() for i in range(args.ntips - 3)]
        trees.append(vec2tree(vec))

    with open(args.out, "w") as f:
        for t in trees:
            f.write(t.write(format=9) + "\n")   # format=9: leaf names only, no branch lengths
    print(f"wrote {len(trees)} sampled topologies -> {args.out}")


if __name__ == "__main__":
    main()
