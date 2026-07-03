#esm2 embedding extraction
from ast import Load

import torch


class ESM2Embedder:
    #Load ESM2-t6 and extract mean-pooled embeddings.

    def __init__(self, device: str = "cuda" if torch.cuda.is_available() else "cpu"):
        self.device = device

        try:
            from transformers import AutoTokenizer, AutoModel
        except ImportError:
            raise ImportError("Install: pip install transformers")

        model_id = "facebook/esm2_t6_8M_UR50D"
        self.tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
        self.model = AutoModel.from_pretrained(model_id, trust_remote_code=True)

        self.model = self.model.to(device)
        self.model.eval()

        # Freeze weights
        for param in self.model.parameters():
            param.requires_grad = False

    def embed_sequence(self, sequence: str) -> torch.Tensor:
        """
        Returns embedding: [320] mean-pooled representation
        """
        with torch.no_grad():
            tokens = self.tokenizer(
                sequence, return_tensors="pt", padding=False, truncation=False
            ).to(self.device)
            output = self.model(**tokens)
            embeddings = output.last_hidden_state[0, 1:-1, :]  # [seq_len, 320]
            pooled = embeddings.mean(dim=0)  # [320]

        return pooled.detach()

    def embed_sequences(self, sequences: list[str], batch_size: int = 32) -> torch.Tensor:
        """
        Embed a list of sequences in mini-batches for GPU efficiency.
        Returns [N, 320] mean-pooled embeddings.
        """
        all_embeddings: list[torch.Tensor] = []
        for start in range(0, len(sequences), batch_size):
            batch_seqs = sequences[start : start + batch_size]
            with torch.no_grad():
                tokens = self.tokenizer(
                    batch_seqs,
                    return_tensors="pt",
                    padding=True,
                    truncation=False,
                ).to(self.device)
                output = self.model(**tokens)
                # last_hidden_state: [B, L_padded, 320]
                # attention_mask: [B, L_padded]
                hidden = output.last_hidden_state  # [B, L, 320]
                mask = tokens["attention_mask"]    # [B, L]  (1 = real, 0 = pad)

                # Exclude CLS (pos 0) and EOS (last real token); use mask to ignore padding
                # CLS is always position 0; EOS is at position (sum(mask)-1) per sequence
                # Simplest correct approach: zero out CLS and pad, then mean over real positions
                # Positions: [1 .. L-1] are real residues + EOS; mask tells us which are real
                # We exclude position 0 (CLS) and last real position (EOS):
                seq_lens = mask.sum(dim=1)  # [B] total tokens including CLS and EOS
                pooled_list = []
                for i in range(len(batch_seqs)):
                    L = seq_lens[i].item()
                    # residue tokens are positions 1 .. L-2 (exclude CLS at 0 and EOS at L-1)
                    residue_emb = hidden[i, 1:L - 1, :]  # [seq_len, 320]
                    if residue_emb.shape[0] == 0:
                        residue_emb = hidden[i, 1:2, :]  # fallback: CLS only
                    pooled_list.append(residue_emb.mean(dim=0))  # [320]
                all_embeddings.append(torch.stack(pooled_list))  # [B, 320]

        return torch.cat(all_embeddings, dim=0)  # [N, 320]


def get_tree_embeddings(tree_state, device: str = "cuda" if torch.cuda.is_available() else "cpu") -> torch.Tensor:
    embedder = ESM2Embedder(device=device)
    sequences = [tree_state.node_seqs[node_id] for node_id in tree_state.node_ids]
    return embedder.embed_sequences(sequences)
