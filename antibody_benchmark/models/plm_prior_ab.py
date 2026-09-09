"""pLM mutation prior on observed Ab topology+BL (ESM-2 per-branch edits)."""

from __future__ import annotations

from antibody_benchmark.models.base import EvolutionModel


class PLMPriorAbModel(EvolutionModel):
    """ESM-2 mutation prior only — sequences evolve on the *given* branch length.

    Matches Track C fairness (observed topology+BL); not a free BD topology draw.
    """

    name = "plm_prior"
    alphabet = "aa"
    AA = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(
        self,
        device: str | None = None,
        esm_id: str = "facebook/esm2_t6_8M_UR50D",
        subs_per_site_scale: float = 1.0,
    ):
        self.esm_id = esm_id
        self.scale = float(subs_per_site_scale)
        if device is None:
            try:
                import torch

                device = "cuda" if torch.cuda.is_available() else "cpu"
            except Exception:
                device = "cpu"
        self.device = device
        self._torch = None
        self._tokenizer = None
        self._model = None
        self._aa_to_idx = {a: i for i, a in enumerate(self.AA)}
        self._load()

    def _load(self) -> None:
        try:
            import torch
            from transformers import AutoTokenizer, EsmForMaskedLM
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: pLM prior needs torch+transformers. Error: {e}"
            ) from e
        self._torch = torch
        self._tokenizer = AutoTokenizer.from_pretrained(self.esm_id)
        self._model = EsmForMaskedLM.from_pretrained(self.esm_id).to(self.device).eval()
        for p in self._model.parameters():
            p.requires_grad = False

    def _lm_probs(self, seq: str):
        torch = self._torch
        assert self._tokenizer is not None and self._model is not None
        toks = self._tokenizer(seq, return_tensors="pt", add_special_tokens=True)
        toks = {k: v.to(self.device) for k, v in toks.items()}
        with torch.no_grad():
            logits = self._model(**toks).logits[0]  # [L+2, V]
        # strip BOS/EOS
        logits = logits[1 : 1 + len(seq)]
        aa_ids = [self._tokenizer.convert_tokens_to_ids(a) for a in self.AA]
        aa_logits = logits[:, aa_ids]  # [L, 20]
        return torch.softmax(aa_logits, dim=-1)

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        import numpy as np

        torch = self._torch
        assert torch is not None
        rng = np.random.default_rng(int(rng_seed) % (2**31 - 1))
        seq = list(parent_seq)
        L = len(seq)
        if L == 0:
            return parent_seq
        bl = max(float(branch_length), 0.0)
        n_sub = int(rng.poisson(bl * L * self.scale))
        if n_sub <= 0:
            return parent_seq
        probs = self._lm_probs(parent_seq)  # [L,20]
        positions = rng.choice(L, size=min(n_sub, L), replace=False)
        for pos in positions:
            pv = probs[int(pos)].detach().cpu().numpy().astype(float).copy()
            cur = seq[int(pos)]
            ci = self._aa_to_idx.get(cur, -1)
            if 0 <= ci < 20:
                pv[ci] = 0.0
            s = float(pv.sum())
            if s <= 0:
                continue
            pv /= s
            seq[int(pos)] = self.AA[int(rng.choice(20, p=pv))]
        return "".join(seq)
