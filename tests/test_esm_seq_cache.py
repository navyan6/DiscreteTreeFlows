"""Smoke tests for incremental seq-keyed ESM emb / R0 caches."""

from __future__ import annotations

import torch

from scripts.eval_single_tree import (
    _seq_keyed_cached_batch,
    embed_sequences_cached,
    get_lm_logits_cached,
)


class _FakeEmbedder:
    def __init__(self):
        self.calls: list[list[str]] = []

    def embed_sequences(self, sequences: list[str]) -> torch.Tensor:
        self.calls.append(list(sequences))
        # Deterministic [N, 4] stand-in for ESM mean-pool.
        rows = []
        for s in sequences:
            h = float(sum(ord(c) for c in s) % 97)
            rows.append(torch.tensor([h, h + 1.0, h + 2.0, len(s)], dtype=torch.float32))
        return torch.stack(rows, dim=0)


def test_seq_keyed_cache_second_pass_hits():
    cache: dict[str, torch.Tensor] = {}
    stats = {"hits": 0, "misses": 0, "unique": 0}
    calls: list[list[str]] = []

    def compute(miss):
        calls.append(list(miss))
        return torch.stack(
            [torch.tensor([float(len(s)), 1.0], dtype=torch.float32) for s in miss],
            dim=0,
        )

    seqs = ["AAA", "CCC", "AAA"]
    out1 = _seq_keyed_cached_batch(seqs, cache, compute, "cpu", stats=stats)
    assert out1.shape == (3, 2)
    assert len(calls) == 1
    assert set(calls[0]) == {"AAA", "CCC"}
    # Within-batch duplicates are deduped before the forward; first pass is all miss.
    assert stats["misses"] == 2
    assert stats["hits"] == 0
    assert stats["unique"] == 2
    assert torch.equal(out1[0], out1[2])

    out2 = _seq_keyed_cached_batch(seqs, cache, compute, "cpu", stats=stats)
    assert torch.equal(out1, out2)
    assert len(calls) == 1  # no new compute — second embed of same seqs hits cache
    assert stats["hits"] == 3
    assert stats["misses"] == 2


def test_embed_sequences_cached_miss_then_hit():
    emb = _FakeEmbedder()
    cache: dict[str, torch.Tensor] = {}
    stats = {"hits": 0, "misses": 0, "unique": 0}

    a = embed_sequences_cached(emb, ["ACDE", "FGHI"], cache, "cpu", stats=stats)
    assert a.shape == (2, 4)
    assert len(emb.calls) == 1
    assert emb.calls[0] == ["ACDE", "FGHI"]
    assert stats["misses"] == 2 and stats["hits"] == 0

    b = embed_sequences_cached(emb, ["FGHI", "ACDE", "ACDE"], cache, "cpu", stats=stats)
    assert len(emb.calls) == 1  # full hit
    assert stats["hits"] == 3
    assert torch.allclose(b[0], a[1])
    assert torch.allclose(b[1], a[0])
    assert torch.allclose(b[2], a[0])


def test_get_lm_logits_cached_uses_exact_recompute_on_miss(monkeypatch):
    cache: dict[str, torch.Tensor] = {}
    stats = {"hits": 0, "misses": 0, "unique": 0}
    n_calls = {"n": 0}

    def fake_get_lm_logits(tokenizer, esm_model, aa_token_ids, sequences, max_seq_len, device):
        n_calls["n"] += 1
        # Distinct rows per sequence length so we can verify identity.
        return torch.stack(
            [
                torch.full((max_seq_len, 20), float(len(s)), dtype=torch.float32)
                for s in sequences
            ],
            dim=0,
        )

    monkeypatch.setattr(
        "scripts.eval_single_tree.get_lm_logits", fake_get_lm_logits
    )

    out1 = get_lm_logits_cached(
        None, None, None, ["AAAA", "CCCC"], 8, "cpu", cache, stats=stats
    )
    assert out1.shape == (2, 8, 20)
    assert n_calls["n"] == 1
    assert stats["misses"] == 2

    out2 = get_lm_logits_cached(
        None, None, None, ["CCCC", "AAAA"], 8, "cpu", cache, stats=stats
    )
    assert n_calls["n"] == 1
    assert stats["hits"] == 2
    assert torch.equal(out2[0], out1[1])
    assert torch.equal(out2[1], out1[0])

    # New sequence forces exact recompute for that miss only.
    out3 = get_lm_logits_cached(
        None, None, None, ["AAAA", "GGGG"], 8, "cpu", cache, stats=stats
    )
    assert n_calls["n"] == 2
    assert stats["misses"] == 3
    assert torch.equal(out3[0], out1[0])
    assert float(out3[1, 0, 0]) == 4.0
