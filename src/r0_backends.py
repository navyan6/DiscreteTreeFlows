"""
Multi-pLM / substitution R0 backends for TreeSBM reference mutation priors.

Paper Table 7 / appendix D.1 compare mutation priors as frozen R0 sources:
  JTT/WAG/LG substitution | ESM-2-650M ± fitness | ESM-C ± fitness | ProGen2 | …
TreeSBM keeps ``log R_θ = log R0(+tilt) + c_θ``; these wrappers only swap how
``log R0`` is produced. They are adapters, not new networks.

Fitness tilting stays in ``src/bridge/fitness_tilt.py`` (orthogonal to backend).
Poisson branching λ(x) stays in ``src/reference_process.py``.
"""

from __future__ import annotations

import os
import sys
import types
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Sequence

import torch
import torch.nn.functional as F

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}

# Canonical names used by CLI / cache tags / paper rows.
BACKEND_ESM2 = "esm2"
BACKEND_ESM2_650M = "esm2_650m"
BACKEND_ESMC = "esmc"
BACKEND_PROGEN2 = "progen2"
BACKEND_EVO2 = "evo2"
BACKEND_JTT = "jtt"
BACKEND_WAG = "wag"
BACKEND_LG = "lg"
BACKEND_NEUTRAL = "neutral"
BACKEND_THRIFTY_AA = "thrifty_aa"
BACKEND_CODON_GY94 = "codon_gy94"

SUBSTITUTION_BACKENDS = frozenset(
    {BACKEND_JTT, BACKEND_WAG, BACKEND_LG, BACKEND_NEUTRAL}
)
STUB_BACKENDS = frozenset({BACKEND_EVO2})

DEFAULT_ESM2_MODEL = "facebook/esm2_t6_8M_UR50D"
DEFAULT_ESM2_650M_MODEL = "facebook/esm2_t33_650M_UR50D"
# Official biohub/ESMC-300M needs transformers with native `esmc` (unreleased as of
# 5.15). Synthyra/ESMplusplus_small is the working HF AutoModelForMaskedLM port of
# ESMC-300M (trust_remote_code, attn sdpa/eager — no flash_attn).
DEFAULT_ESMC_MODEL = "Synthyra/ESMplusplus_small"
ESMC_HF_CANDIDATES = (
    "biohub/ESMC-300M",
    "Synthyra/ESMplusplus_small",
)
# Official enijkamp/progen2 clone + GCS checkpoint (not hugohrban HF AutoModel).
DEFAULT_PROGEN2_HOME = os.environ.get(
    "PROGEN2_HOME",
    "/vast/projects/pranam/lab/nnori/progen2",
)
DEFAULT_PROGEN2_MODEL = "progen2-small"
# progen2-small / medium absolute position limit (config.n_positions).
# With leading terminal ``1``, one forward pass scores at most this many AA − 1.
PROGEN2_MAX_CTX = 1024
# Token ids from enijkamp/progen2 tokenizer.json (A..Z block starts at 5).
PROGEN2_AA_TOKEN_IDS = {
    "A": 5,
    "C": 7,
    "D": 8,
    "E": 9,
    "F": 10,
    "G": 11,
    "H": 12,
    "I": 13,
    "K": 14,
    "L": 15,
    "M": 16,
    "N": 17,
    "P": 19,
    "Q": 20,
    "R": 21,
    "S": 22,
    "T": 23,
    "V": 25,
    "W": 26,
    "Y": 28,
}

# Cache file tag → ``group_{g:03d}_ref_rates{tag}.pt``
# Empty tag preserves legacy ``group_*_ref_rates.pt`` (current ESM-2-8M caches).
BACKEND_CACHE_TAG = {
    BACKEND_ESM2: "",
    BACKEND_ESM2_650M: "_esm2_650m",
    BACKEND_ESMC: "_esmc",
    BACKEND_PROGEN2: "_progen2",
    BACKEND_EVO2: "_evo2",
    BACKEND_JTT: "_jtt",
    BACKEND_WAG: "_wag",
    BACKEND_LG: "_lg",
    BACKEND_NEUTRAL: "_neutral",
    # Antibody OAS Recipe B only. Viral/default caches stay untagged ESM.
    BACKEND_THRIFTY_AA: "_thrifty",
    # Viral codon Q0 (GY94→AA). Does not replace untagged ESM caches.
    BACKEND_CODON_GY94: "_codon_gy94",
}


def normalize_backend_name(name: str) -> str:
    key = name.strip().lower().replace("-", "_")
    aliases = {
        "esm2_8m": BACKEND_ESM2,
        "esm_2": BACKEND_ESM2,
        "esm2_t6": BACKEND_ESM2,
        "esm2_650": BACKEND_ESM2_650M,
        "esm2_t33": BACKEND_ESM2_650M,
        "esm_c": BACKEND_ESMC,
        "esm_cambrian": BACKEND_ESMC,
        "progen": BACKEND_PROGEN2,
        "progen_2": BACKEND_PROGEN2,
        "evo_2": BACKEND_EVO2,
        "substitution": BACKEND_JTT,
        "substitution_only": BACKEND_JTT,
        # Antibody Recipe B only — do not alias a short "thrifty" name that
        # viral train/eval might pass by accident.
        "thrifty_q0": BACKEND_THRIFTY_AA,
        "shm_thrifty": BACKEND_THRIFTY_AA,
        "codon_q0": BACKEND_CODON_GY94,
        "gy94": BACKEND_CODON_GY94,
        "mg94": BACKEND_CODON_GY94,
        "codon_gy94": BACKEND_CODON_GY94,
    }
    return aliases.get(key, key)


def cache_tag_for_backend(backend: str, override: Optional[str] = None) -> str:
    """Return filename tag ('' or '_esmc', …). ``override`` wins if not None."""
    if override is not None:
        tag = override.strip()
        if tag and not tag.startswith("_"):
            tag = "_" + tag
        return tag
    name = normalize_backend_name(backend)
    if name not in BACKEND_CACHE_TAG:
        raise ValueError(
            f"Unknown R0 backend {backend!r}. "
            f"Known: {sorted(BACKEND_CACHE_TAG)}"
        )
    return BACKEND_CACHE_TAG[name]


def ref_rates_filename(group: int, tag: str = "") -> str:
    return f"group_{group:03d}_ref_rates{tag}.pt"


class R0Backend(ABC):
    """Frozen mutation prior → sitewise log-probs over 20 AAs."""

    name: str = "base"

    @abstractmethod
    def log_mutation_rates(
        self,
        sequences: Sequence[str],
        max_seq_len: int,
        device: Optional[torch.device] = None,
    ) -> torch.Tensor:
        """
        Args:
            sequences: AA strings (may contain gaps '-').
            max_seq_len: truncate / pad length L.
            device: optional torch device for compute; output may stay on CPU.

        Returns:
            log_R0: float tensor [N, L, 20] (typically log-softmax over AA).
        """

    def close(self) -> None:
        """Optional resource cleanup."""
        return None


class StubR0Backend(R0Backend):
    """Explicit stub for heavy / unavailable priors (ProGen2, Evo2)."""

    def __init__(self, name: str, install_hint: str):
        self.name = name
        self.install_hint = install_hint

    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        raise NotImplementedError(
            f"R0 backend {self.name!r} is stubbed (too heavy / not wired). "
            f"{self.install_hint}"
        )


class ESM2R0Backend(R0Backend):
    """HuggingFace ESM-2 masked LM → per-position AA log-probs (legacy default)."""

    def __init__(
        self,
        model_id: str = DEFAULT_ESM2_MODEL,
        device: Optional[str] = None,
        name: str = BACKEND_ESM2,
    ):
        from transformers import AutoTokenizer, EsmForMaskedLM

        self.name = name
        self.model_id = model_id
        self.device = torch.device(
            device or ("cuda" if torch.cuda.is_available() else "cpu")
        )
        self.tokenizer = AutoTokenizer.from_pretrained(model_id)
        self.model = EsmForMaskedLM.from_pretrained(model_id).to(self.device)
        self.model.eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.aa_token_ids = torch.tensor(
            [self.tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB],
            dtype=torch.long,
            device=self.device,
        )

    @torch.no_grad()
    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        dev = torch.device(device) if device is not None else self.device
        if dev != self.device:
            # Move model if caller asks for a different device (rare).
            self.model.to(dev)
            self.aa_token_ids = self.aa_token_ids.to(dev)
            self.device = dev

        N, L = len(sequences), max_seq_len
        log_rates = torch.zeros(N, L, 20, dtype=torch.float32)
        # Batch in one shot when N is small; callers may chunk externally.
        tokens = self.tokenizer(
            list(sequences), return_tensors="pt", padding=True, truncation=False
        ).to(self.device)
        logits = self.model(**tokens).logits
        seq_lens = tokens["attention_mask"].sum(dim=1)
        for i in range(N):
            actual_L = int(seq_lens[i].item()) - 2  # strip BOS/EOS
            if actual_L <= 0:
                continue
            aa_logits = logits[i, 1 : actual_L + 1, :][:, self.aa_token_ids]
            log_probs = F.log_softmax(aa_logits, dim=-1)
            clip = min(actual_L, L)
            log_rates[i, :clip, :] = log_probs[:clip].cpu()
        return log_rates


def _load_hf_mlm(model_id: str, device: torch.device):
    """Load a HF masked LM with sdpa/eager attention (avoid flash_attn)."""
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    last_err: Optional[BaseException] = None
    model = None
    for attn in ("sdpa", "eager", None):
        kwargs = {"trust_remote_code": True}
        if attn is not None:
            kwargs["attn_implementation"] = attn
        try:
            model = AutoModelForMaskedLM.from_pretrained(model_id, **kwargs)
            break
        except (TypeError, ValueError, OSError) as e:
            last_err = e
            continue
    if model is None:
        raise RuntimeError(
            f"Failed to load HF MLM {model_id!r} (tried sdpa/eager). Last error: {last_err}"
        )
    model = model.to(device)
    model.eval()
    for p in model.parameters():
        p.requires_grad = False
    return tokenizer, model


class ESMCR0Backend(R0Backend):
    """
    ESM-C via HuggingFace ``AutoModelForMaskedLM`` (no fair-esm, no PYTHONPATH hack).

    Default: ``Synthyra/ESMplusplus_small`` — HF port of ESMC-300M that loads on
    transformers 5.x with ``trust_remote_code`` and ``attn_implementation=sdpa``.
    Also accepts ``biohub/ESMC-300M`` when the installed transformers has native
    ``esmc`` support (PR landing; not in 5.13–5.15 yet).
    """

    def __init__(
        self,
        model_id: str = DEFAULT_ESMC_MODEL,
        device: Optional[str] = None,
    ):
        self.name = BACKEND_ESMC
        self.device = torch.device(
            device or ("cuda" if torch.cuda.is_available() else "cpu")
        )
        candidates = [model_id]
        # Only fall back across HF IDs on architecture/registry errors — not I/O.
        arch_fallback = [
            c for c in ESMC_HF_CANDIDATES if c != model_id
        ]

        last_err: Optional[BaseException] = None
        self.tokenizer = None
        self.model = None
        self.model_id = model_id
        try_list = [model_id]
        for mid in try_list:
            try:
                self.tokenizer, self.model = _load_hf_mlm(mid, self.device)
                self.model_id = mid
                break
            except Exception as e:
                last_err = e
                msg = str(e).lower()
                if any(
                    k in msg
                    for k in (
                        "does not recognize this architecture",
                        "model type `esmc`",
                        "unrecognized configuration",
                        "not a valid model identifier",
                    )
                ):
                    for alt in arch_fallback:
                        if alt not in try_list:
                            try_list.append(alt)
                continue
        if self.model is None or self.tokenizer is None:
            raise RuntimeError(
                "ESM-C HF backend failed. Tried: "
                + ", ".join(try_list)
                + f". Last error: {last_err}. "
                "Default id is Synthyra/ESMplusplus_small (HF port of ESMC-300M). "
                "biohub/ESMC-300M needs transformers with native esmc. "
                "Do NOT use fair-esm or PYTHONPATH=python_esmc."
            ) from last_err

        self.aa_token_ids = torch.tensor(
            [self.tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB],
            dtype=torch.long,
            device=self.device,
        )
        if (self.aa_token_ids < 0).any():
            # Some tokenizers return unk for single letters; encode one-by-one.
            ids = []
            for aa in AA_VOCAB:
                enc = self.tokenizer.encode(aa, add_special_tokens=False)
                if not enc:
                    raise RuntimeError(f"ESM-C tokenizer cannot map amino acid {aa!r}")
                ids.append(int(enc[0]))
            self.aa_token_ids = torch.tensor(ids, dtype=torch.long, device=self.device)

    @torch.no_grad()
    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        dev = torch.device(device) if device is not None else self.device
        if dev != self.device:
            self.model.to(dev)
            self.aa_token_ids = self.aa_token_ids.to(dev)
            self.device = dev

        N, L = len(sequences), max_seq_len
        log_rates = torch.zeros(N, L, 20, dtype=torch.float32)
        clean_seqs = [seq.replace("-", "X")[:L] for seq in sequences]
        tokens = self.tokenizer(
            clean_seqs, return_tensors="pt", padding=True, truncation=False
        ).to(self.device)
        logits = self.model(**tokens).logits  # [N, T, V]
        seq_lens = tokens["attention_mask"].sum(dim=1)
        for i in range(N):
            actual_L = int(seq_lens[i].item()) - 2  # strip BOS/EOS
            if actual_L <= 0:
                continue
            aa_logits = logits[i, 1 : actual_L + 1, :][:, self.aa_token_ids]
            log_probs = F.log_softmax(aa_logits.float(), dim=-1)
            clip = min(actual_L, L, len(clean_seqs[i]))
            log_rates[i, :clip, :] = log_probs[:clip].cpu()
        return log_rates


def _resolve_progen2_paths(
    model_id: Optional[str] = None,
    progen2_home: Optional[str] = None,
) -> tuple[Path, Path, Path]:
    """Return (repo_root, ckpt_dir, tokenizer_json)."""
    home = Path(
        progen2_home
        or os.environ.get("PROGEN2_HOME", DEFAULT_PROGEN2_HOME)
    ).expanduser()
    mid = (model_id or DEFAULT_PROGEN2_MODEL).strip()
    ckpt = Path(mid).expanduser()
    if not ckpt.is_dir():
        ckpt = home / "checkpoints" / mid
    tok = home / "tokenizer.json"
    if not tok.is_file():
        # Some layouts keep tokenizer next to modeling code only.
        alt = home / "progen" / "tokenizer.json"
        tok = alt if alt.is_file() else tok
    return home, ckpt, tok


def _patch_progen_get_head_mask(ProGenModel) -> None:
    """transformers≥5 removed PreTrainedModel.get_head_mask; ProGen still calls it."""
    if hasattr(ProGenModel, "get_head_mask"):
        return

    def get_head_mask(self, head_mask, num_hidden_layers, is_attention_chunked=False):
        if head_mask is not None:
            raise ValueError(
                "ProGen2 head_mask is unsupported under transformers>=5; pass None."
            )
        return [None] * int(num_hidden_layers)

    ProGenModel.get_head_mask = get_head_mask  # type: ignore[attr-defined]


def _ensure_progen_transformers_compat() -> None:
    """
    enijkamp/progen2 still imports transformers.utils.model_parallel_utils and
    calls get_head_mask — both removed in transformers≥5.
    """
    import transformers.utils as tf_utils

    if not hasattr(tf_utils, "model_parallel_utils"):
        stub = types.ModuleType("transformers.utils.model_parallel_utils")

        def get_device_map(*_a, **_k):
            return {}

        def assert_device_map(*_a, **_k):
            return None

        stub.get_device_map = get_device_map  # type: ignore[attr-defined]
        stub.assert_device_map = assert_device_map  # type: ignore[attr-defined]
        sys.modules["transformers.utils.model_parallel_utils"] = stub
        tf_utils.model_parallel_utils = stub  # type: ignore[attr-defined]


def _patch_progen_v5_attrs(*model_classes) -> None:
    """Add attrs removed/renamed in transformers≥5 that ProGen still needs."""
    for cls in model_classes:
        if cls is None:
            continue
        _patch_progen_get_head_mask(cls)
        if not hasattr(cls, "all_tied_weights_keys"):

            def _all_tied_weights_keys(self):
                return getattr(self, "_tied_weights_keys", None) or {}

            cls.all_tied_weights_keys = property(_all_tied_weights_keys)  # type: ignore[attr-defined]


def _import_progen_classes(repo_root: Path):
    """Import ProGen classes from an enijkamp/progen2 checkout."""
    root = str(repo_root.resolve())
    if root not in sys.path:
        sys.path.insert(0, root)
    _ensure_progen_transformers_compat()
    try:
        from progen.configuration_progen import ProGenConfig
        from progen.modeling_progen import ProGenForCausalLM, ProGenModel
    except ImportError:
        try:
            from models.progen.configuration_progen import ProGenConfig  # type: ignore
            from models.progen.modeling_progen import (  # type: ignore
                ProGenForCausalLM,
                ProGenModel,
            )
        except ImportError as e:
            raise ImportError(
                f"Could not import ProGenForCausalLM from {root}. "
                "Clone https://github.com/enijkamp/progen2.git to "
                f"{DEFAULT_PROGEN2_HOME} (or set PROGEN2_HOME)."
            ) from e
    _patch_progen_v5_attrs(ProGenModel, ProGenForCausalLM)
    return ProGenConfig, ProGenForCausalLM


def _load_progen_checkpoint(
    ProGenConfig,
    ProGenForCausalLM,
    ckpt: Path,
    device: torch.device,
    fp16: bool,
):
    """
    Manual load — avoids transformers≥5 ``from_pretrained`` breakage on the
    2022-era ProGen modeling code (``all_tied_weights_keys``, etc.).
    """
    config = ProGenConfig.from_pretrained(str(ckpt))
    model = ProGenForCausalLM(config)
    weight_path = ckpt / "pytorch_model.bin"
    if not weight_path.is_file():
        # Rare alternate layout after untar.
        cands = list(ckpt.glob("**/pytorch_model.bin"))
        if not cands:
            raise FileNotFoundError(f"No pytorch_model.bin under {ckpt}")
        weight_path = cands[0]
    state = torch.load(str(weight_path), map_location="cpu")
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    incompatible = model.load_state_dict(state, strict=False)
    # Prefer a quiet load; surface only if nothing matched.
    missing = getattr(incompatible, "missing_keys", [])
    if missing and len(missing) > 50:
        raise RuntimeError(
            f"ProGen2 state_dict mismatch at {weight_path}: "
            f"{len(missing)} missing keys (e.g. {missing[:3]})"
        )
    model = model.to(device)
    if fp16 and device.type == "cuda":
        model = model.half()
    return model


class ProGen2R0Backend(R0Backend):
    """
    Official enijkamp/progen2 causal LM → sitewise AA log-probs (Table 7).

    Loads ``ProGenForCausalLM`` + ``tokenizer.json`` from ``PROGEN2_HOME``
    (default ``/vast/projects/pranam/lab/nnori/progen2``) and the GCS checkpoint
    ``checkpoints/progen2-small`` (prefer small unless overridden). Scoring follows
    ``likelihood.py``: wrap with terminal ``1``, teacher-forced next-token logits
    over the 20 standard AAs. Left-context only — documented vs ESM masked R0.

    Long sequences (Spike L≈1273 > ``PROGEN2_MAX_CTX``=1024): **chunked** causal
    scoring with sliding left-context windows — never silent truncate. Each forward
    pass is ≤1024 tokens (``1`` + ≤1023 AA); later windows keep max left context and
    only write newly covered sites. Prefer ``progen2-small``.

    Shims transformers≥5 removals (``get_head_mask``, ``model_parallel_utils``)
    and loads weights via ``load_state_dict`` (not broken HF AutoModel path).
    """

    def __init__(
        self,
        model_id: str = DEFAULT_PROGEN2_MODEL,
        device: Optional[str] = None,
        progen2_home: Optional[str] = None,
        fp16: Optional[bool] = None,
        max_ctx: int = PROGEN2_MAX_CTX,
    ):
        from tokenizers import Tokenizer

        self.name = BACKEND_PROGEN2
        self.device = torch.device(
            device or ("cuda" if torch.cuda.is_available() else "cpu")
        )
        home, ckpt, tok_path = _resolve_progen2_paths(model_id, progen2_home)
        self.model_id = str(ckpt)
        self.progen2_home = str(home)

        if not ckpt.is_dir():
            raise RuntimeError(
                f"ProGen2 checkpoint not found at {ckpt}. "
                f"Clone enijkamp/progen2 to {home} and download e.g.\n"
                "  wget -P checkpoints/progen2-small "
                "https://storage.googleapis.com/sfr-progen-research/checkpoints/"
                "progen2-small.tar.gz && tar -xvf "
                "checkpoints/progen2-small/progen2-small.tar.gz "
                "-C checkpoints/progen2-small/"
            )
        if not tok_path.is_file():
            raise RuntimeError(
                f"ProGen2 tokenizer.json missing at {tok_path} "
                f"(expected under {home})."
            )

        ProGenConfig, ProGenForCausalLM = _import_progen_classes(home)
        use_fp16 = (
            bool(fp16) if fp16 is not None else (self.device.type == "cuda")
        )
        try:
            self.model = _load_progen_checkpoint(
                ProGenConfig, ProGenForCausalLM, ckpt, self.device, use_fp16
            )
        except Exception as e:
            raise RuntimeError(
                f"Failed to load ProGen2 checkpoint from {ckpt}: {e}"
            ) from e

        self._fp16 = use_fp16
        # Prefer checkpoint config when present; fall back to documented 1024.
        cfg_ctx = getattr(getattr(self.model, "config", None), "n_positions", None)
        self.max_ctx = int(cfg_ctx) if cfg_ctx else int(max_ctx)
        self.model.eval()
        for p in self.model.parameters():
            p.requires_grad = False

        with open(tok_path, "r") as f:
            self.tokenizer = Tokenizer.from_str(f.read())
        self._aa_token_ids = torch.tensor(
            [PROGEN2_AA_TOKEN_IDS[aa] for aa in AA_VOCAB],
            dtype=torch.long,
            device=self.device,
        )

    def _score_window(self, chunk: str, aa_ids: torch.Tensor) -> torch.Tensor:
        """Teacher-forced AA log-probs for one ≤(max_ctx−1) window. Returns [len(chunk), 20]."""
        context = "1" + chunk
        ids = self.tokenizer.encode(context).ids
        if len(ids) > self.max_ctx:
            # Defensive: tokenizer should be 1 token/char; never exceed n_positions.
            ids = ids[: self.max_ctx]
            chunk = chunk[: len(ids) - 1]
        if not chunk:
            return torch.zeros(0, 20, dtype=torch.float32)
        target = torch.tensor(ids, dtype=torch.long, device=self.device)
        if target.dim() == 1:
            target = target.unsqueeze(0)
        with torch.cuda.amp.autocast(
            enabled=self._fp16 and self.device.type == "cuda"
        ):
            logits = self.model(target).logits
        if logits.dim() == 3:
            logits = logits[0]
        # logits[t] predicts token t+1 → logits[0] predicts first AA after "1".
        n_pred = min(len(chunk), max(0, logits.size(0) - 1))
        out = torch.zeros(n_pred, 20, dtype=torch.float32)
        for j in range(n_pred):
            aa_logits = logits[j, aa_ids]
            out[j] = F.log_softmax(aa_logits.float(), dim=-1).cpu()
        return out

    @torch.no_grad()
    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        """
        Sitewise AA log-probs. Sequences longer than ``max_ctx−1`` are scored with
        sliding left-context chunks (aggregated over the full length) — not truncated.
        """
        N, L = len(sequences), max_seq_len
        log_rates = torch.zeros(N, L, 20, dtype=torch.float32)
        aa_ids = self._aa_token_ids
        # Reserve one token for leading terminal "1".
        max_aa = max(1, self.max_ctx - 1)
        for i, seq in enumerate(sequences):
            clean = seq.replace("-", "X")[:L]
            if not clean:
                continue
            filled = 0
            L_seq = len(clean)
            while filled < L_seq:
                next_end = min(L_seq, filled + max_aa)
                win_start = max(0, next_end - max_aa)
                chunk = clean[win_start:next_end]
                window_lp = self._score_window(chunk, aa_ids)
                for j in range(window_lp.size(0)):
                    abs_j = win_start + j
                    if abs_j < filled or abs_j >= L:
                        continue
                    log_rates[i, abs_j, :] = window_lp[j]
                filled = next_end
        return log_rates


class SubstitutionMatrixR0Backend(R0Backend):
    """
    Sitewise destination log-probs from an AA substitution rate matrix (JTT/WAG/LG).

    For current residue a at each site, build a destination distribution from the
    CTMC jump chain of Q, with residual stay mass on a. Context-independent —
    paper Table 7 / D.1 "Substitution-only" / JTT/WAG/LG row.
    """

    def __init__(self, model: str = "JTT", stay_mass: float = 0.5):
        self.model_name = model.upper() if model.lower() != "neutral" else "NEUTRAL"
        self.name = self.model_name.lower()
        self.stay_mass = float(stay_mass)
        self._Q = self._load_rate_matrix(self.model_name)  # [20, 20] numpy

    @staticmethod
    def _load_rate_matrix(model_name: str):
        import numpy as np

        if model_name == "NEUTRAL":
            Q = np.full((20, 20), 1.0, dtype=np.float64)
            np.fill_diagonal(Q, -19.0)
            return Q
        try:
            import pyvolve
        except ImportError as e:
            raise ImportError(
                f"Substitution backend {model_name} needs pyvolve "
                "(see requirements.txt)."
            ) from e
        m = pyvolve.Model(model_name)
        # pyvolve stores instantaneous rate matrix as .matrix or similar.
        Q = getattr(m, "matrix", None)
        if Q is None and hasattr(m, "params"):
            Q = m.params.get("matrix")
        if Q is None:
            raise RuntimeError(
                f"Could not extract rate matrix from pyvolve.Model({model_name!r})"
            )
        Q = np.asarray(Q, dtype=np.float64)
        if Q.shape != (20, 20):
            raise RuntimeError(f"Expected 20x20 rate matrix, got {Q.shape}")
        return Q

    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        import numpy as np

        N, L = len(sequences), max_seq_len
        log_rates = torch.zeros(N, L, 20, dtype=torch.float32)
        stay = min(max(self.stay_mass, 1e-6), 1.0 - 1e-6)
        Q = self._Q
        for i, seq in enumerate(sequences):
            for pos, aa in enumerate(seq[:L]):
                a = AA_TO_IDX.get(aa)
                if a is None:
                    # gap / unknown → uniform
                    log_rates[i, pos, :] = -np.log(20.0)
                    continue
                row = Q[a].copy()
                off = np.maximum(row, 0.0)
                off[a] = 0.0
                s = off.sum()
                probs = np.full(20, (1.0 - stay) / 19.0, dtype=np.float64)
                if s > 1e-12:
                    probs = off / s * (1.0 - stay)
                probs[a] = stay
                probs = np.clip(probs, 1e-12, None)
                probs /= probs.sum()
                log_rates[i, pos, :] = torch.from_numpy(np.log(probs)).float()
        return log_rates


def build_r0_backend(
    backend: str,
    model_id: Optional[str] = None,
    device: Optional[str] = None,
) -> R0Backend:
    """
    Factory for Table 7 / D.1 R0 backends.

    ``model_id`` overrides the default HF / ProGen2 checkpoint when applicable.
    For ProGen2, ``model_id`` is a checkpoint folder name under ``PROGEN2_HOME``
    (default ``progen2-small``) or an absolute checkpoint path.
    """
    name = normalize_backend_name(backend)
    if name in STUB_BACKENDS:
        hints = {
            BACKEND_EVO2: (
                "Evo2 is heavy (OOM risk on mig GPUs); use API/distill or skip "
                "the appendix D.1 Evo2 row for now."
            ),
        }
        return StubR0Backend(name, hints[name])

    if name == BACKEND_ESM2:
        return ESM2R0Backend(
            model_id=model_id or DEFAULT_ESM2_MODEL,
            device=device,
            name=BACKEND_ESM2,
        )
    if name == BACKEND_ESM2_650M:
        return ESM2R0Backend(
            model_id=model_id or DEFAULT_ESM2_650M_MODEL,
            device=device,
            name=BACKEND_ESM2_650M,
        )
    if name == BACKEND_ESMC:
        return ESMCR0Backend(model_id=model_id or DEFAULT_ESMC_MODEL, device=device)
    if name == BACKEND_PROGEN2:
        return ProGen2R0Backend(
            model_id=model_id or DEFAULT_PROGEN2_MODEL, device=device
        )
    if name in SUBSTITUTION_BACKENDS:
        return SubstitutionMatrixR0Backend(model=name)
    if name == BACKEND_THRIFTY_AA:
        from src.bridge.thrifty_aa_q0 import ThriftyAAR0Backend

        return ThriftyAAR0Backend(
            model_name=model_id or "ThriftyHumV0.2-59",
            device=device,
        )
    if name == BACKEND_CODON_GY94:
        from src.bridge.codon_q0 import CodonGY94R0Backend

        return CodonGY94R0Backend()

    raise ValueError(
        f"Unknown R0 backend {backend!r}. "
        f"Choose from: esm2, esm2_650m, esmc, progen2, evo2, jtt, wag, lg, "
        f"neutral, thrifty_aa, codon_gy94."
    )


def list_backends() -> list[str]:
    return sorted(BACKEND_CACHE_TAG.keys())
