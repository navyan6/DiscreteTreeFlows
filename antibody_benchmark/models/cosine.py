"""CoSiNE adapter: unconditional Gillespie sampling (NOT Guided)."""

from __future__ import annotations

import sys
from pathlib import Path

from antibody_benchmark.models.base import EvolutionModel


class CosineModel(EvolutionModel):
    name = "cosine"
    alphabet = "aa"

    def __init__(self, ckpt_path: str | None = None, device: str | None = None):
        self.ckpt_path = Path(
            ckpt_path
            or "antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt"
        )
        if not self.ckpt_path.is_absolute():
            repo_root = Path(__file__).resolve().parents[2]
            cand = repo_root / self.ckpt_path
            if cand.is_file():
                self.ckpt_path = cand
        if device is None:
            try:
                import torch
                device = "cuda" if torch.cuda.is_available() else "cpu"
            except Exception:
                device = "cpu"
        self.device = device
        self._generator = None
        self._vocab = None
        self._load()

    def _ensure_cosine_on_path(self) -> None:
        repo = Path(__file__).resolve().parents[1] / "data" / "raw" / "repos" / "cosine"
        # Workspace layout: cosine/ + evo/ (evo package lives at evo/evo)
        for p in (repo, repo / "evo"):
            if p.is_dir() and str(p) not in sys.path:
                sys.path.insert(0, str(p))
        try:
            import cosine  # noqa: F401
        except ImportError as e:
            raise RuntimeError(
                "BLOCKER: CoSiNE package not importable. Ensure "
                "antibody_benchmark/data/raw/repos/cosine is present and its "
                f"evo submodule is populated. Error: {e}"
            ) from e
        # HF checkpoint was pickled under the historical package name ``peint``.
        self._install_peint_alias()

    def _install_peint_alias(self) -> None:
        """Map historical ``peint`` package name onto current ``cosine`` package."""
        import cosine  # noqa: F401

        sys.modules.setdefault("peint", sys.modules["cosine"])

        class _PeintAliasFinder:
            def find_module(self, fullname, path=None):  # noqa: ANN001
                if fullname == "peint" or fullname.startswith("peint."):
                    return self
                return None

            def load_module(self, fullname):  # noqa: ANN001
                if fullname in sys.modules:
                    return sys.modules[fullname]
                target = "cosine" if fullname == "peint" else "cosine." + fullname[len("peint.") :]
                mod = __import__(target, fromlist=["*"])
                # Ensure leaf module object
                for part in target.split(".")[1:]:
                    mod = getattr(mod, part)
                sys.modules[fullname] = mod
                return mod

        # Prefer PathFinder-style if present; keep a simple meta_path hook.
        if not any(type(x).__name__ == "_PeintAliasFinder" for x in sys.meta_path):
            sys.meta_path.insert(0, _PeintAliasFinder())
        for name, mod in list(sys.modules.items()):
            if name.startswith("cosine."):
                sys.modules.setdefault("peint." + name[len("cosine.") :], mod)

    def _load(self) -> None:
        if not self.ckpt_path.is_file():
            raise RuntimeError(
                f"BLOCKER: CoSiNE checkpoint not found at {self.ckpt_path}. "
                "Download with: hf download thematrixmaster/cosine cosine_dasm.ckpt "
                "--local-dir antibody_benchmark/data/raw/cosine/checkpoints"
            )
        self._ensure_cosine_on_path()
        try:
            import torch
            from cosine.models.nets.ctmc import NeuralCTMC, NeuralCTMCGenerator  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "BLOCKER: CoSiNE CTMC imports failed. Need cosine + evo on PYTHONPATH "
                f"and lightning installed. Error: {e}"
            ) from e
        self._install_peint_alias()

        try:
            # Official path: Lightning CTMCModule.load_from_checkpoint
            CTMCModule = None
            for modpath in (
                "cosine.models.modules.ctmc_module",
                "cosine.models.ctmc_module",
                "cosine.models.module",
            ):
                try:
                    CTMCModule = __import__(modpath, fromlist=["CTMCModule"]).CTMCModule
                    break
                except Exception:
                    continue

            self._install_peint_alias()
            if CTMCModule is not None:
                # PyTorch 2.6+ defaults weights_only=True; ckpt contains full objects.
                try:
                    module = CTMCModule.load_from_checkpoint(
                        str(self.ckpt_path),
                        map_location=self.device,
                        strict=False,
                        weights_only=False,
                    )
                except TypeError:
                    # Older Lightning without weights_only kw
                    module = CTMCModule.load_from_checkpoint(
                        str(self.ckpt_path), map_location=self.device, strict=False
                    )
                net: NeuralCTMC = module.net
            else:
                ckpt = torch.load(self.ckpt_path, map_location=self.device, weights_only=False)
                if isinstance(ckpt, NeuralCTMC):
                    net = ckpt
                elif isinstance(ckpt, dict) and "net" in ckpt:
                    net = ckpt["net"]
                else:
                    raise RuntimeError(
                        "Could not locate CTMCModule; install cosine with Lightning "
                        "module or provide a checkpoint loadable as NeuralCTMC."
                    )
            net = net.eval().to(self.device)
            self._generator = NeuralCTMCGenerator(neural_ctmc=net)
            self._vocab = net.vocab
            self._torch = torch
        except RuntimeError:
            raise
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: failed to instantiate CoSiNE Gillespie generator from "
                f"{self.ckpt_path}. Error: {e}"
            ) from e

    def _encode(self, seq: str):
        import torch

        vocab = self._vocab
        # Match cosine tokenization: optional BOS/EOS around AA string
        toks = []
        if getattr(vocab, "prepend_bos", False):
            toks.append(vocab.bos_idx)
        for ch in seq:
            toks.append(vocab.tokens_to_idx.get(ch, vocab.unk_idx))
        if getattr(vocab, "append_eos", False):
            toks.append(vocab.eos_idx)
        x = torch.tensor([toks], dtype=torch.long, device=self.device)
        x_sizes = torch.tensor([len(toks)], dtype=torch.long, device=self.device)
        return x, x_sizes

    def _decode(self, y) -> str:
        vocab = self._vocab
        # Prefer official CoSiNE helper when available.
        try:
            from cosine.models.frameworks.ctmc import decode_sequence_from_toks

            if hasattr(y, "dim") and y.dim() >= 2:
                return decode_sequence_from_toks(y[0], vocab)
            return decode_sequence_from_toks(y, vocab)
        except Exception:
            pass
        ids = y[0].tolist() if hasattr(y, "dim") and y.dim() >= 2 else list(y)
        chars = []
        skip = {
            getattr(vocab, "bos_idx", -1),
            getattr(vocab, "eos_idx", -1),
            getattr(vocab, "pad_idx", -1),
            getattr(vocab, "unk_idx", -1),
        }
        for i in ids:
            if i in skip:
                continue
            if i == getattr(vocab, "eos_idx", -2) or i == getattr(vocab, "pad_idx", -2):
                break
            if hasattr(vocab, "token"):
                tok = vocab.token(i)
            elif hasattr(vocab, "idx_to_tokens"):
                tok = vocab.idx_to_tokens[i]
            else:
                raise RuntimeError("BLOCKER: Vocab lacks token()/idx_to_tokens for decode")
            if tok in {".", "-", "<pad>", "<bos>", "<eos>", "<unk>", "<mask>", "<null_1>"}:
                continue
            chars.append(tok)
        return "".join(chars)

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        """Unguided Gillespie only — never set use_guidance / oracle."""
        if self._generator is None:
            raise RuntimeError("BLOCKER: CosineModel generator not loaded.")
        torch = self._torch
        torch.manual_seed(int(rng_seed))
        try:
            if hasattr(self._generator, "encode_sequences"):
                x, x_sizes = self._generator.encode_sequences([parent_seq])  # type: ignore[attr-defined]
            else:
                x, x_sizes = self._encode(parent_seq)
            t = torch.tensor([float(branch_length)], device=self.device)
            # CoSiNE ESM2/FlashAttention requires fp16/bf16 (official scripts use autocast).
            if str(self.device).startswith("cuda"):
                amp = torch.autocast(device_type="cuda", dtype=torch.bfloat16)
            else:
                from contextlib import nullcontext
                amp = nullcontext()
            with amp:
                y = self._generator.generate_with_gillespie(
                    x=x,
                    t=t,
                    x_sizes=x_sizes,
                    temperature=1.0,
                    no_special_toks=True,
                    max_decode_steps=1000,
                    use_scalar_steps=False,
                    verbose=False,
                    use_guidance=False,
                    oracle=None,
                )
            if hasattr(self._generator, "decode_sequences"):
                return self._generator.decode_sequences(y)[0]  # type: ignore[attr-defined]
            return self._decode(y)
        except Exception as e:
            raise RuntimeError(
                "BLOCKER: CoSiNE unguided Gillespie sample failed. Ensure you are "
                "calling generate_with_gillespie WITHOUT oracle / use_guidance. "
                f"Error: {e}"
            ) from e
