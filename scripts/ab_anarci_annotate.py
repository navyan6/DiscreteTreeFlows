#!/usr/bin/env python3
"""
Table 5 — ANARCI annotation on AA antibody sequences.

Prefers the Python `anarci` API (extracts IMGT CDR1–3 + V/J). Falls back to
CLI `--csv` and rebuilds CDRs from numbered IMGT position columns.

Requires ANARCI + HMMER on PATH / in the active env.

Usage:
  python scripts/ab_anarci_annotate.py \\
    --fasta data/ab_t5_500k/sequences.fasta \\
    --out data/ab_t5_500k/anarci.tsv \\
    --batch-size 500 --n-jobs 8
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# IMGT CDR ranges (inclusive position labels as used by ANARCI numbering)
_IMGT_CDR = {
    "cdr1": (27, 38),
    "cdr2": (56, 65),
    "cdr3": (105, 117),
}


def iter_fasta(path: Path):
    header = None
    parts: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts)


def _pos_key_sort(k: str):
    # ANARCI CSV columns: "105", "105A", ...
    k = str(k)
    num = "".join(ch for ch in k if ch.isdigit())
    suf = "".join(ch for ch in k if ch.isalpha())
    return (int(num) if num else 10**9, suf)


def cdr_from_numbered_cols(row: dict) -> dict:
    """Pull CDR AA strings from IMGT-numbered CSV columns."""
    out = {"cdr1": "", "cdr2": "", "cdr3": ""}
    # Collect position -> AA
    pos_aa = {}
    for k, v in row.items():
        if v is None:
            continue
        ks = str(k).strip()
        if not ks or not ks[0].isdigit():
            continue
        aa = str(v).strip()
        if not aa or aa in {"-", ".", "X", "x"}:
            continue
        pos_aa[ks] = aa.upper()
    if not pos_aa:
        return out
    keys = sorted(pos_aa.keys(), key=_pos_key_sort)
    for name, (lo, hi) in _IMGT_CDR.items():
        chars = []
        for k in keys:
            num = int("".join(ch for ch in k if ch.isdigit()) or "0")
            if lo <= num <= hi:
                chars.append(pos_aa[k])
        out[name] = "".join(chars)
    return out


def parse_anarci_csv_row(row: dict) -> dict:
    sid = row.get("Id") or row.get("id") or row.get("sequence_id") or row.get("seq_id") or ""
    v = row.get("v_gene") or row.get("v-gene") or row.get("hmm_species") or ""
    j = row.get("j_gene") or row.get("j-gene") or ""
    # Prefer explicit CDR columns if present
    cdrs = {
        "cdr1": (row.get("cdr1") or row.get("CDR1") or "").upper(),
        "cdr2": (row.get("cdr2") or row.get("CDR2") or "").upper(),
        "cdr3": (row.get("cdr3") or row.get("CDR3") or "").upper(),
    }
    if not cdrs["cdr3"]:
        cdrs = cdr_from_numbered_cols(row)
    return {
        "sequence_id": sid,
        "v_gene": str(v),
        "j_gene": str(j),
        **cdrs,
    }


def run_anarci_python(records: list[tuple[str, str]], scheme: str = "imgt") -> list[dict]:
    """Use anarci Python API; extract CDRs from numbering."""
    from anarci import anarci as anarci_run  # type: ignore

    results = anarci_run(
        records,
        scheme=scheme,
        output=False,
        assign_germline=True,
        allow=set(["H", "K", "L", "A", "B", "G", "D"]),
    )
    # anarci returns (numbering, alignment_details, hit_tables) — version-dependent
    numbering, details, _hits = results
    rows: list[dict] = []
    for (sid, _seq), num_list, det_list in zip(records, numbering, details):
        if not num_list:
            rows.append(
                {
                    "sequence_id": sid,
                    "v_gene": "",
                    "j_gene": "",
                    "cdr1": "",
                    "cdr2": "",
                    "cdr3": "",
                }
            )
            continue
        # Take top hit
        numbered, start, end = num_list[0]
        det = det_list[0] if det_list else {}
        pos_aa = {}
        for (pos, ins), aa in numbered:
            if aa in {"-", "."}:
                continue
            key = f"{pos}{ins or ''}"
            pos_aa[key] = aa
        fake_row = dict(pos_aa)
        cdrs = cdr_from_numbered_cols(fake_row)
        v_gene = ""
        j_gene = ""
        if isinstance(det, dict):
            germ = det.get("germlines") or {}
            # germlines often {"v_gene": [(id, score)], "j_gene": [...]}
            v_list = germ.get("v_gene") or germ.get("V") or []
            j_list = germ.get("j_gene") or germ.get("J") or []
            if v_list:
                v_gene = v_list[0][0] if isinstance(v_list[0], (list, tuple)) else str(v_list[0])
            if j_list:
                j_gene = j_list[0][0] if isinstance(j_list[0], (list, tuple)) else str(j_list[0])
            v_gene = v_gene or str(det.get("v_gene") or det.get("id") or "")
            j_gene = j_gene or str(det.get("j_gene") or "")
        rows.append(
            {
                "sequence_id": sid,
                "v_gene": v_gene,
                "j_gene": j_gene,
                **cdrs,
            }
        )
    return rows


def run_anarci_cli(records: list[tuple[str, str]], scheme: str = "imgt") -> list[dict]:
    if not shutil.which("ANARCI") and not shutil.which("anarci"):
        raise SystemExit(
            "ANARCI not found on PATH. Install (conda/pip) before Table 5 annotate.\n"
            "  e.g. conda install -c bioconda anarci hmmer"
        )
    bin_name = "ANARCI" if shutil.which("ANARCI") else "anarci"
    rows: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        fa = tdir / "batch.fasta"
        with fa.open("w") as f:
            for sid, seq in records:
                f.write(f">{sid}\n{seq}\n")
        out_prefix = tdir / "out"
        cmd = [bin_name, "-i", str(fa), "-o", str(out_prefix), "-s", scheme, "--csv"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write((proc.stderr or "")[-2000:] + "\n")
            raise RuntimeError(f"ANARCI failed rc={proc.returncode}")
        cands = list(tdir.glob("out*.csv"))
        if not cands:
            raise RuntimeError(f"ANARCI produced no CSV under {tdir}")
        for csv_path in cands:
            with csv_path.open() as f:
                for r in csv.DictReader(f):
                    rows.append(parse_anarci_csv_row(dict(r)))
    return rows


def _annotate_batch(args_tuple):
    records, scheme, prefer_python = args_tuple
    if prefer_python:
        try:
            return run_anarci_python(records, scheme=scheme)
        except Exception as e:
            # fall through to CLI
            err = str(e)
            if "No module named" in err or "anarci" in err.lower():
                pass
            else:
                # API present but batch failed — try CLI
                pass
    return run_anarci_cli(records, scheme=scheme)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--scheme", default="imgt")
    ap.add_argument("--max-seqs", type=int, default=None)
    ap.add_argument("--cli-only", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="Count only; do not call ANARCI")
    args = ap.parse_args()

    records = []
    for sid, seq in iter_fasta(args.fasta):
        records.append((sid, seq))
        if args.max_seqs is not None and len(records) >= args.max_seqs:
            break

    print(f"Loaded {len(records):,} sequences from {args.fasta}", flush=True)
    if args.dry_run:
        print("dry-run: ANARCI not invoked")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("sequence_id\tv_gene\tj_gene\tcdr1\tcdr2\tcdr3\tok\n")
        return

    prefer_python = not args.cli_only
    if prefer_python:
        try:
            import anarci  # noqa: F401
            print("Using anarci Python API", flush=True)
        except ImportError:
            print("anarci Python API unavailable; using CLI", flush=True)
            prefer_python = False
            if not shutil.which("ANARCI") and not shutil.which("anarci"):
                raise SystemExit(
                    "ANARCI not found. Install: conda install -c bioconda anarci hmmer"
                )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sequence_id", "v_gene", "j_gene", "cdr1", "cdr2", "cdr3", "ok"]
    n_ok = 0
    batches = []
    for i in range(0, len(records), args.batch_size):
        batches.append(records[i : i + args.batch_size])

    # Map results back in order
    batch_results: list[list[dict] | None] = [None] * len(batches)
    n_jobs = max(1, args.n_jobs)
    if n_jobs == 1 or len(batches) == 1:
        for bi, batch in enumerate(batches):
            print(f"ANARCI batch {bi + 1}/{len(batches)} size={len(batch)}", flush=True)
            try:
                batch_results[bi] = _annotate_batch((batch, args.scheme, prefer_python))
            except Exception as e:
                print(f"WARN: batch {bi + 1} failed: {e}", flush=True)
                batch_results[bi] = []
    else:
        print(f"ANARCI parallel n_jobs={n_jobs} n_batches={len(batches)}", flush=True)
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futs = {
                ex.submit(_annotate_batch, (batch, args.scheme, prefer_python)): bi
                for bi, batch in enumerate(batches)
            }
            done = 0
            for fut in as_completed(futs):
                bi = futs[fut]
                done += 1
                try:
                    batch_results[bi] = fut.result()
                except Exception as e:
                    print(f"WARN: batch {bi + 1} failed: {e}", flush=True)
                    batch_results[bi] = []
                if done % 5 == 0 or done == len(batches):
                    print(f"  … completed {done}/{len(batches)} batches", flush=True)

    with args.out.open("w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for bi, batch in enumerate(batches):
            raw_rows = batch_results[bi] or []
            by_id = {}
            for r in raw_rows:
                sid = r.get("sequence_id") or ""
                if sid:
                    by_id[sid] = r
            for sid, _seq in batch:
                p = by_id.get(sid, {})
                ok = bool(p.get("cdr3"))
                n_ok += int(ok)
                w.writerow(
                    {
                        "sequence_id": sid,
                        "v_gene": p.get("v_gene", ""),
                        "j_gene": p.get("j_gene", ""),
                        "cdr1": p.get("cdr1", ""),
                        "cdr2": p.get("cdr2", ""),
                        "cdr3": p.get("cdr3", ""),
                        "ok": int(ok),
                    }
                )

    summary = {
        "fasta": str(args.fasta),
        "out": str(args.out),
        "n_input": len(records),
        "n_ok_cdr3": n_ok,
        "scheme": args.scheme,
        "n_jobs": n_jobs,
        "batch_size": args.batch_size,
        "backend": "python_api" if prefer_python else "cli",
    }
    args.out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
