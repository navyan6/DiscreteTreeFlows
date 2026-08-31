"""Shared helpers for epidemic-aware TreeSBM data prep (flu / COVID)."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from Bio import SeqIO


def slug_label(raw: str, max_len: int = 64) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "_", (raw or "").strip()).strip("_").lower()
    return s[:max_len] or "unknown"


def nh_flu_season(year: int, month: int) -> str:
    """Northern-hemisphere flu season label Oct(Y-1)–Sep(Y) as Y-1_Y."""
    if month >= 10:
        return f"{year}-{year + 1}"
    return f"{year - 1}-{year}"


def parse_flu_header(desc: str) -> tuple[str, str, int | None, int | None, str | None]:
    """>EPI_ISL,YYYY-MM-DD -> (id, date_str, year, month, season)."""
    parts = desc.split(",", 1)
    rid = parts[0].strip()
    date_str = parts[1].strip() if len(parts) > 1 else "2000-01-01"
    year = month = None
    season = None
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", date_str)
    if m:
        year, month = int(m.group(1)), int(m.group(2))
        season = nh_flu_season(year, month)
    elif re.match(r"(\d{4})", date_str):
        year = int(date_str[:4])
        month = 7
        season = nh_flu_season(year, month)
    return rid, date_str, year, month, season


def assign_temporal_split(
    year: int,
    train_end: int,
    val_year: int | None,
    test_start: int,
    test_end: int,
) -> str | None:
    if year <= train_end:
        return "train"
    if val_year is not None and year == val_year:
        return "val"
    if test_start <= year <= test_end:
        return "test"
    return None


def chunk_records(
    recs: list[tuple],
    group_size: int,
    min_group: int,
    sort_key_idx: int = 2,
) -> list[list[tuple]]:
    """Date-sort and chunk; drop tails smaller than min_group."""
    recs = sorted(recs, key=lambda r: r[sort_key_idx])
    if len(recs) <= group_size:
        return [recs] if len(recs) >= min_group else []
    chunks: list[list[tuple]] = []
    for i in range(0, len(recs), group_size):
        chunk = recs[i : i + group_size]
        if len(chunk) >= min_group:
            chunks.append(chunk)
        elif not chunks and len(chunk) >= max(3, min_group // 2):
            chunks.append(chunk)
    return chunks


def write_fasta_csv_groups(
    group_lists: list[list[tuple]],
    out_dir: Path,
    prefix: str,
    start_group: int,
    extra_csv_cols: list[str] | None = None,
    extra_row_fn=None,
) -> tuple[int, int]:
    """
    Write group FASTA/CSV files.
    Each record tuple: (acc, date_str, sort_key, seq, *extras matching extra_csv_cols).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    g = start_group
    n_seq = 0
    cols = ["name", "date"] + (extra_csv_cols or [])
    for chunk in group_lists:
        gf = out_dir / f"{prefix}_group_{g:03d}.fasta"
        gcsv = out_dir / f"{prefix}_group_{g:03d}.csv"
        with open(gf, "w") as ff, open(gcsv, "w", newline="") as cf:
            w = csv.writer(cf)
            w.writerow(cols)
            for row in chunk:
                acc, date_str = row[0], row[1]
                seq = row[3]
                ff.write(f">{acc},{date_str}\n{seq}\n")
                if extra_row_fn:
                    w.writerow(extra_row_fn(row))
                elif extra_csv_cols:
                    w.writerow([acc, date_str] + list(row[4 : 4 + len(extra_csv_cols)]))
                else:
                    w.writerow([acc, date_str])
        n_seq += len(chunk)
        g += 1
    return g, n_seq


def write_split_protocol(base: Path, protocol: dict[str, Any]) -> None:
    base.mkdir(parents=True, exist_ok=True)
    (base / "SPLIT_PROTOCOL.json").write_text(
        __import__("json").dumps(protocol, indent=2) + "\n"
    )


def count_groups(split_dir: Path) -> tuple[int, int]:
    if not split_dir.is_dir():
        return 0, 0
    fastas = list(split_dir.glob("*_group_*.fasta"))
    n_seq = 0
    for f in fastas:
        n_seq += sum(1 for _ in SeqIO.parse(f, "fasta"))
    return len(fastas), n_seq
