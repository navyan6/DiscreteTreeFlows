#!/usr/bin/env python3
"""
Time the real tree pipeline so the pan-viral schedule is measured, not guessed.

Runs the exact commands from scripts/run_all_groups.py -- mafft --auto,
FastTree -gtr -nt, augur refine --timetree --date-inference marginal, and
augur ancestral --inference joint -- over a sweep of leaf counts, using real
Spike sequences pooled from the existing COVID groups.

augur refine is the stage that decides the schedule: TreeTime's marginal date
inference scales far worse than the alignment or tree search, so extrapolating
from a small run underestimates it badly.

    python scripts/panviral/benchmark_pipeline.py --sizes 75 150 300 600
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ENV = Path("/vast/home/n/nnori/.conda/envs/treesbm/bin")


def sh(cmd: list[str], stdout_to: Path | None = None,
       timeout: int = 14400) -> tuple[float, bool, str]:
    """Run a stage, returning (seconds, ok, tail-of-stderr)."""
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return time.time() - t0, False, f"TIMEOUT after {timeout}s"
    el = time.time() - t0
    if r.returncode != 0:
        return el, False, r.stderr[-400:]
    if stdout_to:
        stdout_to.write_text(r.stdout)
    return el, True, ""


def pool(src_dirs: list[Path], need: int) -> tuple[list[SeqRecord], dict[str, str]]:
    """Gather de-gapped sequences plus their collection dates."""
    recs: list[SeqRecord] = []
    dates: dict[str, str] = {}
    seen: set[str] = set()

    for d in src_dirs:
        for csv_path in sorted(d.glob("*_group_*.csv")) + sorted(d.glob("group_*_meta.csv")):
            try:
                with csv_path.open() as fh:
                    for row in csv.DictReader(fh):
                        if row.get("name") and row.get("date"):
                            dates[row["name"]] = row["date"]
            except Exception:  # noqa: BLE001 - skip unreadable metadata
                continue

    for d in src_dirs:
        for fa in sorted(d.glob("group_*_aligned.fasta")):
            for rec in SeqIO.parse(fa, "fasta"):
                rid = rec.id.split(",")[0]
                if rid in seen or rid not in dates:
                    continue
                # Dates must be fully resolved; TreeTime cannot use "2021-XX-XX".
                if "X" in dates[rid] or len(dates[rid]) < 10:
                    continue
                seq = str(rec.seq).replace("-", "").upper()
                if len(seq) < 3000:
                    continue
                seen.add(rid)
                recs.append(SeqRecord(Seq(seq), id=rid, description=""))
            if len(recs) >= need:
                return recs, dates
    return recs, dates


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", type=int, nargs="+", default=[75, 150, 300, 600])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path,
                    default=Path("/vast/projects/pranam/lab/nnori/panviral_bench"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/panviral/pipeline_timings.json"))
    ap.add_argument("--src", type=Path, nargs="+", default=[
        Path("/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data/covid/test"),
        Path("/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data/covid/train"),
    ])
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    need = max(args.sizes)
    src_dirs = [d for d in args.src if d.exists()]
    print(f"pooling up to {need} sequences from {len(src_dirs)} dirs ...", flush=True)
    recs, dates = pool(src_dirs, need)
    print(f"pooled {len(recs)} dated sequences "
          f"(median len {sorted(len(r.seq) for r in recs)[len(recs)//2]} nt)\n",
          flush=True)
    if len(recs) < min(args.sizes):
        sys.exit("not enough dated sequences to benchmark")

    random.seed(0)
    results = []
    hdr = (f"{'N':>6}{'mafft':>10}{'fasttree':>10}{'refine':>10}"
           f"{'ancestral':>11}{'total':>10}   notes")
    print(hdr)
    print("-" * len(hdr))

    for n in args.sizes:
        if n > len(recs):
            print(f"{n:>6}   skipped (only {len(recs)} available)")
            continue
        d = args.work / f"n{n}"
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

        sub = random.sample(recs, n)
        raw = d / "raw.fasta"
        SeqIO.write(sub, raw, "fasta")
        meta = d / "meta.csv"
        with meta.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["name", "date"])
            w.writeheader()
            for r in sub:
                w.writerow({"name": r.id, "date": dates[r.id]})

        aligned, tree = d / "aligned.fasta", d / "tree.nwk"
        row: dict = {"n": n}
        note = ""

        t_a, ok, err = sh([str(ENV / "mafft"), "--auto", "--thread",
                           str(args.threads), str(raw)], stdout_to=aligned)
        row["mafft"] = round(t_a, 1)
        if not ok:
            print(f"{n:>6}   mafft FAILED: {err[:80]}")
            continue

        t_f, ok, err = sh([str(ENV / "fasttree"), "-gtr", "-nt", str(aligned)],
                          stdout_to=tree)
        row["fasttree"] = round(t_f, 1)
        if not ok:
            print(f"{n:>6}   fasttree FAILED: {err[:80]}")
            continue

        t_r, ok, err = sh([
            str(ENV / "augur"), "refine",
            "--tree", str(tree), "--alignment", str(aligned),
            "--metadata", str(meta), "--metadata-id-columns", "name",
            "--output-tree", str(d / "rooted.nwk"),
            "--output-node-data", str(d / "bl.json"),
            "--timetree", "--coalescent", "opt",
            "--date-inference", "marginal", "--clock-filter-iqd", "4",
        ])
        row["refine"] = round(t_r, 1)
        row["refine_ok"] = ok
        if not ok:
            note = f"refine failed: {err.strip().splitlines()[-1][:60] if err.strip() else '?'}"

        t_n = 0.0
        if ok:
            t_n, ok2, err2 = sh([
                str(ENV / "augur"), "ancestral",
                "--tree", str(d / "rooted.nwk"), "--alignment", str(aligned),
                "--output-sequences", str(d / "anc_nt.fasta"),
                "--inference", "joint",
            ])
            row["ancestral"] = round(t_n, 1)
            if not ok2:
                note = f"ancestral failed: {err2[:60]}"

        row["total"] = round(t_a + t_f + t_r + t_n, 1)
        row["note"] = note
        results.append(row)
        print(f"{n:>6}{t_a:>9.1f}s{t_f:>9.1f}s{t_r:>9.1f}s{t_n:>10.1f}s"
              f"{row['total']:>9.1f}s   {note}", flush=True)

        args.out.write_text(json.dumps({"threads": args.threads,
                                        "rows": results}, indent=2) + "\n")

    # Empirical scaling exponent: t ~ N^k, fitted on the two extreme sizes.
    ok_rows = [r for r in results if r.get("refine_ok")]
    if len(ok_rows) >= 2:
        import math
        a, b = ok_rows[0], ok_rows[-1]
        print("\nscaling exponent k in t ~ N^k (from "
              f"N={a['n']} to N={b['n']}):")
        for stage in ("mafft", "fasttree", "refine", "total"):
            if a.get(stage, 0) > 0 and b.get(stage, 0) > 0:
                k = math.log(b[stage] / a[stage]) / math.log(b["n"] / a["n"])
                print(f"   {stage:<10} k = {k:.2f}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
