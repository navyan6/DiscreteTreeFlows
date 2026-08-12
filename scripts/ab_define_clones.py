#!/usr/bin/env python3
"""
Table 5 — define clonal lineages from donor + V gene + CDR3.

Primary key (locked): clonal lineage holdout unit =
  (donor_id, v_gene_family, cdr3_aa)

Joins subsample metadata.tsv with ANARCI TSV. Filters clone sizes to
[--min-size, --max-size] (default 16–64).

Usage:
  python scripts/ab_define_clones.py \\
    --metadata data/ab_t5_500k/metadata.tsv \\
    --anarci data/ab_t5_500k/anarci.tsv \\
    --out-dir data/ab_t5_500k/clones
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def v_family(v_gene: str) -> str:
    m = re.match(r"(IGH[VDJ]\d+|IGK[VDJ]\d+|IGL[VDJ]\d+)", (v_gene or "").upper())
    return m.group(1) if m else (v_gene or "UNK")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--anarci", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, default=None, help="Optional; write per-clone FASTA if set")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-size", type=int, default=16)
    ap.add_argument("--max-size", type=int, default=64)
    ap.add_argument("--require-donor", action="store_true", default=True)
    ap.add_argument("--allow-unknown-donor", action="store_true")
    args = ap.parse_args()

    anarci = {}
    with args.anarci.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            anarci[row["sequence_id"]] = row

    groups: dict[str, list[dict]] = defaultdict(list)
    n_meta = n_join = n_no_cdr3 = n_no_donor = 0
    with args.metadata.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n_meta += 1
            uid = row["uniq_id"]
            a = anarci.get(uid)
            if not a or not a.get("cdr3"):
                n_no_cdr3 += 1
                continue
            donor = row.get("donor_id") or "unknown"
            if donor == "unknown" and not args.allow_unknown_donor:
                n_no_donor += 1
                continue
            v = a.get("v_gene") or row.get("v_gene") or ""
            fam = v_family(v)
            cdr3 = a["cdr3"].upper().replace("-", "")
            clone_id = f"{donor}|{fam}|{cdr3}"
            n_join += 1
            groups[clone_id].append(
                {
                    "uniq_id": uid,
                    "donor_id": donor,
                    "v_gene": v,
                    "v_family": fam,
                    "cdr3": cdr3,
                    "j_gene": a.get("j_gene", ""),
                    "chain": row.get("chain", ""),
                    "len": row.get("len", ""),
                }
            )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    clones_dir = args.out_dir / "by_clone"
    clones_dir.mkdir(exist_ok=True)

    # Optional seq lookup
    seqs = {}
    if args.fasta and args.fasta.is_file():
        header = None
        parts: list[str] = []
        with args.fasta.open() as f:
            for line in f:
                if line.startswith(">"):
                    if header is not None:
                        seqs[header] = "".join(parts)
                    header = line[1:].strip().split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            if header is not None:
                seqs[header] = "".join(parts)

    kept = []
    size_hist = defaultdict(int)
    for cid, members in groups.items():
        size_hist[len(members)] += 1
        if not (args.min_size <= len(members) <= args.max_size):
            continue
        # dedup identical AA within clone
        seen = set()
        uniq_members = []
        for m in members:
            s = seqs.get(m["uniq_id"], "")
            if s and s in seen:
                continue
            if s:
                seen.add(s)
            uniq_members.append(m)
        if not (args.min_size <= len(uniq_members) <= args.max_size):
            continue
        safe = re.sub(r"[^\w.\-|]", "_", cid)[:180]
        rec = {
            "clone_id": cid,
            "safe_id": safe,
            "n_leaves": len(uniq_members),
            "donor_id": uniq_members[0]["donor_id"],
            "v_family": uniq_members[0]["v_family"],
            "cdr3": uniq_members[0]["cdr3"],
            "members": uniq_members,
        }
        kept.append(rec)
        mem_path = clones_dir / f"{safe}.members.tsv"
        with mem_path.open("w", newline="") as mf:
            w = csv.DictWriter(mf, fieldnames=list(uniq_members[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(uniq_members)
        if seqs:
            fa = clones_dir / f"{safe}.fasta"
            with fa.open("w") as ff:
                for m in uniq_members:
                    s = seqs.get(m["uniq_id"])
                    if not s:
                        continue
                    ff.write(f">{m['uniq_id']}\n")
                    for i in range(0, len(s), 80):
                        ff.write(s[i : i + 80] + "\n")

    index_path = args.out_dir / "clones_index.tsv"
    with index_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["clone_id", "safe_id", "n_leaves", "donor_id", "v_family", "cdr3"],
            delimiter="\t",
        )
        w.writeheader()
        for r in sorted(kept, key=lambda x: -x["n_leaves"]):
            w.writerow({k: r[k] for k in w.fieldnames})

    summary = {
        "n_metadata": n_meta,
        "n_joined_cdr3": n_join,
        "n_no_cdr3": n_no_cdr3,
        "n_no_donor": n_no_donor,
        "n_raw_clone_keys": len(groups),
        "n_kept_clones": len(kept),
        "min_size": args.min_size,
        "max_size": args.max_size,
        "size_hist_top": dict(sorted(size_hist.items(), key=lambda kv: -kv[1])[:30]),
        "holdout_unit": "clonal_lineage = donor_id|v_family|cdr3",
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if len(kept) == 0:
        print(
            "WARN: zero clones in size band — ANARCI CDR3 join may have failed, "
            "or 500K subsample lacks deep lineages. Loosen --min-size or expand subsample.",
            flush=True,
        )


if __name__ == "__main__":
    main()
