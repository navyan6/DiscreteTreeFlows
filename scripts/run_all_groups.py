#!/usr/bin/env python3
"""
Full pipeline for all H3N2 groups (resumable + parallel workers).

Stages per group:
  1. clean      strip dates from FASTA headers  →  group_NNN_clean.fasta
  2. align      MAFFT                           →  group_NNN_aligned.fasta
  3. fasttree   FastTree GTR                    →  group_NNN_tree.nwk
  4. refine     augur refine (temporal root)    →  group_NNN_rooted.nwk + group_NNN_bl.json
  5. anc        augur ancestral                 →  group_NNN_anc_nt.fasta
  6. translate  NT → AA                         →  group_NNN_anc_aa.fasta
  7. embed      NodeEncoder → .pt               →  group_NNN_embeddings.pt

Usage:
    python scripts/run_all_groups.py                    # all 48 groups, 4 workers
    python scripts/run_all_groups.py --workers 8        # more parallelism
    python scripts/run_all_groups.py --groups 2 3 4     # specific groups
    python scripts/run_all_groups.py --stop-after refine
"""

import argparse
import json
import os
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data" / "train"
AUGUR_BIN = os.environ.get("AUGUR_BIN", "augur")  # override with full path if needed
FASTTREE_BIN = "fasttree"

# ── helpers ───────────────────────────────────────────────────────────────────

LOG_DIR = ROOT / "logs"

def log(g: int, msg: str):
    LOG_DIR.mkdir(exist_ok=True)
    line = f"[{g:03d}] {msg}"
    print(line, flush=True)
    with open(LOG_DIR / f"group_{g:03d}.log", "a") as f:
        f.write(line + "\n")


def run(cmd: list[str], desc: str, capture_stdout: bool = False) -> str:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{desc} failed:\n{result.stderr[-2000:]}")
    return result.stdout if capture_stdout else ""


def nextstrain_run(cmd: list[str], desc: str):
    """Run an augur command using the Nextstrain runtime binary directly."""
    # cmd[0] is 'augur' — replace with full path
    full = [AUGUR_BIN] + cmd[1:]
    run(full, desc)


def done(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


# ── pipeline stages ───────────────────────────────────────────────────────────

def stage_clean(g: int) -> tuple[Path, Path]:
    """
    Returns (clean_fasta, meta_csv) — both deduplicated on EPI_ISL ID.
    The meta_csv is used downstream by augur refine so IDs are consistent
    across FASTA, tree leaves, metadata, and all augur outputs.
    """
    src_fasta = DATA / f"master_h3n2_group_{g:03d}.fasta"
    src_csv   = DATA / f"master_h3n2_group_{g:03d}.csv"
    out_fasta = DATA / f"group_{g:03d}_clean.fasta"
    out_csv   = DATA / f"group_{g:03d}_meta.csv"

    if done(out_fasta) and done(out_csv):
        return out_fasta, out_csv

    # Build deduplicated FASTA — keep first occurrence of each ID
    seen = set()
    records = []
    for rec in SeqIO.parse(src_fasta, "fasta"):
        clean_id = rec.id.split(",")[0]
        if clean_id in seen:
            continue
        seen.add(clean_id)
        records.append(SeqRecord(rec.seq.upper(), id=clean_id, description=""))
    SeqIO.write(records, out_fasta, "fasta")

    # Write filtered metadata CSV containing only the kept IDs
    import csv
    kept_ids = {r.id for r in records}
    with open(src_csv) as fin, open(out_csv, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=["name", "date"])
        writer.writeheader()
        for row in reader:
            if row["name"] in kept_ids:
                writer.writerow({"name": row["name"], "date": row["date"]})

    log(g, f"clean → {len(records)} seqs (deduped), meta CSV written")
    return out_fasta, out_csv


def stage_align(g: int, clean: Path, meta_csv: Path) -> Path:
    out = DATA / f"group_{g:03d}_aligned.fasta"
    if done(out):
        return out
    stdout = run(["mafft", "--auto", "--thread", "2", str(clean)], "mafft", capture_stdout=True)
    out.write_text(stdout)
    log(g, "aligned")
    return out


def stage_fasttree(g: int, aligned: Path) -> Path:
    out = DATA / f"group_{g:03d}_tree.nwk"
    if done(out):
        return out
    # FastTree writes tree to stdout, progress to stderr
    result = subprocess.run(
        [FASTTREE_BIN, "-gtr", "-nt", str(aligned)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"fasttree failed:\n{result.stderr[-1000:]}")
    out.write_text(result.stdout)
    log(g, "fasttree done")
    return out


def stage_refine(g: int, treefile: Path, aligned: Path, meta_csv: Path) -> tuple[Path, Path]:
    rooted = DATA / f"group_{g:03d}_rooted.nwk"
    bl_json = DATA / f"group_{g:03d}_bl.json"
    if done(rooted) and done(bl_json):
        return rooted, bl_json
    nextstrain_run([
        "augur", "refine",
        "--tree", str(treefile),
        "--alignment", str(aligned),
        "--metadata", str(meta_csv),
        "--metadata-id-columns", "name",
        "--output-tree", str(rooted),
        "--output-node-data", str(bl_json),
        "--timetree",
        "--coalescent", "opt",
        "--date-inference", "marginal",
        "--clock-filter-iqd", "4",
    ], "augur refine")
    log(g, "refine done")
    return rooted, bl_json


def stage_ancestral(g: int, rooted: Path, aligned: Path) -> Path:
    out = DATA / f"group_{g:03d}_anc_nt.fasta"
    if done(out):
        return out
    nextstrain_run([
        "augur", "ancestral",
        "--tree", str(rooted),
        "--alignment", str(aligned),
        "--output-sequences", str(out),
        "--inference", "joint",
    ], "augur ancestral")
    log(g, "ancestral done")
    return out


def stage_translate(g: int, anc_nt: Path) -> Path:
    out = DATA / f"group_{g:03d}_anc_aa.fasta"
    if done(out):
        return out
    sys.path.insert(0, str(ROOT))
    from src.treeencoder.seq_utils import nt_to_aa
    records = [
        SeqRecord(Seq(nt_to_aa(str(rec.seq))), id=rec.id, description="")
        for rec in SeqIO.parse(anc_nt, "fasta")
    ]
    SeqIO.write(records, out, "fasta")
    log(g, f"translated {len(records)} seqs")
    return out


def stage_embed(g: int, rooted: Path, anc_aa: Path, bl_json: Path) -> Path:
    out = DATA / f"group_{g:03d}_embeddings.pt"
    if done(out):
        return out

    import torch
    sys.path.insert(0, str(ROOT))
    from scripts.embed_tree import parse_newick
    from src.treeencoder.plm_embeddings import ESM2Embedder
    from src.treeencoder.tree_adapter import tree_state_to_encoder_input
    from src.treeencoder.node_encoder import NodeEncoder
    from src.treeencoder.attention_mask import build_temporal_attention_mask
    from src.tree_state import TreeState

    device = "cuda" if torch.cuda.is_available() else "cpu"

    root_id, node_ids, edges, branch_lengths = parse_newick(str(rooted))
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(anc_aa, "fasta")}

    missing = [nid for nid in node_ids if nid not in seqs]
    if missing:
        ref_len = len(next(iter(seqs.values())))
        for nid in missing:
            seqs[nid] = "-" * ref_len

    with open(bl_json) as f:
        node_data = json.load(f)["nodes"]
    node_times = {nid: node_data.get(nid, {}).get("numdate", 0.0) for nid in node_ids}

    has_children = {p for p, _ in edges}
    tree_state = TreeState(
        node_ids=node_ids, root_id=root_id, edges=edges,
        branch_lengths=branch_lengths, node_seqs=seqs,
        active_leaves=[nid for nid in node_ids if nid not in has_children],
    )

    embedder = ESM2Embedder(device=device)
    plm_embeddings = embedder.embed_sequences([seqs[nid] for nid in node_ids])

    enc_input = tree_state_to_encoder_input(tree_state, plm_embeddings, laplacian_dim=8)

    node_enc = NodeEncoder(d_plm=320, d_struct=3, d_laplacian=8, d_node=128).to(device)
    node_enc.eval()
    with torch.no_grad():
        h = node_enc(
            enc_input.x.to(device),
            enc_input.structural_features.to(device),
            enc_input.lap_pe.to(device),
        )

    attn_mask = build_temporal_attention_mask(node_ids, node_times).to(device)

    torch.save({
        "node_ids": node_ids,
        "node_embeddings": h.cpu(),
        "edge_index": enc_input.edge_index,
        "edge_attr": enc_input.edge_attr,
        "edge_type": enc_input.edge_type,
        "attn_mask": attn_mask.cpu(),
        "node_times": torch.tensor([node_times[nid] for nid in node_ids]),
        "root_index": enc_input.root_index,
    }, out)
    log(g, f"embedded {h.shape[0]} nodes → {out.name}")
    return out


# ── per-group entry point (runs in worker process) ────────────────────────────

STAGES = ["clean", "align", "fasttree", "refine", "anc", "translate", "embed"]


def run_group(g: int, stop_after: str | None) -> str:
    try:
        def stop(stage):
            return stop_after is not None and stage == stop_after

        clean, meta_csv = stage_clean(g);
        if stop("clean"):     return f"{g:03d} stopped after clean"
        aligned = stage_align(g, clean, meta_csv);
        if stop("align"):     return f"{g:03d} stopped after align"
        treefile = stage_fasttree(g, aligned);
        if stop("fasttree"):  return f"{g:03d} stopped after fasttree"
        rooted, bl_json = stage_refine(g, treefile, aligned, meta_csv);
        if stop("refine"):    return f"{g:03d} stopped after refine"
        anc_nt  = stage_ancestral(g, rooted, aligned);
        if stop("anc"):       return f"{g:03d} stopped after anc"
        anc_aa  = stage_translate(g, anc_nt);
        if stop("translate"): return f"{g:03d} stopped after translate"
        stage_embed(g, rooted, anc_aa, bl_json)

        return f"{g:03d} done"
    except Exception as e:
        return f"{g:03d} FAILED: {e}\n{traceback.format_exc()}"


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", nargs="+", type=int,
                        help="group numbers (default: 1–48)")
    parser.add_argument("--workers", type=int, default=4,
                        help="parallel workers (default: 4)")
    parser.add_argument("--stop-after", choices=STAGES,
                        help="stop each group after this stage")
    args = parser.parse_args()

    groups = args.groups or list(range(1, 49))
    stop_after = args.stop_after

    print(f"Processing {len(groups)} groups with {args.workers} workers")
    print(f"Stop after: {stop_after or 'all stages'}\n")

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_group, g, stop_after): g for g in groups}
        for fut in as_completed(futures):
            print(fut.result(), flush=True)

    print("\nAll done.")


if __name__ == "__main__":
    main()
