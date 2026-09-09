#!/usr/bin/env python3
"""Post-process ab_clones_* split dirs into TreeDataset layout.

For each group_XXX.{nwk,fasta,anc_aa.fasta} under {train,val,test}:
  1. Midpoint-rooted Newick → group_XXX_rooted.nwk with unique names
     (root=NODE_ROOT; unnamed/bootstrap internals get NODE_XXXXXXX).
     FastTree support values are stripped so BioPython cannot concatenate
     them onto labels (NODE_00000040.92).
  2. Expand group_XXX_anc_aa.fasta: root + all leaf AAs (IDs = parse_newick names);
     missing internals filled by parent-sequence copy (not gaps — gaps poison val).
  3. Write group_XXX_bl.json with numdate = cumulative Newick path length from
     the SAME names src.dataset.parse_newick will use.

Does NOT re-run ANARCI / MAFFT / FastTree.

Usage:
  python scripts/ab_postprocess_clones_for_treesbm.py --data data/ab_clones_1m
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.bridge.sample_bridge_state import sample_bridge_state  # noqa: E402
from src.dataset import (  # noqa: E402
    _cumulative_numdates,
    fill_missing_node_seqs,
    parse_newick,
)
from src.tree_state import TreeState  # noqa: E402
from src.treeencoder.structural_features import compute_depth  # noqa: E402


def _looks_like_support(name: str | None) -> bool:
    if not name:
        return False
    try:
        float(name)
        return True
    except ValueError:
        return False


def _strip_support(tree: Phylo.BaseTree.Tree) -> None:
    """Drop FastTree bootstrap so Phylo.write cannot glue it onto labels."""
    for clade in tree.find_clades():
        clade.confidence = None
        if _looks_like_support(clade.name):
            clade.name = None


def _name_tree_nodes(tree: Phylo.BaseTree.Tree) -> str:
    """Assign unique names matching src.dataset.parse_newick; return root id."""
    counter = [0]

    def name(clade):
        if clade.name and not _looks_like_support(clade.name):
            return clade.name
        n = f"NODE_{counter[0]:07d}"
        clade.name = n
        counter[0] += 1
        return n

    tree.root.name = "NODE_ROOT"
    seen: set[str] = set()

    def walk(parent):
        pn = name(parent)
        if pn in seen and parent is not tree.root:
            # duplicate (e.g. two internals both named after write/read)
            n = f"NODE_{counter[0]:07d}"
            counter[0] += 1
            parent.name = n
            pn = n
        seen.add(pn)
        for child in parent.clades:
            walk(child)

    walk(tree.root)
    tree.rooted = True
    return tree.root.name


def _read_fasta_dict(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).replace("-", "").replace(".", "") for rec in SeqIO.parse(path, "fasta")}


def _verify_group(split_dir: Path, g: int) -> None:
    """TreeDataset-equivalent load + compute_depth + bridge smoke."""
    d = split_dir
    root_id, node_ids, edges, branch_lengths = parse_newick(str(d / f"group_{g:03d}_rooted.nwk"))
    seqs = {
        rec.id: str(rec.seq)
        for rec in SeqIO.parse(d / f"group_{g:03d}_anc_aa.fasta", "fasta")
    }
    seqs = fill_missing_node_seqs(root_id, edges, seqs)
    seqs = {nid: seqs[nid] for nid in node_ids}
    with open(d / f"group_{g:03d}_bl.json") as f:
        node_data = json.load(f)["nodes"]
    path_times = _cumulative_numdates(root_id, edges, branch_lengths)
    times = {
        nid: float(node_data[nid]["numdate"]) if nid in node_data and "numdate" in node_data[nid]
        else float(path_times.get(nid, 0.0))
        for nid in node_ids
    }
    tree = TreeState(
        node_ids=node_ids,
        root_id=root_id,
        edges=edges,
        branch_lengths=branch_lengths,
        node_seqs=seqs,
        active_leaves=[nid for nid in node_ids],
    )
    compute_depth(tree)
    missing = [nid for nid in node_ids if nid not in node_data]
    if missing:
        raise RuntimeError(f"bl.json missing {len(missing)} parse_newick ids e.g. {missing[:3]}")
    for t in (0.0, 0.25, 0.5, 1.0):
        T = sample_bridge_state(t, node_ids, times, edges, branch_lengths, seqs, root_id)
        tree_t = TreeState(
            node_ids=T["node_ids_t"],
            root_id=root_id,
            edges=T["edges_t"],
            branch_lengths=T["branch_lengths_t"],
            node_seqs=T["seqs_t"],
            active_leaves=T["active_leaves_t"],
        )
        compute_depth(tree_t)


def process_group(split_dir: Path, g: int) -> dict:
    stem = f"group_{g:03d}"
    nwk_in = split_dir / f"{stem}.nwk"
    leaf_fa = split_dir / f"{stem}.fasta"
    anc_in = split_dir / f"{stem}_anc_aa.fasta"
    rooted = split_dir / f"{stem}_rooted.nwk"
    anc_out = split_dir / f"{stem}_anc_aa.fasta"
    bl_out = split_dir / f"{stem}_bl.json"

    if not nwk_in.exists():
        raise FileNotFoundError(nwk_in)
    if not leaf_fa.exists():
        raise FileNotFoundError(leaf_fa)

    tree = Phylo.read(str(nwk_in), "newick")
    _strip_support(tree)
    _name_tree_nodes(tree)
    _strip_support(tree)
    Phylo.write(tree, str(rooted), "newick")

    # Re-parse with the same function TreeDataset uses, then write bl/anc
    # from those names so roundtrip cannot desync labels.
    root_id, node_ids, edges, branch_lengths = parse_newick(str(rooted))

    root_seq = ""
    if anc_in.exists():
        for rec in SeqIO.parse(anc_in, "fasta"):
            if rec.id in ("NODE_ROOT", "ROOT", root_id) or root_seq == "":
                root_seq = str(rec.seq).replace("-", "").replace(".", "")
                if rec.id in ("NODE_ROOT", "ROOT", root_id):
                    break
    root_fa = split_dir / f"{stem}_root.fasta"
    if not root_seq and root_fa.exists():
        root_seq = str(next(SeqIO.parse(root_fa, "fasta")).seq).replace("-", "").replace(".", "")
    if not root_seq:
        raise RuntimeError(f"{stem}: missing root AA")

    leaves = _read_fasta_dict(leaf_fa)
    has_children = {p for p, _ in edges}
    tips = {nid for nid in node_ids if nid not in has_children and nid != root_id}
    missing = sorted(tips - set(leaves))
    extra = sorted(set(leaves) - tips)
    if missing:
        raise RuntimeError(f"{stem}: tip IDs missing from leaf fasta: {missing[:5]}")

    records = [SeqRecord(Seq(root_seq), id=root_id, description="")]
    for tip in sorted(tips):
        records.append(SeqRecord(Seq(leaves[tip]), id=tip, description=""))
    # Parent-propagate so anc_aa has every parse_newick node (no gap placeholders).
    tip_and_root = {rec.id: str(rec.seq) for rec in records}
    filled = fill_missing_node_seqs(root_id, edges, tip_and_root)
    records = [
        SeqRecord(Seq(filled[nid]), id=nid, description="")
        for nid in node_ids
    ]
    SeqIO.write(records, str(anc_out), "fasta")

    path_lens = _cumulative_numdates(root_id, edges, branch_lengths)
    nodes = {}
    parent = {c: p for p, c in edges}
    for nid in node_ids:
        bl = 0.0 if nid == root_id else max(0.0, float(branch_lengths.get((parent[nid], nid), 0.0)))
        nodes[nid] = {"numdate": float(path_lens.get(nid, 0.0)), "branch_length": bl}
    payload = {
        "generated_by": "scripts/ab_postprocess_clones_for_treesbm.py",
        "input_tree": str(nwk_in),
        "note": "numdate = cumulative Newick path length from root (genetic distance, not calendar time); keys match parse_newick",
        "root_id": root_id,
        "nodes": nodes,
    }
    bl_out.write_text(json.dumps(payload, indent=2) + "\n")

    _verify_group(split_dir, g)

    nums = [v["numdate"] for v in nodes.values()]
    return {
        "group": g,
        "n_tips": len(tips),
        "n_nodes": len(nodes),
        "n_anc": len(records),
        "extra_leaf_fasta_ids": len(extra),
        "numdate_max": max(nums) if nums else 0.0,
        "root_len": len(root_seq),
        "root_id": root_id,
    }


def discover_groups(split_dir: Path) -> list[int]:
    groups = []
    for p in sorted(split_dir.glob("group_*.nwk")):
        name = p.name
        if name.endswith("_unrooted.nwk") or name.endswith("_rooted.nwk"):
            continue
        stem = p.stem
        parts = stem.split("_")
        if len(parts) == 2 and parts[0] == "group" and parts[1].isdigit():
            groups.append(int(parts[1]))
    return groups


def _quarantine(split_dir: Path, g: int, broken_dir: Path) -> None:
    broken_dir.mkdir(parents=True, exist_ok=True)
    stem = f"group_{g:03d}"
    for suffix in ("_rooted.nwk", "_anc_aa.fasta", "_bl.json"):
        src = split_dir / f"{stem}{suffix}"
        if src.exists():
            shutil.move(str(src), str(broken_dir / src.name))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, required=True, help="e.g. data/ab_clones_1m")
    ap.add_argument("--splits", nargs="+", default=["train", "val", "test"])
    args = ap.parse_args()

    summary: dict = {"data": str(args.data), "splits": {}}
    for split in args.splits:
        d = args.data / split
        if not d.is_dir():
            raise SystemExit(f"missing split dir: {d}")
        broken_dir = d / "broken"
        groups = discover_groups(d)
        ok, fail = [], []
        for g in groups:
            try:
                ok.append(process_group(d, g))
            except Exception as e:
                fail.append({"group": g, "error": str(e)})
                print(f"FAIL {split}/group_{g:03d}: {e}", flush=True)
                try:
                    _quarantine(d, g, broken_dir)
                except Exception as qe:
                    print(f"  quarantine failed: {qe}", flush=True)
        n_rooted = len(list(d.glob("group_*_rooted.nwk")))
        n_anc = len(list(d.glob("group_*_anc_aa.fasta")))
        n_bl = len(list(d.glob("group_*_bl.json")))
        complete = sorted(
            int(p.stem.split("_")[1])
            for p in d.glob("group_*_rooted.nwk")
            if (d / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists()
            and (d / p.name.replace("_rooted.nwk", "_bl.json")).exists()
        )
        summary["splits"][split] = {
            "n_input_nwk": len(groups),
            "n_ok": len(ok),
            "n_fail": len(fail),
            "n_rooted": n_rooted,
            "n_anc": n_anc,
            "n_bl": n_bl,
            "n_complete": len(complete),
            "n_dropped": len(fail),
            "fails": fail[:20],
            "numdate_max_median": (
                sorted(r["numdate_max"] for r in ok)[len(ok) // 2] if ok else None
            ),
        }
        print(
            f"{split}: ok={len(ok)} fail={len(fail)} complete={len(complete)} "
            f"(rooted={n_rooted} anc={n_anc} bl={n_bl})",
            flush=True,
        )

    out = args.data / "postprocess_treesbm_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
