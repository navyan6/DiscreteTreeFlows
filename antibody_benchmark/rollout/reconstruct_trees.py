"""Reconstruct AntibodyTree objects from DASM/Thrifty-style PCP tables.

Root policy (no fabricated edges / no similarity ancestry):
1. Prefer an explicit inferred naive/germline node (parent_is_naive) as root.
2. Else, if the published PCP already has a unique topological root, keep it
   (Zenodo v3 files are upstream ``no-naive``; we do not strip further).
3. If multi-root / no usable single root: exclude the family.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Iterable, Optional

import pandas as pd

from antibody_benchmark.codon import looks_like_nt, translate_nt
from antibody_benchmark.trees import AntibodyTree, Edge, validate_tree


ATTRITION_STAGES = [
    "raw_families",
    "after_root_resolution",
    "after_topology_validation",
    "after_min_leaves",
    "after_min_depth",
    "after_productive_root",
    "final_eligible",
]


def _seq_cols(df: pd.DataFrame, chain: str) -> tuple[str, str]:
    if chain == "paired":
        if "parent_heavy" not in df.columns or "parent_light" not in df.columns:
            raise ValueError("paired chain requested but HL columns missing")
        return "parent_heavy", "child_heavy"  # heavy used for topology seqs; light in metadata
    if "parent_heavy" in df.columns:
        return "parent_heavy", "child_heavy"
    if "parent_light" in df.columns:
        return "parent_light", "child_light"
    if "parent" in df.columns:
        return "parent", "child"
    raise ValueError(f"No parent/child sequence columns in {list(df.columns)}")


def _truthy(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    return str(v).strip().lower() in {"true", "1", "t", "yes"}


def _resolve_root(
    nodes: set[str],
    edges: list[Edge],
    naive_nodes: set[str],
) -> tuple[Optional[str], str, str]:
    """Return (root_id or None, exclusion_reason or '', root_source).

    root_source in {naive, topological, ""}.
    """
    indeg: dict[str, int] = defaultdict(int)
    for _p, c in edges:
        indeg[c] += 1
    topo_roots = [n for n in nodes if indeg[n] == 0]
    naive_roots = [n for n in topo_roots if n in naive_nodes]

    if naive_roots:
        if len(naive_roots) == 1:
            return naive_roots[0], "", "naive"
        return None, "ambiguous_naive_roots", ""
    if len(topo_roots) == 1:
        # Published Zenodo PCPs are already no-naive; unique ancestral endpoint
        # is retained for root-conditioned rollout (not fabricated).
        return topo_roots[0], "", "topological"
    if len(topo_roots) == 0:
        return None, "zero_root", ""
    return None, "multi_root_no_usable_naive", ""


def reconstruct_trees_from_pcp_df(
    df: pd.DataFrame,
    *,
    chain: str = "heavy",
    require_single_root: bool = True,
    min_leaves: int = 1,
    min_depth: int = 1,
    translate_aa: bool = True,
    max_families: Optional[int] = None,
    return_attrition: bool = False,
):
    """Group PCP rows by (sample_id, family) and build validated trees.

    If ``return_attrition=True``, returns ``(trees, attrition_df, reason_counts)``.
    """
    parent_col, child_col = _seq_cols(df, chain)
    trees: list[AntibodyTree] = []
    attrition_rows: list[dict] = []
    reason_counts: Counter = Counter()
    stage_counts: Counter = Counter()

    grouped = df.groupby(["sample_id", "family"], sort=False)

    for (sample_id, family), g in grouped:
        family_id = f"{sample_id}::{family}"
        stage_counts["raw_families"] += 1
        row = {
            "family_id": family_id,
            "donor_id": str(sample_id),
            "family": str(family),
            "n_pcp_rows": int(len(g)),
            "excluded": False,
            "exclusion_reason": "",
            "root_source": "",
            "n_topo_roots": 0,
            "n_naive_flagged": 0,
            "n_leaves": 0,
            "max_depth": 0,
            "passed_stage": "raw_families",
        }

        edges: list[Edge] = []
        bl: dict[Edge, float] = {}
        seqs: dict[str, str] = {}
        is_leaf: dict[str, bool] = {}
        naive_nodes: set[str] = set()
        max_depth = 0
        light_seqs: dict[str, str] = {}

        for _, r in g.iterrows():
            p = str(r["parent_name"])
            c = str(r["child_name"])
            e = (p, c)
            edges.append(e)
            t = float(r["branch_length"])
            bl[e] = t
            seqs[p] = str(r[parent_col])
            seqs[c] = str(r[child_col])
            if "parent_light" in g.columns and chain == "paired":
                light_seqs[p] = str(r["parent_light"])
                light_seqs[c] = str(r["child_light"])
            is_leaf[c] = _truthy(r.get("child_is_leaf", False))
            if _truthy(r.get("parent_is_naive", False)):
                naive_nodes.add(p)
            if "depth" in r and pd.notna(r["depth"]):
                max_depth = max(max_depth, int(r["depth"]))

        nodes = set(seqs)
        indeg: dict[str, int] = defaultdict(int)
        for p, c in edges:
            indeg[c] += 1
        topo_roots = [n for n in nodes if indeg[n] == 0]
        row["n_topo_roots"] = len(topo_roots)
        row["n_naive_flagged"] = len(naive_nodes)
        row["max_depth"] = max_depth

        root_id, excl, root_source = _resolve_root(nodes, edges, naive_nodes)
        row["root_source"] = root_source
        if excl:
            if require_single_root or excl != "multi_root_no_usable_naive":
                row["excluded"] = True
                row["exclusion_reason"] = excl
                reason_counts[excl] += 1
                attrition_rows.append(row)
                continue
            # require_single_root=False: keep first topo root (debug only)
            if not topo_roots:
                row["excluded"] = True
                row["exclusion_reason"] = excl
                reason_counts[excl] += 1
                attrition_rows.append(row)
                continue
            root_id = sorted(topo_roots)[0]
            root_source = "topological_unvalidated"
            row["root_source"] = root_source

        stage_counts["after_root_resolution"] += 1
        row["passed_stage"] = "after_root_resolution"

        parents = {p for p, _ in edges}
        for n in nodes:
            if n not in is_leaf:
                is_leaf[n] = n not in parents

        aa = {}
        if translate_aa:
            for n, s in seqs.items():
                aa[n] = translate_nt(s) if looks_like_nt(s) else s

        tree = AntibodyTree(
            family_id=family_id,
            root_id=root_id,  # type: ignore[arg-type]
            node_ids=sorted(nodes),
            edges=edges,
            branch_lengths=bl,
            true_sequences=seqs,
            is_leaf=is_leaf,
            donor_id=str(sample_id),
            true_aa_sequences=aa,
            metadata={
                "family": str(family),
                "max_depth": max_depth,
                "n_edges": len(edges),
                "root_marked_naive": root_id in naive_nodes,
                "root_source": root_source,
                "chain": chain,
                "light_seqs": light_seqs,
            },
        )
        ok, reason = validate_tree(tree)
        if not ok:
            row["excluded"] = True
            row["exclusion_reason"] = reason
            reason_counts[reason] += 1
            attrition_rows.append(row)
            continue

        stage_counts["after_topology_validation"] += 1
        row["passed_stage"] = "after_topology_validation"

        n_leaves = sum(1 for n in nodes if is_leaf[n])
        row["n_leaves"] = n_leaves
        if n_leaves < min_leaves:
            row["excluded"] = True
            row["exclusion_reason"] = "lt_min_leaves"
            reason_counts["lt_min_leaves"] += 1
            attrition_rows.append(row)
            continue
        stage_counts["after_min_leaves"] += 1
        row["passed_stage"] = "after_min_leaves"

        if max_depth < min_depth:
            row["excluded"] = True
            row["exclusion_reason"] = "lt_min_depth"
            reason_counts["lt_min_depth"] += 1
            attrition_rows.append(row)
            continue
        stage_counts["after_min_depth"] += 1
        row["passed_stage"] = "after_min_depth"

        if translate_aa and "*" in translate_nt(seqs[root_id], stop_as="*"):
            row["excluded"] = True
            row["exclusion_reason"] = "stop_in_root"
            reason_counts["stop_in_root"] += 1
            attrition_rows.append(row)
            continue
        stage_counts["after_productive_root"] += 1
        row["passed_stage"] = "after_productive_root"

        trees.append(tree)
        stage_counts["final_eligible"] += 1
        row["passed_stage"] = "final_eligible"
        attrition_rows.append(row)

        if max_families is not None and len(trees) >= max_families:
            # Still finish attrition for remaining families as truncated? Prefer full audit.
            # Continue scanning for attrition completeness.
            pass

    # If max_families set, truncate returned trees but keep full attrition.
    if max_families is not None:
        trees = trees[: int(max_families)]

    attrition_df = pd.DataFrame(attrition_rows)
    if return_attrition:
        return trees, attrition_df, dict(reason_counts), dict(stage_counts)
    return trees


def summarize_trees(trees: Iterable[AntibodyTree]) -> dict:
    trees = list(trees)
    if not trees:
        return {
            "n_donors": 0,
            "n_families": 0,
            "n_nodes": 0,
            "n_leaves": 0,
            "n_edges": 0,
        }
    import statistics as stats

    donors = {t.donor_id for t in trees}
    n_nodes = sum(len(t.node_ids) for t in trees)
    n_leaves = sum(len(t.leaves()) for t in trees)
    n_edges = sum(len(t.edges) for t in trees)
    depths = [int(t.metadata.get("max_depth", 0)) for t in trees]
    bls = [bl for t in trees for bl in t.branch_lengths.values()]
    rtl = []
    root_sources = Counter(str(t.metadata.get("root_source", "")) for t in trees)
    for t in trees:
        root_aa = t.true_aa_sequences.get(t.root_id) or t.true_sequences[t.root_id]
        for leaf in t.leaves():
            leaf_aa = t.true_aa_sequences.get(leaf) or t.true_sequences[leaf]
            if len(root_aa) == len(leaf_aa):
                rtl.append(sum(a != b for a, b in zip(root_aa, leaf_aa)))
    return {
        "n_donors": len(donors),
        "n_families": len(trees),
        "n_nodes": n_nodes,
        "n_leaves": n_leaves,
        "n_edges": n_edges,
        "mean_depth": stats.mean(depths) if depths else None,
        "median_depth": stats.median(depths) if depths else None,
        "mean_branch_length": stats.mean(bls) if bls else None,
        "median_branch_length": stats.median(bls) if bls else None,
        "mean_root_to_leaf_aa_hamming": stats.mean(rtl) if rtl else None,
        "median_root_to_leaf_aa_hamming": stats.median(rtl) if rtl else None,
        "root_source_counts": dict(root_sources),
    }


def format_attrition_summary(
    stage_counts: dict,
    reason_counts: dict,
    *,
    source_pcp: str,
    n_final: int,
) -> str:
    lines = [
        "# Antibody benchmark filter summary",
        "",
        f"**Source PCP:** `{source_pcp}`",
        "",
        "## Attrition by stage",
        "",
        "| Stage | Count |",
        "|---|---:|",
    ]
    for s in ATTRITION_STAGES:
        if s in stage_counts:
            lines.append(f"| {s} | {stage_counts[s]} |")
    lines += [
        "",
        "## Exclusion reason counts",
        "",
        "| Reason | Count |",
        "|---|---:|",
    ]
    for reason, n in sorted(reason_counts.items(), key=lambda x: (-x[1], x[0])):
        lines.append(f"| {reason} | {n} |")
    lines += [
        "",
        f"**Final eligible families:** {n_final}",
        "",
        "## Root policy",
        "",
        "- Prefer explicit `parent_is_naive` germline/naive node when present.",
        "- Published DASM Zenodo v3 PCPs are upstream `no-naive` (naive edges already removed); "
        "we do **not** strip further and do **not** fabricate naive edges.",
        "- Multi-root forests with no usable single naive/root are **excluded**.",
        "",
    ]
    return "\n".join(lines) + "\n"
