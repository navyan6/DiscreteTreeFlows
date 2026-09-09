#!/usr/bin/env python3
"""
Regression test for the alignment-aware ``positional_recovery``.

The fix is only trustworthy if it is inert where it should be. Three checks:

  1. Synthetic. A leaf with a known deletion is scored against a root and a
     generated leaf built to match it at homologous residues. Legacy scoring
     must fail; aligned scoring must be perfect. This is the only check that
     proves the direction of the fix rather than just its stability.

  2. Gamma (Brazil g003) — the Figure 4/5 case study. P.1 carries no Spike
     deletion and the root is 1273 aa, so aligned and legacy scoring must agree.
     **If Gamma moves, the fix is wrong.**

  3. A deletion-carrying tree, for contrast: the same comparison must move.

Usage:
    python scripts/test_positional_recovery_frame.py                # 1 only
    python scripts/test_positional_recovery_frame.py --with-data    # 1+2+3
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmarks.metrics.sequences import (  # noqa: E402
    align_index_map,
    positional_recovery,
)

TOL = 1e-9


def iter_fasta(path: Path):
    name, buf = None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def is_leaf(name: str) -> bool:
    """Handles both conventions: augur's NODE_* internals and the generator's
    ``<path>|<leaf|internal|root>`` suffix."""
    if "|" in name:
        return name.rsplit("|", 1)[1] == "leaf"
    return not name.startswith("NODE_") and name != "root"


# ── 1. synthetic ────────────────────────────────────────────────────────────

def test_synthetic() -> bool:
    """Deletion at 30-34; the true mutation sits downstream of it.

    The root is pseudo-random rather than a homopolymer: with low-complexity
    sequence the aligner can park the gap at a free terminus instead of over the
    real deletion, which says nothing about real protein data.
    """
    import random

    rng = random.Random(0)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    root = "".join(rng.choice(aa) for _ in range(200))

    # leaf drops root[30:35] and carries a substitution at root index 95
    sub = next(c for c in aa if c != root[95])
    leaf = root[:30] + root[35:95] + sub + root[96:]
    # a generated leaf that got it exactly right, in the root's frame
    gen = root[:95] + sub + root[96:]

    assert len(leaf) == len(root) - 5
    assert len(gen) == len(root), "generated leaf must inherit the root's length"

    ok = True

    m = align_index_map(leaf, root)
    if not all(m[i] == -1 for i in range(30, 35)):
        print("  FAIL: deleted block not mapped to -1")
        ok = False
    if m[95] != 90:
        print(f"  FAIL: root index 95 should map to leaf index 90, got {m[95]}")
        ok = False
    if any(m[i] != i for i in range(30)):
        print("  FAIL: positions upstream of the deletion should be unshifted")
        ok = False

    new = positional_recovery(root, leaf, gen, align=True)
    old = positional_recovery(root, leaf, gen, align=False)

    if abs(new["mut_recovery"] - 1.0) > TOL or abs(new["cons_retention"] - 1.0) > TOL:
        print(f"  FAIL: aligned scoring should be perfect, got "
              f"mut_recovery={new['mut_recovery']:.4f} "
              f"cons_retention={new['cons_retention']:.4f}")
        ok = False
    if new["n_skipped_indel"] != 5:
        print(f"  FAIL: expected 5 skipped positions, got {new['n_skipped_indel']}")
        ok = False
    if old["mut_recovery"] > 0.5:
        print(f"  FAIL: legacy scoring should be badly wrong, got "
              f"mut_recovery={old['mut_recovery']:.4f}")
        ok = False

    print(f"  legacy : mut_recovery={old['mut_recovery']:.4f} "
          f"cons_retention={old['cons_retention']:.4f}")
    print(f"  aligned: mut_recovery={new['mut_recovery']:.4f} "
          f"cons_retention={new['cons_retention']:.4f} "
          f"skipped={new['n_skipped_indel']}")
    return ok


# ── 2/3. real trees ─────────────────────────────────────────────────────────

METRICS = ("mut_recovery", "cons_retention", "aa_acc_given_hit")


def load_tree(obs_fasta: Path, gen_fasta: Path, max_gt: int, max_gen: int):
    obs = list(iter_fasta(obs_fasta))
    internals = {n: s for n, s in obs if not is_leaf(n)}
    gt_leaves = [s for n, s in obs if is_leaf(n)][:max_gt]
    root = internals.get("NODE_0000000") or next(iter(internals.values()))
    gen_leaves = [s for n, s in iter_fasta(gen_fasta) if is_leaf(n)][:max_gen]
    if not gt_leaves or not gen_leaves:
        raise SystemExit(f"no leaves in {obs_fasta} / {gen_fasta}")
    return root, gt_leaves, gen_leaves


def score_tree(obs_fasta: Path, gen_fasta: Path, max_gt: int, max_gen: int) -> dict:
    """Pair each GT leaf with its closest generated leaf, scored both ways.

    Reports the aggregate movement, and separately the movement restricted to
    in-frame GT leaves. The in-frame subset is the real regression invariant:
    those leaves take the identity map by construction, so *any* difference there
    means the fix changed something it had no business touching.
    """
    root, gt_leaves, gen_leaves = load_tree(obs_fasta, gen_fasta, max_gt, max_gen)

    vals: dict[str, dict[str, list[float]]] = {
        m: {"old": [], "new": []} for m in METRICS
    }
    inframe_mismatches = 0
    n_inframe = 0
    skipped = []
    for gt in gt_leaves:
        # closest generated leaf, chosen by a frame-free distance so the choice
        # itself cannot bias the comparison
        best = min(gen_leaves, key=lambda g: sum(
            1 for a, b in zip(g, gt) if a != b) + abs(len(g) - len(gt)))
        old = positional_recovery(root, gt, best, align=False)
        new = positional_recovery(root, gt, best, align=True)
        for m in METRICS:
            for tag, r in (("old", old), ("new", new)):
                if r[m] == r[m]:  # drop nan
                    vals[m][tag].append(r[m])
        skipped.append(new["n_skipped_indel"])

        if len(gt) == len(root) and len(best) == len(root):
            n_inframe += 1
            for m in METRICS:
                o, n = old[m], new[m]
                same = (o != o and n != n) or (o == o and n == n and abs(o - n) < TOL)
                if not same:
                    inframe_mismatches += 1

    out = {
        "root_len": len(root),
        "n_gt": len(gt_leaves),
        "n_gen": len(gen_leaves),
        "frac_gt_same_len_as_root": round(
            sum(len(s) == len(root) for s in gt_leaves) / len(gt_leaves), 4),
        "n_inframe_pairs": n_inframe,
        "inframe_mismatches": inframe_mismatches,
        "mean_skipped_indel": round(mean(skipped), 2) if skipped else 0.0,
    }
    for m in METRICS:
        out[m] = {
            k: (round(mean(v), 6) if v else None) for k, v in vals[m].items()
        }
        out[f"{m}_n"] = {k: len(v) for k, v in vals[m].items()}
    return out


def report(label: str, res: dict, expect_stable: bool) -> bool:
    print(f"\n  {label}")
    print(f"    root {res['root_len']} aa | {res['n_gt']} GT x {res['n_gen']} gen | "
          f"{res['frac_gt_same_len_as_root']:.1%} of GT leaves share the root's length")
    print(f"    {res['mean_skipped_indel']} positions skipped per pair on average")

    for metric in METRICS:
        o, n = res[metric]["old"], res[metric]["new"]
        no, nn = res[f"{metric}_n"]["old"], res[f"{metric}_n"]["new"]
        if o is None and n is None:
            continue
        os_ = "  n/a " if o is None else f"{o:.4f}"
        ns_ = "  n/a " if n is None else f"{n:.4f}"
        note = "" if no == nn else f"   [scored on {no} -> {nn} leaves]"
        print(f"    {metric:<18} legacy {os_} -> aligned {ns_}{note}")

    ok = res["inframe_mismatches"] == 0
    print(f"    in-frame invariant: {res['n_inframe_pairs']} pairs, "
          f"{res['inframe_mismatches']} mismatches -> {'PASS' if ok else 'FAIL'}")

    if not expect_stable:
        moved = any(
            res[m]["old"] is not None and res[m]["new"] is not None
            and abs(res[m]["new"] - res[m]["old"]) > 5e-3
            for m in METRICS
        ) or res["mean_skipped_indel"] > 0
        print(f"    movement on off-frame tree: {'PASS' if moved else 'FAIL'}")
        ok &= moved
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--with-data", action="store_true",
                    help="also run the Gamma and deletion-lineage checks")
    ap.add_argument("--obs-dir", type=Path, default=ROOT / "data" / "covid" / "test")
    ap.add_argument("--gen-dir", type=Path,
                    default=ROOT / "data" / "generated" / "covid" / "matched")
    ap.add_argument("--gamma-group", default="003")
    ap.add_argument("--contrast-group", default=None,
                    help="a deletion-carrying group; auto-picked if omitted")
    ap.add_argument("--max-gt", type=int, default=150)
    ap.add_argument("--max-gen", type=int, default=150)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "voc_threat_panel"
                    / "positional_recovery_regression.json")
    args = ap.parse_args()

    print("1. synthetic deletion")
    results = {"synthetic_pass": test_synthetic()}
    all_ok = results["synthetic_pass"]

    if args.with_data:
        def paths(group: str) -> tuple[Path, Path]:
            return (args.obs_dir / f"group_{group}_anc_aa.fasta",
                    args.gen_dir / f"group_{group}_generated_matched.fasta")

        print("\n2. Gamma g003 (Figure 4/5 control — must not move)")
        obs, gen = paths(args.gamma_group)
        gamma = score_tree(obs, gen, args.max_gt, args.max_gen)
        ok = report(f"Gamma group_{args.gamma_group}", gamma, expect_stable=True)
        results["gamma"] = gamma
        all_ok &= ok

        # contrast: the available group with the most off-frame leaves. Skip
        # groups whose root is not even close to the leaf length -- those are a
        # data defect, not a frame shift, and are reported separately below.
        group = args.contrast_group
        anomalies = []
        if group is None:
            best, best_frac = None, 1.0
            for g in sorted(args.gen_dir.glob("group_*_generated_matched.fasta")):
                gid = g.name.split("_")[1]
                o = args.obs_dir / f"group_{gid}_anc_aa.fasta"
                if not o.exists():
                    continue
                recs = list(iter_fasta(o))
                internals = {n: s for n, s in recs if not is_leaf(n)}
                lv = [s for n, s in recs if is_leaf(n)]
                if not lv or not internals:
                    continue
                rt = internals.get("NODE_0000000") or next(iter(internals.values()))
                modal = Counter(len(s) for s in lv).most_common(1)[0][0]
                if abs(len(rt) - modal) > 0.1 * modal:
                    anomalies.append({"group": gid, "root_len": len(rt),
                                      "modal_leaf_len": modal})
                    continue
                if gid == args.gamma_group:
                    continue
                frac = sum(len(s) == len(rt) for s in lv) / len(lv)
                if frac < best_frac:
                    best, best_frac = gid, frac
            group = best

        if anomalies:
            results["root_length_anomalies"] = anomalies
            print("\n  ! root-length anomalies (unrelated to the frame fix, "
                  "but these trees cannot support any sequence metric):")
            for a in anomalies:
                print(f"      group_{a['group']}: root is {a['root_len']} aa "
                      f"while its leaves are {a['modal_leaf_len']} aa")

        if group:
            print("\n3. deletion-carrying contrast (must move)")
            obs, gen = paths(group)
            contrast = score_tree(obs, gen, args.max_gt, args.max_gen)
            ok = report(f"group_{group}", contrast, expect_stable=False)
            results["contrast"] = {"group": group, **contrast}
            all_ok &= ok
        else:
            print("\n3. skipped — no contrast group available")

        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwrote {args.out}")

    print("\n" + ("ALL CHECKS PASSED" if all_ok else "FAILURES ABOVE"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
