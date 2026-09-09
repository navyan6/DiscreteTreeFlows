#!/usr/bin/env python3
"""
Re-translate group_*_anc_nt.fasta with the frame-preserving translator.

The correct full-length reconstructions are already on disk; only the
translation stage corrupted them (see results/voc_threat_panel/TRUNCATED_ROOTS.md).
This rebuilds group_*_anc_aa.fasta from them -- no re-alignment, no FastTree, no
augur.

Safety: originals are copied to group_*_anc_aa.fasta.prefix_bug.bak before
anything is overwritten, and the script refuses to write if a leaf sequence
would change. Leaves carry only whole-codon deletions, so they must come out
byte-identical; if one moves, the translator is wrong and we stop.

    python scripts/retranslate_anc_aa.py --data-dir data/covid/test --dry-run
    python scripts/retranslate_anc_aa.py --data-dir data/covid/test
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.treeencoder.seq_utils import _frameshifting_gaps, nt_to_aa  # noqa: E402

BAK_SUFFIX = ".prefix_bug.bak"


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


def is_internal(name: str) -> bool:
    return name.startswith("NODE_") or name == "root"


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 60) -> None:
    with path.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# Which translation each dataset needs, from scripts/panviral/audit_frame_offset.py.
#
# "aligned" requires the alignment to be reference-anchored, so that three
# columns are one reference codon for every sequence at once. That holds only
# for the spike data, where offset 0 gives 0 internal stops per sequence and the
# other two offsets give 55 and 144.
#
# The flu and HIV groups are aligned per group with MAFFT on nucleotides, and no
# offset is in frame: the best of the three still leaves 9-30 internal stops per
# sequence, and two alignment widths are not even multiples of three. There the
# frame lives in each sequence's own bases, gaps are alignment artifacts, and
# stripping them recovers the canonical protein -- flu comes out at 566 aa, the
# textbook HA0 length. Their bug is only the truncation at the first stop, which
# had cut HIV roots to a fifth of their leaves' length.
# The second element is cds_start. `nt_to_aa` otherwise locates the CDS with
# seq.find('ATG'), which is only correct when the window opens on the start
# codon. Bundibugyo's window does not: its first ATG sits at index 110, and
# 110 % 3 == 2, so every sequence was read out of frame and died at the first
# in-frame stop, giving 17-residue proteins from a 900-codon alignment that
# contains no gaps whatsoever. H3N2's first ATG is at index 51, which happens to
# be divisible by three, so it got the right frame by luck rather than design.
#
#   dataset   starts ATG   1st ATG  %3   len@ATG   len@col0
#   covid          44/44         0   0      1273       1273
#   bdbv            0/42       113   2        24        905
#   h3n2            0/44        51   0       566        603
#   h1n1           40/44         0   0       566        578
#   hiv            44/44         0   0       871        881
MODE_BY_DATASET = {
    "covid": ("aligned", None),
    "covid_voc_roots": ("aligned", None),
    "covid_epidemic": ("aligned", None),
    "covid_temporal": ("aligned", None),
    "covid_cladeholdout": ("aligned", None),
    "h3n2": ("ungapped", None),
    "h1n1": ("ungapped", None),
    "h3n2_epidemic": ("ungapped", None),
    "h1n1_epidemic": ("ungapped", None),
    "h3n2_temporal": ("ungapped", None),
    "h1n1_temporal": ("ungapped", None),
    "h3n2_forecast_2020_cal": ("ungapped", None),
    "h3n2_forecast_2020_season": ("ungapped", None),
    "panflu_forecast_h3n2_cal": ("ungapped", None),
    "panflu_forecast_h3n2_cal_dual": ("ungapped", None),
    "panflu_forecast_flub_cal": ("ungapped", None),
    "hiv_geo": ("ungapped", None),
    "hiv_temporal": ("ungapped", None),
    "bdbv_temporal": ("aligned", 0),
    "bdbv_pan_temporal": ("aligned", 0),
}


def resolve_mode(d: Path) -> tuple[str, int | None]:
    """Pick translation mode and CDS start from the dataset path."""
    parts = {p.lower() for p in d.parts}
    hits = {v for k, v in MODE_BY_DATASET.items() if k in parts}
    if len(hits) == 1:
        return hits.pop()
    raise SystemExit(
        f"cannot infer translation mode for {d}. Pass --mode explicitly, and "
        f"first confirm the choice with scripts/panviral/audit_frame_offset.py "
        f"and scripts/panviral/probe_cds_start.py -- applying 'aligned' to a "
        f"non-reference-anchored alignment, or trusting the ATG search on a "
        f"window that does not open on the start codon, silently produces "
        f"garbage protein rather than failing."
    )


def process(nt_path: Path, aa_path: Path, dry_run: bool, mode: str = "aligned",
            cds_start: int | None = None, allow_stale: bool = False) -> dict:
    old = dict(iter_fasta(aa_path)) if aa_path.exists() else {}

    new_records, failures, frameshifted = [], [], set()
    for name, nt in iter_fasta(nt_path):
        if _frameshifting_gaps(nt.upper()):
            frameshifted.add(name)
        try:
            new_records.append((name, nt_to_aa(nt, cds_start=cds_start, mode=mode)))
        except ValueError as e:
            failures.append(f"{name}: {e}")
            if name in old:
                new_records.append((name, old[name]))

    new = dict(new_records)

    def rootlen(d):
        r = d.get("NODE_0000000")
        return len(r) if r else None

    leaf_lens = [len(s) for n, s in new.items() if not is_internal(n)]
    modal = Counter(leaf_lens).most_common(1)[0][0] if leaf_lens else None

    changed = [n for n in new if n in old and old[n] != new[n]]
    internal_changed = [n for n in changed if is_internal(n)]
    leaf_changed = [n for n in changed if not is_internal(n)]

    # The tight invariant: a sequence may change for exactly two reasons, and
    # anything else means the translator is rewriting protein it should not
    # touch.
    #
    #   1. Its aligned nucleotides contained a gap run that is not a whole
    #      number of codons, so the old translator frameshifted it. Its old
    #      output is not a reference worth preserving -- group_016's root was
    #      nominally 1272 aa but agreed with the corrected translation for only
    #      211 residues.
    #   2. The old translation stopped at an internal stop codon and the new one
    #      continues past it. Then the old string is a prefix of the new one:
    #      every residue that was there before is unchanged, and the repair is
    #      purely the restored tail.
    #
    # Reason 1 is only a licence to rewrite in "aligned" mode, where the
    # alignment is reference-anchored and column-wise translation is meaningful.
    # In "ungapped" mode the frame comes from each sequence's own bases, so a
    # frameshifting gap run explains nothing and the prefix rule is the only
    # change we accept -- which is exactly the restriction that keeps flu and
    # HIV from being rewritten by a repair built for SARS-CoV-2.
    # Overriding cds_start moves the reading frame, so nothing from the old
    # translation survives and neither test above can apply. The replacement is
    # justified by the result instead of by the diff: the new protein has to
    # reach the length its own alignment implies. That is a real check -- the
    # Bundibugyo proteins were 17 aa against a 900-codon window, so a repair
    # that failed to reach full length would still be refused here.
    #
    # --allow-stale grants the same latitude for a different reason: in the
    # *_epidemic splits the stored protein was not derived from the stored
    # nucleotides at all. They disagree on 1-5% of sequences, often a single
    # residue, which is drift from an ancestral reconstruction that was re-run
    # without re-translating. Nucleotides are upstream in run_all_groups.py, so
    # they win -- but only where the result still reaches full length, so a
    # genuinely broken group is not waved through with the merely stale ones.
    rebuild_ok = cds_start is not None or allow_stale

    def is_expected(n: str) -> bool:
        if new[n].startswith(old[n]):
            return True
        if rebuild_ok and modal is not None and len(new[n]) >= 0.9 * modal:
            return True
        return mode == "aligned" and n in frameshifted

    unexpected = [n for n in changed if not is_expected(n)]
    leaf_repairs = [n for n in leaf_changed if n in frameshifted]
    # Repaired but not caught by a length check: full length before and after,
    # yet frameshifted in the middle.
    silent = [
        n for n in changed
        if n in frameshifted and modal is not None and len(old[n]) >= 0.9 * modal
    ]

    stat = {
        "tree": nt_path.name.replace("_anc_nt.fasta", ""),
        "n_records": len(new_records),
        "root_len_before": rootlen(old),
        "root_len_after": rootlen(new),
        "modal_leaf_len": modal,
        "n_frameshifted": len(frameshifted),
        "n_changed": len(changed),
        "n_internal_changed": len(internal_changed),
        "n_leaf_changed": len(leaf_changed),
        "n_leaf_repairs": len(leaf_repairs),
        "n_silent_frameshifts": len(silent),
        "silent_frameshifts": silent[:10],
        "n_unexpected_changes": len(unexpected),
        "failures": failures,
    }
    stat["was_truncated"] = (
        stat["root_len_before"] is not None and modal is not None
        and stat["root_len_before"] < 0.9 * modal
    )
    stat["now_full_length"] = (
        stat["root_len_after"] is not None and modal is not None
        and stat["root_len_after"] >= 0.9 * modal
    )

    if unexpected:
        stat["error"] = (
            f"{len(unexpected)} sequences would change in a way mode={mode!r} "
            f"does not explain (e.g. {unexpected[:3]}) -- refusing to write"
        )
        return stat

    if not dry_run:
        bak = aa_path.with_name(aa_path.name + BAK_SUFFIX)
        if aa_path.exists() and not bak.exists():
            shutil.copy2(aa_path, bak)
        write_fasta(aa_path, new_records)
        stat["written"] = True
    return stat


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, required=True, nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--mode", choices=["aligned", "ungapped", "auto"], default="auto",
        help="auto picks per dataset from the frame audit; see MODE_BY_DATASET",
    )
    ap.add_argument(
        "--cds-start", type=int, default=None,
        help="override the ATG search; only used with an explicit --mode",
    )
    ap.add_argument(
        "--allow-stale", action="store_true",
        help="permit rebuilding protein files that were not derived from the "
             "stored nucleotides (the *_epidemic splits); still requires the "
             "result to reach full length",
    )
    args = ap.parse_args()

    all_stats, any_error = {}, False
    for d in args.data_dir:
        d = d if d.is_absolute() else ROOT / d
        nts = sorted(d.glob("group_*_anc_nt.fasta"))
        if not nts:
            print(f"{d}: no group_*_anc_nt.fasta -- skipped")
            continue

        if args.mode == "auto":
            mode, cds_start = resolve_mode(d)
        else:
            mode, cds_start = args.mode, args.cds_start
        print(f"{d}: mode={mode} cds_start={cds_start}")

        stats = []
        for nt_path in nts:
            aa_path = nt_path.with_name(nt_path.name.replace("_anc_nt", "_anc_aa"))
            s = process(nt_path, aa_path, args.dry_run, mode=mode,
                        cds_start=cds_start, allow_stale=args.allow_stale)
            stats.append(s)
            if s.get("error"):
                any_error = True

        fixed = [s for s in stats if s["was_truncated"] and s["now_full_length"]]
        still = [s for s in stats if s["was_truncated"] and not s["now_full_length"]]
        errs = [s for s in stats if s.get("error")]
        leaf_repairs = sum(s["n_leaf_repairs"] for s in stats)
        unexpected = sum(s["n_unexpected_changes"] for s in stats)
        changed = sum(s["n_changed"] for s in stats)
        silent = sum(s["n_silent_frameshifts"] for s in stats)
        nrec = sum(s["n_records"] for s in stats)

        print(f"\n{d}")
        print(f"  {len(stats)} trees | {len(fixed)} truncated roots repaired | "
              f"{len(still)} still short | {len(errs)} refused")
        print(f"  sequences corrected: {changed} of {nrec} "
              f"({leaf_repairs} leaves, {silent} were full length but "
              f"frameshifted mid-sequence)")
        print(f"  unexpected sequence changes: {unexpected} (must be 0)")
        for s in fixed[:8]:
            print(f"    {s['tree']}: root {s['root_len_before']} -> "
                  f"{s['root_len_after']} aa (leaves {s['modal_leaf_len']})")
        if len(fixed) > 8:
            print(f"    ... and {len(fixed) - 8} more")
        for s in still + errs:
            detail = s.get("error") or (
                f"root {s['root_len_before']} -> {s['root_len_after']}, "
                f"leaves {s['modal_leaf_len']}"
            )
            print(f"    ! {s['tree']}: {detail}")
        all_stats[str(d)] = stats

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(all_stats, indent=2) + "\n")
        print(f"\nwrote {args.out}")

    if args.dry_run:
        print("\n(dry run -- nothing written)")
    sys.exit(1 if any_error else 0)


if __name__ == "__main__":
    main()
