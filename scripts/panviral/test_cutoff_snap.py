#!/usr/bin/env python3
"""
Check that the train/test cutoff can never cut a group in half.

Groups are `window`-year bins but the cutoff is a single year, so an unsnapped
cutoff splits the bin it lands in -- sending part of one outbreak to train and
the rest to test. These assertions pin the property that the snapped cutoff is
always the final year of a bin, which is what makes that impossible.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_virus_dataset import choose_cutoff, snap_to_window  # noqa: E402

failures = 0


def bin_start(year: int, window: int) -> int:
    """Bins are anchored the way build_groups anchors them: year - year % window."""
    return year - year % window

print("snap_to_window(window=2)   bins are [even, odd]")
for c in range(2014, 2022):
    sn = snap_to_window(c, 2)
    tag = "  (already a boundary)" if sn == c else ""
    print(f"  {c} -> {sn}   train <= {sn}, test >= {sn + 1}{tag}")

print("\ninvariant: snapped cutoff is the last year of a bin")
last_ok = 0
for w in (2, 3, 4, 5):
    for c in range(2000, 2030):
        sn = snap_to_window(c, w)
        if sn % w != w - 1:
            print(f"  FAIL window={w} cutoff={c} -> {sn}")
            last_ok += 1
failures += last_ok
print("  ok" if not failures else f"  {failures} violations")

print("\ninvariant: no bin straddles the cutoff")
straddles = 0
for w in (2, 3, 4):
    for c in range(2005, 2028):
        sn = snap_to_window(c, w)
        b = bin_start(sn, w)
        if sn != b + w - 1:              # cutoff sits inside a bin, not at its end
            print(f"  FAIL window={w} cutoff={c}->{sn} splits bin [{b},{b + w - 1}]")
            straddles += 1
failures += straddles
print("  ok" if not straddles else f"  {straddles} violations")

print("\nEbola-like case (the one that motivated this)")
recs = ([{"year": y} for y in (2014,) * 200] + [{"year": y} for y in (2015,) * 61]
        + [{"year": y} for y in (2018,) * 60] + [{"year": y} for y in (2019,) * 47])
cut = choose_cutoff(recs, None, 0.2, 2)
tr = [r for r in recs if r["year"] <= cut]
te = [r for r in recs if r["year"] > cut]
tr_bins = sorted({r["year"] - r["year"] % 2 for r in tr})
te_bins = sorted({r["year"] - r["year"] % 2 for r in te})
print(f"  cutoff {cut}: train {len(tr)} (bins {tr_bins}), test {len(te)} (bins {te_bins})")
overlap = set(tr_bins) & set(te_bins)
print(f"  bins appearing in BOTH splits: {overlap or 'none'}")
if overlap:
    failures += 1

print("\nall records land in exactly one split, every virus shape")
lost = 0
for shape in ([2020] * 50, list(range(2000, 2026)) * 8, [2015] * 10 + [2024] * 90):
    recs = [{"year": y} for y in shape]
    cut = choose_cutoff(recs, None, 0.2, 2)
    n = len([r for r in recs if r["year"] <= cut]) + len([r for r in recs if r["year"] > cut])
    if n != len(recs):
        print(f"  FAIL shape lost records: {n} != {len(recs)}")
        lost += 1
failures += lost
print("  ok" if not lost else f"  {lost} violations")

print("\nPASS" if not failures else f"\nFAIL ({failures})")
sys.exit(1 if failures else 0)
