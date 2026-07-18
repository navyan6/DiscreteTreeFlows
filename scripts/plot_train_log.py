#!/usr/bin/env python3
"""
Parse a TreeSBM training log into per-epoch loss curves + an overfitting report.

Reads the `Epoch NNN  train=... (rate=.. mut=.. cons=.. top=.. br=.. stop=.. pll=..)
val=...` lines that train.py prints, writes a CSV, prints a train/val overfitting
analysis, and (if matplotlib is available) saves a loss-curve PNG.

    python scripts/plot_train_log.py logs/h3n2_train_12345.log
    python scripts/plot_train_log.py logs/h3n2_train_12345.log --out curves/h3n2
"""

import argparse
import csv
import re
from pathlib import Path

LINE = re.compile(
    r"Epoch\s+(\d+)\s+train=([-\d.]+)\s+"
    r"\((?:rate|seq)=([-\d.]+)\s+mut=([-\d.]+)\s+cons=([-\d.]+)\s+"
    r"top=([-\d.]+)\s+br=([-\d.]+)\s+stop=([-\d.]+)\s+pll=([-\d.]+)\)\s+"
    r"val=([-\d.]+)"
)
COLS = ["epoch", "train", "rate", "mut", "cons", "top", "br", "stop", "pll", "val"]


def parse(log_path: str) -> list[dict]:
    rows = []
    for line in Path(log_path).read_text(errors="ignore").splitlines():
        m = LINE.search(line)
        if m:
            vals = [int(m.group(1))] + [float(x) for x in m.groups()[1:]]
            rows.append(dict(zip(COLS, vals)))
    return rows


def report(rows: list[dict]):
    if not rows:
        print("No epoch lines parsed — check the log format / path.")
        return
    epochs = [r["epoch"] for r in rows]
    train = [r["train"] for r in rows]
    val = [r["val"] for r in rows]

    best_i = min(range(len(val)), key=lambda i: val[i])
    best_ep, best_val = epochs[best_i], val[best_i]
    final_ep, final_val, final_train = epochs[-1], val[-1], train[-1]

    print(f"\nEpochs logged: {epochs[0]}..{final_ep}  ({len(rows)} points)")
    print(f"Best val: {best_val:.4f} @ epoch {best_ep}   "
          f"(train there = {train[best_i]:.4f}, gap = {val[best_i]-train[best_i]:+.4f})")
    print(f"Final:    train={final_train:.4f}  val={final_val:.4f}  "
          f"gap={final_val-final_train:+.4f}")

    # overfitting: val rose past best while train kept falling
    val_rise = final_val - best_val
    train_fell = train[best_i] - final_train
    epochs_past_best = final_ep - best_ep
    print(f"\nAfter best-val epoch: val moved {val_rise:+.4f}, train moved {-train_fell:+.4f}, "
          f"{epochs_past_best} epochs later.")
    if val_rise > 0.02 and train_fell > 0.0:
        print("VERDICT: OVERFITTING — val degrades while train improves past the best epoch.")
    elif epochs_past_best >= 1 and val_rise <= 0.0:
        print("VERDICT: still improving / not overfitting at stop.")
    else:
        print("VERDICT: roughly flat post-best (small val val set — treat early-stopping cautiously).")

    # component trend (first vs last)
    print("\nComponent losses (first -> last epoch):")
    for c in ["rate", "mut", "cons", "top", "br", "stop", "pll"]:
        print(f"  {c:5s} {rows[0][c]:.3f} -> {rows[-1][c]:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--out", default=None, help="output prefix (writes .csv and .png)")
    args = ap.parse_args()

    rows = parse(args.log)
    report(rows)
    if not rows:
        return

    out = args.out or str(Path(args.log).with_suffix(""))
    Path(out).parent.mkdir(parents=True, exist_ok=True)

    with open(out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out}.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ep = [r["epoch"] for r in rows]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        ax1.plot(ep, [r["train"] for r in rows], label="train total", lw=2)
        ax1.plot(ep, [r["val"] for r in rows], label="val total", lw=2)
        bi = min(range(len(rows)), key=lambda i: rows[i]["val"])
        ax1.axvline(rows[bi]["epoch"], color="k", ls="--", alpha=0.5, label=f"best val @ {rows[bi]['epoch']}")
        ax1.set_xlabel("epoch"); ax1.set_ylabel("loss"); ax1.set_title("Train vs Val"); ax1.legend()
        for c in ["rate", "mut", "cons", "top", "br", "stop", "pll"]:
            ax2.plot(ep, [r[c] for r in rows], label=c)
        ax2.set_xlabel("epoch"); ax2.set_ylabel("loss"); ax2.set_title("Components (train)"); ax2.legend()
        fig.tight_layout(); fig.savefig(out + ".png", dpi=120)
        print(f"Wrote {out}.png")
    except ImportError:
        print("(matplotlib not available — CSV written; plot elsewhere)")


if __name__ == "__main__":
    main()
