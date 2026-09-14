#!/usr/bin/env python3
"""Donor composition using only cells >= the UMI cutoff, both methods, one colour scheme.

Reads combined_bar_counts.tsv (written by combine_batches.py), keeps only the
">= cutoff UMI" rows, and draws one stacked bar per (batch, method) - no hatching,
since UMI is now a hard filter rather than a split. Same donor colours as
combine_batches.py's figures.

Usage
-----
  python scripts/filter_umi500.py --results-dir results
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

METHODS = ["Vireo", "All informative SNPs"]


def load_ge_cutoff(src):
    counts, panels = {}, {}
    rows = list(csv.DictReader(open(src), delimiter="\t"))
    batches = list(dict.fromkeys(r["batch"] for r in rows))
    for b in batches:
        counts[b] = {m: Counter() for m in METHODS}
    for r in rows:
        if not (r["UMI_group"].startswith("≥")):
            continue
        counts[r["batch"]][r["method"]][r["assignment"]] += int(r["cell_count"])
    for b in batches:
        # union across both methods: a donor Vireo never confidently calls (or vice
        # versa) must still appear in the panel/legend.
        seen = set(counts[b]["Vireo"]) | set(counts[b]["All informative SNPs"])
        panels[b] = sorted(a for a in seen if a not in ("Unresolved", "Doublet"))
    return batches, counts, panels


def draw(batches, counts, panels, donor_color, normalize, out_png):
    fig, ax = plt.subplots(figsize=(max(11, 1.8 * len(batches)), 9))
    bar_w = 0.38
    offsets = {"Vireo": -0.21, "All informative SNPs": 0.21}
    xc = np.arange(len(batches))

    for bi, batch in enumerate(batches):
        stack = panels[batch] + ["Unresolved", "Doublet"]
        vtot = sum(counts[batch]["Vireo"].values())
        for method in METHODS:
            x = xc[bi] + offsets[method]
            denom = vtot if normalize else 1
            bottom = 0.0
            for cat in stack:
                n = counts[batch][method][cat]
                if not n:
                    continue
                h = 100 * n / denom if normalize else n
                ax.bar(x, h, bar_w, bottom=bottom, color=donor_color[cat],
                       edgecolor="white", linewidth=0.3)
                bottom += h
        ax.annotate("\n".join(panels[batch]), (xc[bi], -0.045),
                    xycoords=("data", "axes fraction"), ha="center", va="top",
                    fontsize=7, color="#555555")

    ax.set_xticks(xc)
    ax.set_xticklabels(batches, fontsize=11)
    ax.set_xlabel("batch   —   left bar = Vireo,  right bar = all informative SNPs "
                  "(>= cutoff UMI cells only)", labelpad=42)
    ax.set_ylabel("Percent of >= cutoff UMI cells" if normalize else "Cell count (>= cutoff UMI)")
    ax.set_title("Donor assignment of >= cutoff UMI cells: Vireo vs all informative SNPs "
                 f"({'per-batch %' if normalize else 'cell counts'})")
    ax.margins(x=0.01)

    all_donors = list(dict.fromkeys(d for p in panels.values() for d in p))
    handles = ([Patch(facecolor=donor_color[d], label=d) for d in all_donors]
               + [Patch(facecolor=donor_color["Unresolved"], label="Unresolved"),
                  Patch(facecolor=donor_color["Doublet"], label="Doublet (Vireo only)")])
    leg = ax.legend(handles=handles, ncol=1, fontsize=9, title="Donor / status",
                    loc="center left", bbox_to_anchor=(1.005, 0.5))
    fig.subplots_adjust(bottom=0.20)
    fig.savefig(out_png, dpi=170, bbox_extra_artists=(leg,), bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_png)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results",
                   help="directory containing combined_bar_counts.tsv (default: results)")
    ap.add_argument("--out-dir", default=None, help="default: same as --results-dir")
    args = ap.parse_args()
    out_dir = Path(args.out_dir or args.results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    batches, counts, panels = load_ge_cutoff(Path(args.results_dir) / "combined_bar_counts.tsv")

    with open(out_dir / "UMI_filtered_assignment_counts.tsv", "w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["batch", "method", "assignment", "cell_count_ge_cutoff_UMI"])
        for b in batches:
            for m in METHODS:
                for cat in panels[b] + ["Unresolved", "Doublet"]:
                    w.writerow([b, m, cat, counts[b][m][cat]])

    all_donors = list(dict.fromkeys(d for p in panels.values() for d in p))
    donor_color = {d: c for d, c in zip(all_donors, plt.get_cmap("tab10").colors)}
    donor_color["Unresolved"] = "#BDBDBD"
    donor_color["Doublet"] = "#7B3294"

    draw(batches, counts, panels, donor_color, normalize=False,
        out_png=out_dir / "UMI_filtered_assignments_counts.png")
    draw(batches, counts, panels, donor_color, normalize=True,
        out_png=out_dir / "UMI_filtered_assignments_percent.png")


if __name__ == "__main__":
    main()
