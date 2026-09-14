#!/usr/bin/env python3
"""One figure across every batch: Vireo vs all-informative-SNPs, coloured by donor.

Reads the per-batch output of assign_donors_from_snps.py
(<results-dir>/<batch>_all_SNPs_vs_vireo/all_SNP_cell_assignments_UMI_flagged.tsv)
and produces, per batch, two grouped bars (left = Vireo, right = all-SNP), each
stacked by donor (+ Unresolved + Doublet), split by UMI (solid = >=cutoff, faded
+ hatched = <cutoff). Donor panels and colours are inferred from the data, not
hard-coded, so this works for any set of batches/donors.

Usage
-----
  python scripts/combine_batches.py --results-dir results
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

STATUS = {"Unassigned", "Unresolved", "Doublet", "Ambiguous", "No evidence"}
METHODS = ["Vireo", "All informative SNPs"]


def discover_batches(results_dir):
    return sorted(p.name[:-len("_all_SNPs_vs_vireo")]
                 for p in Path(results_dir).glob("*_all_SNPs_vs_vireo"))


def load(results_dir, batches):
    """counts[batch][method][(category, umi_group)] = n ; donor panel per batch; tidy rows."""
    counts, panels, tidy = {}, {}, []
    for batch in batches:
        f = Path(results_dir) / f"{batch}_all_SNPs_vs_vireo" / "all_SNP_cell_assignments_UMI_flagged.tsv"
        rows = list(csv.DictReader(open(f), delimiter="\t"))
        donors = sorted({r["Vireo_assignment"] for r in rows if r["Vireo_assignment"] not in STATUS}
                        | {r["all_SNP_assignment"] for r in rows if r["all_SNP_assignment"] not in STATUS})
        panels[batch] = donors
        cb = {m: Counter() for m in METHODS}
        for r in rows:
            grp = "≥cutoff UMI" if r["UMI_group"] in ("≥cutoff UMI", "≥500 UMI") else "<cutoff UMI"
            va = r["Vireo_assignment"]
            vcat = va if va in donors else ("Doublet" if va == "Doublet" else "Unresolved")
            sa = r["all_SNP_assignment"]
            scat = sa if sa in donors else "Unresolved"
            cb["Vireo"][(vcat, grp)] += 1
            cb["All informative SNPs"][(scat, grp)] += 1
        counts[batch] = cb
        for m in METHODS:
            for (cat, grp), n in sorted(cb[m].items()):
                tidy.append({"batch": batch, "method": m, "assignment": cat,
                             "UMI_group": grp, "cell_count": n})
    return counts, panels, tidy


def draw(counts, panels, batches, donor_color, all_donors, normalize, out_png):
    x = np.arange(len(batches))
    bar_w = 0.38
    offsets = {"Vireo": -0.21, "All informative SNPs": 0.21}
    fig, ax = plt.subplots(figsize=(max(11, 1.8 * len(batches)), 9))

    for bi, batch in enumerate(batches):
        stack = panels[batch] + ["Unresolved", "Doublet"]
        total = sum(counts[batch]["Vireo"].values()) if normalize else 1
        for method in METHODS:
            xpos = x[bi] + offsets[method]
            bottom = 0.0
            for cat in stack:
                for grp in ("≥cutoff UMI", "<cutoff UMI"):
                    n = counts[batch][method][(cat, grp)]
                    if not n:
                        continue
                    h = 100 * n / total if normalize else n
                    faded = grp == "<cutoff UMI"
                    ax.bar(xpos, h, bar_w, bottom=bottom, color=donor_color[cat],
                           alpha=0.45 if faded else 1.0, hatch="////" if faded else None,
                           edgecolor="white", linewidth=0.3)
                    bottom += h
        ax.annotate("\n".join(panels[batch]), (x[bi], -0.045),
                    xycoords=("data", "axes fraction"), ha="center", va="top",
                    fontsize=7, color="#555555")

    ax.set_xticks(x)
    ax.set_xticklabels(batches, fontsize=11)
    ax.set_xlabel("batch   —   left bar = Vireo,  right bar = all informative SNPs", labelpad=42)
    ax.set_ylabel("Percent of batch cells" if normalize else "Cell count")
    ax.set_title("Vireo vs all-informative-SNP assignment "
                 f"({'per-batch %' if normalize else 'cell counts'}; faded + hatched = <cutoff UMI)")
    ax.margins(x=0.01)

    handles = ([Patch(facecolor=donor_color[d], label=d) for d in all_donors]
               + [Patch(facecolor=donor_color["Unresolved"], label="Unresolved"),
                  Patch(facecolor=donor_color["Doublet"], label="Doublet"),
                  Patch(facecolor="#888888", label="≥cutoff UMI"),
                  Patch(facecolor="#888888", alpha=0.45, hatch="////", label="<cutoff UMI")])
    leg = ax.legend(handles=handles, ncol=1, fontsize=9, title="Donor / status / UMI",
                    loc="center left", bbox_to_anchor=(1.005, 0.5))
    fig.subplots_adjust(bottom=0.20)
    fig.savefig(out_png, dpi=170, bbox_extra_artists=(leg,), bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_png)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results",
                   help="directory containing <batch>_all_SNPs_vs_vireo/ subdirs (default: results)")
    ap.add_argument("--out-dir", default=None, help="default: same as --results-dir")
    ap.add_argument("--batches", default=None, help="comma-separated; default: autodetect")
    args = ap.parse_args()
    out_dir = Path(args.out_dir or args.results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    batches = args.batches.split(",") if args.batches else discover_batches(args.results_dir)
    if not batches:
        raise SystemExit(f"No *_all_SNPs_vs_vireo dirs found in {args.results_dir}")

    counts, panels, tidy = load(args.results_dir, batches)
    with open(out_dir / "combined_bar_counts.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["batch", "method", "assignment", "UMI_group", "cell_count"],
                          delimiter="\t")
        w.writeheader(); w.writerows(tidy)

    all_donors = list(dict.fromkeys(d for p in panels.values() for d in p))
    donor_color = {d: c for d, c in zip(all_donors, plt.get_cmap("tab10").colors)}
    donor_color["Unresolved"] = "#BDBDBD"
    donor_color["Doublet"] = "#7B3294"

    draw(counts, panels, batches, donor_color, all_donors, normalize=False,
        out_png=out_dir / "combined_vireo_vs_allSNP_UMI_counts.png")
    draw(counts, panels, batches, donor_color, all_donors, normalize=True,
        out_png=out_dir / "combined_vireo_vs_allSNP_UMI_percent.png")


if __name__ == "__main__":
    main()
