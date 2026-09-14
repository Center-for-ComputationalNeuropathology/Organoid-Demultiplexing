#!/usr/bin/env python3
"""One headline figure + table: how Vireo and the all-SNP approach reconcile,
across every cell in every batch.

Partitions every cell into exactly one of:
  agree            - both methods assign the same donor
  rescued_ge       - Vireo did not assign (Unassigned/Doublet); all-SNP did, cell is >= cutoff UMI
  rescued_lt       - same, but < cutoff UMI
  discordant       - both assign, but to different donors; or Vireo assigns and all-SNP can't
  unresolved_both  - neither method resolves the cell (Ambiguous / No evidence)

Reads every <results-dir>/<batch>_all_SNPs_vs_vireo/all_SNP_cell_assignments_UMI_flagged.tsv
(the output of assign_donors_from_snps.py).

Usage
-----
  python scripts/plot_reconciliation_summary.py --results-dir results
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATUS = {"Unassigned", "Doublet", "Ambiguous", "No evidence"}
COLORS = {
    "agree": "#2E7D32",
    "rescued_ge": "#1565C0",
    "rescued_lt": "#90CAF9",
    "discordant": "#D32F2F",
    "unresolved_both": "#BDBDBD",
}
LABELS = {
    "agree": "Agree (same donor)",
    "rescued_ge": "Rescued, ≥cutoff UMI",
    "rescued_lt": "Rescued, <cutoff UMI",
    "discordant": "Discordant",
    "unresolved_both": "Unresolved by both",
}
ORDER = ["agree", "rescued_ge", "rescued_lt", "discordant", "unresolved_both"]


def discover_batches(results_dir):
    return sorted(p.name[:-len("_all_SNPs_vs_vireo")]
                 for p in Path(results_dir).glob("*_all_SNPs_vs_vireo"))


def classify_and_count(results_dir, batches):
    c = Counter()
    n_ambig = n_noev = n_snp_assigned = n_vireo_assigned = n_vireo_unassigned = n_vireo_doublet = 0
    for b in batches:
        f = Path(results_dir) / f"{b}_all_SNPs_vs_vireo" / "all_SNP_cell_assignments_UMI_flagged.tsv"
        for r in csv.DictReader(open(f), delimiter="\t"):
            va, sa = r["Vireo_assignment"], r["all_SNP_assignment"]
            ge = r["UMI_group"] in ("≥cutoff UMI", "≥500 UMI")
            v_assigned = va not in STATUS
            s_assigned = sa not in STATUS
            n_vireo_assigned += v_assigned
            n_vireo_unassigned += va == "Unassigned"
            n_vireo_doublet += va == "Doublet"
            n_snp_assigned += s_assigned
            n_ambig += sa == "Ambiguous"
            n_noev += sa == "No evidence"
            if v_assigned and s_assigned and va == sa:
                c["agree"] += 1
            elif s_assigned and not v_assigned:
                c["rescued_ge" if ge else "rescued_lt"] += 1
            elif (v_assigned and s_assigned and va != sa) or (v_assigned and not s_assigned):
                c["discordant"] += 1
            else:
                c["unresolved_both"] += 1
    totals = {"total": sum(c.values()), "vireo_assigned": n_vireo_assigned,
             "vireo_unassigned": n_vireo_unassigned, "vireo_doublet": n_vireo_doublet,
             "snp_assigned": n_snp_assigned, "snp_ambiguous": n_ambig, "snp_no_evidence": n_noev}
    return c, totals


def draw(counts, total, out_png):
    fig, ax = plt.subplots(figsize=(12, 3.4))
    left = 0.0
    label_y_above = []
    for key in ORDER:
        n = counts[key]
        width = 100 * n / total
        ax.barh(0, width, left=left, color=COLORS[key], edgecolor="white", linewidth=1.2,
               height=0.6, zorder=2)
        center = left + width / 2
        text = f"{LABELS[key]}\n{n:,} ({width:.1f}%)"
        if width >= 11:
            ax.text(center, 0, text, ha="center", va="center", fontsize=9.5,
                    color="white" if key in ("agree", "rescued_ge", "discordant") else "#222222",
                    fontweight="bold", linespacing=1.4, zorder=3)
        else:
            label_y_above.append((center, text, key))
        left += width

    # external callouts for slivers too thin to label in place
    for i, (center, text, key) in enumerate(label_y_above):
        y = 0.62 + 0.42 * (i + 1)
        ax.plot([center, center], [0.31, y - 0.06], color="#666666", linewidth=1, zorder=1)
        ax.text(center, y, text, ha="center", va="bottom", fontsize=8.8, color="#333333",
                linespacing=1.3)

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.55, 0.55 + 0.42 * (len(label_y_above) + 1))
    ax.set_yticks([])
    ax.set_xlabel("Percent of all cells")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title(f"Vireo vs all-informative-SNP assignment: how {total:,} cells reconcile",
                fontsize=13, fontweight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    print("wrote", out_png)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default=None, help="default: same as --results-dir")
    ap.add_argument("--batches", default=None)
    args = ap.parse_args()
    out_dir = Path(args.out_dir or args.results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    batches = args.batches.split(",") if args.batches else discover_batches(args.results_dir)
    if not batches:
        raise SystemExit(f"No *_all_SNPs_vs_vireo dirs found in {args.results_dir}")

    counts, totals = classify_and_count(args.results_dir, batches)
    with open(out_dir / "reconciliation_summary.tsv", "w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["category", "cell_count", "percent_of_all_cells"])
        for key in ORDER:
            w.writerow([LABELS[key], counts[key], round(100 * counts[key] / totals["total"], 3)])

    draw(counts, totals["total"], out_dir / "reconciliation_summary.png")

    print(f"\ntotal cells: {totals['total']:,}")
    print(f"Vireo:  assigned {totals['vireo_assigned']:,} | unassigned {totals['vireo_unassigned']:,} "
          f"| doublet {totals['vireo_doublet']:,}")
    print(f"all-SNP: assigned {totals['snp_assigned']:,} | ambiguous {totals['snp_ambiguous']:,} "
          f"| no evidence {totals['snp_no_evidence']:,}")
    for key in ORDER:
        print(f"  {LABELS[key]:24}: {counts[key]:>8,}  ({100*counts[key]/totals['total']:.2f}%)")


if __name__ == "__main__":
    main()
