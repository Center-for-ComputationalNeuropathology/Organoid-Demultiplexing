#!/usr/bin/env python3
"""Tabular summary of the combine_batches.py figures.

Reshapes combined_bar_counts.tsv (written by combine_batches.py) into, per batch,
one row per assignment category with Vireo vs all-informative-SNP cell counts split
by UMI, plus the net Vireo -> all-SNP change in total.

Usage
-----
  python scripts/summarize_combined_figure.py --results-dir results
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results",
                   help="directory containing combined_bar_counts.tsv (default: results)")
    ap.add_argument("--out-dir", default=None, help="default: same as --results-dir")
    args = ap.parse_args()
    out_dir = Path(args.out_dir or args.results_dir)
    src = Path(args.results_dir) / "combined_bar_counts.tsv"
    out = out_dir / "combined_bar_figure_summary.tsv"

    rows = list(csv.DictReader(open(src), delimiter="\t"))
    batches = list(dict.fromkeys(r["batch"] for r in rows))
    c = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    panel = {}
    for r in rows:
        grp = "ge" if r["UMI_group"] in ("≥cutoff UMI", "≥500 UMI") else "lt"
        c[r["batch"]][r["method"]][r["assignment"]][grp] += int(r["cell_count"])
    for b in batches:
        # union across both methods: a donor either one never confidently calls
        # must still appear in the panel.
        seen = set(c[b]["Vireo"]) | set(c[b]["All informative SNPs"])
        panel[b] = sorted(a for a in seen if a not in ("Unresolved", "Doublet"))

    fields = ["batch", "assignment", "Vireo_ge", "Vireo_lt", "Vireo_total",
              "allSNP_ge", "allSNP_lt", "allSNP_total", "total_change_Vireo_to_allSNP"]
    out_rows = []
    for b in batches:
        order = panel[b] + ["Unresolved", "Doublet"]
        for a in order:
            v, s = c[b]["Vireo"][a], c[b]["All informative SNPs"][a]
            vt, st = v["ge"] + v["lt"], s["ge"] + s["lt"]
            out_rows.append({"batch": b, "assignment": a,
                             "Vireo_ge": v["ge"], "Vireo_lt": v["lt"], "Vireo_total": vt,
                             "allSNP_ge": s["ge"], "allSNP_lt": s["lt"], "allSNP_total": st,
                             "total_change_Vireo_to_allSNP": st - vt})

    with open(out, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(out_rows)

    for b in batches:
        br = [r for r in out_rows if r["batch"] == b]
        tot = sum(r["Vireo_total"] for r in br)
        print(f"\n### {b}  ({tot:,} cells)  — donors: {', '.join(panel[b])}")
        print("| assignment | Vireo ge | Vireo lt | Vireo total | allSNP ge | allSNP lt | allSNP total | Δ total |")
        print("|---|--:|--:|--:|--:|--:|--:|--:|")
        for r in br:
            print(f"| {r['assignment']} | {r['Vireo_ge']:,} | {r['Vireo_lt']:,} | {r['Vireo_total']:,} | "
                  f"{r['allSNP_ge']:,} | {r['allSNP_lt']:,} | {r['allSNP_total']:,} | "
                  f"{r['total_change_Vireo_to_allSNP']:+,} |")

    print("\n\n### ALL batches combined")
    agg = defaultdict(lambda: defaultdict(int))
    for r in out_rows:
        key = "Donor" if r["assignment"] not in ("Unresolved", "Doublet") else r["assignment"]
        for k in ("Vireo_ge", "Vireo_lt", "Vireo_total", "allSNP_ge", "allSNP_lt", "allSNP_total"):
            agg[key][k] += r[k]
    print("| assignment | Vireo ge | Vireo lt | Vireo total | allSNP ge | allSNP lt | allSNP total | Δ total |")
    print("|---|--:|--:|--:|--:|--:|--:|--:|")
    for key in ("Donor", "Unresolved", "Doublet"):
        a = agg[key]
        print(f"| {key} | {a['Vireo_ge']:,} | {a['Vireo_lt']:,} | {a['Vireo_total']:,} | "
              f"{a['allSNP_ge']:,} | {a['allSNP_lt']:,} | {a['allSNP_total']:,} | "
              f"{a['allSNP_total'] - a['Vireo_total']:+,} |")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
