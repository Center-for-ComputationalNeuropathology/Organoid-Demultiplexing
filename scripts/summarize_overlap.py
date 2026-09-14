#!/usr/bin/env python3
"""Per-batch: total cells, Vireo assigned/unassigned/doublet, manual all-SNP
assigned/unresolved, and the overlap of the two approaches - all split by UMI.

Reads <results-dir>/<batch>_all_SNPs_vs_vireo/all_SNP_cell_assignments_UMI_flagged.tsv
(the output of assign_donors_from_snps.py). Prints four markdown tables and writes
two wide TSVs.

Usage
-----
  python scripts/summarize_overlap.py --results-dir results
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

STATUS = {"Unassigned", "Doublet", "Ambiguous", "No evidence", "Unresolved"}


def discover_batches(results_dir):
    return sorted(p.name[:-len("_all_SNPs_vs_vireo")]
                 for p in Path(results_dir).glob("*_all_SNPs_vs_vireo"))


def load(results_dir, batch):
    f = Path(results_dir) / f"{batch}_all_SNPs_vs_vireo" / "all_SNP_cell_assignments_UMI_flagged.tsv"
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    for r in rows:
        r["_g"] = "ge" if r["UMI_group"] in ("≥cutoff UMI", "≥500 UMI") else "lt"
    return rows


def classify(rows):
    for r in rows:
        va, sa, g = r["Vireo_assignment"], r["all_SNP_assignment"], r["_g"]
        v_assigned = va not in STATUS
        s_assigned = sa not in STATUS
        vireo_kind = "assigned" if v_assigned else ("doublet" if va == "Doublet" else "unassigned")
        snp_kind = "assigned" if s_assigned else "unresolved"
        if v_assigned and s_assigned:
            ov = "both_same_donor" if va == sa else "both_diff_donor"
        elif v_assigned:
            ov = "vireo_only"
        elif s_assigned:
            ov = "snp_only_rescued"
        else:
            ov = "neither"
        yield vireo_kind, snp_kind, ov, g


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

    marg_fields = ["batch", "total_cells",
                   "vireo_assigned", "vireo_assigned_ge", "vireo_assigned_lt",
                   "vireo_unassigned", "vireo_unassigned_ge", "vireo_unassigned_lt",
                   "vireo_doublet",
                   "snp_assigned", "snp_assigned_ge", "snp_assigned_lt",
                   "snp_unresolved", "snp_unresolved_ge", "snp_unresolved_lt"]
    ovl_fields = ["batch", "total_cells",
                  "both_same_donor", "both_same_donor_ge", "both_same_donor_lt",
                  "both_diff_donor", "both_diff_donor_ge", "both_diff_donor_lt",
                  "vireo_only", "vireo_only_ge", "vireo_only_lt",
                  "snp_only_rescued", "snp_only_rescued_ge", "snp_only_rescued_lt",
                  "neither", "neither_ge", "neither_lt"]
    mrows, orows = [], []
    tot_m, tot_o = Counter(), Counter()

    for b in batches:
        rows = load(args.results_dir, b)
        n = len(rows)
        vk, vk_u, sk, sk_u, ov = Counter(), Counter(), Counter(), Counter(), Counter()
        for vireo_kind, snp_kind, ovk, g in classify(rows):
            vk[vireo_kind] += 1; vk_u[(vireo_kind, g)] += 1
            sk[snp_kind] += 1; sk_u[(snp_kind, g)] += 1
            ov[ovk] += 1; ov[(ovk, g)] += 1

        m = {"batch": b, "total_cells": n,
             "vireo_assigned": vk["assigned"], "vireo_assigned_ge": vk_u[("assigned", "ge")],
             "vireo_assigned_lt": vk_u[("assigned", "lt")],
             "vireo_unassigned": vk["unassigned"], "vireo_unassigned_ge": vk_u[("unassigned", "ge")],
             "vireo_unassigned_lt": vk_u[("unassigned", "lt")], "vireo_doublet": vk["doublet"],
             "snp_assigned": sk["assigned"], "snp_assigned_ge": sk_u[("assigned", "ge")],
             "snp_assigned_lt": sk_u[("assigned", "lt")], "snp_unresolved": sk["unresolved"],
             "snp_unresolved_ge": sk_u[("unresolved", "ge")], "snp_unresolved_lt": sk_u[("unresolved", "lt")]}
        mrows.append(m)
        for k, v in m.items():
            if k != "batch":
                tot_m[k] += v

        o = {"batch": b, "total_cells": n}
        for kind in ("both_same_donor", "both_diff_donor", "vireo_only", "snp_only_rescued", "neither"):
            o[kind] = ov[kind]; o[f"{kind}_ge"] = ov[(kind, "ge")]; o[f"{kind}_lt"] = ov[(kind, "lt")]
        orows.append(o)
        for k, v in o.items():
            if k != "batch":
                tot_o[k] += v

    mrows.append({"batch": "ALL", **{k: tot_m[k] for k in marg_fields if k != "batch"}})
    orows.append({"batch": "ALL", **{k: tot_o[k] for k in ovl_fields if k != "batch"}})

    with open(out_dir / "vireo_vs_snp_method_totals.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=marg_fields, delimiter="\t")
        w.writeheader(); w.writerows(mrows)
    with open(out_dir / "vireo_vs_snp_overlap.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=ovl_fields, delimiter="\t")
        w.writeheader(); w.writerows(orows)

    def show(rows, cols, title):
        print(f"\n## {title}\n")
        print("| " + " | ".join(cols) + " |")
        print("|" + "|".join(["---"] + ["--:"] * (len(cols) - 1)) + "|")
        for r in rows:
            print("| " + " | ".join(f"{r[c]:,}" if isinstance(r[c], int) else str(r[c]) for c in cols) + " |")

    show(mrows, ["batch", "total_cells", "vireo_assigned", "vireo_unassigned", "vireo_doublet",
                 "snp_assigned", "snp_unresolved"], "Method totals (confident assignments)")
    show(orows, ["batch", "total_cells", "both_same_donor", "both_diff_donor",
                 "vireo_only", "snp_only_rescued", "neither"], "Overlap of the two approaches")
    print(f"\nwrote vireo_vs_snp_method_totals.tsv and vireo_vs_snp_overlap.tsv in {out_dir}")


if __name__ == "__main__":
    main()
