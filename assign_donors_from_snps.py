#!/usr/bin/env python3
"""Donor assignment from raw allele evidence, independent of Vireo.

For each batch this expects two gzipped TSVs in --data-dir:

  <batch>_SNP_evidence_part01.tsv.gz   one row per (cell, SNP locus):
      cell, chrom, pos, REF, ALT, REF_reads, ALT_reads, GT_<donor> (one per donor)
  <batch>_cell_scores.tsv.gz           one row per cell:
      cell, vireo_status, vireo_donor, UMI, vireo_probability

Method
------
1. A SNP is "donor-informative" if every donor has a valid genotype
   (0/0, 0/1, 1/1) and at least two donors differ there.
2. For each cell and each candidate donor, sum over informative SNPs:
       loglik += ALT_reads * ln(p) + REF_reads * ln(1 - p)
   where p = P(observe an ALT read | genotype) = {0/0: .01, 0/1: .5, 1/1: .99}.
3. The donor with the highest total wins. An exact tie -> "Ambiguous".
   Zero informative reads -> "No evidence".
4. Cells are flagged <cutoff> UMI (default 500) vs. >= cutoff, and the
   manual call is compared against Vireo's own status/donor.

Output (per batch, under --out-dir/<batch>_all_SNPs_vs_vireo/)
----------------------------------------------------------------
  all_SNP_cell_assignments_UMI_flagged.tsv      one row per cell
  vireo_vs_all_SNPs_counts_by_UMI_flag.tsv      counts x (method, category, UMI group)
  vireo_vs_all_SNPs_harmonized_counts_by_UMI_flag.tsv   same, unresolved states folded
  vireo_vs_all_informative_SNPs_UMI_flagged_bars.png
  vireo_vs_all_informative_SNPs_harmonized_UMI_bars.png

Usage
-----
  python assign_donors_from_snps.py --data-dir data/example --out-dir results
  python assign_donors_from_snps.py --data-dir /path/to/batches --batches BATCH-A,BATCH-B
"""

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from donor_mapping import map_donor

GT_ALT_PROB = {"0/0": 0.01, "0/1": 0.5, "1/0": 0.5, "1/1": 0.99}
METHOD_COLORS = {"Vireo": ("#4C78A8", "#9ECAE9"),
                 "All informative SNPs": ("#E45756", "#F2A7A5")}
WIDTH = 0.34
EVIDENCE_SUFFIX = "_SNP_evidence_part01.tsv.gz"
SCORES_SUFFIX = "_cell_scores.tsv.gz"


def discover_batches(data_dir):
    return sorted(p.name[:-len(EVIDENCE_SUFFIX)]
                 for p in Path(data_dir).glob(f"*{EVIDENCE_SUFFIX}"))


def donors_from_header(evidence_path):
    with gzip.open(evidence_path, "rt", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    return [col[3:] for col in header if col.startswith("GT_")]


def assign_from_all_snps(evidence_path, cells, donors):
    """Log-likelihood donor call per cell using only donor-informative loci."""
    all_loci, informative_loci = set(), set()
    with gzip.open(evidence_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], row["pos"], row["REF"], row["ALT"])
            all_loci.add(key)
            gts = tuple(row[f"GT_{d}"] for d in donors)
            if all(gt in GT_ALT_PROB for gt in gts) and len(set(gts)) >= 2:
                informative_loci.add(key)

    ll = defaultdict(lambda: {d: 0.0 for d in donors})
    sites, reads = Counter(), Counter()
    with gzip.open(evidence_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], row["pos"], row["REF"], row["ALT"])
            if key not in informative_loci:
                continue
            cell = row["cell"]
            rr, ar = int(row["REF_reads"]), int(row["ALT_reads"])
            sites[cell] += 1
            reads[cell] += rr + ar
            for donor in donors:
                p = GT_ALT_PROB[row[f"GT_{donor}"]]
                ll[cell][donor] += ar * math.log(p) + rr * math.log(1 - p)

    assignments, margins = {}, {}
    for cell in cells:
        if sites[cell] == 0:
            assignments[cell], margins[cell] = "No evidence", "NA"
            continue
        ordered = sorted(ll[cell].items(), key=lambda item: item[1], reverse=True)
        best = ordered[0][1]
        winners = [d for d, v in ordered if abs(v - best) < 1e-12]
        assignments[cell] = winners[0] if len(winners) == 1 else "Ambiguous"
        margins[cell] = best - ordered[1][1]
    return assignments, margins, ll, sites, reads, all_loci, informative_loci


def grouped_stack_plot(out_png, categories, counts, methods, title, xlabel):
    x = np.arange(len(categories))
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for i, method in enumerate(methods):
        xpos = x + (-WIDTH / 2 if i == 0 else WIDTH / 2)
        high = [counts[method][(c, "≥cutoff UMI")] for c in categories]
        low = [counts[method][(c, "<cutoff UMI")] for c in categories]
        dark, light = METHOD_COLORS[method]
        ax.bar(xpos, high, WIDTH, color=dark, label=f"{method}: ≥cutoff UMI")
        ax.bar(xpos, low, WIDTH, bottom=high, color=light, hatch="//",
               label=f"{method}: <cutoff UMI")
    ax.set_xticks(x, categories, rotation=20, ha="right")
    ax.set_ylabel("Cell count")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def process_batch(batch, data_dir, out_dir, umi_cutoff):
    evidence = Path(data_dir) / f"{batch}{EVIDENCE_SUFFIX}"
    scores = Path(data_dir) / f"{batch}{SCORES_SUFFIX}"
    if not (evidence.exists() and scores.exists()):
        print(f"{batch}: missing evidence or cell-scores file, skipped")
        return
    out = Path(out_dir) / f"{batch}_all_SNPs_vs_vireo"
    out.mkdir(parents=True, exist_ok=True)

    donors = donors_from_header(evidence)          # raw labels (match GT_/loglik_ columns)
    dmap = {d: map_donor(d) for d in donors}       # raw -> harmonised
    disp_donors = [dmap[d] for d in donors]        # harmonised, panel order
    meta = {}
    with gzip.open(scores, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            meta[row["cell"]] = row

    assignments, margins, ll, sites, reads, all_loci, informative_loci = \
        assign_from_all_snps(evidence, meta, donors)

    def vireo_assignment(row):
        if row["vireo_status"] == "Singlet" and row["vireo_donor"] in donors:
            return dmap[row["vireo_donor"]]
        return "Doublet" if row["vireo_status"] == "Doublet" else "Unassigned"

    flagged_fields = ["Cell_barcode", "Vireo_status", "Vireo_donor", "all_informative_SNPs_observed",
                      "informative_reads", "all_SNP_assignment", "loglik_margin",
                      *[f"loglik_{dmap[d]}" for d in donors],
                      "Vireo_assignment", "UMI_count", "Vireo_probability",
                      "UMI_below_cutoff", "UMI_group"]
    flagged = []
    for cell, row in meta.items():
        low = int(float(row["UMI"])) < umi_cutoff
        flagged.append({
            "Cell_barcode": cell, "Vireo_status": row["vireo_status"],
            "Vireo_donor": map_donor(row["vireo_donor"]),
            "all_informative_SNPs_observed": sites[cell], "informative_reads": reads[cell],
            "all_SNP_assignment": map_donor(assignments[cell]), "loglik_margin": margins[cell],
            **{f"loglik_{dmap[d]}": ll[cell][d] for d in donors},
            "Vireo_assignment": vireo_assignment(row),
            "UMI_count": row["UMI"], "Vireo_probability": row["vireo_probability"],
            "UMI_below_cutoff": low, "UMI_group": "<cutoff UMI" if low else "≥cutoff UMI",
        })

    with open(out / "all_SNP_cell_assignments_UMI_flagged.tsv", "w", newline="") as handle:
        wr = csv.DictWriter(handle, fieldnames=flagged_fields, delimiter="\t")
        wr.writeheader(); wr.writerows(flagged)

    methods = {"Vireo": "Vireo_assignment", "All informative SNPs": "all_SNP_assignment"}
    groups = ["≥cutoff UMI", "<cutoff UMI"]
    n = len(flagged)

    categories = [*disp_donors, "Unassigned", "Doublet", "Ambiguous", "No evidence"]
    counts, count_rows = {}, []
    for method, field in methods.items():
        counts[method] = Counter((r[field], r["UMI_group"]) for r in flagged)
        for c in categories:
            for g in groups:
                v = counts[method][(c, g)]
                count_rows.append({"method": method, "assignment": c, "UMI_group": g,
                                   "cell_count": v, "percent_of_all_cells": 100 * v / n})
    with open(out / "vireo_vs_all_SNPs_counts_by_UMI_flag.tsv", "w", newline="") as handle:
        wr = csv.DictWriter(handle, fieldnames=list(count_rows[0]), delimiter="\t")
        wr.writeheader(); wr.writerows(count_rows)
    grouped_stack_plot(out / "vireo_vs_all_informative_SNPs_UMI_flagged_bars.png",
                       categories, counts, methods,
                       f"{batch}: Vireo versus all informative SNPs (all {n:,} cells)",
                       "Assignment")

    harmonized = [*disp_donors, "Unresolved", "Doublet"]

    def harmonize(method, call):
        if method == "Vireo":
            return "Unresolved" if call == "Unassigned" else call
        return "Unresolved" if call in {"Ambiguous", "No evidence"} else call

    harm_counts, harm_rows = {}, []
    for method, field in methods.items():
        counter = Counter((harmonize(method, r[field]), r["UMI_group"]) for r in flagged)
        harm_counts[method] = counter
        for c in harmonized:
            for g in groups:
                v = counter[(c, g)]
                harm_rows.append({"method": method, "harmonized_assignment": c,
                                  "UMI_group": g, "cell_count": v,
                                  "percent_of_all_cells": 100 * v / n})
    with open(out / "vireo_vs_all_SNPs_harmonized_counts_by_UMI_flag.tsv", "w", newline="") as handle:
        wr = csv.DictWriter(handle, fieldnames=list(harm_rows[0]), delimiter="\t")
        wr.writeheader(); wr.writerows(harm_rows)
    grouped_stack_plot(out / "vireo_vs_all_informative_SNPs_harmonized_UMI_bars.png",
                       harmonized, harm_counts, methods,
                       f"{batch}: Vireo versus all informative SNPs (all {n:,} cells)",
                       "Harmonized assignment")

    below = sum(r["UMI_below_cutoff"] for r in flagged)
    print(f"{batch}: donors={[f'{d}->{dmap[d]}' for d in donors]} cells={n} "
          f"<{umi_cutoff} UMI={below} loci={len(all_loci)} informative={len(informative_loci)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="data",
                   help="directory containing <batch>_SNP_evidence_part01.tsv.gz and "
                        "<batch>_cell_scores.tsv.gz (default: data)")
    ap.add_argument("--out-dir", default="results",
                   help="directory to write <batch>_all_SNPs_vs_vireo/ into (default: results)")
    ap.add_argument("--batches", default=None,
                   help="comma-separated batch names; default: autodetect every "
                        f"*{EVIDENCE_SUFFIX} in --data-dir")
    ap.add_argument("--umi-cutoff", type=int, default=500,
                   help="UMI count threshold for the low/high split (default: 500)")
    args = ap.parse_args()

    batches = args.batches.split(",") if args.batches else discover_batches(args.data_dir)
    if not batches:
        raise SystemExit(f"No {EVIDENCE_SUFFIX} files found in {args.data_dir}")
    for batch in batches:
        process_batch(batch, args.data_dir, args.out_dir, args.umi_cutoff)


if __name__ == "__main__":
    main()
