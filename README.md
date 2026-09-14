# Donor assignment from all informative SNPs

A from-scratch, Vireo-independent way to assign pooled single-cell RNA-seq
cells back to their donor of origin using every donor-informative SNP a cell
has read coverage at, plus scripts to compare that call against Vireo's own
assignment across any number of batches.

## Method

1. **Find donor-informative SNPs.** A locus is informative for a batch if
   every candidate donor has a valid genotype there (`0/0`, `0/1`, `1/1`) and
   at least two donors differ.
2. **Score each cell against each donor.** Summing over every informative
   locus the cell has reads at:

   ```
   loglik[donor] += ALT_reads * ln(p) + REF_reads * ln(1 - p)
   ```

   where `p = P(observe an ALT read | genotype)` = `{0/0: 0.01, 0/1: 0.5, 1/1: 0.99}`.
3. **Call the donor with the highest total.** An exact tie -> `Ambiguous`.
   Zero informative reads -> `No evidence`.
4. **Flag by UMI** (`<cutoff` vs `>= cutoff`, default 500) and compare against
   Vireo's own `vireo_status` / `vireo_donor`.

This is deliberately simple - no doublet model, no prior, no genotype
uncertainty - so where it disagrees with Vireo, the disagreement is
informative: it usually means the confidence is being spent on very few
reads either way.

## Install

```
pip install -r requirements.txt
```

Everything else is the standard library.

## Quickstart (synthetic example)

```
python assign_donors_from_snps.py --data-dir data/example --out-dir results
python scripts/combine_batches.py         --results-dir results
python scripts/filter_umi500.py           --results-dir results
python scripts/summarize_overlap.py       --results-dir results
python scripts/summarize_combined_figure.py --results-dir results
```

`data/example/` is a 5-cell, fully-synthetic dataset (see `data/README.md`)
that exercises every outcome the pipeline can produce - run the five
commands above and check `results/` to see what each script contributes.

## On your own data

```
python assign_donors_from_snps.py --data-dir /path/to/batches --out-dir results
```

See `data/README.md` for the exact file naming and columns each batch needs.
Batches are auto-discovered from `*_SNP_evidence_part01.tsv.gz` in
`--data-dir`, or pass `--batches BATCH-A,BATCH-B` explicitly. Real donor
genotype data should stay out of version control - `.gitignore` already
excludes everything under `data/` except `data/example/`.

## Scripts

| Script | Input | Output |
|---|---|---|
| `assign_donors_from_snps.py` | `<batch>_SNP_evidence_part01.tsv.gz`, `<batch>_cell_scores.tsv.gz` | per-batch `<batch>_all_SNPs_vs_vireo/`: per-cell assignments TSV, count TSVs, two bar-chart PNGs |
| `scripts/combine_batches.py` | the per-batch dirs above | one figure across all batches (counts + %), coloured by donor, split by UMI |
| `scripts/filter_umi500.py` | `combine_batches.py`'s tidy TSV | the same comparison restricted to `>= cutoff` UMI cells, both methods on one colour scheme |
| `scripts/summarize_overlap.py` | the per-batch dirs | per-batch tables: Vireo assigned/unassigned/doublet, manual assigned/unresolved, and the cell-level overlap of the two (agree / disagree / rescued / neither) |
| `scripts/summarize_combined_figure.py` | `combine_batches.py`'s tidy TSV | the combined figure as a table, with the net Vireo -> manual change per category |

`donor_mapping.py` harmonises donor aliases (e.g. `"VAMD05-C"` and `"VAMD05"`
are the same donor) - edit `DONOR_MAPPING` for your own donor IDs; anything
not listed passes through unchanged.

## A note on interpreting disagreements

Wherever a cell has reasonable coverage (dozens of informative SNPs), Vireo
and this method agree essentially every time. Most of what this method
"rescues" beyond Vireo are low-UMI cells resolved on a handful of reads -
useful to know, but treat those calls as provisional, not equivalent in
confidence to a Vireo singlet call.
