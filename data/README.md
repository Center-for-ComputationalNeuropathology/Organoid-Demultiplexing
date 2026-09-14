# Data

Real batch data is **not** tracked in this repo (see `.gitignore`) - only the
synthetic `example/` fixture is. Drop your own batches anywhere and point
`--data-dir` at that directory.

## Required input files, per batch

For a batch named `<batch>` (any string, e.g. `RUN-3`), the pipeline expects
two gzipped, tab-separated files with these exact names in `--data-dir`:

### `<batch>_SNP_evidence_part01.tsv.gz`
One row per (cell, SNP locus). Required columns:

| column | meaning |
|---|---|
| `cell` | cell barcode |
| `chrom`, `pos` | locus coordinates |
| `REF`, `ALT` | reference / alternate base |
| `REF_reads`, `ALT_reads` | that cell's read counts at this locus |
| `GT_<donor>` (one column per candidate donor in the pool) | donor's genotype: `0/0`, `0/1`, `1/1`, or `./.` for missing |

The donor panel for a batch is read straight from these `GT_` column names -
there's nothing to configure per batch.

### `<batch>_cell_scores.tsv.gz`
One row per cell. Required columns:

| column | meaning |
|---|---|
| `cell` | cell barcode (must match the evidence file) |
| `vireo_status` | `Singlet`, `Doublet`, or `Unassigned` |
| `vireo_donor` | donor name if `Singlet`, else any placeholder |
| `UMI` | UMI count for the low/high split |
| `vireo_probability` | Vireo's posterior for its call (carried through to outputs, not thresholded here) |

## Example

`example/make_example_data.py` builds a 5-cell, 3-donor, fully-synthetic
`BATCH-A` dataset that exercises every outcome the pipeline can produce:
a clean singlet, a Vireo-unassigned cell the manual method rescues, an
exact tie (`Ambiguous`), a cell with zero informative reads (`No evidence`),
and a Vireo doublet the manual method still assigns to one donor. Regenerate
it any time with:

```
python data/example/make_example_data.py
```
