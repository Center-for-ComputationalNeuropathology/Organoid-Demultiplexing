#!/usr/bin/env python3
"""Build a tiny, fully-synthetic BATCH-A dataset that exercises every outcome
the pipeline can produce. No real donor or patient data - every ID, read
count, and genotype below is fabricated for demonstration.

Run from this directory: python make_example_data.py
Produces BATCH-A_SNP_evidence_part01.tsv.gz and BATCH-A_cell_scores.tsv.gz.
"""
import csv
import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Raw donor labels chosen to demo donor_mapping.py: two have known aliases
# ("VAMD05-C" -> "VAMD05", "051064-A" -> "051064"), one is novel and passes
# through unchanged ("DonorZ").
DONORS = ["VAMD05-C", "051064-A", "DonorZ"]

CELLS = {
    # cell_barcode: (vireo_status, vireo_donor, UMI, vireo_probability)
    "CELL_CONCORDANT": ("Singlet", "DonorZ", 5000, 0.95),
    "CELL_RESCUED":    ("Unassigned", "unassigned", 800, 0.55),
    "CELL_TIE":        ("Unassigned", "unassigned", 600, 0.30),
    "CELL_NO_EVIDENCE": ("Unassigned", "unassigned", 40, 0.20),
    "CELL_DOUBLET":    ("Doublet", "doublet", 9000, 0.03),
}

# Each entry: (chrom, pos, REF, ALT, REF_reads, ALT_reads, {donor: GT})
# CELL_NO_EVIDENCE intentionally has zero rows - it should have no allele data at all.
EVIDENCE = {
    "CELL_CONCORDANT": [
        (1, 1001, "A", "G", 0, 5, {"VAMD05-C": "0/0", "051064-A": "0/0", "DonorZ": "1/1"}),
        (1, 2002, "C", "T", 4, 0, {"VAMD05-C": "0/1", "051064-A": "0/1", "DonorZ": "0/0"}),
        (2, 3003, "G", "A", 0, 3, {"VAMD05-C": "0/0", "051064-A": "0/1", "DonorZ": "1/1"}),
        (2, 4004, "T", "C", 3, 0, {"VAMD05-C": "0/1", "051064-A": "1/1", "DonorZ": "0/0"}),
        (3, 5005, "A", "T", 0, 2, {"VAMD05-C": "0/1", "051064-A": "0/0", "DonorZ": "1/1"}),
        (3, 6006, "C", "G", 2, 0, {"VAMD05-C": "1/1", "051064-A": "0/1", "DonorZ": "0/0"}),
    ],
    "CELL_RESCUED": [
        (1, 1001, "A", "G", 0, 3, {"VAMD05-C": "1/1", "051064-A": "0/0", "DonorZ": "0/0"}),
        (1, 2002, "G", "C", 2, 0, {"VAMD05-C": "0/0", "051064-A": "0/1", "DonorZ": "0/1"}),
        (2, 3003, "T", "A", 0, 2, {"VAMD05-C": "0/1", "051064-A": "0/0", "DonorZ": "0/0"}),
        (2, 4004, "C", "T", 0, 1, {"VAMD05-C": "1/1", "051064-A": "0/0", "DonorZ": "0/1"}),
    ],
    "CELL_TIE": [
        # VAMD05-C and 051064-A share the same genotype at both loci -> exact tie.
        (1, 1001, "A", "G", 0, 2, {"VAMD05-C": "0/1", "051064-A": "0/1", "DonorZ": "0/0"}),
        (1, 2002, "C", "T", 1, 0, {"VAMD05-C": "0/1", "051064-A": "0/1", "DonorZ": "1/1"}),
    ],
    "CELL_DOUBLET": [
        (1, 1001, "A", "G", 0, 4, {"VAMD05-C": "0/0", "051064-A": "1/1", "DonorZ": "0/0"}),
        (1, 2002, "T", "C", 3, 0, {"VAMD05-C": "0/1", "051064-A": "0/0", "DonorZ": "0/1"}),
        (2, 3003, "G", "A", 0, 2, {"VAMD05-C": "0/1", "051064-A": "1/1", "DonorZ": "0/0"}),
    ],
}


def write_scores():
    fields = ["cell", "vireo_status", "vireo_donor", "UMI", "vireo_probability"]
    with gzip.open(HERE / "BATCH-A_cell_scores.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for cell, (status, donor, umi, prob) in CELLS.items():
            w.writerow({"cell": cell, "vireo_status": status, "vireo_donor": donor,
                       "UMI": umi, "vireo_probability": prob})


def write_evidence():
    fields = ["cell", "chrom", "pos", "REF", "ALT", "REF_reads", "ALT_reads"] + \
             [f"GT_{d}" for d in DONORS]
    with gzip.open(HERE / "BATCH-A_SNP_evidence_part01.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for cell, loci in EVIDENCE.items():
            for chrom, pos, ref, alt, rr, ar, gts in loci:
                row = {"cell": cell, "chrom": chrom, "pos": pos, "REF": ref, "ALT": alt,
                      "REF_reads": rr, "ALT_reads": ar}
                row.update({f"GT_{d}": gts[d] for d in DONORS})
                w.writerow(row)


if __name__ == "__main__":
    write_scores()
    write_evidence()
    print("wrote BATCH-A_cell_scores.tsv.gz and BATCH-A_SNP_evidence_part01.tsv.gz")
