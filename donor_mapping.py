#!/usr/bin/env python3
"""Donor-name harmonisation.

Different sequencing runs sometimes label the same donor differently
(e.g. "VAMD05-C" vs "VAMD05", or a sample-sheet ID vs the SRA accession).
DONOR_MAPPING collapses every known alias to one canonical name.

map_donor() is the only thing other scripts should import: it looks a raw
label up in DONOR_MAPPING and falls back to returning it unchanged (so an
unmapped donor ID still works, it just isn't harmonised with any alias).
Status labels (Unassigned, Doublet, Ambiguous, ...) are never touched.
"""

DONOR_MAPPING = {
    "VAMD05-C": "VAMD05",
    "VAMD05": "VAMD05",
    "VAMD08-B": "VAMD08",
    "VAMD08": "VAMD08",
    "MSN25-B": "MSN25",
    "MSN25": "MSN25",
    "F13505-B": "F13505",
    "F13505": "F13505",
    "051064-A": "051064",
    "051064": "051064",
    "ik208": "051064",
    "NPBB219-A": "NPBB219",
    "NPBB219": "NPBB219",
    "RAJBrain_NPBB219": "NPBB219",
    "RAJBrain_NDBB060": "NPBB60",
    "NDBB060": "NPBB60",
    "MSN04-A": "MSN04",
    "MSN04": "MSN04",
    "SRR13291835": "MSN04",
    "MSN08-C": "MSN08",
    "MSN08": "MSN08",
    "VAMD04-A": "VAMD04",
    "VAMD04-C": "VAMD04",
    "VAMD04": "VAMD04",
    "NPBB60-C": "NPBB60",
    "NPBB60-A": "NPBB60",
    "NPBB60": "NPBB60",
    "GSA8_0_F13505.1.N-NEA11": "F13505",
}

# Labels that are call statuses, not donors - never remapped.
NON_DONOR = {"Unassigned", "Unresolved", "Doublet", "Ambiguous", "No evidence",
             "unassigned", "doublet", "NA", ""}


def map_donor(name):
    if name in NON_DONOR:
        return name
    return DONOR_MAPPING.get(name, name)
