#!/bin/bash

# ========================================
# SRA Data Download Script (Generic Template)
# ========================================
# This script downloads raw FASTQ files from the Sequence Read Archive (SRA)
# using fasterq-dump.
# Requires: SRA Toolkit (version ≥ 3.0.0)

module load bio/SRA-Toolkit/3.0.0-centos_linux64
export OMP_NUM_THREADS=39

# List of SRA accession numbers (replace with your own)
ACCESSION_LIST=(
    SRRXXXXXXX1
    SRRXXXXXXX2
    SRRXXXXXXX3
    # Add more accessions as needed
)

for ACCESSION in "${ACCESSION_LIST[@]}"; do
    fasterq-dump --split-files "$ACCESSION"
done