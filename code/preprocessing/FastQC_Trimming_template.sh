#!/bin/bash
# RNA-seq Preprocessing Pipeline: FastQC, Trimmomatic
# Compatible with high-performance computing environments.

# =========================================
# Adapter Trimming Setup (IMPORTANT):
# =========================================
# Provide an adapter FASTA file combining all adapter sequences required.
# Example format (Adapters.fa):
#
# >NexteraPE-PE
# [adapter_sequence]
# >TruSeq3-PE
# [adapter_sequence]
#
# You may concatenate adapters from vendor-provided files (e.g., NexteraPE-PE.fa, TruSeq3-PE.fa)
# into a single file, separated by headers (>).
# Specify this file in the Trimmomatic command below:
# ILLUMINACLIP:Adapters.fa:2:30:10:8:TRUE
# =========================================

# Load FastQC module (adjust to your system's module name)
module load bio/FastQC/0.11.8-Java-1.8.0_191

# Run FastQC on raw reads
mkdir -p fastqc_raw
export OMP_NUM_THREADS=39
fastqc -o fastqc_raw *_1.fastq *_2.fastq

# Load Trimmomatic environment (adjust as needed)
module load lang/Anaconda3
source activate trimmomatic

# Adapter file (provide your own as described above)
adapter_file="Adapters.fa"

# Trim reads
for forward_read in *_1.fastq; do
  reverse_read="${forward_read%_1.fastq}_2.fastq"

  output_forward_paired="${forward_read%_1.fastq}_forward_paired.fastq"
  output_forward_unpaired="${forward_read%_1.fastq}_forward_unpaired.fastq"
  output_reverse_paired="${forward_read%_1.fastq}_reverse_paired.fastq"
  output_reverse_unpaired="${forward_read%_1.fastq}_reverse_unpaired.fastq"

  ~/.conda/envs/trimmomatic/bin/trimmomatic PE -threads 39 -phred33 \
    "$forward_read" "$reverse_read" \
    "$output_forward_paired" "$output_forward_unpaired" \
    "$output_reverse_paired" "$output_reverse_unpaired" \
    ILLUMINACLIP:"$adapter_file":2:30:10:8:TRUE \
    LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:20
done

# Run FastQC on trimmed reads
mkdir -p fastqc_trimmed
fastqc -o fastqc_trimmed *_paired.fastq