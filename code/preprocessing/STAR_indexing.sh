#!/bin/bash

# ========================================
# STAR Genome Indexing Script
# For Human (GRCh38.p13) and SARS-CoV-2 (wuhCor1)
# ========================================

# Human Genome Indexing
# Downloaded from Ensembl: GRCh38.p13 (Release 109)
# Reference Genome FASTA and GTF annotation
human_genome_fasta="STAR_INDEX/GRCh38.primary_assembly.genome.fa"
human_gtf_file="STAR_INDEX/gencode.v38.annotation.gtf"
human_output_dir="STAR_INDEX/HUMAN"

mkdir -p "${human_output_dir}"

STAR --runMode genomeGenerate \
     --runThreadN 39 \
     --genomeDir "${human_output_dir}" \
     --genomeFastaFiles "${human_genome_fasta}" \
     --sjdbGTFfile "${human_gtf_file}" \
     --sjdbOverhang 100


# SARS-CoV-2 Genome Indexing
# Reference: wuhCor1 from Ensembl
sars_genome_fasta="STAR_INDEX/wuhCor1.fa"
sars_gtf_file="STAR_INDEX/ncbiGenes.gtf"
sars_output_dir="STAR_INDEX/sars"

mkdir -p "${sars_output_dir}"

STAR --runMode genomeGenerate \
     --runThreadN 39 \
     --genomeDir "${sars_output_dir}" \
     --genomeFastaFiles "${sars_genome_fasta}" \
     --sjdbGTFfile "${sars_gtf_file}" \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 6