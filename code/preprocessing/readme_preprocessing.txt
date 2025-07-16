README.txt  
Preprocessing Pipeline for RNA-seq Data (LongCOVID_RNAseq_meta)
=================================================================

This folder contains scripts for preprocessing bulk RNA-seq data for COVID-19 studies.
All scripts assume paired-end RNA-seq data.

Folder Contents:
-----------------
- STAR_indexing.sh            → Builds genome indices for Human and SARS-CoV-2 (required for STAR alignment)
- SRA_download_template.sh    → Downloads FASTQ files from SRA (edit SRR IDs inside the script)
- FastQC_Trimming_template.sh → Runs FastQC, trims adapters, and re-runs FastQC
- STAR_alignment.sh           → Aligns reads to Human and SARS-CoV-2 genomes using STAR
- Adapters.fa                 → Example adapter file for Trimmomatic (modify as needed)
- .gitkeep                    → Placeholder file

Processing Steps (In Order):
----------------------------

1. **Genome Indexing**
   - Run `STAR_indexing.sh` to generate indices for Human and SARS-CoV-2 genomes.
   - Output directories:  
     - `STAR_INDEX/human/`  
     - `STAR_INDEX/sars/`

2. **Download RNA-seq Data**
   - Edit `SRA_download_template.sh` to include your list of SRR IDs.
   - Run the script to download paired-end FASTQ files.

3. **Quality Control and Trimming**
   - Run `FastQC_Trimming_template.sh` to:
     - Check quality of raw reads
     - Trim adapters and low-quality bases
     - Check quality of trimmed reads
   - You must provide an adapter file in Trimmomatic format (example included: `Adapters.fa`).

4. **RNA-seq Alignment**
   - Run `STAR_alignment.sh` to align trimmed reads to both genomes.
   - Produces sorted BAM files and gene counts.
   - Output folders:
     - `STAR/human_quants/`
     - `STAR/sars_quants/`

Directory Structure:
--------------------
STAR_INDEX/
├── human/         → Human STAR index (created in Step 1)
├── sars/          → SARS-CoV-2 STAR index (created in Step 1)
├── GRCh38.primary_assembly.genome.fa
├── gencode.v38.annotation.gtf
├── wuhCor1.fa
└── ncbiGenes.gtf

Adapters.fa        → Adapter file for Trimmomatic (edit as needed)

fastqc_raw/        → FastQC results for raw reads  
fastqc_trimmed/    → FastQC results for trimmed reads  

STAR/
├── human_quants/  → STAR results for Human genome  
└── sars_quants/   → STAR results for SARS-CoV-2 genome  

Notes:
------
- Run each script in the order shown above.
- All scripts must be edited before use to adjust file paths, SRR IDs, and reference files.
- Designed for HPC environments (multi-threaded processing).
- Minimal console output (no `echo` commands).

If using these scripts, please cite:
Maison et al., Peripheral Immune Progression to Long COVID is Associated with Mitochondrial Gene Transcription: a Meta-Analysis (2025).