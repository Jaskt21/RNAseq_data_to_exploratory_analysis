#!/bin/bash

#SBATCH --job-name=multiQC_salmon
#SBATCH --output=multiQC_salmon-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00

#####################################################

module load py-multiqc/1.28-3ouxuso

#####################################################

# FastQC directory
FASTQC_DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastQC

# Salmon directory
SALMON_DIR=/home/jkt21/rnaseq_workflow/GSE81089/salmon_practice/salmon_quant2

# Star directory
STAR_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_featureCounts2

# Output directory for MultiQC report
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089/multiQC_full

# Run MultiQC on both directories
multiqc ${FASTQC_DIR} ${SALMON_DIR} ${STAR_DIR} -o ${OUT_DIR}