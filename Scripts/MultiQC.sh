#!/bin/bash

#####################################################

#SBATCH --job-name=MultiQC
#SBATCH --output=MultiQC-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00

#####################################################

module load py-multiqc/1.28-3ouxuso

#####################################################

DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastQC
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089/multiQC2

for fastqc in ${DIR}; 
do 
multiqc "$fastqc" -o ${OUT_DIR}
done


























