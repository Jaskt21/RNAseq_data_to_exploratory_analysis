#!/bin/bash

#####################################################

#SBATCH --job-name=sra_prefetch
#SBATCH --output=sra_prefetch_%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00

#####################################################

module load  sra-tools/3.0.3-mrfsga6

#####################################################

DIR=/home/jkt21/rnaseq_workflow/GSE81089
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089


cat ${DIR}/srr_accessions.txt | parallel -j 10 prefetch {} 





