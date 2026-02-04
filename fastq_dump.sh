#!/bin/bash

#####################################################

#SBATCH --job-name=fastq_dump
#SBATCH --output=fastq-dump_%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00

######################################################


module load sra-tools/3.0.3-mrfsga6

######################################################

DIR=/home/jkt21/rnaseq_workflow/GSE81089
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastq

for i in ${DIR}/*/*.sra; 
do 
fasterq-dump "$i" --split-3 --outdir ${OUT_DIR}
done


