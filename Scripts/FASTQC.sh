#!/bin/bash

#####################################################

#SBATCH --job-name=Fastqc
#SBATCH --output=Fastqc-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00b

#####################################################

module load fastqc/0.12.1-u6klb7y

#####################################################

DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastq
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastQC

for fastq in ${DIR}/*.fastq; 
do 
fastqc "$fastq" -o ${OUT_DIR}
done
