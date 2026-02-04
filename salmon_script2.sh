#!/bin/bash

#####################################################

#SBATCH --job-name=Salmon2
#SBATCH --output=Salmon2-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-24:00:00

#####################################################

module load salmon/1.10.2-ffgkyki

#####################################################

DIR=/home/jkt21/rnaseq_workflow/GSE81089/fastq
INDEX_DIR=/home/jkt21/rnaseq_workflow/GSE81089/salmon_practice/salmon_index
OUT_DIR=/home/jkt21/rnaseq_workflow/GSE81089/salmon_practice/salmon_quant2

for r1 in ${DIR}/*_1.fastq
do
    samp=$(basename "$r1" _1.fastq)
    r2=${DIR}/${samp}_2.fastq

    salmon quant \
        -i ${INDEX_DIR} \
        -l A \
        -1 ${r1} \
        -2 ${r2} \
        -p 2 \
        --validateMappings \
        -o ${OUT_DIR}/${samp}_quant2
done
