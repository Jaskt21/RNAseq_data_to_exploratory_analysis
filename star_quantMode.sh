#!/bin/bash

#####################################################

#SBATCH --job-name=star_quantMode
#SBATCH --output=star_quantMode-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0-24:00:00

#####################################################

module load star/2.7.11a-pgsk3s4

#####################################################

FASTQ_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/rnaseq_rawdata
GENOME_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_index
OUT_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_featureCounts2


for r1 in ${FASTQ_DIR}/*_1.fastq
do
    samp=$(basename "$r1" _1.fastq)
    r2=${FASTQ_DIR}/${samp}_2.fastq

    STAR \
      --runThreadN 4 \
      --genomeDir ${GENOME_DIR} \
      --readFilesIn ${r1} ${r2} \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${OUT_DIR}/${samp}_ \
      --quantMode GeneCounts
done