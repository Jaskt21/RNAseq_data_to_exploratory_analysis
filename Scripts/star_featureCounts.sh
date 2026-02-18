#!/bin/bash

#####################################################

#SBATCH --job-name=star_featureCounts
#SBATCH --output=featureCounts-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0-24:00:00

#####################################################

module load subread/2.0.6-plshecg

#####################################################

BAM_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_aligned
OUT_DIR=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_featureCounts
GTF=/home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/Homo_sapiens.GRCh38.115.gtf


for bam in ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam
do
    samp=$(basename "$bam" _Aligned.sortedByCoord.out.bam)


    featureCounts -p --countReadPairs -B \
  -T 4 -t exon -g gene_id -s 0 \
  -a ${GTF} \
  -o ${OUT_DIR}/${samp}_featureCounts_exon.txt \
  ${bam}

done