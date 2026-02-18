#!/bin/bash

#####################################################

#SBATCH --job-name=star_index
#SBATCH --output=star_index-%j.out
#SBATCH --partition=cpu
#SBATCH --qos=qos_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0-24:00:00

#####################################################

module load star/2.7.11a-pgsk3s4

#####################################################

STAR \
--runMode genomeGenerate \
--runThreadN 4 \
--genomeDir /home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/star_index \
--genomeFastaFiles /home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /home/jkt21/mnt/network/bioinformatics/users/jkt21/rnaseq_GSE81089_workflow/star/Homo_sapiens.GRCh38.115.gtf \
--sjdbOverhang 100
