# define working dir paths
datadir = "/home/ubuntu/workspace/rnaseq/expression/htseq_counts"
outdir = "/home/ubuntu/workspace/rnaseq/de/htseq_counts/deseq2"

# load R libraries we will use in this section
library(DESeq2)
library(data.table)
library(ggplot2)

# set working directory to data dir
setwd(datadir)

# read in the RNAseq read counts for each gene (produced by htseq-count)
htseqCounts = fread("gene_read_counts_table_all_final.tsv")

# set working directory to the output dir
setwd(outdir)
