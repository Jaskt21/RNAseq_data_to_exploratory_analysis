###################################################################################################################

# Script: 01_import_tximport.R
# Purpose: Template to import salmon quant.sf files to build a raw DESeq2 object- reference script
# Note: Update the placeholder paths and design formula in this framework script to match your data before running

###################################################################################################################

library(tximport)
library(DESeq2)

# Paths (placeholders - update to match project)

data_dir     <- "data/salmon_quants"  
meta_path    <- "metadata/samples.csv"
tx2gene_path <- "reference/tx2gene.csv"

# Load Sample metadata
metadata <- read.csv(meta_path, stringsAsFactors = FALSE)
metadata$quant_file <- file.path(data_dir, metadata$sample_id, "quant.sf")

# Named vector required by tximport: file paths named by sample ID
files <- setNames(metadata$quant_file, metadata$sample_id)

# Load transcript to gene reference table
tx2gene <- readRDS(tx2gene_path)
 
# Import Salmon quantification
txi <- tximport(
  files           = files,
  type            = "salmon",
  tx2gene         = tx2gene,
  ignoreTxVersion = TRUE
)
 
# Build DESeq2 object 
metadata$condition <- factor(metadata$condition) 
 
dds <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = metadata,
  design  = ~ condition
)
 
# Save output
saveRDS(dds, "results/dds_raw.rds")
