#####################################################################################################

# Script: 01_import_star.R
# Purpose: TEMPLATE- import STAR ReadsPerGene.out.tab and build a DESeq2 object
#
# NOTE: This is a reference/framework script. in this repository. The paths below are placeholders,
# update paths and the design formula to match your own data and metadata before running.
#####################################################################################################

library(data.table)
library(DESeq2)

# Paths
quant_dir <- "data/star_output"        
meta_path <- "metadata/samples.csv"  

# Locate STAR count files
files <- Sys.glob(file.path(quant_dir, "*ReadsPerGene.out.tab"))

# Recover sample names, and name the file vector
sample_names <- gsub("_ReadsPerGene.out.tab", "", basename(files))
names(files) <- sample_names
 
# Merge count files into a single count matrix 
# Column 4 = reverse-stranded counts
countData <- data.frame(fread(files[1]))[c(1, 4)]
 
for (i in 2:length(files)) {
  countData <- cbind(countData, data.frame(fread(files[i]))[4])
}
 
# Remove STAR's first 4 summary rows
countData <- countData[5:nrow(countData), ]
 
colnames(countData) <- c("GeneID", sample_names)
rownames(countData) <- countData$GeneID
countData <- countData[, -1]
 
# Load sample metadata 
metadata <- read.csv(meta_path, row.names = 1)  
 
# Ensure count matrix columns match metadata row order (required by DESeq2)
countData <- countData[, rownames(metadata)]

# Build DESeq2 object 
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = metadata,
  design    = ~ condition   # update to match your metadata column name
)
 
# Save output
saveRDS(dds, "results/dds_raw.rds")
