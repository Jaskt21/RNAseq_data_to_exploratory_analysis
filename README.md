# RNAseq_data_to_exploratory_analysis

## 🧬 Analysis of GSE81089 

### Project Overview 
This repository contains a pipeline for processing of raw RNA-seq data from the NCBI Gene Expression Omnibus (GEO).  
**Note:** For the purpose of demonstrating this pipeline, a subset of 10 samples was selected from the primary study GSE81089.

### Biological Context
The original study aimed to investigate the expression of cancer-testis antigens (CTAs) in non-small cell lung cancer (NSCLC) using RNA-seq. CTAs are genes that are normally restricted to testicular tissue but can be aberrantly expressed in cancer. 
In this project, the dataset is used only to demonstrate RNA-seq data processing and analysis steps, without reproducing the biological analyses of the study.

### Dataset Information
**GEO Accession:** GSE81089  
**Organism:** _Homo sapiens_  
**Experiment type:** Expression profiling by high throughput sequencing (RNA-seq)  
**Library layout:** Paired-end

## ⚒️ Upstream Processing Pipeline   
The following steps outline the pipeline used to process raw SRA data into count matrices, in preparation for downstream exploratory data analysis and differential expression in R.

**1. Data Retrival and Extraction**  
- Tools: sra-tools/3.0.3 (prefetch, fasterq-dump)
- Context: Data on GEO is stored in a highly compressed .sra format to save space.
    - **prefetch:** prefetch is used to download the data securely from the NCBI servers using the unique Run IDs  
      [View Script- prefetch](sra_script.sh)
    - **fasterq-dump:** then converts these files (.sra) into raw .fastq format  
      [View Script- fasterq-dump](fastq_dump.sh)

**2. Quality Control**  
- Tools: fastqc/0.12.1, py-multiqc/1.28 (FastQC, MultiQC)
- Context: The quality control process specifically involves checking for adequate read quality and verifying the overall GC content distribution of the dataset in a html output.
    - **FastQC:** produces an individual report for each sample you run it on  
      [View Script- FastQC](FASTQC.sh)
    - **MultiQC:** then takes all those individual reports and combines them into one summary dashboard, allowing you to compare all 10 samples simultaneously  
      [View Script- MultiQC](MULTIQC.sh)

**3. Alignment, Mapping & Quantification**  
- Tools: salmon/1.10.2, star/2.7.11, subread/2.0.6 (Salmon, STAR)
- Context: The goal at this stage is to transform raw unordered sequencing reads (FASTQ files) into a structured count matrix. Whether utilising a traditional aligner like STAR or a pseudo-aligner like Salmon, the objective remains the same: to identify the genomic origin of each read and convert these digital signals into numerical values (counts) per gene. This "Count Matrix" serves as the essential entry point for all downstream statistical analysis, providing the raw data required for identifying differentially expressed genes (DGE) in R.
    - **Salmon (Pseudoalignment):**  




