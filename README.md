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
      [View Script FastQC](FASTQC.sh)
    - **MultiQC:** then takes all those individual reports and combines them into one summary dashboard, allowing you to compare all 10 samples simultaneously  
      [View Script MultiQC](MULTIQC.sh)

**3. Alignment, Mapping & Quantification**  
- Tools: salmon/1.10.2, star/2.7.11, subread/2.0.6 (Salmon, STAR)
- Context: The goal at this stage is to transform raw unordered sequencing reads (FASTQ files) into a structured count matrix. Whether utilising a traditional aligner like STAR or a pseudo-aligner like Salmon, the objective remains the same: to identify the genomic origin of each read and convert these digital signals into numerical values (counts) per gene. This "Count Matrix" serves as the essential entry point for all downstream statistical analysis, providing the raw data required for identifying differentially expressed genes (DGE) in R.
    - **Salmon (Pseudoalignment):** Salmon serves as a pseudoaligner, it involves breaking the cDNA transcriptome into small, fixed-length fragments called k-mers to create a searchable index.  Salmon then simply checks which transcripts are "compatible" with a read's k-mer signature, rather than the base-by-base alignment. This results in .sf files that contain both "Estimated Counts" for statistical testing and "TPM" values, which are counts adjusted for gene length and sequencing depth to allow for easy comparison between samples.
 
      **Indexing:** We build a Salmon Index using the transcriptome (mRNA sequences) to create a searchable map of all known transcripts   
      `salmon index -t Homo_sapiens.GRCh38.cdna.all.fa -i salmon_index`

      **Quantification:** Using the salmon quant command, the software performs quasi-mapping. Instead of finding the exact base to base alignment, it quickly determines which transcript a read likely belongs to  
      [View Script Salmon Quantification](salmon_script2.sh)

    - **STAR (Traditional Alignment):** STAR workflow provides a high-resolution, spatial alignment by mapping reads to the entire reference genome. Using the genomic DNA (FASTA) and a GTF annotation file. STAR aligns raw reads to the reference genome to produce BAM files, which record the exact coordinates of every sequence fragment. Since BAM files only contain locations, featureCounts must then cross-reference these coordinates against a GTF annotation to generate the final gene-count matrix.
 
       **Indexing:** STAR uses the genomic DNA (FASTA) and GTF annotation to build a searchable "map" of the genome. This index allows the software to quickly find the coordinates of millions of short reads  
     [View Script Star Index](star_indexing.sh)

      **Alignment:** The STAR aligner maps raw reads to the genome, producing a BAM file    
      [View Script Star Alignment](star_mapping2.sh)

       **Quantification via featureCounts:** Since a BAM file is just a list of locations, featureCounts is used to tally them up. It cross-references the BAM coordinates with the gene boundaries in the GTF file to determine which gene each read belongs to. Creating a txt file containing Gene IDs and their corresponding raw integer counts, which can be used for downstream analysis  
      [View Script Star Quantifictaion](star_featureCounts.sh)

      
      

  

      
      




