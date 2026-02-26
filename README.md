# RNAseq_data_to_exploratory_analysis

## Analysis of GSE81089 

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

## Upstream Processing Pipeline   
The following steps outline the pipeline used to process raw SRA data into count matrices, in preparation for downstream exploratory data analysis and differential expression in R.

**1. Data Retrival and Extraction**  
- Tools: sra-tools/3.0.3 (prefetch, fasterq-dump)
- Context: Data on GEO is stored in a highly compressed .sra format to save space.
    - **prefetch:** prefetch is used to download the data securely from the NCBI servers using the unique Run IDs  
      [View Script: prefetch](Scripts/sra_script.sh)
    - **fasterq-dump:** then converts these files (.sra) into raw .fastq format  
      [View Script: fasterq-dump](Scripts/fastq_dump.sh)

      What do we see in fastq files:
      - Line 1: The sequence identifier  
        ID, Instrument info (HWI-) , Floecell Details, Coordinates (lane and tile location), Read length
      - Line 2: The raw sequences  
        This is the actual biological data a string of nucleotide bases (A, C, G, T)
      - Line 3: The Separator (+)  
        A simple plus sign used to separate the sequence from the quality data
      - Line 4: Quality Scores  
        The string of symbols (like B7BFFBF) represents the Phred Quality Score for every single base in Line 2    

**2. Quality Control**  
- Tools: fastqc/0.12.1, py-multiqc/1.28 (FastQC, MultiQC)
- Context: The quality control process specifically involves checking for adequate read quality and verifying the overall GC content distribution of the dataset in a html output.
    - **FastQC:** produces an individual report for each sample you run it on  
      [View Script: FastQC](Scripts/FASTQC.sh)
    - **MultiQC:** then takes all those individual reports and combines them into one summary dashboard, allowing you to compare all 10 samples simultaneously  
      [View Script: MultiQC](Scripts/MULTIQC.sh)

      <ins>What does our multiQC output html contain?</ins>    
      General stats:   
      Dups (sequence duplications): 20-60%, is normal however higher percentage of duplication of 80-100% could be due to a number of reasons:
      - PRC over amplification
      - Low library complexity  
      
      GC Content:
      - Percentage of bases that are G or C
      - This is usually around 40-45%
      - All samples should have a similar percentage
      - A big shift in the GC content could be due to cross contamination or biases
      
      Seq:
      - Total number of reads sequenced
      - Paired end _1 and _2 should have identical counts

      **We can also use MultiQC to aggregate and visualise alignment statistics from STAR and quantification metrics from Salmon**
      [View Script: Full MultiQC](Scripts/multiQC_full.sh)    
      
      
**3. Alignment, Mapping & Quantification**  
- Tools: salmon/1.10.2, star/2.7.11, subread/2.0.6 (Salmon, STAR)
- Context: The goal at this stage is to transform raw unordered sequencing reads (FASTQ files) into a structured count matrix. Whether utilising a traditional aligner like STAR or a pseudo-aligner like Salmon, the objective remains the same: to identify the genomic origin of each read and convert these digital signals into numerical values (counts) per gene. This "Count Matrix" serves as the essential entry point for all downstream statistical analysis, providing the raw data required for identifying differentially expressed genes R.
  
    - **Salmon (Pseudoalignment):** Salmon serves as a pseudoaligner, it involves breaking the cDNA transcriptome into small, fixed-length fragments called k-mers to create a searchable index.  Salmon then determines which transcripts are compatible with a read’s k-mer signature rather than performing base-by-base alignment.This results in .sf files that contain both Estimated Counts for statistical testing and TPM values, which are counts normalised for gene length and sequencing depth to allow for comparison between samples.
      
      **Indexing:** We build a Salmon Index using the transcriptome (mRNA sequences) to create a searchable map of all known transcripts   
      `salmon index -t Homo_sapiens.GRCh38.cdna.all.fa -i salmon_index`

      **Quantification:** Using the salmon quant command, the software performs quasi-mapping. Instead of finding the exact base to base alignment, it quickly determines which transcript a read likely belongs to  
      [View Script: Salmon Quantification](Scripts/salmon_script2.sh)  

      **Output files:**
      - quant.sf: This is the primary results file. It is a tab-separated text file containing the quantification estimates for each transcript, including columns like Name, Length, EffectiveLength, TPM (Transcripts Per Million), and NumReads.
      - cmd_info.json: A JSON file that records the exact command-line parameters and options you used to run Salmon. It’s very useful for reproducibility.
      - lib_format_counts.json: Contains details about the sequencing library format. It lists how many fragments were compatible with the expected library type (e.g., stranded vs. unstranded) and helps check if the data matches your assumptions.
      - logs/: Contains the log files for the run, which record the progress, warnings, and any errors that occurred during the execution.  

      
    - **STAR (Traditional Alignment):** The STAR workflow provides high-resolution, spatial alignment by mapping reads to the reference genome using genomic DNA (FASTA) and a GTF annotation file. STAR aligns raw reads to the genome to produce BAM files, which record the exact genomic coordinates of each sequence fragment. Since BAM files only contain alignment locations, featureCounts must then cross reference these coordinates against a GTF annotation to generate the final gene count matrix.  
 
       **Indexing:** STAR uses the genomic DNA (FASTA) and GTF annotation to build a searchable "map" of the genome. This index allows the software to quickly find the coordinates of millions of short reads  
     [View Script: Star Index](Scripts/star_indexing.sh)

      **Alignment:** The STAR aligner maps raw reads to the genome, producing coordinate sorted BAM files      
      [View Script: Star Alignment](Scripts/star_mapping2.sh)

       **1) Quantification via featureCounts:** Since a BAM file is a list of aligned read coordinates, featureCounts is used to quantify gene expression. It cross references BAM alignments with gene features defined in the GTF file to determine which gene each read overlaps. This produces a .txt file containing gene IDs and their corresponding raw integer counts, which are suitable for downstream analysis  
      [View Script: Star Quantifictaion (featureCounts)](Scripts/star_featureCounts.sh)

       **2) Quantification via --quantMode GeneCounts:** Instead of using tool like featureCounts, STAR can perform quantification during the alignment step itself using the --quantMode GeneCounts. This requires the same STAR Index and FASTQ files.  

        -Output: This process generates a specific output file named ReadsPerGene.out.tab  
      -Structure: The file contains four columns: Gene ID, unstranded counts, counts for the 1st read strand, and counts for the 2nd read strand  
      -R input: The ReadsPerGene.out.tab file is the primary input for downstream differential expression analysis in RStudio  
       [View Script: Star Quantifictaion (quantMode](Scripts/star_quantMode.sh)    
      
      
      
      
      
      

      
      

  

      
      




