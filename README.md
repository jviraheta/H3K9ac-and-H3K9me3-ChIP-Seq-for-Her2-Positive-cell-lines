# H3K9ac-and-H3K9me3-ChIP-Seq-for-Her2-Positive-cell-lines

# ChIP-seq Analysis Pipeline for GEO-GSE85158

This repository provides a complete pipeline for processing and visualizing ChIP-seq data, including:

- **SRA Download**
- **Bowtie2 Indexing**
- **SAM to BAM Conversion**
- **GC Bias Calculation and Correction**
- **BAM to BigWig Conversion**
- **Visualization on IGV**
- **Coverage Analysis of Promoters**
- **Heatmap Generation for H3K9Ac and H3K9me3 Marks**
---

## Table of Contents

1. [SRA Download](#1-sra-download)
2. [Bowtie2 Indexing](#2-bowtie2-indexing)
3. [SAM to BAM Covnversion](#3-sam-to-bam-conversion)
4. [GC Bias Calculation and Correction](#4-gc-bias-calculation-and-correction)
    - [GC Bias Calculation](#gc-bias-calculation)
    - [GC Bias Correction](#gc-bias-correction)
5. [BAM to BigWig Conversion](#5-bam-to-bigwig-conversion)
6. [Visualization on IGV](#6-visualization-on-IGV)
7. [Coverage Analysis of Promoters](#7-coverage-analysis-of-promoters)
8. [Heatmap Generation for H3K9Ac and H3K9me3 Marks](#8-heatmap-generation-for-h3k9ac-and-h3k9me3-marks)

--- 
## 1. SRA Download
Download raw sequencing data (SRA files) and converts them to FASTQ format.

```bash
for ACC in "SRRYYYYY" ; do
    echo "processing $ACC"
    prefetch $ACC && fastq-dump --gzip $ACC
done
```

---

## 2. Bowtie2 Indexing 

Building a Bowtie2 index for the hg19 reference genome.

```bash
# Build the Bowtie2 index
bowtie2-build $REFERENCE_FILE $INDEX_PREFIX
```
FASTQ files were aligned to the hg19 genome using Bowtie2. Alignments were performed in local mode to allow soft clipping (--local) and base quality encoding was set to --phred33.

```bash
for FASTQ_FILE in $FASTQ_DIR/*.fastq.gz; do
    BASENAME=$(basename $FASTQ_FILE .fastq.gz)
    
     bowtie2 -p 8 --local -x $INDEX_PREFIX -U $FASTQ_FILE -S $OUTPUT_DIR/$BASENAME.sam \
        --phred33 --local
done
```

---
## 3. SAM to BAM conversion

Aligned reads were converted to BAM format using Samtools (version 1.9)

```bash
# Convert each SAM file to BAM
for SAM_FILE in $SAM_DIR/*.sam; do
    BASENAME=$(basename $SAM_FILE .sam)
    
    samtools view -@ 4 -b $SAM_FILE > $BAM_DIR/$BASENAME.bam
    samtools sort -@ 4 -o $BAM_DIR/$BASENAME.sorted.bam $BAM_DIR/$BASENAME.bam
    samtools index $BAM_DIR/$BASENAME.sorted.bam
done

```
---

## 4. GC Bias Calculation and Correction

Calculate and correct for GC bias in ChIP-seq data using DeepTools(version 3.5.1)

## GC Bias Calculation
The computeGCBias tool was used to calculate the GC content of the reads compared to the expected distribution. An effective genome size of 2,864,785,220 bases was set, with a fragment length of 35 bp.

```bash
# Compute GC bias
for BAM_FILE in $IN_DIR/*.sorted.bam; do
    BASENAME=$(basename $BAM_FILE .sorted.bam)
    echo "processing $BAM_FILE -> $BASENAME"
    computeGCBias -b $BAM_FILE -g $GENOME --effectiveGenomeSize $GENOME_SIZE --GCbiasFrequenciesFile $OUT_DIR/$BASENAME.freq.txt
done
```

---

## GC Bias Correction
The correctGCBias tool was used to adjust the GC bias in the BAM files.

```bash
# Correct GC bias
for BAM_FILE in $IN_DIR/*.sorted.bam; do
    BASENAME=$(basename $BAM_FILE .sorted.bam)
    echo "processing $BAM_FILE -> $BASENAME"
    correctGCBias -p 10 -b $BAM_FILE --effectiveGenomeSize $GENOME_SIZE -g $GENOME --GCbiasFrequenciesFile $FREQ_DIR/$BASENAME.freq.txt -o $OUT_DIR/$BASENAME.gcCorrected.bam
done
```

---

## 5. BAM to BigWig Conversion

After GC bias correction, we normalized the ChIP-seq data using bamCoverage from Deeptools. The normalization was done using Counts Per Million mapped reads (CPM) with a bin size of 10 bp and a smoothing length of 30 bp.

```bash
# Convert each sorted BAM file to BigWig
for BAM_FILE in $BAM_DIR/*.sorted.bam; do
    BASENAME=$(basename $BAM_FILE .sorted.bam)
    
    bamCoverage -b $BAM_FILE -o $BIGWIG_DIR/$BASENAME.bw --binSize 10 --normalizeUsing CPM --smoothLength 30 --effectiveGenomeSize 2700000000
done
```

--- 

##6. Visualization on IGV

Sorted BAM files were converted to BigWig format using the bamCoverage tool from Deeptools and loaded into IGV for visualization.

---
## 7. Coverage Analysis of Promoters

The bedtools map command (-o max) for each region was used to compute the maximum cpm value within +/- 5 kb of the transcription start site (TSS) for each gene.

```bash
# Analyze coverage over gene promoters
for FILE_PATH in compare_bam_cpm/*.bg.gz; do
    FILE="${FILE_PATH##*/}"
    BASE="${FILE%.bg.gz}"
    OUT=compare_bam_cpm_1kb_promoter/$BASE.geneList_promoter_5kb_max
    zcat $FILE_PATH | sort -k1,1 -k2,2n | gzip > compare_bam_cpm/$BASE.sort.bg.gz
    bedtools map -c 4 -o max -a ../gene_list.5k_promoter.bed -b compare_bam_cpm/$BASE.sort.bg.gz | gzip > $OUT.gz
done
```

---
## 8. Heatmap Generation for H3K9Ac and H3K9me3 Marks

R script generates a heatmap visualizing the maximum CPM values across selected genes for H3K9Ac and H3K9me3 histone marks.
Look at H3K9ac-H3K9me_heatmaps.R file

<img width="911" alt="Screenshot 2024-10-16 at 11 25 58 PM" src="https://github.com/user-attachments/assets/421c2190-f366-4387-aee5-403eb0f8cdff">



---
