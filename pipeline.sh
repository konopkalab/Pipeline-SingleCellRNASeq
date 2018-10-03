#!/bin/bash

# Raw sequencing data was acquired from the McDermott center next generation sequencing core at UT Southwestern Medical Center in the form of binary base call (BCL) files. 
# BCL files were then de-multiplexed with the 10X Genomics i7 index (used during library preparation) using Illumina’s bcl2fastq v2.17.1.14 [1] and mkfastq command from 10X Genomics CellRanger v2.1.1 [2] tools.
cellranger mkfastq --run=<path-to-bcl-files> --csv=<samplesheet.csv> --use-bases-mask Y26n,I8n,Y124n

# Example 'samplesheet.csv' as below
# Lane,Sample,Index
# 1,SAMPLE1,SI-GA-G1
# 1,SAMPLE2,SI-GA-G2
# 1,SAMPLE3,SI-GA-G3
# 2,SAMPLE1,SI-GA-G1
# 2,SAMPLE2,SI-GA-G2
# 2,SAMPLE3,SI-GA-G3
# 3,SAMPLE1,SI-GA-G1
# 3,SAMPLE2,SI-GA-G2
# 3,SAMPLE3,SI-GA-G3
# 4,SAMPLE1,SI-GA-G1
# 4,SAMPLE2,SI-GA-G2
# 4,SAMPLE3,SI-GA-G3


# Extracted paired-end fastq files (26 bp long R1 - cell barcode and UMI sequence information, 124 bp long R2 - transcript sequence information) were checked for read quality using FASTQC v0.11.5 [3].
fastqc <input-R1.fastq.gz>
fastqc <input-R2.fastq.gz>


# R1 reads were then used to estimate and identify real cells using whitelist command from UMI-tools v0.5.4 [4] program. 
umi_tools whitelist --stdin=<input-R1.fastq.gz> --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --expect-cells=<number-of-expected-cells> --plot-prefix=<prefix> --log=<prefix-whitelist.out> --error=<prefix-whitelist.err> --stdout=<prefix-whitelist.txt>


# A whitelist of cell-barcodes (putative real cells) and R2 fastq files were later used to extract reads corresponding to real cells only (excluding sequence information representing empty beads, doublets, low quality/degrading cells, etc.) using extract command from UMI-tools v0.5.4 [4]. 
# This step also appends the cell-barcode and UMI sequence information from R1 to read names in R2 fastq file. 
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin=<input-R1.fastq.gz> --stdout=<extracted-R1.fastq.gz> --read2-in=<input-R2.fastq.gz> --read2-out=<extracted-R2.fastq.gz> --filter-cell-barcode --whitelist=<prefix-whitelist.txt> --log=<prefix-extract.out> --error=<prefix-extract.err>


# Extracted R2 reads were then aligned to reference genome (HG38/GRCh38p12 for Human, MM10/GRCm38p6 for Mouse) from UCSC genome browser [5] and reference annotation (Gencode v28 for human [6], Gencode vM17 for mouse [7]) using STAR aligner v2.5.2b [8] allowing upto 5 mismatches. 
STAR --runThreadN 48 \
	 --genomeDir <star-genome-index> \
     --readFilesIn <extracted-R2.fastq.gz> \
     --readFilesCommand zcat \
     --sjdbGTFfile <gtf-annotation> \
     --outFilterType BySJout  \
     --outFilterMismatchNoverReadLmax 0.04 \
     --outFilterMultimapNmax 10 \
     --alignSJoverhangMin 10 \
     --alignSJDBoverhangMin 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMismatchNmax 5 \
     --twopassMode Basic \
     --outFileNamePrefix <output-prefix>


# Uniquely mapped reads were then assigned to exons (for single cell datasets) or genes (for single nuclei datasets) using featureCounts program from Subread package (v1.6.2) [9]. 
featureCounts \
	--primary \
	-R BAM \
	-T 48 \
	-s 1 \
	-t <exon-or-gene> \
	-g <gene-name-or-gene-id> \
	-a <gtf-annotation> \
	-o <prefix-primary-assigned> \
	<star-out-bam>


# Assigned reads sorted and indexed using Samtools v1.6 [10] 
samtools sort <feature-counts-out-bam> -o <aligned-assigned-sorted>
samtools index <aligned-assigned-sorted>


# Assigned, sorted and indexed reads were then used to generate raw expression UMI count tables using count command from UMI-tools v0.5.4 [4] program. 
umi_tools count --per-gene --gene-tag=XT --per-cell --stdin=<aligned-assigned-sorted.bam> --stdout=<prefix-counts.tsv.gz> --log=<prefix-counts.log> --error=<prefix-counts.err> --wide-format-cell-counts


# This raw expression matrix contains cells as rows and genes as columns and can be further used for downstream analysis such as normalization, clustering, differentially expressed genes, etc.

# References
#1.	bcl2fastq. Illumina Inc.
#2.	CellRanger. 10X Genomics.
#3.	Andrews, S., FastQC - A quality control tool for high throughput sequence data. Babraham Bioinformatics.
#4.	Smith, T., A. Heger, and I. Sudbery, UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res, 2017. 27(3): p. 491-499.
#5.	Kent, W.J., et al., The human genome browser at UCSC. Genome Res, 2002. 12(6): p. 996-1006.
#6.	Gencode Genes v28. 2017.
#7.	Gencode Genes vM17. 2017.
#8.	Dobin, A., et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013. 29(1): p. 15-21.
#9.	Liao, Y., G.K. Smyth, and W. Shi, featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 2014. 30(7): p. 923-30.
#10.	Li, H., et al., The Sequence Alignment/Map format and SAMtools. Bioinformatics, 2009. 25(16): p. 2078-9.



