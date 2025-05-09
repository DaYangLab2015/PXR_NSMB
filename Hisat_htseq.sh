#!/bin/bash


# hisat2 mapping / HTSeq quantification pipeline
# usage: from an empty working directory, run
# ./Hisat_htseq.sh (read1) (read2) (output_dir)

# input: gzipped fastq file read1 & read2 for paired-end
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read2
output_dir=$3 ## output directory

# output: all in the working directory, fixed names
# Aligned.sam					# alignments, strangeness RF sam
# Aligned.srt.bam				# BAM sorted by name
# Aligned.srt.filtered.bam			# BAM sorted by name filter by q 30
# Aligned.Introns.counts.txt			# Counting reads that map to intronic segments of each gene
# Aligned.consExons.counts.txt			# Counting reads that map to exonic segments of each gene

#### setting output directory
mkdir $output_dir
cd $output_dir

#### alignment using hisat2 (mapping reads) for ### pair-end strand-specific seq data
## stranded RF

echo hisat2 -t -p 4 -x /bgfs/dyang/RNA-seq/XiaofeiWang/grch38/genome -1 $read1 -2 $read2 -S Aligned.sam --rna-strandness RF &> stdout.log

hisat2 -t -p 4 -x /bgfs/dyang/RNA-seq/XiaofeiWang/grch38/genome -1 $read1 -2 $read2 -S Aligned.sam --rna-strandness RF &>> stdout.log

###### strand specific detection (package rseqc)

echo infer_experiment.py -r /bgfs/dyang/RNA-seq/XiaofeiWang/Ensembl/Homo_sapiens.GRCh38.103.bed -i Aligned.sam &>> stdout.log

infer_experiment.py -r /bgfs/dyang/RNA-seq/XiaofeiWang/Ensembl/Homo_sapiens.GRCh38.103.bed -i Aligned.sam &>> stdout.log


#### .sam file to .bam file, sorting by names for pair-end strand-specific seq data, filter -q score 30

samtools sort -n -o Aligned.srt.bam Aligned.sam

samtools view -b -F 1548 -q 30 -o Aligned.srt.filtered.bam Aligned.srt.bam

#### Counting reads that map to intronic or exonic segments of each gene

htseq-count -m union -f bam -t intron -r pos -s reverse Aligned.srt.filtered.bam /bgfs/dyang/RNA-seq/XiaofeiWang/Ensembl/Homo_sapiens.GRCh38.103.chr.Introns.gtf > Aligned.Introns.counts.txt

htseq-count -m intersection-strict -f bam -t exon -r pos -s reverse Aligned.srt.filtered.bam /bgfs/dyang/RNA-seq/XiaofeiWang/Ensembl/Homo_sapiens.GRCh38.103.chr.consExons.gtf > Aligned.consExons.counts.txt




