#trimming the adaptors

for f in /XW1/XW1_CKDL240034851-1A_22HVYCLT4_L4_2; do
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.trim.fastq.gz /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/1rawdata/01.RawData/$f.fq.gz > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.cutadpt.log
done

for f in /XW2/XW2_CKDL240034852-1A_22HVYCLT4_L4_2; do
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 \
-a ATTGCTTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a ACAAGCCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.trim.fastq.gz /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/1rawdata/01.RawData/$f.fq.gz > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.cutadpt.log
done

for f in /XW3/XW3_CKDL240034853-1A_22HVYCLT4_L4_2; do
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 \
-a AACTTGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a AGGACCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.trim.fastq.gz /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/1rawdata/01.RawData/$f.fq.gz > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.cutadpt.log
done

#collapse duplicates

for f in mpxr.IP1 mpxr.IP2
do
perl /usr/local/ctk/fastq2collapse.pl /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/2trimming/$f.fastq.gz - | gzip -c > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/3collapse_duplicates/$f.trim.c.fastq.gz
done

#strip randomers

for f in mpxr.IP1 mpxr.IP2; do
perl /usr/local/ctk/stripBarcode.pl -format fastq -len 10 /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/3collapse_duplicates/$f.trim.c.fastq.gz - | gzip -c > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/4strip_randomer/$f.trim.c.tag.fastq.gz
done

#mapping

#!/bin/bash
#
#BATCH -N 4
#SBATCH -t 0-24:00
#SBATCH --cpus-perp-task=4
#SBATCH --mem=100g

module load gcc/8.2.0
module load samtools/1.9
module load bamtools/2.5.1
module load bwa/0.7.17

for f in mpxr.input mpxr.IP1 mpxr.IP2; do
bwa aln -t 4 -n 0.06 -q 20 /bgfs/dyang/RNA-seq/reference/human/hg38/hg38_ucsc /bgfs/dyang/xiaofei/mPXR_LS174T/4strip_randomer/$f.trim.c.tag.fastq.gz > /bgfs/dyang/xiaofei/mPXR_LS174T/5mapping/$f.sai
bwa samse /bgfs/dyang/RNA-seq/reference/human/hg38/hg38_ucsc /bgfs/dyang/xiaofei/mPXR_LS174T/5mapping/$f.sai /bgfs/dyang/xiaofei/mPXR_LS174T/4strip_randomer/$f.trim.c.tag.fastq.gz | gzip -c > /bgfs/dyang/xiaofei/mPXR_LS174T/5mapping/$f.sam.gz


done

#parsing

for f in mpxr.input mpxr.IP1 mpxr.IP2 ; do
perl /usr/local/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/6parsing/$f.mutation.txt /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/5mapping/$f.sam.gz -| gzip -c > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/6parsing/$f.tag.bed.gz
done

#removing rRNAs

for f in mpxr.input mpxr.IP1 mpxr.IP2 ; do
perl /usr/local/ctk/tagoverlap.pl -big -region /usr/local/ctk/annotation/genomes/hg38/annotation/rmsk.RNA.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/6parsing/$f.tag.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/7removing_rRNA/$f.tag.norRNA.bed
done

#removing snoRNAs

for f in mpxr.input mpxr.IP1 mpxr.IP2; do
perl /usr/local/ctk/tagoverlap.pl -big -region /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/8removing_snoRNA/snoRNA_snRNA_scaRNA_gencodev22.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/7removing_rRNA/$f.tag.norRNA.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/8removing_snoRNA/$f.tag.norRNA.nosnoRNA.bed
done

#collapse duplicates 2

for f in mpxr.input mpxr.IP1 mpxr.IP2; do
perl /usr/local/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment\
 -weight --weight-in-name --keep-max-score --keep-tag-name /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/8removing_snoRNA/$f.tag.norRNA.nosnoRNA.bed\
  /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/$f.R2.fastq.tag.uniq.bed
done

#peakcalling

awk -vFS="\t" -vOFS="\t" '{if($6=="+") print $0}' /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP1.R2.fastq.tag.uniq.bed | sort -k1,1 -k2,2n > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP1.fastq.tag.uniq.pos.bed
awk -vFS="\t" -vOFS="\t" '{if($6=="-") print $0}' /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP1.R2.fastq.tag.uniq.bed | sort -k1,1 -k2,2n > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP1.fastq.tag.uniq.neg.bed

awk -vFS="\t" -vOFS="\t" '{if($6=="+") print $0}' /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP2.R2.fastq.tag.uniq.bed | sort -k1,1 -k2,2n > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP2.fastq.tag.uniq.pos.bed
awk -vFS="\t" -vOFS="\t" '{if($6=="-") print $0}' /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP2.R2.fastq.tag.uniq.bed | sort -k1,1 -k2,2n > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP2.fastq.tag.uniq.neg.bed

macs2 callpeak --keep-dup all --nomodel -n mpxr.IP.pos -t /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP1.fastq.tag.uniq.pos.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP2.fastq.tag.uniq.pos.bed 
macs2 callpeak --keep-dup all --nomodel -n mpxr.IP.neg -t /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP1.fastq.tag.uniq.neg.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.IP2.fastq.tag.uniq.neg.bed

#annotation 

/home/xiaofei/miniconda3/share/homer-4.9.1-6/bin/annotatePeaks.pl /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.LS174T.IP.bed hg38 > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/11annotation/mpxr.LS174T.IP.annotation.xls

#badgarph

perl /usr/local/ctk/tag2profile.pl -v -ss -exact -of bedgraph -n "mPXR.IP" /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP.R2.fastq.tag.uniq.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/12bedgraph/LS174T.mPXR.IP.R2.tag.uniq.bedgraph

perl /usr/local/ctk/tag2profile.pl -v -ss -exact -of bedgraph -n "mPXR.input" /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.input.R2.fastq.tag.uniq.bed /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/12bedgraph/LS174T.mPXR.input.R2.tag.uniq.bedgraph


#distribution

PATH=$PATH:/home/xiaofei/miniconda3/share/homer-4.9.1-6/bin/

annotatePeaks.pl /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.LS174T.IP.bed hg38 -ann /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/new_analysis/LS174T_WT_Ctrl/14distribution_homer/annotations_1000TSSTTS.3utrfirst.txt > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/distribution_homer/mpxr.LS174T.annot.txt

#counting reads in peaks

bedtools intersect -a /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/10peakcalling/mpxr.LS174T.IP.sorted.bed -b /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/9collapse_duplicates2/mpxr.IP2.R2.fastq.tag.uniq.bed -c > /home/xiaofei/Desktop/Documents/xiaofei/pxr_eclip/mPXR_LS174T/readcount/mpxr.IP2.count.bedgraph
