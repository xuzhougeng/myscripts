#!/bin/bash

set -e
set -u
set -o pipefail

bam=$1
width=20
depth=10
threads=10

# program not in PATH
spade=~/opt/biosoft/SPAdes-3.11.1-Linux/bin/spades.py

# find the unalignment position
if [ ! -f tmp.bed ]; then
    bedtools merge -i $bam -bed > tmp.bed
fi

# filter unalignment that caused by gap(N) or other genomic problem
if [ ! -f tmp_flt1.bed ]; then
    gawk -v width="$width" 'NR==1{chr=$1;a=$3}; NR>1{if ($2-a<width && chr==$1) print $1"\t"a"\t"$2"\t"$2-a;a=$3;chr=$1}' tmp.bed | grep  -v "Chr\(M\|Ch\)" > gap.bed
    gawk '{print $1"\t"$2-2"\t"$2}' gap.bed > left_gap.bed
    gawk '{print $1"\t"$3"\t"$3+2}' gap.bed > right_gap.bed
    samtools bedcov right_gap.bed $bam > right_gap_cov.bed
    samtools bedcov left_gap.bed $bam > left_gap_cov.bed
    paste left_gap_cov.bed right_gap_cov.bed  | \
        gawk -v depth="$depth" '{if ($4>depth && $8 >depth) print $0}' > insertion.bed

fi

samtools view -@ ${threads} -b $bam -L insertion.bed  | samtools sort -n |\
    samtools fasta -1 query_1.fa -2 query_2.fa -

# assembly contig
if [ ! -f assembly/contigs.fasta ];then
    mkdir -p assembly
    samtools view -@ ${threads} -b -F 0x2 $bam | samtools sort -n | \
    samtools fastq -1 reads_1.fq -2 reads_2.fq -
    ${spade} -1 reads_1.fq -2 reads_2.fq -o assembly
fi

# build blastdb
if [ ! -f assembly/contigs.fasta.nhr ];then
    makeblastdb -dbtype nucl -in assembly/contigs.fasta
fi

blastn -query query_1.fa -db assembly/contigs.fasta -outfmt 6 -qcov_hsp_perc 100 > query_1.blastn.out

