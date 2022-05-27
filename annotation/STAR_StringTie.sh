#!/bin/bash

set -e
set -u
set -o pipefail

# Usage

# Load STAR 
#STAR=/opt/biosoft/STAR-2.6.1c/bin/Linux_x86_64/STAR
module load STAR/2.6.1.d
STAR=`which STAR`

reference=$1
sample=$2

# building index for alignment
ulimit -n 4096
index=STAR
mkdir -p ${index}

if [ ! -f ${index}/SA ]
then
    ${STAR} \
    --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeDir $index \
    --genomeFastaFiles ${reference}
fi
# create outdir
mkdir -p 01-clean-data
mkdir -p 02-read-align

# align RNA-seq to reference
exec 0< $2
while read prefix
do
	fq1=00-raw-data/${prefix}_1.fq.gz
	fq2=00-raw-data/${prefix}_2.fq.gz
	fastp -w 8 -i $fq1 -I $fq2 -o 01-clean-data/${prefix}_1.fq.gz -O 01-clean-data/${prefix}_2.fq.gz \
		-j 01-clean-data/${prefix}.json -h 01-clean-data/${prefix}.html
    STAR \
    	--genomeDir $index \
    	--runThreadN 20 \
    	--readFilesIn 01-clean-data/${prefix}_1.fq.gz 01-clean-data/${prefix}_2.fq.gz \
    	--readFilesCommand zcat \
    	--outFileNamePrefix 02-read-align/${prefix}_ \
    	--outSAMtype BAM SortedByCoordinate \
    	--outBAMsortingThreadN 40 \
    	--outSAMstrandField intronMotif \
    	--outFilterIntronMotifs RemoveNoncanonical \
		--alignIntronMax 10000

done
#if [ ! -f merged.bam ] ;then
#	samtools merge merged.bam *sortedByCoord.out.bam
#fi
