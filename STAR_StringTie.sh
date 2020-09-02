#!/bin/bash

set -e
set -u
set -o pipefail

# Load STAR 
#STAR=/opt/biosoft/STAR-2.6.1c/bin/Linux_x86_64/STAR
module load STAR/2.6.1.d

reference=$1
sample=$2

THREADS=30

# building index for alignment
index=STAR
mkdir -p ${index}

if [ ! -f ${index}/SA ]
then
    ${STAR} \
    --runThreadN $threads \
    --runMode genomeGenerate \
    --genomeDir $index \
    --genomeFastaFiles ${reference}
fi

# align RNA-seq to reference
exec 0< $2
while read prefix
do
	fq1=${prefix}_R1.fq.gz
	fq2=${prefix}_R2.fq.gz
    STAR \
    	--genomeDir $index \
    	--runThreadN $THREADS \
    	--readFilesIn $fq1 $fq2 \
    	--readFilesCommand zcat \
    	--outFileNamePrefix ${prefix}_ \
    	--outSAMtype BAM SortedByCoordinate \
    	--outBAMsortingThreadN $THREADS \
    	--outSAMstrandField intronMotif \
    	--outFilterIntronMotifs RemoveNoncanonical
done
#if [ ! -f merged.bam ] ;then
#	samtools merge merged.bam *sortedByCoord.out.bam
#fi
