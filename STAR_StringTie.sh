#!/bin/bash

set -e
set -u
set -o pipefail

# Load STAR 
#STAR=/opt/biosoft/STAR-2.6.1c/bin/Linux_x86_64/STAR
module load STAR/2.6.1.d

index=$1
fq1=$2
fq2=$3
prefix=$4

reference=$(realpath $1)
outdir=$(realpath $2)
threads=20

THREADS=30

mkdir -p ${outdir}

if [ ! -f ${outdir}/SA ]
then
    ${STAR} \
    --runThreadN $threads \
    --runMode genomeGenerate \
    --genomeDir $outdir \
    --genomeFastaFiles ${reference}
fi

STAR \
	--genomeDir $index \
	--runThreadN $THREADS \
	--readFilesIn $fq1 $fq2 \
	--readFilesCommand zcat \
	--outFileNamePrefix $prefix\
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN $THREADS \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical

if [ ! -f merged.bam ] ;then
	samtools merge merged.bam *sortedByCoord.out.bam
fi
