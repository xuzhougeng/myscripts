#!/bin/bash

# For more information, please visit
# https://hicexplorer.readthedocs.io/en/latest/content/example_usage.html

set -e
set -u
set -o pipefail

if [ $# -lt 1 ]; then
	echo "$0 input.fa read_R1.fq.gz read_R2.fq.gz"
	exit 1
fi


FASTA=$1
R1=$2
R2=$3
ENZYME=GATC
THREADS=50

# requirement: hicexplorer-3.5.1
# conda activate /opt/sysoft/miniconda3/envs/hic
which hicexplorer
if [ $? -eq 1 ];then
	echo "hicexplorer is not installed"
	exit 1
fi

module load bwa/0.7.17
module load samtools/1.10

# build index
if [ ! -f reference/genome.fa.bwt ]; then
    mkdir -p reference
    seqkit seq -w 80 $FASTA > reference/genome.fa
    cd reference
    bwa index genome.fa &> /dev/null &
    hicFindRestSite --fasta genome.fa --searchPattern $ENZYME -o ${ENZYME}.bed
    cd ..
fi

wait

# align
bwa mem -t $THREADS -A1 -B4 -E50 -L0 reference/genome.fa $R1 2> bwa_mem_R1.log | \
    samtools view -Shb - > bwa_mem_R1.bam
bwa mem -t $THREADS -A1 -B4 -E50 -L0 reference/genome.fa $R2 2> bwa_mem_R2.log | \
    samtools view -Shb - > bwa_mem_R2.bam

# build matrix from independently mated read pairs
# the restriction sequence GATC is recognized by the DpnII restriction enzyme

hicBuildMatrix --samFiles bwa_mem_R1.bam bwa_mem_R2.bam \
                 --binSize 10000 \
                 --restrictionSequence $ENZYME \
                 --danglingSequence $ENZYME \
                 --restrictionCutFile reference/${ENZYME}.bed \
                 --threads $THREADS \
                 --inputBufferSize 100000 \
                 --outBam hic.bam \
                 -o hic_matrix.h5 \
                 --QCfolder ./hicQC

#correct Hi-C matrix
hicCorrectMatrix diagnostic_plot -m hic_matrix.h5 -o hic_corrected.png
hicCorrectMatrix correct -m hic_matrix.h5 --filterThreshold -1.5 5 -o hic_corrected.h5
