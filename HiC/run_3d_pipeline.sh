#!/bin/bash

set -e
set -u
set -o pipefail

if [ $# -lt 1 ]; then
	echo "$0 input.fa"
	exit 1
fi


FASTA=$1
module load bwa/0.7.17

if [ ! -f reference/genome.fa.bwt ]; then
    mkdir -p reference
    seqkit seq -w 80 $FASTA > reference/genome.fa
    cd reference
    bwa index genome.fa &> /dev/null &
    python /opt/biosoft/juicer/misc/generate_site_positions.py DpnII genome genome.fa
    awk 'BEGIN{OFS="\t"}{print $1, $NF}' genome_DpnII.txt > genome.chrom.size
    cd ..
fi

wait

/opt/biosoft/juicer/scripts/juicer.sh \
	-g genome \
	-s MboI \
	-z reference/genome.fa \
	-y reference/genome_DpnII.txt \
	-p reference/genome.chrom.size \
	-D /opt/biosoft/juicer \
	-t 100 &> juicer.log 

#bash /opt/biosoft/3d-dna/run-asm-pipeline.sh -r 2 reference/genome.fa aligned/merged_nodups.txt &> 3d.log 
