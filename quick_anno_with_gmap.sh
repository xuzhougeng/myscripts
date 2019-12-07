#!/bin/bash

set -e
set -u
set -o pipefail

# Required
# - GMAP
# - gffread

function usage {
	echo -e "usage: $0 ref cds"
}

if [ $# -lt 1 ];then
	usage
	exit
fi

ref=$1
cds=$2

threads=40

module load GMAP/20180704

gmap_build -D tmp -d tmp $ref

gmap -t $threads -D tmp -d tmp -f gff3_gene $cds > cds_genes.gff3

gffread cds_genes.gff3 -g $ref -x cds_genes.cds
