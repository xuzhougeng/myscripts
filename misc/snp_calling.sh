#!/usr/bin/bash
# 

set -e
set -u
set -o pipefail

REF=$1
FQ1=$2
FQ2=$3


THREADS=30
source /etc/profile.d/modules.sh
# read align
module load bwa/0.7.17
module load samtools/1.10

mkdir -p tmp
mkdir -p index

# build index
if [ ! -f index/ref.bwt ]; then
    bwa index -p index/ref ${REF}
fi

if [ ! -f align.bam.bai ]; then
    bwa mem -v 1 -t $THREADS index/ref $FQ1 $FQ2 2> bwa_error.log |\
        samtools sort -T tmp/align -@ $(($THREADS / 2)) -O bam -o align.bam && \
    samtools index -@ $THREADS align.bam
fi

# mark duplication
if [ ! -f align_mkdup.bam ]; then
    sambamba markdup -t 10 align.bam align_mkdup.bam
fi

module load bcftools/1.10

bcftools mpileup -Ou -f ${REF}  align_mkdup.bam | bcftools call -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
