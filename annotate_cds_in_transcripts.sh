#!/bin/bash

set -e
set -u
set -o pipefail

ref=`realpath $1`
gtf=`realpath $2`
db=`realpath $3`

echo
echo "Processing"
echo reference: $ref
echo gtf: $gtf
echo db: $db
echo

# set the environment
unset PERL5LIB
module load TransDecoder/5.5
export PATH=$PATH:/opt/biosoft/diamond

# get the fasta from GTF
mkdir -p tmp

if [ ! -f tmp/transcripts.fasta ]; then
    gtf_genome_to_cdna_fasta.pl $gtf $ref > tmp/transcripts.fasta
fi

# convert gtf to gff3
if [ ! -f tmp/transcripts.gff3 ]; then
    gtf_to_alignment_gff3.pl $gtf > tmp/transcripts.gff3
fi

# predict the ORF in transcripts
cd tmp
TransDecoder.LongOrfs -t transcripts.fasta

# search orf in Uniport database
diamond blastp -d $db -q transcripts.fasta.transdecoder_dir/longest_orfs.pep \
	--evalue 1e-5 --max-target-seqs 1 > blastp.outfmt6

# predict the coding region
TransDecoder.Predict \
    -t transcripts.fasta \
    --retain_blastp_hits blastp.outfmt6 

# generate the annotation file
cdna_alignment_orf_to_genome_orf.pl \
    transcripts.fasta.transdecoder.gff3 \
    transcripts.gff3 \
    transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

cd ..
mv tmp/transcripts.fasta.transdecoder.genome.gff3 .
mv tmp/transcripts.fasta.transdecoder.cds . 
rm -rf tmp
