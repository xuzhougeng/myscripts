#!/usr/bin/bash

set -e
set -u
set -o pipefail

# usage:
# ./maker2apollo.sh prefix genome.fa annotation.gff3

ds_index=$1
genome=$2
gff=$3

# get the annotation from datastore index of maker
# - est alignment (stringtie is est_gff:stringtie)
# - protein alignments
# - repeat alignments

gff3_merge -s -n -d $ds_index > tmp.gff3
awk '{ if ($2 ~ "est") print $0 }' tmp.gff3 | sort -k1,1V -k4,4n >  maker.est2genome.gff3
awk '{ if ($2 ~ "protein2genome") print $0 }' tmp.gff3 | sort -k1,1V -k4,4n > maker.protein2genome.gff3
awk '{ if ($2 ~ "repeat") print $0 }' tmp.gff3 | sort -k1,1V -k4,4n >  maker.repeats.gff3
# delete the tmp file
rm tmp.gff3 

# create GMOD output directory
mkdir -p $PWD/GMOD/jbrowse_data
mkdir -p $PWD/GMOD/postgres_data
mkdir -p $PWD/GMOD/apollo_data 

# generate the output
# module load jbrowse
prepare-refseqs.pl --fasta $genome --out $PWD/GMOD/jbrowse_data 

# gene
flatfile-to-json.pl --gff $gff \
    --trackLabel "gene" --key "Gene spans" \
    --className "feature5" --subfeatureClasses '{"mRNA": "mRNA"}'\
    --type "gene" --out $PWD/GMOD/jbrowse_data &

# MAKER
flatfile-to-json.pl --gff $gff \
    --trackLabel "maker" --key "Transcripts" \
    --className "transcript" \
    --subfeatureClasses '{"exon": "exon", "CDS": "CDS", "five_prime_UTR": "five_prime_UTR", "three_prime_UTR": "three_prime_UTR"}' \
    --type "mRNA" --out $PWD/GMOD/jbrowse_data &

# EST
flatfile-to-json.pl --gff maker.est2genome.gff3 \
    --trackLabel "est2genome" --key "ESTs" \
    --className "generic_parent" --subfeatureClasses '{"match_part": "est2genome_part"}' \
    --type "expressed_sequence_match" --out $PWD/GMOD/jbrowse_data &

# Protein
flatfile-to-json.pl --gff maker.protein2genome.gff3 \
    --trackLabel "protein2genome" --key "proteins" \
    --className "generic_parent" --subfeatureClasses '{"match_part": "protein2genome_part"}' \
    --type "protein_match" --out $PWD/GMOD/jbrowse_data  &

# RepeatMasker and RepeatRunner
flatfile-to-json.pl --gff maker.repeats.gff3 \
    --trackLabel "RepeatMaskr" --key "RepeatMaskr" \
    --className "generic_parent" --subfeatureClasses '{"match_part": "repeat_part"}' \
    --type "match" --out $PWD/GMOD/jbrowse_data  &
    
wait;
echo "done"
