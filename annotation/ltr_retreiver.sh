#!/bin/bash
set -e
set -u
set -o pipefail

REF=$1
THREADS=80
RATE=7e-9

if [ ! -f ltr.finder.scn ]; then

/opt/biosoft/LTR_Finder/source/ltr_finder \
	-D 15000 \
	-d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 \
	$REF > ltr.finder.scn &

fi

if [ ! -f ltr.harvest.scn ] ;then

/opt/biosoft/genometools-1.5.10/bin/gt suffixerator \
	-db $REF -indexname $REF \
	-tis -suf -lcp -des -ssp -sds -dna

/opt/biosoft/genometools-1.5.10/bin/gt ltrharvest \
	-index $REF \
	-similar 90 -vic 10 -seed 20 -seqids yes \
	-minlenltr 100 -maxlenltr 7000 -mintsd 4 \
	-maxtsd 6 -motif TCGA -motifmis 1 > ltr.harvest.scn &

fi

wait

num1=$(grep '^\[' ltr.finder.scn | wc -l)
num2=$(grep -v '^#' ltr.harvest.scn | wc -l)

echo ""
echo "LTR of LTR_FINDER: $num1"
echo "LTR of LTR_harvest: $num2"

/opt/biosoft/LTR_retriever/LTR_retriever \
	-genome $REF  \
	-inharvest ltr.harvest.scn \
	-infinder ltr.finder.scn \
	-threads $THREADS -u $RATE
