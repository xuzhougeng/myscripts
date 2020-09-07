#!/bin/bash

if [ $# -lt 2 ];then
	echo "$0 pb_asm.fa input.fofn"
	exit 1
fi

asm=$1
fofn=$2

if [ ! -f cutoffs ]; then
	cat $fofn | while read file;
    do
	    fn=$(basename $file)
		    minimap2 -t 100 -x map-pb $asm $file | pigz -c - > ${fn}_aln.paf.gz 
    done

	/opt/biosoft/purge_dups/bin/pbcstat *aln.paf.gz
	/opt/biosoft/purge_dups/bin/calcuts PB.stat > cutoffs 2> calcults.log
fi


# Split an assembly and do a self-self alignment
/opt/biosoft/purge_dups/bin/split_fa $asm > asm.split
minimap2 -t 20 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz

/opt/biosoft/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log

/opt/biosoft/purge_dups/bin/get_seqs dups.bed $asm
