#!/bin/bash

set -e
set -u
set -o pipefail

HIC_PRO=~/HiC-Pro_2.11.1

function usage {
    echo -e "usage: hic_pro_pre.sh -e MboI -i input.fa -p prefix"
	echo -e "Use option -h|--help for more information"
}

function help {
    usage;
	echo
	echo "Generate HiC-Pro input with fasta and enzyme"
	echo "---------------"
	echo "OPTIONS"
	echo
	echo "	-i|--input input.fa"
	echo "	-e|--enzyme enzyme(MboI, HindI) or enzyme string"
	echo "	-p|--prefix output file prefix" 
	echo "	-t|--threads bowtie-build threads " 
	echo " [-h|--help]: help"
	exit;

}

if [ $# -lt 1 ]
then
	usage
	exit
fi

# Transform long options to short ones
for arg in "$@";do
	shift
	case "$arg" in
		"--input") set -- "$@" "-i" ;;
		"--enzyme") set -- "$@" "-e" ;;
		"--prefix") set -- "$@" "-p" ;;
		"--threads") set -- "$@" "-t" ;;
		"--help") set -- "$@" "-h" ;;
		*) set -- "$@" "$arg"
	esac
done

fa=""
prefix=""
enzyme=""
threads=20

while getopts ":i:e:p:ch" OPT
do
	case $OPT in
		i) fa=$OPTARG;;
		e) enzyme=$OPTARG;;
		p) prefix=$OPTARG;;
		t) threads=$OPTARG;;
		h) help ;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			usage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			usage
			exit 1
			;;
	esac
done

echo 
echo "-------------------------------"
echo Fasta is $fa 
echo Preifx is $prefix
echo Enzyme is $enzyme
echo "------------------------------"
echo


if [[ -z $fa || -z $prefix || -z $enzyme ]];then
	usage
	exit
fi

# build enzyme bed
$HIC_PRO/bin/utils/digest_genome.py -r $enzyme -o ${prefix}_${enzyme}.bed $fa &

# chromosome size
seqkit fx2tab -nl $fa | awk '{print $1"\t"$2}' > ${prefix}.chrom.size &

# bowtie2 index
bowtie2-build --threads $threads $fa $prefix
