#!/usr/bin/bash

set -e
set -u
set -o pipefail

function usage {
    echo -e "Usage: $0 -t 20 -f 6 -w 21 target.fa query.fa"
	echo -e "Use option -h | --help for more information"
}

function help {
	usage;
	echo
	echo "Simple BLASTN Wrapper"
	echo "---------------"
	echo "OPTIONS"
	echo
	echo " -t|--thread thread"
	echo " -f|--outfmt blastn output format"
	echo " -w|--wordsize minimum word size to initialize a alignment"
	echo "[-h|--help]: help"

}


if [ $# -lt 1 ]
then
	usage
	exit
fi

for arg in "$@"; do
	shift
	case "$arg" in
		"--thread") set -- "$@" "-t" ;;
		"--outfmt") set -- "$@" "-f" ;;
		"--wordsize") set -- "$@" "-w" ;;
		"--help") set -- "$@" "-h" ;;
		*) set -- "$@" $arg
	esac
done


THREAD=20
OUTFMT=6
WORDSIZE=21

while getopts ":t:w:f:ch" OPT
do
    case $OPT in
        t) THREAD=$OPTARG;;
        f) OUTFMT=$OPTARG;;
		w) WORDSIZE=$OPTARG;;
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

if [ $# -lt 2 ]
then
	usage
	exit
fi

TARGET=$1
QUERY=$2

if [ ! -f $TARGET.nsq ]; then
    makeblastdb -in $TARGET -dbtype nucl
fi

blastn -query $QUERY -db $TARGET -outfmt $OUTFMT -num_threads $THREAD -word_size $WORDSIZE

