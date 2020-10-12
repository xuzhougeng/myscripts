#!/usr/bin/bash

set -e
set -u
set -o pipefail

if [ $# -lt 4 ]; then
	echo "Usage: "
	echo "  $0 genome lgsreads R1.fq.gz R2.fq.gz"
	exit 1
fi


genome=$1
lgsreads=$2
read1=$3
read2=$4
threads=100
round=3


NextPolish=/opt/biosoft/NextPolish-1.2.2
module load samtools/1.10
module load bwa/0.7.17
module load racon

if [ ! -f genome.lgspolish.fa ]; then
    minimap2 -x map-pb -t ${threads} ${genome} ${lgsreads} -o tmp.paf
	racon -t ${threads} ${lgsreads} tmp.paf $genome > tmp.fa
	minimap2 -x map-pb -t ${threads} tmp.fa ${lgsreads} > tmp.paf
	racon -t ${threads} ${lgsreads} tmp.paf tmp.fa > genome.lgspolish.fa
fi 

#Set input and parameters
input=genome.lgspolish.fa
for ((i=1; i<=${round};i++)); do
#step 1:
	#index the genome file and do alignment
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${read1} ${read2} | samtools view --threads 10 -F 0x4 -b - | samtools sort - -m 2g --threads 20 -o sgs.sort.bam;
	#index bam and genome files
	samtools index -@ 20 sgs.sort.bam;
	samtools faidx ${input};
	#polish genome file
	python $NextPolish/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
	input=genome.polishtemp.fa;
#step2:
	#index genome file and do alignment
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${read1} ${read2} | samtools view --threads 10 -F 0x4 -b - |samtools sort - -m 2g --threads 20 -o sgs.sort.bam;
	#index bam and genome files
	samtools index -@ 20 sgs.sort.bam;
	samtools faidx ${input};
	#polish genome file
	python $NextPolish/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
	input=genome.nextpolish.fa;
done;
#Finally polished genome file: genome.nextpolish.fa
