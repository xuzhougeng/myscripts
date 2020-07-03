#!/bin/bash
set -e
set -u
set -o pipefail

REF=$1
INDEX=$2
SAMPLES=$3
SPECIES=Arabidopsis_thaliana
THREADS=80

module load bwa/0.7.17
module load samtools/1.9

mkdir -p 1-align
mkdir -p 2-snp-calling
mkdir -p tmp
mkdir -p log

outdir=1-align
snpdir=2-snp-calling

exec 0< $SAMPLES
# FQ to BAM
while read SM
do
	echo "Processing ${SM}"
	if [ ! -f ${outdir}/${SM}_mkdup.bam ]; then
	bwa mem -v 1 -t ${THREADS} -R "@RG\\tID:${SM}\\tSM:${SM}\\tPL:illumina" \
		$INDEX ${SM}_R1.fq.gz ${SM}_R2.fq.gz 2> log/${SM}.log | \
		samtools sort -@  $((THREADS/2))  > ${outdir}/${SM}_sort.bam &&
		sambamba markdup -t 10 ${outdir}/${SM}_sort.bam ${outdir}/${SM}_mkdup.bam
    fi 
done

# calling variants per chromosome
module load GATK/4.0.9.0 

if [ ! -f ${REF%%.fa}.dict ]; then
    gatk CreateSequenceDictionary -R ${REF}
fi

if [ ! -f ${REF}.fai ]; then
    samtools faidx ${REF}	 
fi
chroms=($(grep '>' $REF |sed 's/>//' | tr '\n' ' '))

for chr in ${chroms[@]}
do
    if [ ! -f ${snpdir}/gatk-hc.${chr}.vcf.gz ]; then
        gatk HaplotypeCaller -R $REF  \
            $(for bam in 1-align/*_mkdup.bam  ; do echo "-I $bam";done) \
            --genotyping-mode DISCOVERY \
            --intervals ${chr} -stand-call-conf 30 --sample-ploidy 2 \
            -O ${snpdir}/gatk-hc.${chr}.vcf.gz &
    fi
done && wait

## merge chr.vcf
if [ ! -f ${snpdir}/gatk.hc.vcf.gz ]; then
        merge_vcfs=""
        for chr in ${chroms[@]}; do
            merge_vcfs=${merge_vcfs}" -I ${snpdir}/gatk-hc.${chr}.vcf.gz"
        done 
        gatk MergeVcfs ${merge_vcfs} -O ${snpdir}/gatk-hc.vcf.gz && echo "Vcfs haved been successfully merged" 
fi

## Hard filter for SNP
gatk SelectVariants \
    -select-type SNP \
    -V ${snpdir}/gatk-hc.vcf.gz \
    -O ${snpdir}/gatk-hc.snp.vcf.gz

gatk VariantFiltration \
    -V ${snpdir}/gatk-hc.snp.vcf.gz \
    -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "LowQual" -O ${snpdir}/gatk-hc.snp.filter.vcf.gz

gatk SelectVariants \
    --exclude-filtered true \
    -V ${snpdir}/gatk-hc.snp.filter.vcf.gz \
    -O ${snpdir}/gatk-hc.final.snp.vcf.gz

# Annotation Using snpEff
if [ 1 -lt 2 ] ;then
java -Xmx4g -jar /opt/biosoft/snpEff/snpEff.jar ann \
    -o gatk ${SPECIES} \
    ${snpdir}/gatk-hc.final.snp.vcf.gz > ${snpdir}/gatk-hc.final.snp.ann.vcf

fi
