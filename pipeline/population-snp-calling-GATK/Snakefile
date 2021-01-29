# multiple sample alignment and variant calling script shell.prefix("set -eo pipefail; echo BEGIN at $(date);")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

# import python package
from os.path import join

def build_interval(chr_len_file, threads = 50, debug=False):
    """Build interval based on threads
    """

    from collections import defaultdict 
    # get the contig length and contig size
    chr_lens = open(chr_len_file, "r").readlines()
    chr_num = len(chr_lens)
    # initize the list for store 
    chr_name = [ '' for i in range(chr_num) ] 
    chr_len  = [ '' for i in range(chr_num) ]

    for i in range(chr_num):
        cn,cl = chr_lens[i].strip().split()
        chr_name[i] = cn
        chr_len[i] = int(cl)

    total_len = sum(chr_len)
    interval_len = total_len // threads 
    interval_dict = defaultdict(list)
    interval_name_dict = defaultdict(list)
    j = 0
    for i in range(chr_num):
        cur_len = interval_dict[j]
        if cur_len == []:
            interval_dict[j].append(chr_len[i])
            interval_name_dict[j].append(chr_name[i])
        elif sum(cur_len) < interval_len:
            interval_dict[j].append(chr_len[i])
            interval_name_dict[j].append(chr_name[i])
        else:
            j += 1
            interval_dict[j].append(chr_len[i])
            interval_name_dict[j].append(chr_name[i])
    if debug == True:
        print("Debug Info")
        for i in range(j):
            print("{}: {} \t {} ".format(i, sum(interval_dict[i]), " ".join(interval_name_dict[i] )))

    return interval_name_dict


# configuration
#configfile: "config.yaml"
FASTA    = config['FASTA']
INDEX    = config['INDEX']
PLOIDY   = config['PLOIDY']
CHRFILE  = config['CHRFILE']
VCFJOBS  = config['VCFJOBS']
SAMPLES  = [line.strip() for line in open(config['samples'])] # samples
ADAPTER  = ""

# get the interval name dict
interval_name_dict = build_interval(CHRFILE, VCFJOBS)
# How to deal with duplication
DEDUPLICATE = config['DEDUPLICATE']
readfilter = "--disable-read-filter NotDuplicateReadFilter"

if DEDUPLICATE:
    readfilter = "--read-filter NotDuplicateReadFilter"

# make workding directory
RAW_DIR   = join("analysis", "00-raw-data")
CLEAN_DIR = join("analysis", "01-clean-data")
ALIGN_DIR = join("analysis", "02-read-alignment")
GVCF_DIR  = join("analysis", "03-gvcf-calling")
GENOMEDB  = join("analysis", "04-genomeDB")
VCF_DIR   = join("analysis", "05-variant-calling")
QC_DIR    = join("report", "fastp")
TMPDIR    = join("analysis", "temp")
LOGDIR    = join("analysis", "log")

shell("mkdir -p {CLEAN_DIR} {ALIGN_DIR} {GVCF_DIR} {GENOMEDB} {VCF_DIR}")
shell("mkdir -p {QC_DIR} {TMPDIR} {LOGDIR}")

ALL_BAM  = [join(ALIGN_DIR, "{}_markdup_sorted.bam".format(sample)) for sample in SAMPLES]
ALL_DB   = [join(GENOMEDB, "{}".format(chr)) for chr in list(interval_name_dict.keys()) ]
ALL_VCF  = [join(VCF_DIR, "{}.vcf.gz".format(chr)) for chr in list(interval_name_dict.keys()) ]
FINAL_VCF = join(VCF_DIR, "merged.vcf.gz")
#VCFS     = join(VCF_DIR, "freebayes.var.vcf") # FreeBayes


rule all:
    input:
        FINAL_VCF
        #ALL_BAM, VCFS

rule data_clean:
    input:
        r1 = join(RAW_DIR, "{sample}_1.fq.gz"),
        r2 = join(RAW_DIR, "{sample}_2.fq.gz")
    params:
        adapter = ADAPTER,
        length  = 50,
        quality = 20,
        prefix  = lambda wildcards: join(QC_DIR, wildcards.sample)
    output:
        r1 = join(CLEAN_DIR, "{sample}_1.fq.gz"),
        r2 = join(CLEAN_DIR, "{sample}_2.fq.gz")
    log: join(LOGDIR, "{sample}_fastp.log")
    threads: 4
    message: "----- processing {input.r1} and {input.r2} with fastp ------"
    shell: """
        fastp $(if [ -n "{params.adapter}" ]; then echo '-a {params.adapter}'; fi;) \
        -w {threads} -j {params.prefix}.json -h {params.prefix}.html -l {params.length} -q {params.quality}\
        -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}
    """

# map to genome
rule bwa_mem:
    input:
        r1 = join(CLEAN_DIR, "{sample}_1.fq.gz"),
        r2 = join(CLEAN_DIR, "{sample}_2.fq.gz")
    params:
        index  = INDEX,
        rg     = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        tmpdir = TMPDIR,
        logdir = LOGDIR,
        split  = join(ALIGN_DIR,"{sample}_split.sam"),
        disc   = join(ALIGN_DIR,"{sample}_disc.sam")
    output:
        join(ALIGN_DIR,"{sample}_markdup_sorted.bam"),
    log:
        log1=join(LOGDIR, "{sample}_align_warning.log"),
        log2=join(LOGDIR, "{sample}_samblaster.log")
    threads: 20
    message: "align {input.r1} and {input.r2} to {params.index} with {threads} threads"
    run:
        if DEDUPLICATE:
            shell("""
            module load bwa/0.7.17
            bwa mem -v 2 -t {threads} -R '{params.rg}' {params.index} {input.r1} {input.r2} 2> {log.log1} |\
            samblaster --acceptDupMarks --addMateTags --splitterFile {params.split} --discordantFile {params.disc} 2> {log.log2} |\
            sambamba view -S -f bam -l 0 /dev/stdin |\
            sambamba sort -t {threads} -m 2G --tmpdir {params.tmpdir} -o {output[0]} /dev/stdin
            """)
        else:
            shell("""
            module load bwa/0.7.17
            bwa mem -v 2 -t {threads} -R '{params.rg}' {params.index} {input.r1} {input.r2} 2> {log.log1} |\
            sambamba view -S -f bam -l 0 /dev/stdin |\
            sambamba sort -t {threads} -m 2G --tmpdir {params.tmpdir} -o {output[0]} /dev/stdin
            """)

# Split Reference For GATK to
rule singleGVCF:
    input:
        join(ALIGN_DIR,"{sample}_markdup_sorted.bam"),
    params:
        ref       = FASTA,
        gatkfilter= readfilter,
        interval = lambda wildcards: interval_name_dict[int(wildcards.chr)]
    output:
       join(GVCF_DIR, "{sample}.{chr}.g.vcf.gz")
    threads: 6 
    shell:"""
    module load GATK/4.1.2.0
	gatk HaplotypeCaller -R {params.ref} \
	    --genotyping-mode DISCOVERY \
		--input {input} \
	    $(for val in {params.interval}; do echo "--intervals $val"; done) \
		-stand-call-conf 30 --sample-ploidy 2 -ERC GVCF \
		{params.gatkfilter} \
		--output {output}
	"""

rule buildGenomicsDB:
    input:
        expand(join(GVCF_DIR, "{sample}.{chr}.g.vcf.gz"), sample = SAMPLES, chr = "{chr}")
    params:
         interval = lambda wildcards: interval_name_dict[int(wildcards.chr)]
    output:
        directory(join(GENOMEDB, "{chr}" ))
    threads: 4 
    shell:"""
	module load GATK/4.1.2.0
	gatk GenomicsDBImport --genomicsdb-workspace-path {output} \
	   $(for val in {input}; do echo "--variant $val"; done) \
	   $(for val2 in {params.interval}; do echo "--intervals $val2"; done)
	"""

rule SingleGenotype:
    input:
        directory(join(GENOMEDB, "{chr}" ))
    params:
        ref       = FASTA,
    output:
        join(VCF_DIR, "{chr}.vcf.gz")
    threads: 6 
    shell:"""
	module load GATK/4.1.2.0
	gatk GenotypeGVCFs -R {params.ref} -V gendb://{input} -O {output}
    """


rule MergeVCFS:
    input:
        ALL_VCF
    output:
        FINAL_VCF
    threads: 10
    shell:"""
	module load bcftools/1.9 
	bcftools concat --threads {threads} {ALL_VCF} -Oz -o {FINAL_VCF}
	#module load GATK/4.1.2.0
    #gatk MergeVcfs  \
	#   $(for val in {input}; do echo "-I $val"; done) \
	#      -O {output}
    """

# Using FreeBayes to calling variants
rule freebayes:
    input: ALL_BAM
    params:
        ref    = FASTA,
        ploidy = PLOIDY,
        cov    = 10 
    output: join(VCF_DIR, "freebayes.var.vcf")
    threads: 100
    shell:"""
        module load freebayes
        freebayes-parallel <(fasta_generate_regions.py {params.ref}.fai 100000) {threads} \
            -p {params.ploidy} \
            --standard-filters \
            --min-coverage {params.cov} \
            -f {params.ref} {input} > {output}
    """
# filter the results
