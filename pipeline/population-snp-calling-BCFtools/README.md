
# Requirement

- Snakemake
- BWA
- SAMTools/1.9
- BCFTools/1.9

# How to use it

Step1: copy this pipeline to the work directory 

```bash
cd /path/to/work
# copy
cp -r ~/myscripts/pipeline/population-snp-calling-BCFtools/{Snakefile,config.yaml} .
```

Step2: prepare the reference and its index

```bash
# ref diretory
bwa index ref.fa
samtools faidx ref.fa
```

Step3: prepare the chr_len.txt,  two columns, contig name and its length

```bash
cut -f 1,2 ref.fa.fai > chr_len.txt
#seqkit fx2tab -nli ref.fa | sort -k 2,2nr > chr_len.txt
```

Step4: prepare the samples.tsv

```bash
mkdir -p analysis/00-raw-data
# build the softlink to analysis/00-raw-data
# Fastq name format: xxx_1.fq.gz xxx_2.fq.gz
ls -1 analysis/00-raw-data | sed 's/_[12].fq.gz//' | sort | uniq > samples.tsv
```

Step5: modify the configure file `vim config.yaml`

```bash
# 参考基因组和索引的位置
FASTA: ref/ref.fa
INDEX: ref/ref.fa

#是否去重复,GBS设置False,WGS设置True
DEDUPLICATE: False

# 染色体长度信息
CHRFILE: ref/chr_len.txt
# 染色体拆分
VCFJOBS: 50

# 样本名信息
samples: samples.tsv
depth: 20
PLOIDY: 2
```

Step6: Test and run pipeline

```bash
# test
snakemake -np
# run
# -j cores for parallel
snakemake -j 20
```

