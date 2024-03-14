## 使用方法

第一步: 将项目克隆到本地

```bash
git clone https://github.com/xuzhougeng/myscripts.git
```

第二步: 创建分析环境

```bash
# recommend mamba 
conda create -c bioconda -c conda-forge -n atac-seq fastp sambamba bowtie2 samtools deeptools macs2  htslib snakemake
```

保证能够调用 fastp, sambamba, bowtie2, samtools, deeptools, macs2即可, 例如，如果管理可以允许使用module管理加载

```bash
module load bowtie/2.3.4.3
module load samtools/1.10
module load deeptools/2.0
module load MACS/2.1.2
```

另外还需要编译我写的两个工具

```bash
conda activate atac-seq
gcc -o bam_basic_stats bam_basic_stats.c -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lhts -Wl,-rpath,$CONDA_PREFIX/lib
gcc -o bam_frag_stats bam_frag_stats.c -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lhts -Wl,-rpath,$CONDA_PREFIX/lib
```

第三步: 编辑config.yaml

```yaml
FASTA: 参考基因组fasta的位置, 
GFF: 参考基因组gff所对应的位置，为空时，不绘制gene body相关图
GENOMESIZE: 基因组大小, 例如 1.1e8 表示110M, 不能为空
ADAPTER: "CTGTCTCTTATACACATCT"  接头序列
SAMPLES: samples.txt 记录样本信息
MAP_QUALITY: 10 最低比对质量
REMOVE_PLASTOME: True/False/或者不包含质体的序列编号的txt文件
```

FASTA对应的fasta文件名中间不要有多余的 "."，例如 A.thaliana.fa, 请改成 Athaliana.fa或 A_thaliana.fa

REMOVE_PLASTOME:

如果输入的是False或者空白, 表示不过滤

如果输入的是数字，则表示认为低于 X bp都是潜在的细胞器， 推荐200000

一些质体的测序文章

- https://www.frontiersin.org/articles/10.3389/fpls.2018.00493/full
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7076112

如果输入的是 bed 文件，那么就会只保留在 bed 文件里面的序列，每行一个，例如

```txt
chr1    0   1000000
chr2    0   2000000
```

其中 samples.tsv, 文件格式如下

```text
id,sample,rep,fq1,fq2
r1_1,s1,1,/path/to/s1.fastq.gz,/path/to/s2.fastq.gz
r1_2,s2,2,/path/to/s2.fastq.gz,/path/to/s2.fastq.gz
r2_1,s3,1,/path/to/s3.fastq.gz,/path/to/s3.fastq.gz
r2_2,s4,2,/path/to/s4.fastq.gz,/path/to/s4.fastq.gz
```

id表示样本的编号，必须**唯一**, sample和rep分别表示样本编号和重复信息, fq1和fq2表示的fastq的位置

第四步: 运行snakemake

```bash
# 使用默认的config.yaml 测试
snakemake -np
# 使用默认的config.yaml
snakemake -j 120
# 使用给定的config文件
snakemake -j 120 --configfile config_diy.yaml -s Snakefile
```