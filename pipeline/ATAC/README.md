
使用方法

第一步: 将项目克隆到本地

```bash
git clone https://github.com/xuzhougeng/myscripts.git
```

第二步: 创建分析环境

```bash
conda create -c bioconda -c conda-forge -n atac-seq fastp sambamba bowtie2 samtools deeptools macs2 
```

或者保证能够调用 fastp, sambamba, bowtie2, samtools, deeptools, macs2即可, 利用通过module管理加载

```bash
module load bowtie/2.3.4.3
module load samtools/1.10
module load deeptools/2.0
module load MACS/2.1.2
```

第三步: 编辑config.yaml

```yaml
FASTA: 参考基因组fasta的位置
GFF: 参考基因组gff所对应的位置，为空时，不绘制gene body相关图
GENOMESIZE: 基因组大小, 例如 1.1e8 表示110M, 不能为空
ADAPTER: "CTGTCTCTTATACACATCT"  接头序列
SAMPLES: samples.txt 记录样本信息
```

第四步: 创建 samples.tsv, 文件格式如下

```text
id,sample,rep,fq1,fq2
r1_1,s1,1,/path/to/s1.fastq.gz,/path/to/s2.fastq.gz
r1_2,s2,2,/path/to/s2.fastq.gz,/path/to/s2.fastq.gz
r2_1,s3,1,/path/to/s3.fastq.gz,/path/to/s3.fastq.gz
r2_2,s4,2,/path/to/s4.fastq.gz,/path/to/s4.fastq.gz
```

id表示样本的编号，必须**唯一**, sample和rep分别表示样本编号和重复信息, fq1和fq2表示的fastq的位置

第三步: 运行snakemake

```bash
#测试
snakemake -np
# 运行
snakemake -j 120
```