# Snakemake pipeline for HiC alignment with bwa-aln and bwa-sampe

Edit the config.yaml for specifying the input

- draft.fasta
- allele.table.txt
- R1.fastq.gz, R2.fastq.gz

Usage:

```bash
snakemake 
```

Requirement:

- fastp
- Bwa
- Samtools >= 1.8
- ALLHiC

The output is bwa_aln.bam
