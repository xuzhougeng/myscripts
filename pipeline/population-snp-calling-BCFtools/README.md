To run this pipeline, you should provided the following files

- Refercence fasta
- BWA INDEX

```bash
bwa index ref.fa
```

- Reference Dict

```bash
gatk CreateSequenceDictionary -R ref.fa
```

- chr_len.txt : two columns, contig name and its length

```bash
seqkit fx2tab -nl ref.fa | sort -k 2,2nr > chr_len.txt
```
