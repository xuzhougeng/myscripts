#!/usr/bin/env python3

"""
convert GFF3 to GTF in ENSEMBLE fomrats

It will only use the following features

- gene
- mRNA -> transcript
- CDS
- exon
- five_prime_UTR -> five_prime_utr
- three_prime_UTR -> three_prime_utr

TO DO:add additional two features

- start_codon
- stop_codon

I only test this script in Arabidopsis_thaliana.TAIR10.44.gff3.gz
"""

from sys import argv
import re

for line in open(argv[1]):
    line = line.strip()
    if line.startswith("###"):
        continue
    if line.startswith("#"):
        print(line)
        continue
    fields = line.split("\t")
    feature = fields[2]
    if feature == "gene":
        gene_id = re.findall(r'gene_id=(.*?);',fields[8])[0]
        fields[8] = "gene_id \"{}\";gene_name \"{}\";".format(gene_id, gene_id)
    elif feature == "mRNA" or feature == "transcript":
        tx_id = re.findall(r'transcript_id=(.*?)$',fields[8])[0]
        gene_id = tx_id.split(".")[0]
        fields[2] = "transcript"
        fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
    elif feature == "CDS":
        tx_id = re.findall(r'protein_id=(.*?)$',fields[8])[0]
        gene_id = tx_id.split(".")[0]
        fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
    elif feature == "exon":
        exon_id = re.findall(r'exon_id=(.*?);',fields[8])[0]
        tx_id = ".".join(exon_id.split(".")[0:2])
        gene_id = tx_id.split(".")[0]
        fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
    elif feature == "five_prime_UTR" or feature == "five_prime_utr" or feature == "5UTR":
        fields[2] = "five_prime_utr"
        tx_id = fields[8].split(":")[1]
        gene_id = tx_id.split(".")[0]
        fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
    elif feature == "three_prime_UTR" or feature == "three_prime_utr" or feature == "3UTR":
        fields[2] = "three_prime_utr"
        tx_id = fields[8].split(":")[1]
        gene_id = tx_id.split(".")[0]
        fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
    else:
        continue
    print("\t".join(fields))
