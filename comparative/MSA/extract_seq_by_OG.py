#!/usr/bin/env python3


import sys
from Bio import SeqIO

# 需要CDS编号需要和SCO的编号一样
SCO = sys.argv[1]
CDS = sys.argv[2]

outdir = "./SingleCopyGenes/"

seqs = SeqIO.to_dict(SeqIO.parse(CDS, "fasta"))

with open(SCO, "r") as file:
    for line in file:
        line = line.strip()
        rec = line.split("\t")
        OG_id = rec[0]
        gene_id = rec[1:]
        out_seq = [ seqs[gene] for gene in gene_id ]
        SeqIO.write(out_seq, outdir + OG_id + ".fa", "fasta") 
