#!/usr/bin/env python

import sys
import glob
from os.path import join

from Bio import SeqIO

# find all the multi-alignment sequence  and read into a single object
msa_dir = "MultipleSequenceAlignment"
msa_files = glob.glob( join(msa_dir, "**", "*final_align_NT.aln") )
msa_seqs = []
for msa in msa_files:
    msa_seqs += list(SeqIO.parse(msa, "fasta"))

msa_dict = SeqIO.to_dict(msa_seqs)

# get the orthgroups
orthgroups = open(sys.argv[1])
# get the species name
species_name = next(orthgroups).strip().split("\t")[1:]
# output seq
seqs = [ "" for i in range(len(species_name)) ]

# parse the orth groups
for line in orthgroups:
    items = line.strip().split("\t")
    seqs_temp = [ "" for i in range(len(items) - 1) ]
    OG_name = items[0]
    flag = 1
    for idx, gene_id in enumerate(items[1:]):
        seq = msa_dict.get(gene_id)
        if seq is not None:
            seqs_temp[idx] = seq.seq
        else:
            flag = 0
            print( OG_name + " is not complete", file = sys.stderr)
            break
    if flag != 0:
        for idx, seq in enumerate(seqs_temp):
            seqs[idx] = seqs[idx] + seqs_temp[idx]


# output 
seq_len = int(len(seqs[0]) / 3)
specie_num = len(species_name)

with open("input.seq", "w") as file:
    print(str(specie_num) + ' ' + str(seq_len), file=file)
    for idx, seq in enumerate(seqs):
        print(species_name[idx] + '\n' + seq[0::3], file=file)
    print(str(specie_num) + ' ' + str(seq_len), file=file)
    for idx, seq in enumerate(seqs):
        print(species_name[idx] + '\n' + seq[1::3], file=file)
    print(str(specie_num) + ' ' + str(seq_len), file=file)
    for idx, seq in enumerate(seqs):
        print(species_name[idx] + '\n' + seq[2::3], file=file)



