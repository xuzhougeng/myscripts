#!/usr/bin/env python3

# Convert the 3d-dna assembly format to human readable format

from sys import argv

fn = argv[1]

contig_id = {}
contig_len = {}
scaffolds = []

for line in open(fn):
    line = line.strip()
    if line.startswith(">"):
        seq_name, seq_id, seq_len = line.split(' ')
        contig_id[seq_id] = seq_name[1:]
        contig_len[seq_id] = seq_len
    else:
        scaffolds.append(line)

HiC_order = 1

for line in scaffolds:
    for seq_id in line.split(' '):
        if seq_id.startswith('-'):
            seq_name = contig_id[seq_id[1:]]
            seq_len  = contig_len[seq_id[1:]]
            orient = "-"
        else:
            seq_name = contig_id[seq_id]
            seq_len = contig_len[seq_id]
            orient = "+"
        print("HiC_scaffold_{}\t{}\t{}\t{}".format(str(HiC_order), seq_name, orient, seq_len))
    HiC_order += 1
    
