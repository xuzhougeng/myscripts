#!/usr/bin/env python
# simple2links

from sys import argv

simple_file = argv[1]

ref_bed = simple_file.split(".")[0] + ".bed"
qry_bed = simple_file.split(".")[1] + ".bed"

ref_dict = {line.split("\t")[3]:line.split("\t")[0:3] for line in open(ref_bed)}
qry_dict = {line.split("\t")[3]:line.split("\t")[0:3] for line in open(qry_bed)}

fo = open(simple_file + "_link.txt", "w")

for line in open(simple_file):
    if line.startswith("#"):
        continue
    items = line.strip().split("\t") 
    ref_start_gene = items[0]
    ref_end_gene = items[1]
    qry_start_gene = items[2]
    qry_end_gene = items[3]
    
    ref_chr, ref_start = ref_dict[ref_start_gene][0:2]
    ref_end = ref_dict[ref_end_gene][2]
    qry_chr, qry_start = qry_dict[qry_start_gene][0:2]
    qry_end = qry_dict[qry_end_gene][2]
    
    circos_input = [ref_chr, ref_start, ref_end, qry_chr, qry_start, qry_end]
    fo.writelines('\t'.join(circos_input) + '\n')
    
fo.close()
