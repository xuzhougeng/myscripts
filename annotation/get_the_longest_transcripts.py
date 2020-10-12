#!/usr/bin/env python3

# v0.1:fixing the problem that the feature of exon/mRNA is before gene
# v0.2: Use CDS to decide which transcript is primary
# v0.3: Fix the bug that if the CDS_parent include '\r' , the lst_len will be 0

# GFF must have CDS feature
# GFF must have ID and Parent in column 9

import re
import fileinput
from collections import defaultdict

gene_dict = defaultdict(list)
tx_pos_dict = defaultdict(list)
CDS_dict = defaultdict(list)

for line in fileinput.input():
    if line.startswith("#"):
        continue
    content = line.split("\t")
    if len(content) <= 8:
        continue
    #print(content)
    '''
    if content[2] == 'gene':
        gene_id = re.search(r'ID=(.*?)[;\n]', content[8]).group(1)
        #print(gene_id)
        if gene_id not in gene_dict:
            gene_dict[gene_id] = []
    '''
    if content[2] == "transcript" or content[2] == "mRNA":

        # use regual expression to extract the gene ID
        # match the pattern ID=xxxxx; or ID=xxxxx

        tx_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
        tx_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
        tx_parent = tx_parent.strip() # remove the 'r' and '\n'

        # if the parent of transcript is not in the gene_dict, create it rather than append
        if tx_parent in gene_dict:
            gene_dict[tx_parent].append(tx_id)
        else:
            gene_dict[tx_parent] = [tx_id]
        tx_pos_dict[tx_id] = [content[0],content[3], content[4], content[6] ]
    # GFF must have CDS feature
    if content[2] == 'CDS':
        width = int(content[4]) - int(content[3])
        CDS_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
        CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
        CDS_dict[CDS_parent].append(width)

for gene, txs in gene_dict.items():
    tmp = 0
    for tx in txs:
        tx_len = sum(CDS_dict[tx])
        if tx_len > tmp:
            lst_tx = tx
            tmp = tx_len
    tx_chrom = tx_pos_dict[lst_tx][0]
    tx_start = tx_pos_dict[lst_tx][1]
    tx_end   = tx_pos_dict[lst_tx][2]
    tx_strand = tx_pos_dict[lst_tx][3]
    print("{gene}\t{tx}\t{chrom}\t{start}\t{end}\t{strand}".format(gene=gene,tx=lst_tx,chrom=tx_chrom,start=tx_start,end=tx_end, strand=tx_strand))
