#!/usr/bin/env python3

# fix bug: use position and ID together as the key of CDS
# different transcripts may shared the same CDS

# get the complete CDS in TransDecoder results


import re
import fileinput
from sys import stderr
from sys import argv
from collections import defaultdict

# filter info count
by_homology=0

# set the dictionary for storing record
gene_dict = {}
mRNA_dict = {}
other_dict = {}

CDS_dict = defaultdict(list) # store the mRNA CDS width

# Parent-Child Dictionary
gene_tx_dict = defaultdict(list)
mRNA_other_dict = defaultdict(list)

#lines = [ line for line in fileinput.input()]

# get the Parent ID in the attribute
def get_parent_id(attribute):
    parent = re.search(r'Parent=(.*?)[;\n]', attribute).group(1)
    parent = parent.strip() # remove the 'r' and '\n'
    return parent

# get the ID in the attribute
def get_id(x):
    feature_id = re.search(r'ID=(.*?)[;\n]', attribute).group(1)
    return feature_id

# parse the gff
for line in fileinput.input():

    # filter un-related line
    if line.startswith("#"):
        continue
    content = line.split("\t")
    if len(content) <= 8:
        continue

    attribute = content[8]
    feature_id = get_id(attribute)

    if content[2] == "gene":
        gene_dict[feature_id] = line

    elif content[2] == "mRNA":

        # only keep complete and homology-supported gene
        if not re.search('complete', line):
            continue
        # extrac the homology score
        score = re.findall("score(.*)", line)[0].split("%2C")
        if len(score) != 2:
            by_homology += 1
            continue

        # store the mRNA and its parents
        mRNA_dict[feature_id] = line
        parent = get_parent_id(content[8])

        # if the parent of transcript is not in the gene_dict, create it rather than append
        if parent in gene_tx_dict:
            gene_tx_dict[parent].append(feature_id)
        else:
            gene_tx_dict[parent] = [feature_id]

    else:
        parent = get_parent_id(attribute)
        if content[2] == 'CDS':
            width = int(content[4]) - int(content[3])
            CDS_dict[parent].append(width)
            feature_id = '{}-{}-{}-{}-cds'.format(feature_id,content[0],content[3],content[4])
            other_dict[feature_id] = line
            if parent in mRNA_other_dict:
                mRNA_other_dict[parent].append(feature_id)
            else:
                mRNA_other_dict[parent] = [feature_id]
        else:
            other_dict[feature_id] = line
            if parent in mRNA_other_dict:
                mRNA_other_dict[parent].append(feature_id)
            else:
                mRNA_other_dict[parent] = [feature_id]

stderr.write("Number of entry filtered by homology: {}\n".format(by_homology))

# filter the uncomplete mRNA
#print(mRNA_other_dict.items())
for gene, txs in gene_tx_dict.items():
    new_txs = []
    for tx in txs:
        feature = ';'.join(mRNA_other_dict[tx])
        if re.search('cds',feature) and re.search('exon',feature) and \
           re.search('utr5p',feature ) and re.search('utr3p',feature):
            new_txs.append(tx)
    gene_tx_dict[gene] = new_txs

#print(gene_tx_dict.items())
# get the high quality geneID and mRNA ID 
for gene, txs in gene_tx_dict.items():
    tmp = 0

    if len(txs) == 0:
        continue

    for tx in txs:
        tx_len = sum(CDS_dict[tx])
        if tx_len > tmp:
            lst_tx = tx
            tmp = tx_len

    print(gene_dict[gene].strip())
    print(mRNA_dict[lst_tx].strip())     
    #print(mRNA_other_dict[lst_tx]) # for debug
    for value in set(mRNA_other_dict[lst_tx]):
        print(other_dict[value].strip())
