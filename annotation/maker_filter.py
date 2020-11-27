#!/usr/bin/env python3

"""
filter the maker output gff by:
    - AED
    - eAED
    - QI
"""

import re
import sys
import argparse


def get_opt():
    group = argparse.ArgumentParser()
    group.add_argument("-c", "--ss", type = float, default = -1, 
            help = "The faction of splice sites confirmed by an EST alignment, default -1", required = False)
    group.add_argument("-e", "--exon", type = float, default = -1, 
            help = "The faction of exons that overlap an EST alignment, default -1", required = False)
    group.add_argument("-o", "--exon2", type = float, default = -1, 
            help = "The faction of exons that overlap an EST/Protein alignment, default -1", required = False)
    group.add_argument("-a", "--ss2", type = float, default = -1, 
            help = "The faction of splice sites confirmed by an ab-initio prediction, default -1", required = False)
    group.add_argument("-t", "--exon3", type = float, default = -1, 
            help = "The faction of exon confirmed by an ab-initio prediction, default -1", required = False)
    group.add_argument("-l", "--length", type = int, default = 0, 
            help = "The min length of the protein sequence produced by the mRNA", required = False)
    group.add_argument("-d", "--AED", type = float, default = 1, required = False,
            help = "Max AED  to allow, default is 1")
    group.add_argument("-i", "--geneid", required = False,
            help = "filter by given gene list")
    group.add_argument("gff", help = "input gff file")
    return group.parse_args()

def is_high_confidence(AED, QI, thrAED, thresh):
    flag = True
    if AED > float(thrAED):
        flag = False
    QI = QI.split('|') 
    QI = QI[1:6] + QI[8:]
    for i in range(0,len(thresh)):
        if float(QI[i]) < float(thresh[i]):
            flag = False
            break
    return flag

def parse_anno(col):
    """
    parse the annotation column, and return a dict
    """
    anno = re.split('[;=]', col)
    if anno[-1] == '':
        anno.pop()

    anno_dict = {}
    for i in range(0, len(anno), 2):
        anno_dict[anno[i]] = anno[i+1]
    return anno_dict


def parse_gff(gff, thrAED, thresh, geneid = None):
    """
    parse the gff file 
    """
    good_gene = set()
    good_mRNA = set()
    for line in open(gff, "r"):
        if line.startswith("#"):
            continue

        line = line.strip()
        cols = line.split('\t')
        if len(cols) != 9:
            continue
        if cols[2] != "mRNA":
            continue

        anno_dict = parse_anno(cols[8])
        AED = float(anno_dict['_AED'])
        QI  = anno_dict['_QI']

        if geneid is not None and anno_dict['Parent'] not in geneid:
            continue

        if is_high_confidence(AED, QI, thrAED, thresh): 
            good_mRNA.add(anno_dict['ID'])
            good_gene.add(anno_dict['Parent'])
        #print(good_gene)

    return (good_gene, good_mRNA)

def filter_gff(gff, good_gene, good_mRNA):
    """
    filter the gff by the gene list and mRNA list
    """

    for line in open(gff, "r"):
        if line.startswith("#"):
            continue
        line = line.strip()
        cols = line.split('\t')
        if len(cols) != 9:
            continue

        anno_dict = parse_anno(cols[8])

        if cols[2] == 'gene':
            if anno_dict['ID'] in good_gene:
                print(line)
            else:
                continue
        elif cols[2] == 'mRNA':
            if anno_dict['ID'] in good_mRNA:
                print(line)
            else:
                continue
        else:
            if anno_dict['Parent'] in good_mRNA:
                print(line)
            else:
                continue

if __name__ == "__main__":
    opts = get_opt()
    gff = opts.gff
    thrAED = opts.AED
    thresh = [opts.ss, opts.exon, opts.exon2, opts.ss2, opts.exon3, opts.length]
    if opts.geneid:
        geneid = set( gene.strip() for gene in open(opts.geneid) )
        good_gene,good_mRNA = parse_gff(gff, thrAED, thresh, geneid)
    else:
        good_gene,good_mRNA = parse_gff(gff, thrAED, thresh)

    filter_gff(gff, good_gene, good_mRNA)

