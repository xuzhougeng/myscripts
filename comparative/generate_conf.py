#!/usr/bin/env python3

# GFF must have CDS feature
# GFF must have ID and Parent in column 9

import re
import argparse
from collections import defaultdict
from collections import OrderedDict

def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta',
                        help="fasta file name")
    parser.add_argument('gff',
                        help="gff file name")
    parser.add_argument('-p','--prefix', type=str, default="output",
                        help="prefix for ouput ")

    return parser



# get the fasta  len
def get_fasta_len(fasta):
    fasta_dict = OrderedDict()
    handle = open(fasta, "r")
    active_sequence_name = ""
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"): 
            active_sequence_name = line[1:]
            active_sequence_name = active_sequence_name.split(" ")[0]
        if active_sequence_name not in fasta_dict:
            fasta_dict[active_sequence_name] = 0
            continue
        sequence = line
        fasta_dict[active_sequence_name] += len(sequence)
    handle.close()
    return fasta_dict

# parse the gff 
def parse_gff(gff):

    gene_dict = OrderedDict()
    tx_pos_dict = defaultdict(list)
    CDS_dict = defaultdict(list)

    handle = open(gff, "r")

    for line in handle:
        if line.startswith("#"):
            continue
        content = line.split("\t")
        if len(content) <= 8:
            continue
        #print(content)
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
    handle.close()
    return [gene_dict, tx_pos_dict, CDS_dict]

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    fa_dict = get_fasta_len( args.fasta)
    gene_dict, tx_pos_dict, CDS_dict= parse_gff( args.gff )
    gene_count = {}

    # outfile
    len_file = open(args.prefix + ".len", "w")
    gff_file = open(args.prefix + ".gff", "w")

    for gene, txs in gene_dict.items():
        tmp = 0
        for tx in txs:
            tx_len = sum(CDS_dict[tx])
            if tx_len > tmp:
                lst_tx = tx
                tmp = tx_len
        tx_chrom = tx_pos_dict[lst_tx][0]
        if tx_chrom not in gene_count:
            gene_count[tx_chrom] = 1
        else:
            gene_count[tx_chrom] += 1
        tx_start = tx_pos_dict[lst_tx][1]
        tx_end   = tx_pos_dict[lst_tx][2]
        tx_strand = tx_pos_dict[lst_tx][3]

        print("{chrom}\t{gene}\t{start}\t{end}\t{strand}\t{order}\t{tx}".format(\
                chrom=tx_chrom,gene=gene,start=tx_start,end=tx_end,strand=tx_strand,order=gene_count[tx_chrom],tx=lst_tx), file=gff_file )

    for chrom,lens in fa_dict.items():
        print("{chrom}\t{lens}\t{count}".format(\
                chrom=chrom,lens=lens,count=gene_count.get(chrom,0)), file=len_file)
    len_file.close()
    gff_file.close()
