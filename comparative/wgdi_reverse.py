
# 作用
# 对给定的染色体进行顺序调整
# 输入文件: xxx.len, xxx.gff

import sys
import pandas as pd
import numpy as np

def read_gff(gff_file):
    names = ['chrom','gene','start','end', 'strand','index','transcript']
    gff = pd.read_csv(gff_file, delimiter = '\t', names = names)

    return gff

def read_len(len_file):
    len_dict = {}
    with open(len_file, 'r') as handler:
        for line in handler:
            chrom,length,gene_num = line.strip().split('\t')
            len_dict[chrom] = [length, gene_num]

    return len_dict


def split_df(gff_df):
    
    gff_dict = {}
    chroms = np.unique(gff_df['chrom'])
    for chrom in chroms:
        gff_dict[chrom] = gff_df.loc[gff_df['chrom'] == chrom,:]

    return gff_dict

def reverse_chrom(gff_dict, len_dict, chrom):
    gff_df = gff_dict[chrom]
    length,gene_num = len_dict[chrom]
    gff_df = gff_df.iloc[::-1]

    # reverse the strand
    gff_df['strand'] = ["+" if x == "-" else "-" for x in gff_df['strand'] ]
    # revese the order
    new_end = [ int(length) - x + 1 for x in  gff_df['start'] ]
    new_start  = [ int(length) - x + 1 for x in  gff_df['end'] ]
    
    gff_df['start'] = new_start
    gff_df['end'] = new_end

    # reserve the index
    gff_df['index'] = [ int(gene_num) - x + 1 for x in  gff_df['index'] ]
    
    return gff_df
    

def main(args):

    gff_file = args[1]
    len_file = args[2]
    chroms = args[3].split(",")
    out_gff = args[4]

    gff_df = read_gff(gff_file)
    all_chroms = np.unique(gff_df['chrom'])
    len_dict = read_len(len_file)

    gff_dict = split_df(gff_df)

    new_gff_dict = {}
    for chrom in all_chroms:
        if chrom in chroms:
            new_gff_dict[chrom] = reverse_chrom(gff_dict, len_dict, chrom)
        else:
            # keep order
            new_gff_dict[chrom] = gff_df.loc[gff_df['chrom'] == chrom,:]
    
    new_gff = pd.concat(new_gff_dict)

    new_gff.to_csv(out_gff, sep="\t", header=False, index=False)


if __name__ == "__main__":
    args = sys.argv
    if len(args) != 5:
        print("Usage: wgdi_reverse.py input_gff len chroms out_gff")
        sys.exit(-1)
    main(args)