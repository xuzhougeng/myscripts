import pandas as pd
import sys

# 作用
# 根据WGDI提供共线性信息，整理注释的注释
# 输入是: 两个物种的gff, 以及bi的输出结果

# gff1: query
# gff2: target
# strand: +: same direction, otherwise "-"

def read_gff(gff_file):
    names = ['chrom','gene','start','end', 'strand','index','transcript']
    gff = pd.read_csv(gff_file, delimiter = '\t', names = names) 
    gff['chrom'] = [str(x) for x in gff['chrom']]

    return gff


def get_gff_pos(gff, chrom, index):

    target = gff.loc[ (gff['index'] == int(index)) & (gff['chrom'] == str(chrom)), :]

    gene = str(target.gene.to_list()[0])
    #print(gene)
    return gene

def read_block_info(block_file):
   block_df =  pd.read_csv(block_file)
   return block_df


def extract_pair_gene(block_df, gff1, gff2):
    pair_index = {}

    gene_block_len = {} #记录基因共线性所在的block的长度

    for rows in block_df.iterrows():
        chr1 = rows[1]['chr1']
        chr2 = rows[1]['chr2']
        x = rows[1]['block1']
        y = rows[1]['block2']
        block_len = int(rows[1]['length'])
        
        for i,j in zip(x.split("_"), y.split("_")):

            gene_i = get_gff_pos(gff1, chr1, i )
            gene_j = get_gff_pos(gff2, chr2, j )

            if gene_i in gene_block_len:
                if block_len > gene_block_len[gene_i]:
                    gene_block_len[gene_i] = block_len
                    pair_index[gene_i] = [gene_i, gene_j] 
                    
            else:
                gene_block_len[gene_i] = block_len
                pair_index[gene_i] = [gene_i, gene_j] 

    return pair_index

def write_pair(gene_pair, out_file):
    with open(out_file, 'w') as handler:
        for x in gene_pair.values():
            #x = [str(x) for x in s ]
            handler.writelines('\t'.join(x) + '\n')


def main(args):
    block_file = args[1]
    gff1 = args[2]
    gff2 = args[3]
    out_file = args[4]

    block_df = read_block_info(block_file=block_file)
    #block_df = filter_block(block_df)
    gff1 = read_gff(gff1)
    gff2 = read_gff(gff2)
    #print(block_df)
    gene_pair = extract_pair_gene(block_df, gff1, gff2)

    write_pair(gene_pair, out_file)    


if __name__ == "__main__":
    main(sys.argv)
