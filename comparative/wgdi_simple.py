import pandas as pd
import sys

# 作用
# 将WGDI的block information表(csv)导出为JCVI的simple格式

# gff1: query
# gff2: target
# strand: +: same direction, otherwise "-"

def read_gff(gff_file):
    names = ['chrom','gene','start','end', 'strand','index','transcript']
    gff = pd.read_csv(gff_file, delimiter = '\t', names = names)
    gff['chrom'] = [str(x) for x in gff['chrom']]
    return gff


def get_gff_pos(gff, chrom, index):
    target = gff.loc[ (gff['index'] == int(index)) & (gff['chrom'] == chrom), :]


    gene = str(target.gene.to_list()[0])
    #print(gene)
    return gene

def read_block_info(block_file):
   block_df =  pd.read_csv(block_file)
   block_df['chr1'] = [str(x) for x in block_df['chr1'] ]
   block_df['chr2'] = [str(x) for x in block_df['chr2'] ]
   return block_df


def filter_block(block_df, ):
    pass


def block2simple(block_df, gff1, gff2):

    simple_list = []
    
    for rows in block_df.iterrows():

        # 获取起始基因和结束基因的位置
        chr1 = rows[1]['chr1']
        start1 = get_gff_pos(gff1, chr1, rows[1]['start1'])
        end1 = get_gff_pos(gff1, chr1, rows[1]['end1'] )

        chr2 = rows[1]['chr2']
        start2 = get_gff_pos(gff2, chr2, rows[1]['start2'])
        end2 = get_gff_pos(gff2,chr2, rows[1]['end2'])

        block_length = rows[1]['length']

        direction =  ( rows[1]['end1'] - rows[1]['start1'] ) * (rows[1]['end2'] - rows[1]['start2'])
        strand = "+" if direction >=0 else "-"

        if strand == "-":
            start1,end1 = end1,start1

        #print([ start1, end1, start2, end2, block_length, strand])

        simple_list.append([ start1, end1, start2, end2, block_length, strand])

    return simple_list

def write_simple(simple_list, simple_file):
    with open(simple_file, 'w') as handler:
        for s in simple_list:
            x = [str(x) for x in s ]
            handler.writelines('\t'.join(x) + '\n')



def main(args):
    block_file = args[1]
    gff1 = args[2]
    gff2 = args[3]
    simple_file = args[4]

    block_df = read_block_info(block_file=block_file)
    #block_df = filter_block(block_df)
    gff1 = read_gff(gff1)
    gff2 = read_gff(gff2)
    print(block_df)
    simple_df = block2simple(block_df, gff1, gff2)

    write_simple(simple_df, simple_file)    


if __name__ == "__main__":
    main(sys.argv)
