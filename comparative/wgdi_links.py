import pandas as pd 
import sys

def get_gff_pos(gff, chrom, index):
    target = gff.loc[ (gff['index'] == int(index)) & (gff['chrom'] == chrom), :]
    #print(chrom, index)
    #print(target)
    pos = int(target.start)
    return pos

if len(sys.argv) < 4:
    print("Usage: prog gff_file block_file")
    sys.exit(1)

gff_file   = sys.argv[1]
block_file = sys.argv[2]
link_file  = sys.argv[3]

# load gff
gff = pd.read_csv(gff_file, delimiter = '\t', names = ['chrom','gene','start','end', 'strand','index','transcript'])
# load block info
block_info = pd.read_csv(block_file)

# 获取起始和结束位点
links = []
for rows in block_info.iterrows():

    # 获取起始基因和结束基因的位置
    chr1 = rows[1]['chr1']
    start1 = get_gff_pos(gff, chr1, rows[1]['start1'])
    end1 = get_gff_pos(gff, chr1, rows[1]['end1'] )

    chr2 = rows[1]['chr2']
    start2 = get_gff_pos(gff, chr2, rows[1]['start2'])
    end2 = get_gff_pos(gff,chr2, rows[1]['end2'])

    links.append([ chr1, start1, end1, chr2, start2, end2])


links_df = pd.DataFrame(links)
links_df.to_csv(link_file, header = None, index=None)
