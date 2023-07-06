#!/usr/bin/env python3

# convet GFF formt to BED12
# BED12 ref to: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

# requirment:
## GFF must have CDS feature
## GFF must have ID and Parent in column 9

import re
from collections import defaultdict

def parse_gff(gff_file):

    gene_dict = defaultdict(list) # 记录基因对应的tx_id
    tx_pos_dict = defaultdict(list) #记录转录本中的chr,start,end,strand
    CDS_dict = defaultdict(list)    #记录CDS的width, start

    for line in open(gff_file, 'r'):
        if line.startswith("#"):
            continue
        content = line.split("\t")
        if len(content) <= 8:
            continue

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
            tx_pos_dict[tx_id] = [content[0],int(content[3])-1, int(content[4]), content[6] ]
        
        # GFF must have CDS feature
        if content[2] == 'CDS':
            width = int(content[4]) - int(content[3]) + 1 
            start = int(content[3]) - 1
            CDS_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
            CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
            CDS_dict[CDS_parent].append([width, start])

    return gene_dict,tx_pos_dict,CDS_dict

def get_longest_tx(txs, CDS_dict):
    
    # get the longest transcript id
    tmp = 0

    for tx in txs:
        tx_len = sum( [x[0] for x in  CDS_dict[tx] ])
        if tx_len > tmp:
            lst_tx = tx
            tmp = tx_len
        else:
            lst_tx = tx
    
    return lst_tx

def gff2bed12(gene_dict, tx_pos_dict, CDS_dict):
    bed_list = []
    for gene, txs in gene_dict.items():
        
        lst_tx = get_longest_tx(txs, CDS_dict)    


        tx_chrom = tx_pos_dict[lst_tx][0]
        tx_start = tx_pos_dict[lst_tx][1]
        tx_end   = tx_pos_dict[lst_tx][2]
        tx_strand = tx_pos_dict[lst_tx][3]

        #print(CDS_dict)
        CDS_list = CDS_dict[lst_tx]

        blockCount = 0
        block_size_list = []
        block_start_list = []

        #print(CDS_list)

        for CDS in CDS_list:
            blockCount += 1
            block_size_list.append(str(CDS[0]))
            # start  relative to chromStart
            block_start_list.append(str(CDS[1] - int(tx_start)))
        
        if tx_strand == "-":
            block_size_list = block_size_list[::-1]
            block_start_list = block_start_list[::-1]

        # bed12 format
        # chrom, start,end,name,score,strand,start,end,0,blockCount,blockSize,blockStart
        blockSize = ','.join(block_size_list)
        blockStart = ','.join(block_start_list)

        #print(blockSize, blockStart)

        bed = f'{tx_chrom}\t{tx_start}\t{tx_end}\t{lst_tx}\t0\t{tx_strand}\t{tx_start}\t{tx_end}\t0\t{blockCount}\t{blockSize}\t{blockStart}'
        
        bed_list.append(bed)

    return bed_list


if __name__ == "__main__":
    import sys
    gff_file = sys.argv[1]
    bed_file = sys.argv[2]

    gene_dict,tx_pos_dict,CDS_dict = parse_gff(gff_file)
    bed_list = gff2bed12(gene_dict,tx_pos_dict,CDS_dict)

    with open(bed_file, 'w') as handler:
        for bed in bed_list:
            handler.write(bed + '\n')
    
