#!/usr/bin/env python3

"""
一个非常简单的脚本，主要目的是将maker_map_ids输出的编号，改成类似于AT1G10010这种形式
输入参数为两个, maker_map_ids的输出文件, 还有你的物种名的缩写, 例如Athaliana-> AT, Alyrata->AL
scaffold，则是用U表示
"""

import re
import sys

if len(sys.argv) < 3:
    print("Usage: {} map.id prefix".format(sys.argv[0]))
    sys.exit()

gff = sys.argv[1]
prefix = sys.argv[2]

chr_name = ""
gene_count = 1
count = 0

with open(gff, "r") as file:
    lines = file.readlines()

for line in lines:
    raw_name = line.split("\t")[0]
    if raw_name.find("scaffold") > 0:
        break
    count += 1

    # if it has mRNA, keep the count unchange, and mRNA_count add 1
    if raw_name.find("mRNA") >=0:
        mRNA_count = raw_name.split("-")[6]
        print("{}\t{}{}G".format(raw_name, prefix, chr_id) + str(gene_count).zfill(4) + "0." + str(mRNA_count ) )
    else:
        # check whether chrosome name is same as before
        if chr_name == raw_name.split("-")[1]:
            gene_count += 1
            print("{}\t{}{}G".format(raw_name, prefix, chr_id) + str(gene_count).zfill(4) + "0" )
        else:
            chr_name = raw_name.split("-")[1]
            chr_id = chr_name.lstrip("chr")
            gene_count = 1
            print("{}\t{}{}G".format(raw_name, prefix, chr_id) + str(gene_count).zfill(4) + "0" )

gene_count = 0
for line in lines[count:]:
    raw_name = line.split("\t")[0]
    chr_id = "U"
    if raw_name.find("mRNA") >=0:
        mRNA_count = raw_name.split("-")[6]
        print("{}\t{}{}G".format(raw_name, prefix, chr_id) + str(gene_count).zfill(4) + "0." + str(mRNA_count ) )
    else:
        gene_count += 1
        print("{}\t{}{}G".format(raw_name, prefix, chr_id) + str(gene_count).zfill(4) + "0" )
