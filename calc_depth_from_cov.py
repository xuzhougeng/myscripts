#!/usr/bin/env python3

# calculate the mean depth in PB.base.cov from purge_dups

from sys import argv

if len(argv) <= 2:
    print("")
    print("Program: Calculate the mean depth of each contig in PB.base.cov from purge_dups")
    print("")
    print("Usage:   {} PB.base.cov output".format(argv[0]))
    exit(1)


fi = argv[1]
fo = argv[2]
depth = 0

contig_dict = {}

ifp = open(fi, "r")
ofp = open(fo, "w")

# read the first line
line = ifp.readline()
name, size = line.strip().split('\t')

for line in ifp:
    line = line.strip()
    if line.startswith(">"):
        contig_dict[name[1:]] = [int(size), depth / int(size)]
        #print(size + str(depth))
        name, size = line.split('\t')
        depth = 0
    else:
        start,end,cov= line.split('\t')
        depth = depth + ( int(end) - int(start) + 1 ) * int(cov)
        

# final 
contig_dict[name[1:]] = [int(size), depth / int(size)]

for key, value in contig_dict.items():
    ofp.writelines(key + '\t' + str(value[0]) + '\t' + str(value[1]) + '\n')

ofp.close()
ifp.close()
