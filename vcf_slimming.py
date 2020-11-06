#!/usr/env/bin python3

# A script to filter GBS/RAD data 

import pysam
import argparse
from collections import deque


def get_opt():
    group = argparse.ArgumentParser()
    group.add_argument('vcf', help="VCF file or BCF")
    group.add_argument('-o', '--out', type=str, default='-',help="output file name")
    group.add_argument('-l', '--lib_size', type=int, default=500,help="library size, default is 500")
    group.add_argument('-p1', '--parent_1', type=str, help="parent 1")
    group.add_argument('-p2', '--parent_2', type=str, help="parent 2")
    return group.parse_args()


def select_snp( snp_buffer ):

    # the first record to comparsion
    rec = snp_buffer.pop()

    # iterate the snp to select the rec with highest DP
    while snp_buffer :
        tmp = snp_buffer.pop()
        if tmp.info['DP'] > rec.info['DP']:
            rec = tmp

    return rec


def filter_snp(vcf, lib_size = 500, out_fn = "-"):
    vcf_in  = pysam.VariantFile(vcf)
    vcf_out = pysam.VariantFile(out_fn, 'w', header = vcf_in.header)

    # read the first record
    rec = next(vcf_in)
    snp_chrom = rec.chrom
    snp_start = rec.pos

    # set the buffer queue
    # only a small part of match sequence will be different
    snp_buffer = deque( maxlen = int(lib_size / 10 ))

    for rec in vcf_in.fetch():
        #print(f'rec.pos: {rec.pos}; rec.chrom: {rec.chrom}')
        #print(f'snp_start: {snp_start}; snp_chrom: {snp_chrom}')
        if rec.pos < (snp_start + lib_size) and rec.chrom == snp_chrom:
            #print(f'append {rec.pos} : {rec.chrom}')
            snp_buffer.append(rec)
            continue
        
        snp_start = rec.pos
        snp_chrom = rec.chrom
        if snp_buffer:
            rec = select_snp(snp_buffer)
            vcf_out.write(rec)

    # process the last record
    ## if the buffer queue is empty, it means
    ## - the last record is not at the same chromosome/scaffold/contig as the last but one
    ## - or the last record position is larger than the last but one
    if snp_buffer:
        rec = select_snp(snp_buffer)

    vcf_out.write(rec)
                
if __name__ == "__main__":
    opts = get_opt()
    vcf = opts.vcf
    out_fn = opts.out
    lib_size = opts.lib_size
    filter_snp(vcf, lib_size, out_fn)

