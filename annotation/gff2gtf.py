#!/usr/bin/env python3

"""
convert GFF3 to GTF in ENSEMBLE fomrats

It will only use the following features

- gene
- mRNA -> transcript
- CDS
- exon
- five_prime_UTR -> five_prime_utr
- three_prime_UTR -> three_prime_utr

TO DO:add additional two features

- start_codon
- stop_codon

I only test this script in Arabidopsis_thaliana.TAIR10.44.gff3.gz
"""

import re
import sys

def ensemble_style(input_file):

    for line in open(input_file):
        line = line.strip()
        if line.startswith("###"):
            continue
        if line.startswith("#"):
            print(line)
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "gene":
            gene_id = re.findall(r'gene_id=(.*?);',fields[8])[0]
            fields[8] = "gene_id \"{}\";gene_name \"{}\";".format(gene_id, gene_id)
        elif feature == "mRNA" or feature == "transcript":
            tx_id = re.findall(r'transcript_id=(.*?)$',fields[8])[0]
            gene_id = tx_id.split(".")[0]
            fields[2] = "transcript"
            fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
        elif feature == "CDS":
            tx_id = re.findall(r'protein_id=(.*?)$',fields[8])[0]
            gene_id = tx_id.split(".")[0]
            fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
        elif feature == "exon":
            exon_id = re.findall(r'exon_id=(.*?);',fields[8])[0]
            tx_id = ".".join(exon_id.split(".")[0:2])
            gene_id = tx_id.split(".")[0]
            fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
        elif feature == "five_prime_UTR" or feature == "five_prime_utr" or feature == "5UTR":
            fields[2] = "five_prime_utr"
            tx_id = fields[8].split(":")[1]
            gene_id = tx_id.split(".")[0]
            fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
        elif feature == "three_prime_UTR" or feature == "three_prime_utr" or feature == "3UTR":
            fields[2] = "three_prime_utr"
            tx_id = fields[8].split(":")[1]
            gene_id = tx_id.split(".")[0]
            fields[8] = "gene_id \"{}\";transcript_id \"{}\";gene_name \"{}\";".format(gene_id, tx_id, gene_id)
        else:
            continue
        print("\t".join(fields))

def re_find(name, string):
    pattern = r'{}=(.*?)(;|$)'.format(name)
    matches = re.findall(pattern, string)
    if matches:
        return matches
    else:
        return None

def star_style(input_file, output_file=None):

    gene_attr_dict = {} # gene : gene_name
    gene_dict = {} # gene [ tx1, tx2, ...]
    tx_dict = {}   # tx : [exon1, exon2, ...]
    exon_dict = {} # exon: [chr, start, end, strand]

    for line in open(input_file):
        line = line.strip()
        if line.startswith("##"):
            continue
        
        fields = line.split("\t")
        feature = fields[2]
        if feature == "gene":
            # exctrat gene_id from ID
            matches  = re_find("ID", fields[8])
            if matches is None:
                raise ValueError("No ID found in gene feature.")
            gene_id = matches[0][0]

            # extract gene_name from Name
            gene_name = re_find("Name", fields[8])[0][0]
            gene_attr_dict[gene_id] = gene_name

        elif feature == "mRNA" or feature == "transcript":
            # extract transcript_id from ID
            matches = re_find("ID", fields[8])
            if matches is None:
                raise ValueError("No ID found in mRNA feature.")
            tx_id = matches[0][0]
            # extract gene_id from Parent
            matches = re_find("Parent", fields[8])
            if matches is None:
                raise ValueError("No Parent found in mRNA feature.")
            gene_id = matches[0][0]
            
            if gene_id not in gene_dict:
                gene_dict[gene_id] = [tx_id]
            else:
                gene_dict[gene_id].append(tx_id)
        elif feature == "exon":
            # extract transcript_id from ID
            matches = re_find("ID", fields[8])
            if matches is None:
                raise ValueError("No ID found in exon feature.")
            exon_id = matches[0][0]
            # extract gene_id from Parent
            matches = re_find("Parent", fields[8])
            if matches is None:
                raise ValueError("No Parent found in exon feature.")
            tx_id = matches[0][0]

            if tx_id not in tx_dict:
                tx_dict[tx_id] = [exon_id]
            else:
                tx_dict[tx_id].append(exon_id)
            exon_dict[exon_id] = [fields[0], fields[3], fields[4], fields[6]]

    if output_file is None:
        output = sys.stdout
    else:
        output = open(output_file, "w")

    for gene,tx_ids in gene_dict.items():
        gene_name = gene_attr_dict[gene]

        for tx_id in tx_ids:
            exon_ids = tx_dict[tx_id]
            for exon_id in exon_ids:
                chrom, start, end, strand = exon_dict[exon_id]
                output.write("{}\tSTAR\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\";\n".format(chrom, start, end, strand, gene, tx_id, gene_name))


    if output_file is not None:
        output.close()


def main(args):
    input_file = args.input
    output_file = args.output

    if args.format == "ensemble":
        ensemble_style(input_file)
    elif args.format.lower() == "star":
        star_style(input_file, output_file)
    else:
        print("Please provide the correct format.")
        exit(1)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="This script convert the gff to gtf.")

    parser.add_argument("-i", "--input", type=str, required=True, help="input gff file")
    parser.add_argument("-o", "--output", type=str, required=False, help="output gtf file")
    parser.add_argument("-f", "--format", type=str, required=True, choices=["ensemble", "star"], help="out gff format")
    args = parser.parse_args()

    main(args)

    
