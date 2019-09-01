#!/user/bin/env python
import re
import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("usage: add_annotation_from_dat.py blastp dat")
    sys.exit(1)

blastp_file = sys.argv[1]
dat = sys.argv[2]

uniprot = SeqIO.index(dat, "swiss")

out_file = open("swiss_annotation.tsv", "w")

for line in open(blastp_file, "r"):
    gene,acc,ident = line.strip().split()[0:3]
    if not uniprot.get(acc.strip(";")):
        continue
    record = uniprot.get_raw(acc.strip(";"))
    GO_RECORD = re.findall(r"GO; (GO:\d+); ([F|P|C]):.*?; (.*):",record)

    SPECIES = re.findall(r"OS   (.*)\.", record)
    if len(SPECIES) == 0:
        SPECIES = [""]
    else:
        SPECIES = [ SPECIES[0].replace("\t"," ") ]

    ENSEMBLE_Plant = re.findall(r"EnsemblPlants; (.*?);", record)
    if len(ENSEMBLE_Plant) == 0:
        ENSEMBLE_Plant = [""]
    else:
        ENSEMBLE_Plant = [ENSEMBLE_Plant[0].replace("\t"," ")]

    for GO in GO_RECORD:
        outline = "\t".join([gene, acc,ident] + SPECIES + ENSEMBLE_Plant + list(GO))
        out_file.writelines(outline + "\n")

out_file.close()
