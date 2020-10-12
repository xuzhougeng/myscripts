from sys import argv

f = open(argv[1])
prefix = argv[2]
non_chr = 'scaffold'

prev_chrom = 'flag'
count = 1
sub_count = 1
non_chr_count = 0

for line in f.readlines():
    gene_name = line.strip().split()[0]
    chrom = gene_name.split("-")[1]
    if chrom == prev_chrom:
        if 'mRNA' in gene_name:
            if non_chr in chrom:
                print(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0.{}'.format(sub_count))
            else:
                chrom_id = chrom[3] 
                print(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0.{}'.format(sub_count))
            count += 1
            sub_count += 1
        else:
            sub_count = 1
            if non_chr in chrom:
                print(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0')
            else:
                chrom_id = chrom[3] 
                print(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0')

        prev_chrom = chrom
    else :
        count = 1 
        if non_chr in chrom:
            non_chr_count += 1
            print(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0')
        else:
            chrom_id = chrom[3] 
            print(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0')
        prev_chrom = chrom

f.close()
