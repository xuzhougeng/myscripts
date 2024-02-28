from typing import List, Dict, Tuple
import argparse
import subprocess

# TODO: CDS prediction

def build_gmap_index(genome_file, index_name, gmap_build="gmap_build", dir_name = ".") :
    """
    Build gmap index
    """
    cmd = f'{gmap_build} -D {dir_name} -d {index_name} {genome_file}'

    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(e)
        print("Build gmap index error")

def gmap_alignment( index_name, query_file, gmap="gmap", dir_name = ".", psl_out = None):
    """
    Gmap alignment
    """
    thread = 8
    if not psl_out:
        psl_out = f'{query_file}.psl'
    
    cmd = f'{gmap} -D {dir_name} -d {index_name} -f psl -n 0 -t {thread} {query_file} > {psl_out}'
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(e)
        print("Gmap alignment error")

def read_psl(psl_file) -> List[Dict]:
    """
    Read psl file
    """
    psl_list = []
    with open(psl_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("match"):
                continue
            line_list = line.split("\t")
            psl_dict = {}
            psl_dict["match"] = int(line_list[0])
            psl_dict["mismatch"] = int(line_list[1])
            psl_dict["repmatch"] = int(line_list[2])
            psl_dict["Ns"] = int(line_list[3])
            psl_dict["Qgapcount"] = int(line_list[4])
            psl_dict["Qgapbases"] = int(line_list[5])
            psl_dict["Tgapcount"] = int(line_list[6])
            psl_dict["Tgapbases"] = int(line_list[7])
            psl_dict["strand"] = line_list[8]
            psl_dict["Qname"] = line_list[9]
            psl_dict["Qsize"] = int(line_list[10])
            psl_dict["Qstart"] = int(line_list[11])
            psl_dict["Qend"] = int(line_list[12])
            psl_dict["Tname"] = line_list[13]
            psl_dict["Tsize"] = int(line_list[14])
            psl_dict["Tstart"] = int(line_list[15])
            psl_dict["Tend"] = int(line_list[16])
            psl_dict["blockcount"] = int(line_list[17])
            psl_dict["blocksizes"] = line_list[18]
            psl_dict["qstarts"] = line_list[19]
            psl_dict["tstarts"] = line_list[20]
            psl_list.append(psl_dict)
    return psl_list


def psl_split_by_Tname(psl_list: List[Dict]) -> Dict:
    """
    Split psl list by Tname
    """
    psl_dict = {}
    for psl_record in psl_list:
        Tname = psl_record["Tname"]
        if Tname not in psl_dict:
            psl_dict[Tname] = []
        psl_dict[Tname].append(psl_record)
    return psl_dict


def psl2gtf(gene_name, psl_records: List[Dict]) -> str:
    tx_list = []

    for psl_record in psl_records:
        tx_name = psl_record["Qname"]
        tx_ref = psl_record["Tname"]
        tx_start = psl_record["Tstart"]
        tx_end = psl_record["Tend"]
        tx_strand = psl_record["strand"]
        
        exon_start = [ int(i)  for i in psl_record["tstarts"].split(",")[:-1] ]
        exon_size = [ int(i)  for i in psl_record["blocksizes"].split(",")[:-1] ]
        exon_end = [exon_start[i] + exon_size[i] for i in range(len(exon_start))]

        tx_list.append( (tx_name, tx_ref, tx_start, tx_end, tx_strand, exon_start, exon_end) )
    # sort tx_list by tx_start
    tx_list.sort(key=lambda x: x[2])

    # gene_ref = tx_list[0][1]
    # gene_start = tx_list[0][2]
    # gene_end = tx_list[-1][3]
    # gene_strand = tx_list[0][4]


    # gene_line = f"{gene_ref}\tGMAP\tgene\t{gene_start}\t{gene_end}\t.\t{gene_strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{gene_name}\";"
    tx_lines = []
    for tx in tx_list:
        tx_name, tx_ref, tx_start, tx_end, tx_strand, exon_start, exon_end = tx
        tx_line = f"{tx_ref}\tGMAP\ttranscript\t{tx_start}\t{tx_end}\t.\t{tx_strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{tx_name}\";"
        tx_lines.append(tx_line)
        for i in range(len(exon_start)):
            exon_line = f"{tx_ref}\tGMAP\texon\t{exon_start[i]}\t{exon_end[i]}\t.\t{tx_strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{tx_name}\";"
            tx_lines.append(exon_line)
    
    return "\n".join(tx_lines)


def main( genome_file, query_file, gmap_build="gmap_build", gmap="gmap", dir_name = ".", prefix = None, is_cds = False):
    """
    Main function
    """
    if prefix:
        psl_out = f'{prefix}.psl'
        gff_out = f'{prefix}.gtf'

    index_name = genome_file.split("/")[-1].split(".")[0]

    build_gmap_index(genome_file, index_name, gmap_build, dir_name)
    
    gmap_alignment(index_name, query_file, gmap, dir_name, psl_out)
    
    psl_list = read_psl(psl_out)

    psl_dict = psl_split_by_Tname(psl_list)

    # sort psl_dict by Tstart
    for Tname, psl_records in psl_dict.items():
        psl_dict[Tname] = sorted(psl_records, key=lambda x: x["Tstart"])

    gtf_list = []

    for Tname, psl_records in psl_dict.items():

        gene_dict: Dict[str, List[Dict]] = {}

        for psl_record in psl_records:
            Qname = psl_record["Qname"]
            gene_name = ".".join(Qname.split(".")[:-1])
            if gene_name not in gene_dict:
                gene_dict[gene_name] = []
            gene_dict[gene_name].append(psl_record)
    
    
        for gene_name, psl_records in gene_dict.items():
            gtf_list.append(  psl2gtf(gene_name, psl_records) )

    # write gff file
    with open(gff_out, "w") as f:
        f.write("\n".join(gtf_list))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Using GMAP for alignment and process the output to gff format, with only exon and CDS features (if possible)')
    parser.add_argument('-b', '--build', help='gmap_build path, default is gmap_build', default="gmap_build")
    parser.add_argument('-m', '--map', help='gmap path, default is gmap', default="gmap")
    parser.add_argument('-d', '--dir', help='dir name', default=".")
    parser.add_argument("--is_cds", action="store_true", default=False, help="input file is cds file")
    parser.add_argument('-p', '--prefix', help='prefix os output file', default="out")
    parser.add_argument('genome', help='genome file')
    parser.add_argument('query', help='query file, cds or transcript')
    args = parser.parse_args()
    main(args.genome, args.query, args.build, args.map, args.dir, args.prefix, args.is_cds)
