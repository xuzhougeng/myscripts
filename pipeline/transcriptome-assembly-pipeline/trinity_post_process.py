import Bio
from Bio import SeqIO
import re

def parse_cds_file(cds_file):
    """
    Parses a CDS file and returns three dictionaries containing gene IDs, transcript IDs, and protein IDs.

    Args:
        cds_file (str): The path to the CDS file.

    Returns:
        tuple: A tuple containing three dictionaries:
            - g_id_dict: A dictionary containing gene IDs as keys and a list of transcript IDs as values.
            - t_id_dict: A dictionary containing transcript IDs as keys and a list of protein IDs as values.
            - p_id_dict: A dictionary containing protein IDs as keys and their descriptions as values.
    """
    g_id_dict = {}
    t_id_dict = {}
    p_id_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        p_id = record.id
        t_id = p_id.split(".")[0]
        g_id = "_".join(t_id.split("_")[:-1])
        
        if g_id in g_id_dict:
            t_ids = g_id_dict[g_id]
            if t_id  not in t_ids:
                g_id_dict[g_id].append(t_id)
        else:
            g_id_dict[g_id] = [t_id]
        
        if t_id in t_id_dict:
            p_ids = t_id_dict[t_id]
            if p_id not in p_ids:
                t_id_dict[t_id].append(p_id)
        else:
            t_id_dict[t_id] = [p_id]
        
        p_id_dict[p_id] = record.description

    return g_id_dict, t_id_dict, p_id_dict

def stat(g_id_dict, t_id_dict):
    """
    This function takes in two dictionaries, g_id_dict and t_id_dict, where the keys are gene IDs and the values are lists of
    transcript IDs. It then calculates the average and maximum number of transcripts per gene for both dictionaries and prints
    the results.

    Args:
    - g_id_dict (dict): A dictionary where the keys are gene IDs and the values are lists of transcript IDs.
    - t_id_dict (dict): A dictionary where the keys are gene IDs and the values are lists of predicted transcript IDs.

    Returns:
    - None
    """
    # 统计每个基因的平均转录本数和最多转录本数
    gene_transcript_counts = []
    max_transcript_count = 0

    for gene_id, transcript_list in g_id_dict.items():
        transcript_count = len(transcript_list)
        gene_transcript_counts.append(transcript_count)
        max_transcript_count = max(max_transcript_count, transcript_count)

    average_transcript_count = sum(gene_transcript_counts) / len(gene_transcript_counts)

    print("平均转录本数：", average_transcript_count)
    print("最多转录本数：", max_transcript_count)

    # 统计每个基因的平均转录本数和最多转录本数
    gene_transcript_counts = []
    max_transcript_count = 0

    for gene_id, transcript_list in t_id_dict.items():
        transcript_count = len(transcript_list)
        gene_transcript_counts.append(transcript_count)
        max_transcript_count = max(max_transcript_count, transcript_count)

    average_transcript_count = sum(gene_transcript_counts) / len(gene_transcript_counts)

    print("平均预测：", average_transcript_count)
    print("最多预测：", max_transcript_count)

def filter_proteins(g_id_dict, t_id_dict, p_id_dict):
    """
    Filters out incomplete proteins from a dictionary of gene IDs, transcript IDs, and protein IDs.
    If a transcript has no protein, the program will exit with an error message.
    If a transcript has more than one complete protein, the one with the highest score will be chosen.
    If a transcript has only one complete protein, it will be kept.
    If a transcript has no complete protein, it will be deleted.
    If a gene has no complete transcript, it will be deleted.
    
    Args:
    - g_id_dict: a dictionary of gene IDs as keys and lists of transcript IDs as values
    - t_id_dict: a dictionary of transcript IDs as keys and lists of protein IDs as values
    - p_id_dict: a dictionary of protein IDs as keys and their descriptions as values
    
    Returns:
    - g_id_dict: the filtered gene ID dictionary
    - t_id_dict: the filtered transcript ID dictionary
    """
    g_id_to_delete = []
    t_id_to_delete = []

    for g_id in g_id_dict.keys():

        t_ids = g_id_dict[g_id]
        complete_tx = 0 # flag for whether the gene is complete

        for t_id in t_ids:

            p_ids = t_id_dict[t_id]

            for p_id in p_ids:
                if len(p_id) == 0:
                    # terminate the program if there is no protein
                    print("Error: no protein for transcript", t_id)
                    exit(1)
                # filter out the incomplete protein
                if "complete" not in p_id_dict[p_id]:
                    t_id_dict[t_id].remove(p_id)

            complete_num = len(t_id_dict[t_id])
            # if there is only one complete protein, keep it
            if complete_num == 1:
                complete_tx += 1
                continue

            # if there are more than one complete protein, choose the one with highest score
            if complete_num > 1:
                complete_tx += 1
                best_p_id = ""
                current_best_score = 0
                
                p_ids = t_id_dict[t_id]
                for p_id in p_ids:
                    # extract the score from description with regular expression, the pattern is score=xxx.yyy or score=-xxx.yyy
                    #score = float(re.findall(r"score=(\d+\.\d+)", p_id_dict[p_id])[0])
                    score = float(re.findall(r"score=(-?\d+\.\d+)", p_id_dict[p_id])[0])
                    if score > current_best_score:
                        best_p_id = p_id
                        current_best_score = score
                
                t_id_dict[t_id] = [best_p_id]
    
        
            # if there is no complete protein, delete the transcript
            if complete_num == 0:
                t_id_to_delete.append(t_id)
        
        # if there is no complete protein, delete the gene
        if complete_tx == 0:
            g_id_to_delete.append(g_id)


    for g_id in g_id_to_delete:
        del g_id_dict[g_id]

    for t_id in t_id_to_delete:
        if t_id in t_id_dict:
            del t_id_dict[t_id]
    
    return g_id_dict, t_id_dict

def extract_fasta_to_file(output_file, fasta_file, g_id_dict, t_id_dict, prefix=None):
    """
    Extracts transcripts from a CDS file and writes them to an output file.

    Args:
    output_file (str): The path to the output file.
    cds_file (str): The path to the CDS file.
    g_id_dict (dict): A dictionary containing gene IDs as keys and a boolean value indicating whether the gene should be included as values.
    t_id_dict (dict): A dictionary containing transcript IDs as keys and a set of protein IDs associated with each transcript as values.

    Returns:
    None
    """
    
    gene_count = 0
    new_name_dict = {}

    delim = "."

    with open(output_file, "w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            p_id = record.id
            t_id = p_id.split(".")[0]
            g_id = "_".join(t_id.split("_")[:-1])

            # if gene id deleted , skip
            if not g_id in g_id_dict:
                continue

            # if the transcript is deleted, skip
            if not t_id in t_id_dict:
                continue

            if p_id in t_id_dict[t_id]:
                # update the header with transcript ID
                if prefix is  None:
                    record.id = t_id
                else:
                    if g_id not in new_name_dict:
                        gene_count += 1
                        new_name_dict[g_id] = f"{prefix}{delim}G{gene_count:05d}"
                    record.id = new_name_dict[g_id] + {delim} + t_id.split("_")[-1]

                record.description = ""

                # write the record to output file
                SeqIO.write(record, f, "fasta")

    return new_name_dict


def save_gene_transcript_mapping_to_file(out_file, g_id_dict):
    """
    Save the gene id to transcript id mapping to a file.

    Args:
        out_file (str): The path to the output file.
        g_id_dict (dict): A dictionary containing gene ids as keys and a list of corresponding transcript ids as values.

    Returns:
        None
    """
    # save the gene id to transcript id mapping
    with open(out_file, "w") as f:
        for g_id in g_id_dict:
            t_ids = g_id_dict[g_id]
            for t_id in t_ids:
                f.write(g_id + "\t" + t_id + "\n")


def transdecoder_filter(args):
    cds_in = args.cds_in
    cds_out = args.cds_out
    pep_in = args.pep_in
    pepe_out = args.pep_out
    prefix = args.prefix
    mapping_file = args.mapping

    g_id_dict, t_id_dict, p_id_dict = parse_cds_file(cds_in)

    stat(g_id_dict, t_id_dict)

    g_id_dict, t_id_dict = filter_proteins(g_id_dict, t_id_dict, p_id_dict)

    new_name_dict = extract_fasta_to_file(cds_out, cds_in, g_id_dict, t_id_dict, prefix)
    if len(new_name_dict) > 0:
        with open(mapping_file, "w") as f:
            for k,v in new_name_dict.items():
                 f.write(f"{k}\t{v}\n")

    extract_fasta_to_file(pepe_out, pep_in, g_id_dict, t_id_dict, prefix)

    #save_gene_transcript_mapping_to_file(out_mapping_file, g_id_dict, prefix)

def trinity_extract(args):
    fasta_file = args.input
    output_file = args.output
    mapping_file = args.mapping

    # read mapping file
    mapping_dict = {}
    with open(mapping_file, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            items = line.split("\t")
            mapping_dict[items[0]] = items[1]
    
    # read trinity.fasta file
    with open(output_file, "w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            t_id = record.id
            g_id = "_".join(t_id.split("_")[:-1])

            # if gene id deleted , skip
            if not g_id in mapping_dict:
                continue

            record.id = mapping_dict[g_id] + "_" + t_id.split("_")[-1]  

            record.description = ""

            # write the record to output file
            SeqIO.write(record, f, "fasta")

def longest_extract(args):
    """ extract the longest cds/pep for each gene from input.fasta
    """
    fasta_file = args.input
    output_file = args.output
    
    gene_dict = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            t_id = record.id
            t_size = len(record.seq)
            g_id = "_".join(t_id.split("_")[:-1])
            if g_id in gene_dict:
                if gene_dict[g_id][1] < t_size:
                    gene_dict[g_id] = (t_id, t_size, record.seq)
            else:
                gene_dict[g_id] = (t_id, t_size, record.seq)
    
    # write the longest transcript to output file
    with open(output_file, "w") as f:
        for k,v in gene_dict.items():
            f.write(f">{v[0]}\n{v[2]}\n")


def refine_lace_gff(args):
    """
    refine lace gff file for STAR
    """
    gff_file = args.input

    with open(gff_file, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            items = line.split("\t")
            if items[2] == "gene":
                items[2] = "transcript"
                items[8] = items[8].replace("gene_id", "transcript_id")
                print("\t".join(items))
            elif items[2] == "transcript":
                items[2] = "exon"
                items[8] = items[8].replace("transcript_id", "exon_id")
                print("\t".join(items))
            else:
                print(line) 


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="This script filters out incomplete proteins from a Trinity CDS file.")

    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help')

    # subparser for the filter command
    parser_filter = subparsers.add_parser('filter', help='filter transdecoder cds and pepe file ') 
    parser_filter.add_argument("-c", "--cds_in", type=str, help="The path to the Trinity CDS file.", required=True)
    parser_filter.add_argument("-p", "--pep_in", type=str, help="The path to the Trinity protein file.", required=True)
    parser_filter.add_argument("-o", "--cds_out", type=str, help="The path to the output CDS file.", required=True)
    parser_filter.add_argument("-e", "--pep_out", type=str, help="The path to the output protein file.", required=True)
    parser_filter.add_argument("--prefix", type=str, help="prefix for cds/pep ", required=False)
    parser_filter.add_argument("-m", "--mapping", type=str, help="The path to the output mapping file.", required=False)
    parser_filter.set_defaults(func=transdecoder_filter)  

    # subparser for the extract command
    parser_extract = subparsers.add_parser('extract', help='extract Trinity.fasta from map.txt generated by filter')
    parser_extract.add_argument("-i", "--input", type=str, help="The path to the Trinity.fasta file.", required=True)
    parser_extract.add_argument("-o", "--output", type=str, help="The path to the output file.", required=True)
    parser_extract.add_argument("-m", "--mapping", type=str, help="The path to the mapping file.", required=True)
    parser_extract.set_defaults(func=trinity_extract)

    # subparser for the longest command
    parser_longest = subparsers.add_parser('longest', help='extract longest transcript from Trinity.fasta')
    parser_longest.add_argument("-i", "--input", type=str, help="The path to the input.fasta file.", required=True)
    parser_longest.add_argument("-o", "--output", type=str, help="The path to the output file.", required=True)
    parser_longest.set_defaults(func=longest_extract)
    # subparser for 



    args = parser.parse_args()
    
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
