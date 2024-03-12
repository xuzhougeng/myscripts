import re

def read_fasta(fasta_file):
    """
    Read a FASTA file and return a dictionary with sequence identifiers as keys
    and sequences as values. Sequences are joined at the end to reduce memory usage.
    
    Parameters:
    fasta_file (str): Path to the FASTA file.
    
    Returns:
    dict: Dictionary with sequence identifiers as keys and sequences as values.
    """
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence identifier
                if sequence_id is not None:
                    sequences[sequence_id] = ''.join(sequence)  # Store the previous sequence
                sequence_id = line[1:]  # Remove the '>' character
                sequence = []  # Initialize a new sequence list
            else:
                sequence.append(line)
        if sequence_id is not None:
            sequences[sequence_id] = ''.join(sequence)  # Store the last sequence
    return sequences
    
def get_chromosome_sizes(sequences):
    """
    Calculate the sizes of chromosomes from the given sequences.

    Parameters:
    sequences (dict): Dictionary with chromosome identifiers as keys and sequences as values.

    Returns:
    dict: A dictionary where keys are chromosome identifiers and values are their sizes (length of the sequence).
    """
    chrom_sizes = {seq_id: len(sequence) for seq_id, sequence in sequences.items()}
    return chrom_sizes


def find_continuous_n_regions(sequences):
    """
    Find continuous regions of the base 'N' in each sequence using regular expressions
    for improved performance.
    
    Parameters:
    sequences (dict): Dictionary with sequence identifiers as keys and sequences as values.
    
    Returns:
    dict: A dictionary where keys are sequence identifiers, and values are lists of tuples,
          each tuple representing the start and end positions of a continuous 'N' region.
    """
    n_regions = {}
    for seq_id, sequence in sequences.items():
        # Use regular expression to find all continuous 'N' regions
        regions = [(match.start() + 1, match.end()) for match in re.finditer(r'N+', sequence)]
        if regions:
            n_regions[seq_id] = regions
    return n_regions



def create_genome_split_position(chrom_size, n_regions, max_chrom_size):
    result = []
    for chrom, size in chrom_size.items():
        if size <= max_chrom_size:
            # 染色体长度小于等于 max_chrom_size，直接输出
            result.append((chrom, 0, size, chrom ))
        else:
            # 需要分割
            splits = n_regions.get(chrom, [])
            start_pos = 0
            end_pos = 0
            chrom_counter = 1
            for i, split in enumerate(splits):
                # 尽可能在 N 区域分割，同时确保分割后的长度不超过 max_chrom_size
                next_split_start = size if i + 1 == len(splits) else splits[i + 1][0]
                if split[0] - start_pos <= max_chrom_size and next_split_start - start_pos > max_chrom_size:
                    end_pos = split[0]
                    result.append((chrom, start_pos, end_pos, chrom + f"00{chrom_counter}"))
                    start_pos = split[1]  # 下一个分割开始位置是当前 'N' 区域的结束位置
                    chrom_counter += 1
            # 添加最后一个分割，如果有的话
            if start_pos < size:
                result.append((chrom, start_pos, size, chrom + f"00{chrom_counter}"))
    return result




def create_split_fasta(sequences , split_info, output_file):
    """
    Create a single split FASTA file directly from the original FASTA based on the provided chromosome split information.

    Parameters:
    fasta_file (str): Path to the original FASTA file.
    split_info (list of tuples): Split information as returned by `create_genome_split_position`.
    output_file (str): Path to the output file for the split FASTA sequences.
    """
    with open(output_file, 'w') as outfile:
        for chrom, start, end, new_chrom in split_info:
            if chrom in sequences:
                sequence = sequences[chrom][start:end]
                outfile.write(f'>{new_chrom} {chrom}:{start}-{end}\n')
                for i in range(0, len(sequence), 80):
                    outfile.write(sequence[i:i+80] + '\n')

def create_adjusted_gtf(gtf_file, split_info, output_file):
    """
    Create a single adjusted GTF file directly from the original GTF based on chromosome split information.

    Parameters:
    gtf_file (str): Path to the original GTF file.
    split_info (list of tuples): Split information as returned by `create_genome_split_position`.
    output_file (str): Path to the output file for the adjusted GTF entries.
    """
    # Preprocess split_info for quick access during adjustments
    split_map = {}
    for chrom, start, end, new_chrom in split_info:
        if chrom not in split_map:
            split_map[chrom] = []
        split_map[chrom].append((start, end, new_chrom))
    
    with open(output_file, 'w') as outfile:
        with open(gtf_file, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue  # Ignore header lines
                parts = line.strip().split('\t')
                chrom, feature_start, feature_end = parts[0], int(parts[3]), int(parts[4])
                if chrom in split_map:
                    for start, end, new_chrom in split_map[chrom]:
                        if feature_start >= start and feature_end <= end:
                            adjusted_start = feature_start - start
                            adjusted_end = feature_end - start
                            parts[0] = new_chrom
                            parts[3] = str(adjusted_start)
                            parts[4] = str(adjusted_end)
                            adjusted_line = '\t'.join(parts)
                            outfile.write(adjusted_line + '\n')
                            break


def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Split FASTA and adjust GTF based on chromosome size and N regions.')
    parser.add_argument('-f', '--fasta', required=True, help='Path to the FASTA file.')
    parser.add_argument('-g', '--gtf', required=True, help='Path to the GTF file.')
    parser.add_argument('-m', '--maxsize', type=int, default=500000000, help='Maximum chromosome size for splitting.')
    parser.add_argument('-o', '--output', default='split', help='Prefix for the output FASTA and GTF files. Default is "split".')

    args = parser.parse_args()

    if args.maxsize > 1000000000:
        print("too large. decrease it")
        sys.exit(-1)

    # Read and process the FASTA file
    sequences = read_fasta(args.fasta)
    chrom_size = get_chromosome_sizes(sequences)
    n_regions = find_continuous_n_regions(sequences)
    result = create_genome_split_position(chrom_size, n_regions, args.maxsize)

    # Create split FASTA and adjusted GTF files
    fasta_output = f'{args.output}.fasta'
    gtf_output = f'{args.output}.gtf'
    create_split_fasta(sequences, result, fasta_output)
    create_adjusted_gtf(args.gtf, result, gtf_output)

    print(f'Split FASTA file saved as {fasta_output}')
    print(f'Adjusted GTF file saved as {gtf_output}')

def main2():
    import sys
    fasta = sys.argv[1]
    gff = sys.argv[2]
    maxsize = int(sys.argv[3])
    outdir = sys.argv[4]
    if maxsize > 1000000000:
        print("too large. decrease it")
        sys.exit(-1)

    # Read and process the FASTA file
    sequences = read_fasta(fasta)
    chrom_size = get_chromosome_sizes(sequences)
    n_regions = find_continuous_n_regions(sequences)
    result = create_genome_split_position(chrom_size, n_regions, maxsize)
    
    fasta_output = f'{outdir}/splited.fasta'
    gtf_output = f'{outdir}/splited.gtf'
    create_split_fasta(sequences, result, fasta_output)
    create_adjusted_gtf(gff, result, gtf_output)



if __name__ == '__main__':
    main()