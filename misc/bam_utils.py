import pysam
import os
import argparse
import time
import logging
import multiprocessing

NUM_RECORD = 1_000_000
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_barcodes(barcodes_file):
    with open(barcodes_file, 'r') as f:
        return [line.strip() for line in f]



def process_region(region, input_bam, barcodes):
    barcode_to_reads = {barcode: [] for barcode in barcodes}
    
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        count = 0
        start_time = time.time()

        for read in infile.fetch(region[0], region[1], region[2]):
            cb = read.get_tag("CB")
            if cb in barcode_to_reads:
                barcode_to_reads[cb].append(read)

            count += 1
            if count % NUM_RECORD == 0:
                elapsed_time = time.time() - start_time
                logging.info(f"Processed {count} records in region {region[0]}:{region[1]}-{region[2]}. Time consumed for last {NUM_RECORD:,} records: {elapsed_time:.2f} seconds.")
                start_time = time.time()


def filter_bam_by_cb(input_bam, output_dir, barcodes, num_processes=4):
    # Get the references and their lengths from the BAM file to determine regions
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        refs = infile.references
        lens = infile.lengths
    
    # Create regions based on references and their lengths
    regions = [(ref, 0, length) for ref, length in zip(refs, lens)]
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.starmap(process_region, [(region, input_bam, barcodes) for region in regions])

    
    # Collect results
    barcode_to_reads = {barcode: [] for barcode in barcodes}
    for result in results:
        for barcode, reads in result.items():
            barcode_to_reads[barcode].extend(reads)

    # Write results to output BAM files
    with pysam.AlignmentFile(input_bam, "rb") as template_file:
        for barcode, reads in barcode_to_reads.items():
            if len(reads) == 0:
                continue
            with pysam.AlignmentFile(os.path.join(output_dir, f"{barcode}.bam"), "wb", template=template_file) as outfile:
                for read in reads:
                    outfile.write(read)

    logging.info(f"Processed BAM file into {len(barcodes)} separate files based on barcodes.")


def filter_bam_by_cb_single(input_bam, output_dir, barcodes):

    barcode_to_reads = {barcode: [] for barcode in barcodes}

    with pysam.AlignmentFile(input_bam, "rb") as infile:
        count = 0
        start_time = time.time()  # Capture the start time

        for read in infile:
            cb = read.get_tag("CB")
            if cb in barcode_to_reads:
                barcode_to_reads[cb].append(read)

            count += 1
            if count % NUM_RECORD == 0:
                elapsed_time = time.time() - start_time  # Calculate elapsed time
                logging.info(f"Processed {count} records. Time consumed for last {NUM_RECORD:,} records: {elapsed_time:.2f} seconds.")
                start_time = time.time()  # Reset start time for the next 10,000 records

    # Now write out the reads for each barcode to individual BAM files
    with pysam.AlignmentFile(input_bam, "rb") as template_file:
        for barcode, reads in barcode_to_reads.items():
            with pysam.AlignmentFile(os.path.join(output_dir, f"{barcode}.bam"), "wb", template=template_file) as outfile:
                for read in reads:
                    outfile.write(read)

    logging.info(f"Processed BAM file into {len(barcodes)} separate files based on barcodes.")



# bam compare by read name
import pysam
import concurrent.futures
import pickle

def generate_hierarchy_dict_for_bam(bam_file_path, tag_value):
    hierarchy_dict = {}
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file:
            components = read.query_name.split(":")[3:]
            update_hierarchy_dict(hierarchy_dict, components, tag_value)
    return hierarchy_dict

def update_hierarchy_dict(hierarchy_dict, components, value):
    if len(components) == 1:
        hierarchy_dict[components[0]] = value
    else:
        if components[0] not in hierarchy_dict:
            hierarchy_dict[components[0]] = {}
        update_hierarchy_dict(hierarchy_dict[components[0]], components[1:], value)

def check_existence_in_hierarchy_dict(hierarchy_dict, components):
    current_dict = hierarchy_dict
    for component in components:
        if component in current_dict:
            current_dict = current_dict[component]
        else:
            return False  # Not found
    return True  # Found

def process_bam_and_write_unique(bam_file_path, output_file_path, other_hierarchy_dict):
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file, \
        pysam.AlignmentFile(output_file_path, "wb", template=bam_file) as out_file:
        for read in bam_file:
            components = read.query_name.split(":")[3:]
            if not check_existence_in_hierarchy_dict(other_hierarchy_dict, components):
                out_file.write(read)

def bam_compare_by_read_name(bam1_file_path, bam2_file_path):

    """
    read-name pattern: A00234:1245:HHKYLDSX7:4:2415:5755:7044
    """

    out_bam1_name = bam1_file_path.replace(".bam", "_unique.bam")
    out_bam2_name = bam2_file_path.replace(".bam", "_unique.bam")


    # 使用ThreadPoolExecutor并行生成两个BAM文件的递归字典
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        future_bam1 = executor.submit(generate_hierarchy_dict_for_bam, bam1_file_path, 1)
        future_bam2 = executor.submit(generate_hierarchy_dict_for_bam, bam2_file_path, 2)
        hierarchy_dict_bam1 = future_bam1.result()
        hierarchy_dict_bam2 = future_bam2.result()

    # save the hierarchy_dict_bam1 and hierarchy_dict_bam2 to file as pickle
    out_bam1_pickle = bam1_file_path.replace(".bam", "_hierarchy_dict.pickle")
    out_bam2_pickle = bam2_file_path.replace(".bam", "_hierarchy_dict.pickle")
    with open(out_bam1_pickle, "wb") as f:
        pickle.dump(hierarchy_dict_bam1, f)
    with open(out_bam2_pickle, "wb") as f:
        pickle.dump(hierarchy_dict_bam2, f)


    # 使用ThreadPoolExecutor并行处理两个BAM文件
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(process_bam_and_write_unique, bam1_file_path, out_bam1_name, hierarchy_dict_bam2)
        executor.submit(process_bam_and_write_unique, bam2_file_path, out_bam2_name, hierarchy_dict_bam1)





# bam_compare_by_read_name_set
def read_bam_file(bam_file_path):
    read_names_set = set()
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # 没比过加和不加.fetch()的速度，我猜直接遍历AlignmentFile会快一点
        for read in bam_file:  
            read_names_set.add(read.query_name)
    return read_names_set

def write_unique_reads(bam_file_path, out_bam_name, unique_reads_set):
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file, pysam.AlignmentFile(out_bam_name, "wb", template=bam_file) as out_file:
        for read in bam_file: 
            if read.query_name in unique_reads_set:
                out_file.write(read)

def bam_compare_by_read_name_set(bam1_file_path, bam2_file_path):
    # 生成输出文件名
    out_bam1_name = bam1_file_path.replace(".bam", "_unique.bam")
    out_bam2_name = bam2_file_path.replace(".bam", "_unique.bam")

    # 使用ProcessPoolExecutor并行读取BAM文件
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        future_bam1 = executor.submit(read_bam_file, bam1_file_path)
        future_bam2 = executor.submit(read_bam_file, bam2_file_path)
        read_names_set1, read_names_set2 = future_bam1.result(), future_bam2.result()

    # 使用集合运算找到各自独有的read-names
    unique_to_bam1 = read_names_set1 - read_names_set2
    unique_to_bam2 = read_names_set2 - read_names_set1

    # 并行写入独有的reads到新的BAM文件
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        executor.submit(write_unique_reads, bam1_file_path, out_bam1_name, unique_to_bam1)
        executor.submit(write_unique_reads, bam2_file_path, out_bam2_name, unique_to_bam2)


def parse_args():
    parser = argparse.ArgumentParser(description="Filter a BAM file by cell barcodes and split into separate files based on barcodes.")

    parser.add_argument("-bf", "--barcode-file", type=str, help="Path to the file containing the list of barcodes.")
    parser.add_argument("-b", "--barcode", type=str, help="Barcode to filter on.")

    parser.add_argument("-od", "--output-dir", type=str, help="Path to the output directory where the separate BAM files will be saved.", default="./")
    
    parser.add_argument("input_bam", type=str, help="Input BAM file to filter and split.")
    
    return parser.parse_args()

if __name__ == "__main__":
    
    args = parse_args()

    if not (args.barcode_file or args.barcode):
        print("Error: At least one of --barcodes-file or --barcode must be provided.")
        exit(1)
    
    if args.barcode_file and args.barcode:
        print("Error: Only one of --barcodes-file or --barcode should be provided.")
        exit(1) 
    
    if args.barcode_file:
        barcode = load_barcodes(args.barcode_file)
    else:
        barcode = args.barcode

    # Measure the size of the input file in bytes
    input_size_bytes = os.path.getsize(args.input_bam)

    # Capture the start time
    start_time = time.time()

    filter_bam_by_cb_single(args.input_bam, args.output_dir, barcode)
    
    # Capture the end time
    end_time = time.time()

    # Calculate and print the runtime
    runtime_seconds = end_time - start_time
    logging.info(f"Total runtime: {runtime_seconds:.2f} seconds")

    # Calculate and print the I/O speed in MB/s
    io_speed_mbps = (input_size_bytes / (1024 * 1024)) / runtime_seconds
    logging.info(f"I/O speed: {io_speed_mbps:.2f} MB/s")

    




