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

    




