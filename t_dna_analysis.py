#!/usr/bin/env python3

################################################################################################
#### Author: xuzhougeng
#### Requriements:
####   - Linux
####   - Python3
####   - Python3/pysam
####   - Minimap2
####   - SPades
####   - NCBI-BLAST
####   - seqkit
####   - matplotlib
####   - numpy
####   - sklearn
################################################################################################

import sys
import argparse
import pysam
from os import system
from os import path
from os import mkdir


def has_tool(name):
    """
    if tool existed return True otherwise False
    """

    from shutil import which

    return which(name) is not None

def check_tools():

    for tool in ["bwa", "samtools", "spades.py", "blastn", "seqkit", "minimap2"]:
        if not has_tool(tool):
            print("{} is not avaible".format(tool))
            sys.exit(1)


# align 
def bwa_mem_align(index, prefix, threads):
    """
    Parameters:
    ----------
    index: str
        the index of bwa mem
    prefix: str
        prefix is the sample name, 
        e.g CYJ029-10-1_R1.fq.gz CYJ029-10-1_R2.fq.gz, the prefix is CYJ029-10-1
    threads: int
        the threads for bwa mem

    Returns:
    ---------- 
    bamfile: str
        the bam file from bwa and sort by samtools sort
    """

    threads2 = str(round(int(threads) / 2 ))
    threads  = str(threads)
    fastq1   = "{}_R1.fq.gz".format(prefix)
    fastq2   = "{}_R2.fq.gz".format(prefix)
    bamfile  = "{}_sort.bam".format(prefix)

    if not path.isfile(bamfile):
        print("Processing {} and {} with bwa mem".format(fastq1, fastq2))
        bwa_mem_cmd = """
        bwa mem -v 2 -t {THREADS} -R "@RG\\tID:{SM}\\tSM:{SM}\\tPL:ILLUMINA" {INDEX} {FQ1} {FQ2} 2> bwa.log | \
                samtools sort -@ {THREADS2} > {OUT}  && samtools index -@ {THREADS2} {OUT} || rm {OUT}
        """.format(THREADS=threads, SM=prefix, INDEX=index,
            FQ1=fastq1, FQ2=fastq2, THREADS2=threads2, OUT=bamfile)
        print(bwa_mem_cmd)

        system(bwa_mem_cmd)

    return bamfile

def get_soft_clip(bamfile, prefix, disimilarity=10):
    """get soft clip reads from BAM
    Parameters
    ----------
    bamfile: str
        bam file name 
    disiilartiy: int
        disimilarity of alignment, from 0 to 100, default is 10

    Returns:
    ---------
    fqfile: str
        fastq file of soft clip reads
    """

    disimilarity = args.disimilarity / 100
    bamFile = pysam.AlignmentFile(bamfile)
    fqfile  = prefix + ".fq"

    print("Extract the soft-clip reads from bam")

    if path.isfile(fqfile):
        return fqfile

    file_out = open(fqfile, "w")

    for record in bamFile:
        align_num = 0
        clip_num = 0
        other_num = 0
        clip_percent = 0.0
        if not record.is_unmapped:
            cigar = record.cigartuples
            for op,num in cigar:
                if op == 0:
                    align_num += num
                elif op == 4:
                    clip_num += num
                else:
                    other_num +=num
            clip_percent = clip_num / (align_num + clip_num + other_num)

            if clip_percent >  disimilarity:
                 fq_record = "@{}\n{}\n+\n{}\n".format(record.qname, record.seq, record.qual)
                 file_out.writelines(fq_record)
    file_out.close()

    print("Finished Extract") 
    return fqfile

def genome_assembly(fqfile, prefix, threads):

    print("Assembly the soft-clip reads with SPades")

    contig = "{}/contigs.fasta".format(prefix)

    if not path.isfile(contig):
        spades_cmd = """
        spades.py -k 77 -t {} -s {} -o {}
        """.format(threads, fqfile, prefix)
        system(spades_cmd)
        
    print("Finished Assembly ")

    return contig

def get_depth(bamfile, contig, start, end):
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    depth_dict = {}
    for pileupcolumn in bamfile.pileup(contig = contig, 
            start = start, end=end, 
            stepper="nofilter"):
        pos = pileupcolumn.reference_pos + 1
        depth = pileupcolumn.nsegments
        depth_dict[pos] = depth

    return depth_dict

def depth_plot(data, contig, start, end, file_name):
    """save the pictures
    """
    import matplotlib.pyplot as plt
    x = []
    y = []
    title = '{}:{}-{}'.format(contig, start, end)
    
    for i in range(start, end+1):
        x.append(i)
        y.append(data.get(i, 0))

    fig = plt.figure()
    plt.plot(x, y)
    plt.xlim(left=start, right=end)
    plt.ylim(0)
    plt.title(title)
    fig.savefig(file_name)
    plt.close()

def text_plot(data, contig, start, end, file_name):
    """Plot using text
    """
    f =  open(file_name, "w")
    f.writelines("{}:{}-{}\n".format(contig, str(start), str(end)))
    for i in range(start, end + 1):
        depth = data.get(i,0)
        outlines = "{} {}:{}".format(contig, str(i), "*" * depth)
        f.writelines(outlines + "\n") 

def raw_output(data, contig, start, end, file_name):
    """Plot with raw number
    """
    f =  open(file_name, "w")
    f.writelines("#{}:{}-{}\n".format(contig, str(start), str(end)))
    for i in range(start, end + 1):
        depth = data.get(i,0)
        f.writelines(str(depth) + ",")
    f.writelines("\n")

def run_blastn(ref, prefix, contig):
    """Run blastn
    Parameters:
    ----------
    ref: str
        reference genome file path
    prefix: str
        output prefix
    contig: str
        spades assembly result
    """
    
    blastn_out = "contig_blastn.txt"
    # run blast
    if not path.isfile(blastn_out):
        blast_cmd = """
        makeblastdb -in {REF} -out tmp/{PREFIX} -dbtype nucl
        blastn -query {CONTIG} -db tmp/{PREFIX} -outfmt 6 -num_threads 20 > {OUT}
        """.format(REF=ref, PREFIX=prefix, CONTIG=contig, OUT=blastn_out)
        system(blast_cmd)

    return blastn_out
    
def blastn_filter(blastn_out, plastid, deviation, min_ident, min_cov, max_cov, min_len, max_len):
    """
    Filter BLASTN result based on the contig lengthe, contig coverage and alignment similarity

    Parameters:
    -----------
    blastn_out: str
       blastn output with outfmt 6
    min_len: int
       minimum query length
    max_len: int
       maximum query length
    min_cov: int
       minimum coverage 
    max_cov: int
       maximum coverage
    max_hit: int
       maximum coverage
    min_ident: int
       minimum identity

    Returns:
    ----------
    tmp: list
       a list stores filter records
    """

    black_list = []
    tmp = []

    for line in open(blastn_out, "r"):
        list_from_line = line.strip().split("\t")
        if list_from_line[1] in plastid:
            black_list.append(list_from_line[1])
            continue
        elif list_from_line[1] in black_list:
            continue
        elif float(list_from_line[2]) < min_ident:
            continue
        else:
            qname = list_from_line[0]

            # filter the contig whose coverage too low or too high
            qcov  = float(qname.split("_")[-1])
            if qcov < min_cov or qcov > max_cov:
                continue
            # filter the contig too short or too long
            qlen  = int(qname.split("_")[-3])
            if (qlen < min_len) or (qlen > max_len):
                continue
            # filter the hit fully match to genome
            # this should be strict
            match_len = int(list_from_line[7]) - int(list_from_line[6]) + 1
            if match_len == qlen:
                continue

            # filter the alignments of inner
            # this should be soft
            if (int(list_from_line[6]) > deviation) and \
                    (qlen - int(list_from_line[7]) > deviation):
                continue
            
            tmp.append(list_from_line)
    return tmp

def break_filter(ref, contig, start, end):

    fafile = pysam.FastaFile(ref)
    sequences = fafile.fetch(reference=contig, start=start, end=end)
    pos = sequences.find("N")
    #print("{}:{}-{} {}".format(contig, start, end, pos))
    return pos
    

def get_candidate_with_vector(ref, prefix, contig, vector):
    """
    Parameters:
    ----
    bamfile: str
    """

    mmp2_cmd = """
        minimap2 -x asm5  {CONTIG}  {VECTOR} > align.paf
        cut -f 6 align.paf > id.txt
        seqkit grep -f id.txt {CONTIG} > candidate.fa
    """.format(VECTOR = vector, CONTIG=contig)

    system(mmp2_cmd)

    blast_cmd = """
        makeblastdb -in {REF} -out tmp/{PREFIX} -dbtype nucl
        blastn -query candidate.fa -db tmp/{PREFIX} -outfmt 6 > candidate_blastn.txt
        """.format(REF=ref, PREFIX=prefix)
    system(blast_cmd)

def get_candidate_without_vector(fa, prefix, blastn_list, bamfile,
        deviation = 20, max_hit = 5, max_depth = 50, datatype = "png"):
    """
    Parameters:
    ----
    blastn_list: list
        blastn out filtered list
    prefix: str
        output directory name
    bamfile: str
        bwa-mem alignment output
    max_hit: int
        maximum coverage
    max_depth: int
        maximum depth 
    datatype: str
        output data type: png, text, pdf
    """

    # filter the query with many hits
    # count the hits of each query
    query_dict = {}
    for l in blastn_list:
        query_dict[l[0]] = query_dict.get(l[0],0) + 1
    qname_list = [key for key, value in query_dict.items() if value < max_hit ]

    prev_qname = []
    version    = 1
    flank      = 100

    # create directory for save data
    outdir = "{}_depth_pattern".format(prefix)
    if not path.isdir(outdir):
        mkdir(outdir)

    for i in range(len(blastn_list)):
        if blastn_list[i][0] in qname_list:

            ref       = blastn_list[i][1] 
            ref_left  = int(blastn_list[i][8])
            ref_right = int(blastn_list[i][9])
            # insertion in the left and reference in the right
            if int(blastn_list[i][6]) < deviation:
                start = ref_right - flank
                end   = ref_right + flank
            else:
                start = ref_left - flank
                end   = ref_left + flank

            # 
            if start < 0:
                continue

            # filter the break region
            pos = break_filter(ref = fa, contig = ref, start = start, end = end)
            #print(blastn_list[i])
            #print(pos)
            
            if pos > 0:
                continue

            # get the depth of each base in the target region    
            depth_dict = get_depth(bamfile, ref, start, end)

            # filter candidater which is to high
            mean_depth = sum(depth_dict.values()) / len(depth_dict)
            if mean_depth > max_depth:
                continue
                
            # file name 
            if blastn_list[i][0] in prev_qname: 
                file_name = "{}/{}_aln{}_depth.{}".format(outdir, blastn_list[i][0], str(version), datatype)
                version += 1
            else:
                version = 1
                file_name = "{}/{}_depth.{}".format(outdir, blastn_list[i][0], datatype)
                prev_qname.append(blastn_list[i][0])     
      
            # output the depth 
            if datatype == "txt":
                text_plot(depth_dict, ref, start, end, file_name)

            elif datatype == "png":
                depth_plot(depth_dict, ref, start, end, file_name)
            
            elif datatype == "raw":
                raw_output(depth_dict, ref, start, end, file_name)
            else:
                print("warnning: No out data type select")


def cluster_by_depth(n_clusters, prefix):
    import glob
    from sklearn.cluster import KMeans
    import numpy as np
    from shutil import rmtree
    import matplotlib.pyplot as plt

    title = []
    data  = []
    depth_files = glob.glob("{}_depth_pattern/*.raw".format(prefix))
    for file in depth_files:
        f = open(file, "r")
        for line in f:
            if line.startswith("#"):
                start,end = line.strip().split(":")[1].split("-")
                if (int(end) - int(start)) < 200:
                    next(f)
                    continue
                title.append(line.strip()[1:])
                
            else:
                data_set = [ int(i) for i in line.strip().split(",")[:-1] ]
                data.append( data_set )       
    data_array = np.array(data, dtype="uint8")
 
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data_array)
    y_kmeans = kmeans.predict(data_array) 

    for i in range(n_clusters):
        dir_name = "{}_depth_pattern/Cluster{}".format(prefix, str(i))
        if path.isdir(dir_name):
            rmtree(dir_name)
        mkdir(dir_name)

    x = np.arange(0,201)
    for i in range(0, len(y_kmeans)):
        cluster = y_kmeans[i]
        file_name = path.basename(depth_files[i])
        png_name = "{}_depth_pattern/Cluster{}/{}.png".format(prefix, cluster, file_name)
        y = data_array[i]
        fig = plt.figure()
        plt.plot(x, y)
        plt.xlim(left=0, right=201)
        plt.ylim(0)
        plt.title(title[i])
        fig.savefig(png_name)
        plt.close()    

def main(args):

    disimilarity = args.disimilarity
    threads    = args.threads
    ref        = args.ref
    index      = args.index
    prefix     = args.prefix
    vector     = args.vector

    # These parameters are used for find insertion sites without vector
    min_len    = args.min_len
    max_len    = args.max_len
    min_cov    = args.min_cov
    max_cov    = args.max_cov
    max_hit    = args.max_hit
    min_ident  = args.min_ident
    max_depth  = args.max_depth
    datatype   = args.data_type
    n_clusters = args.clusters

    bamfile = bwa_mem_align(index, prefix, threads)
    fqfile  = get_soft_clip(bamfile, prefix, disimilarity)
    contig  = genome_assembly(fqfile, prefix, threads)

    # 
    if vector is not None:
        get_candidate_with_vector(ref, prefix, contig, vector)
    else:
        blastn_out = run_blastn(ref, prefix, contig)
        plastid = ["ChrM", "ChrCh"]
        deviation = 20 
        blastn_list = blastn_filter(blastn_out, plastid, deviation, min_ident, min_cov, max_cov, min_len, max_len)
        get_candidate_without_vector(ref, prefix, blastn_list, bamfile, deviation, max_hit, max_depth, datatype)

        if datatype == "raw":
            cluster_by_depth(n_clusters, prefix)


if __name__ == "__main__":

    check_tools()

    parser = argparse.ArgumentParser()
    
    # Optional arguments
    parser.add_argument("--disimilarity", type = int, default=10,
            help="disimilrity of alignments of soft-clip reads, default is 10")
    parser.add_argument("--threads", type = int, default=80,
            help="threads of bwa mem, default is 80")    
    parser.add_argument("--vector", help="vector sequences")
    parser.add_argument("--min_len", type = int, default=200,
            help="minimum query length")
    parser.add_argument("--max_len", type = int, default=500,
            help="maxium query length")
    parser.add_argument("--min_cov", type = int, default=2,
            help="minimum contig coverage")
    parser.add_argument("--max_cov", type = int, default=20,
            help="maxium contig coverage")
    parser.add_argument("--max_hit", type = int, default=5,
            help="maxium query associated hits")
    parser.add_argument("--min_ident", type = int, default=90,
            help="minimum identity between query and reference")
    parser.add_argument("--max_depth", type = int, default=50,
            help="maxium depth of query associated region")
    parser.add_argument("--data_type", default="raw",
            help="data type for depth output: png, text, raw")
    parser.add_argument("-k", "--clusters", type = int, default=21,
            help="cluster number for K-means clustering")

    # Required arguments
    parser.add_argument("ref", help="reference fasta")
    parser.add_argument("index", help="index of species")
    parser.add_argument("prefix", help ="sample name")

    args = parser.parse_args()
    sys.exit(main(args))

