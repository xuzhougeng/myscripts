import pyfastx
import concurrent.futures


def process_batch(reads):
    sequences = ["GA", "AG", "CT", "TC"]
    counts_per_read = []
    
    # 确保read是正确的对象
    for read in reads:
        seq_name = read[0]
        seq = read[1]  # 假设read是一个tuple，其中第二个元素是序列
        seq_len = len(seq)
        counts = {sequence: seq.count(sequence) for sequence in sequences}
        counts['seq_name'] = seq_name
        counts['seq_len'] = seq_len
        counts_per_read.append(counts)

    
    return counts_per_read

def process_fastq_in_batches(fastq_path, batch_size=10000):
    # 使用pyfastx读取FASTQ文件
    fq = pyfastx.Fastq(fastq_path, build_index=False)
    
    all_counts = []
    
    # 使用concurrent.futures并行处理
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        batch_reads = []
        
        # 手动迭代Fastq对象来管理批次
        for read in fq:
            batch_reads.append(read)
            if len(batch_reads) == batch_size:
                futures.append(executor.submit(process_batch, batch_reads))
                batch_reads = []  # 重置批次列表
        
        # 处理最后一个批次（如果它不为空）
        if batch_reads:
            futures.append(executor.submit(process_batch, batch_reads))
        
        # 收集结果
        for future in concurrent.futures.as_completed(futures):
            all_counts.extend(future.result())
    
    return all_counts


if __name__ == "__main__":
    import sys
    fastq_path = sys.argv[1]
    counts = process_fastq_in_batches(fastq_path)

    # save to file
    with open(fastq_path + '.counts', 'w') as f:
        # write the header
        f.write("GA\tAG\tCT\tTC\tseq_len\tseq_name\n")
        for count in counts:
            f.write(f"{count['GA']}\t{count['AG']}\t{count['CT']}\t{count['TC']}\t{count['seq_len']}\t{count['seq_name']} \n")