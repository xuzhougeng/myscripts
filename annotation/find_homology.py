import argparse

def read_blast_results(filename, evalue_thresh=1e-5):
    """
    读取BLAST结果文件，返回每个查询序列的所有匹配结果，按得分排序，应用E-value阈值过滤。
    """
    results = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            query_id, sbj_id, score, e_value = parts[0], parts[1], float(parts[11]), float(parts[10])
            if e_value <= evalue_thresh:  # 应用E-value阈值过滤
                if query_id not in results:
                    results[query_id] = []
                results[query_id].append((sbj_id, score))
        # 对每个查询的结果按得分排序
        for query_id in results:
            results[query_id].sort(key=lambda x: x[1], reverse=True)
    return results

def find_rbh(qry_sbj, sbj_qry):
    rbh_pairs = []
    for query_id, matches in qry_sbj.items():
        if matches:
            best_match = matches[0][0]  # 基于得分排序取得最佳匹配的基因ID
            if best_match in sbj_qry and sbj_qry[best_match]:
                reverse_best_match = sbj_qry[best_match][0][0]
                if query_id == reverse_best_match:
                    rbh_pairs.append((query_id, best_match))
    return rbh_pairs

def find_extended_matches(qry_sbj, sbj_qry, rbh_pairs):
    rbh_dict = {qry: sbj for qry, sbj in rbh_pairs}  
    rbh_rev_dict =  {sbj: qry for qry, sbj in rbh_pairs}  
    extended_matches = {'rank2': [] , 'rank3': [], 'rank4': [], 'rank5':[] }
    for qry_id, matches in qry_sbj.items():
        if qry_id in rbh_dict:  # Skip if already part of RBH
            continue
        best_match = matches[0][0]
        if best_match in rbh_rev_dict:
            if qry_id in [x[0] for x in sbj_qry[best_match] ]:
                 extended_matches['rank2'].append( (qry_id, best_match, rbh_rev_dict[best_match]) )
            else:
                 extended_matches['rank3'].append( (qry_id, best_match, rbh_rev_dict[best_match]) )
        else:
            is_found = False
            for sbj_id,_ in matches:
                if sbj_id in sbj_qry and   qry_id in  [x[0] for x in sbj_qry[sbj_id] ]:
                     extended_matches['rank4'].append( (qry_id, sbj_id, '-') )
                     is_found = True
                     break
            if not is_found:    
                extended_matches['rank5'].append( (qry_id, '-', '-') )
                
    return extended_matches

def main():
    parser = argparse.ArgumentParser(description="Find RBH and extended matches from BLAST results with E-value filtering.")
    parser.add_argument("qry_sbj_file", type=str, help="BLAST result file from query to subject.")
    parser.add_argument("sbj_qry_file", type=str, help="BLAST result file from subject to query.")
    parser.add_argument("outfile", type=str, help="output file.")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for filtering (default: 1e-5).")
    
    args = parser.parse_args()

    qry_sbj = read_blast_results(args.qry_sbj_file, args.evalue)
    sbj_qry = read_blast_results(args.sbj_qry_file, args.evalue)

    rbh_pairs = find_rbh(qry_sbj, sbj_qry)
    print(f"Found {len(rbh_pairs)} RBH pairs.")

    extended_matches = find_extended_matches(qry_sbj, sbj_qry, rbh_pairs)
    
    with open(args.outfile, 'w') as fh:
        for qry_id, sbj_id in rbh_pairs:
            fh.write(f'{qry_id}\t{sbj_id}\trank1\trbh\n')
        for rank, gene_pairs in extended_matches.items():
            for qry_id, sbj_id, sub_rbh in gene_pairs:
                fh.write(f'{qry_id}\t{sbj_id}\t{rank}\t{sub_rbh}\n')
            


if __name__ == "__main__":
    main()
