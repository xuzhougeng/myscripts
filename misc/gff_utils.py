import re
import argparse


## fix the gff3 file
def parse_gff3_line( parts ):
    """解析GFF3文件的一行，返回一个字典"""
    attributes_text = parts[8]
    attributes = {}
    for attr in attributes_text.split(';'):
        key, value = attr.split('=') if '=' in attr else (attr, '')
        key = key.strip()
        value = value.strip()
        if len(value) > 0 and len(key) > 0:
            attributes[key] = value
    # 处理科学计数法
    start = int(float(parts[3])) if re.match(r'\d+\.?\d*e[\+\-]?\d+', parts[3]) else parts[3]
    end = int(float(parts[4])) if re.match(r'\d+\.?\d*e[\+\-]?\d+', parts[4]) else parts[4]

    return {
        'seqid': parts[0],
        'source': parts[1],
        'type': parts[2],
        'start': str(start),
        'end': str(end),
        'score': parts[5],
        'strand': parts[6],
        'phase': parts[7],
        'attributes': attributes
    }

def read_gff3(filename):
    """读取GFF3文件，返回一个包含所有行的字典"""
    gene_dict = {} # 记录基因
    gene2mRAN_dict = {} # 记录基因对应的mRNA
    mRNA2other_dict = {} # 记录mRNA对应的其他特征

    for line in open(filename):
        line = line.strip()
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 9:
            continue
        feature = parse_gff3_line(parts)
        if feature['type'] == 'gene':
            gene_dict[feature['attributes']['ID']] = feature
        elif feature['type'] == 'mRNA' or feature['type'] == 'transcript':
            gene_id = feature['attributes']['Parent']
            if gene_id not in gene2mRAN_dict:
                gene2mRAN_dict[gene_id] = []
            gene2mRAN_dict[gene_id].append(feature)
        else:
            mRNA_id = feature['attributes']['Parent']
            if mRNA_id not in mRNA2other_dict:
                mRNA2other_dict[mRNA_id] = []
            mRNA2other_dict[mRNA_id].append(feature)
    
    return {
        'gene': gene_dict,
        'mRNA': gene2mRAN_dict,
        'other': mRNA2other_dict
    }


def add_ids_to_features(mRNA2other_dict):
    """检查每个特征是否有ID，如果没有，则添加ID"""
    for mRNA_id, features in mRNA2other_dict.items():
        for index, feature in enumerate(features):
            if 'ID' not in feature['attributes']:
                # 构建一个唯一的ID：使用mRNA的ID，特征的类型和索引
                unique_id = f"{mRNA_id}_{feature['type']}_{index}"
                feature['attributes']['ID'] = unique_id

def build_gff3_line(feature):
    """构建GFF3文件的一行"""
    attributes = feature['attributes']
    attributes_text = ';'.join([f"{key}={value}" for key, value in attributes.items()])
    return f"{feature['seqid']}\t{feature['source']}\t{feature['type']}\t{feature['start']}\t{feature['end']}\t{feature['score']}\t{feature['strand']}\t{feature['phase']}\t{attributes_text}"

def build_gff3_lines(gff_dict):
    gene, mRNA, other = gff_dict['gene'], gff_dict['mRNA'], gff_dict['other']
    lines = []
    for gene_id, gene_feature in gene.items():
        lines.append(build_gff3_line(gene_feature))
        for mRNA_feature in mRNA[gene_id]:
            lines.append(build_gff3_line(mRNA_feature))
            for other_feature in other[mRNA_feature['attributes']['ID']]:
                lines.append(build_gff3_line(other_feature))
    return lines


def reformat_gff(gff_filename, output_filename):
    gff_dict = read_gff3(gff_filename)
    add_ids_to_features(gff_dict['other'])
    gff_lines = build_gff3_lines(gff_dict)
    with open(output_filename, 'w') as f:
        for line in gff_lines:
            f.write(line + '\n')


# convert gff3 to gtf
def split(input_string, delim='\t'):
    """Split the string by a delimiter and return a list of strings."""
    return input_string.split(delim)

def merge(items, delim='\t'):
    """Merge a list of strings into a single string, separated by a delimiter."""
    return delim.join(items)

def extract(items, key):
    """Extract a value for a given key from a list of key-value pairs."""
    for item in items:
        if item.startswith(key):
            # the key is ID=, Parent=, etc. so we need to remove the key from the value
            return item[len(key):]
    return ''

def preprocess(item, gene_ids, feature_maps, feature_infos):
    """Process each item and populate gene_ids, feature_maps, and feature_infos."""
    feature = item[2]
    attr = split(item[8], ';')
    feature_id = extract(attr, "ID=")
    feature_infos[feature_id] = item
    if feature == "gene":
        gene_ids.append(feature_id)
    elif feature in ["mRNA", "exon"]:
        parent_id = extract(attr, "Parent=")
        if ',' in parent_id:
            parent_ids = split(parent_id, ',')
        else:
            parent_ids = [parent_id]
        for pid in parent_ids:
            if pid not in feature_maps:
                feature_maps[pid] = []
            feature_maps[pid].append(feature_id)

def generate_output(gene_id, feature_maps, feature_infos):
    """Generate output for each gene_id."""
    results = []
    gene_info = feature_infos[gene_id]
    results.append(merge(modify(gene_info, "", ""), '\t'))
    if gene_id in feature_maps:
        for transcript_id in feature_maps[gene_id]:
            tx_info = feature_infos[transcript_id]
            results.append(merge(modify(tx_info, gene_id, ""), '\t'))
            if transcript_id in feature_maps:
                for exon_id in feature_maps[transcript_id]:
                    exon_info = feature_infos[exon_id]
                    results.append(merge(modify(exon_info, gene_id, transcript_id), '\t'))
    return merge(results, '\n')

def modify(item, gene_id, transcript_id):
    """Modify the item based on the feature type."""
    feature = item[2]
    attr = split(item[8], ';')
    feature_id = extract(attr, "ID=")
    if feature == "gene":
        item[8] = f'gene_id "{feature_id}";'
    elif feature == "mRNA":
        item[2] = "transcript"
        item[8] = f'gene_id "{gene_id}"; transcript_id "{feature_id}";'
    elif feature == "exon":
        item[8] = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_id "{feature_id}";'
    else:
        raise ValueError("Unknown feature")
    return item

def convert_gff_to_gtf(gff_filename, gtf_filename):

    gene_ids = []
    feature_maps = {}
    feature_infos = {}

    with open(gff_filename) as gff_file, open(gtf_filename, 'w') as gtf_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            item = split(line.strip(), '\t')
            if item[2] in ["gene", "mRNA", "exon"]:
                preprocess(item, gene_ids, feature_maps, feature_infos)
        
        for gene_id in gene_ids:
            gtf_file.write(generate_output(gene_id, feature_maps, feature_infos) + '\n')

    

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description='GFF Tools')
    subparsers = parser.add_subparsers(help='sub-command help')

    # 创建reformat子命令
    parser_reformat = subparsers.add_parser('reformat', help='Reformat a GFF file')
    parser_reformat.add_argument('gff_filename', type=str, help='The GFF filename to reformat')
    parser_reformat.add_argument('output_filename', type=str, help='The output filename')
    parser_reformat.set_defaults(func=reformat_gff)

    # 创建convert子命令
    parser_convert = subparsers.add_parser('convert', help='Convert a GFF file to GTF format')
    parser_convert.add_argument('gff_filename', type=str, help='The GFF filename to convert')
    parser_convert.add_argument('gtf_filename', type=str, help='The GTF filename')
    parser_convert.set_defaults(func=convert_gff_to_gtf)

    # 解析命令行参数
    args = parser.parse_args()

    # 根据提供的子命令调用相应的函数
    if 'func' in args:
        args.func(**vars(args))
    else:
        # 如果没有提供子命令，打印帮助信息
        parser.print_help()
