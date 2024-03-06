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

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <gff_filename> <gtf_filename>")
        sys.exit(1)

    gff_filename = sys.argv[1]
    gtf_filename = sys.argv[2]

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
