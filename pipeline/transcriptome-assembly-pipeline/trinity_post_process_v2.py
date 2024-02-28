import os
os.environ["NUM_THREADS"] = "1"
import pandas as pd
from Bio import SeqIO
from enum import Enum
import math
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import gridspec
import re
from typing import Any, List
import pysam
import joblib


# read transdecoder output
class ORF:
    """
    Store the ORF (Open Reading Frame) predicted by TransDecoder.

    Attributes:
        transcript_name (str): Name of the transcript.
        transcript_length (int): Length of the transcript.
        orf_id (str): Identifier for the ORF.
        orf_type (str): Type of the ORF.
        orf_score (float): Score of the ORF prediction.
        strand (str): Strand information, '+' or '-'.
        orf_start (int): Start position of the ORF.
        orf_end (int): End position of the ORF.
        lr_support (bool): Indicator of long read support, default False.
    """
    
    def __init__(self, transcript_name:str, transcript_length:str, 
                orf_id:str, orf_type:str, orf_score:float, strand:str,
                orf_start:int, orf_end:int):
        self.transcript_name = transcript_name
        self.transcript_length = transcript_length
        self.orf_id = orf_id
        self.orf_type = orf_type
        self.orf_score = orf_score
        self.strand = strand
        self.orf_start = orf_start
        self.orf_end = orf_end
        #self.lr_support = lr_support


    def __repr__(self):
        return f"{self.transcript_name}\t{self.transcript_length}\t{self.orf_id}\t{self.orf_type}\t{self.orf_score}\t{self.strand}\t{self.orf_start}\t{self.orf_end}\t{self.lr_support}"
    
    __str__ = __repr__

    @classmethod
    def parse_bed_line(cls, line):
        """
        Parse a BED line and return an ORF object
        """
        fields = line.strip().split('\t')
        
        # Extract information from fields
        transcript_name = fields[0]
        transcript_start = int(fields[1])
        transcript_end = int(fields[2])
        transcript_length = transcript_end - transcript_start + 1
        info = fields[3]
        strand = fields[5]
        orf_start = int(fields[6])
        orf_end = int(fields[7])
        
        # Parse the info field to extract ORF ID and ORF type
        info_parts = info.split(';')
        orf_id = info_parts[0].split('=')[1]
        orf_type = info_parts[2].split(',')[0].split(':')[1].strip()
        orf_score = float( info_parts[2].split(',')[1].split('=')[1].strip())

        return cls(transcript_name, transcript_length, orf_id, orf_type, orf_score, strand, orf_start, orf_end)

class IsoformRank(Enum):
    # isoform 分级
    single_complete = 1
    multi_complete = 2
    multi_complete_partial = 3
    multi_partial = 4
    single_partial = 5

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

# 首先定义一个函数来找到不与当前区间重叠的前一个区间的索引
def find_previous(interval, intervals):
    for i in range(len(intervals) - 1, -1, -1):
        if intervals[i][1] <= interval[0]:
            return i
    return -1

# 动态规划函数来计算最大得分
def max_score(intervals):
    if len(intervals) == 0:
        return 0
    
    if len(intervals) == 1:
        return intervals[0][2]

    # 根据结束时间对区间进行排序
    intervals.sort(key=lambda x: x[1])

    # 创建一个表来存储到每个点为止的最大得分
    dp = [0] * (len(intervals) + 1)

    # 初始化第一个区间的得分
    dp[1] = intervals[0][2]

    # 动态规划填表过程
    for i in range(2, len(intervals) + 1):
        # 找到不与当前区间重叠的前一个区间的索引
        prev_index = find_previous(intervals[i - 1], intervals)
        # 选择不包括当前区间的得分或者包括当前区间及其最大兼容区间得分的最大值
        dp[i] = max(dp[i - 1], intervals[i - 1][2] + (dp[prev_index + 1] if prev_index != -1 else 0))

    return dp[-1]

class Isoform:
    """
    isoform with multiple orfs
    """
    
    def __init__(self, name:str=None, orfs: List[ORF]=None):
        """
        Args:
            name: isoform name
            orfs: orf list
        
            
        """
        self.name = name
        self.orfs = orfs
        self.size = len(orfs)
        self.rank  = self._rank()
        self.score = self._calc_score()

    def _rank(self) -> IsoformRank:
            
        if len(self.orfs) == 1:
            if "complete" in self.orfs[0].orf_type:
                return IsoformRank.single_complete
            else:
                return IsoformRank.single_partial
        else:
            complete_count = 0
            partial_count = 0
            for orf in self.orfs:
                if "complete" in orf.orf_type:
                    complete_count += 1
                else:
                    partial_count += 1
                
            if complete_count == 0 and partial_count > 0:
                return IsoformRank.multi_partial
            
            if complete_count > 0 and partial_count > 0:
                return IsoformRank.multi_complete_partial

            if complete_count > 0 and partial_count == 0:
                return IsoformRank.multi_complete    
    
    def _calc_score(self):
        """
        """
        orf_plus_region = []
        orf_minus_region = []

        if len(self.orfs) == 1:
            return self.orfs[0].orf_score
        
        # 获取每条链上的最高得分
        for orf in self.orfs:
            if orf.strand == "+":
                orf_plus_region.append([ int(orf.orf_start), int(orf.orf_end), orf.orf_score ])
            else:
                orf_minus_region.append([ int(orf.orf_start), int(orf.orf_end), orf.orf_score ])
        
        # 计算正链上的最大得分
        plus_max_score = max_score(orf_plus_region)
        # 计算负链上的最大得分
        minus_max_score = max_score(orf_minus_region)

        # 计算两条链上的最大得分之和
        return plus_max_score + minus_max_score

    def add(self, orf: ORF):
        self.orfs.append(orf)
        self.size += 1

    def get_orf(self, index:int):
        if index >= self.size:
            raise IndexError(f"index out of range, size: {self.size}, index: {index}")
        return self.orfs[index]
    
    def delete_orf(self, index:int):
        del self.orfs[index]
        self.size -= 1

    def __repr__(self):
        return f"{self.name}: rank={self.rank} score={self.score} size={self.size}"

class Node:
    def __init__(self, name, size=None) -> None:
        self.name = name
        self.size = size
    
    def __repr__(self) -> str:
        return f"{self.name}:{self.size}" if self.size is not None else f"{self.name}"

    # __str__ is the same as __repr__ in this case
    __str__ = __repr__

class OrderedNodes:
    def __init__(self, nodes: List, sizes: List=None) -> None:
        self.nodes: List[Node] = []
        for index, name in enumerate(nodes):
            size = None if sizes is None else sizes[index]
            self.add_node(Node(name=name, size=size))  # No need to pass size separately

    def add_node(self, node: Node):  # Changed to accept only a Node object
        self.nodes.append(node)

    def get_node(self, index):
        return self.nodes[index]

    def get_nodes(self):
        return self.nodes
    
    def get_node_size(self):
        return len(self.nodes)

    def __repr__(self) -> str:
        return f"{self.nodes}"
    
    __str__ = __repr__

    @classmethod
    def path_parser(cls, path):
        nodes = []
        node_sizes = []
        for part in path.split(" "):
            node_id, range_part = part.split(":")
            start, end = map(int, range_part.split("-"))
            node_size = end - start + 1
            
            nodes.append(int(node_id))
            node_sizes.append(node_size)
        
        return cls(nodes, node_sizes)
    
    @staticmethod
    def create_node_from_string(description):
        match = re.search(r'path=\[(.*?)\]', description)
        if match:
            result = match.group(1)
            # Assuming path_parser can be called as is, but ideally should be OrderedNodes.path_parser
            return OrderedNodes.path_parser(result)
        else:
            return None

def plot_gene_nodes(gene_nodes:List[OrderedNodes], highlight_start=True, highlight_common=False, subplot=False, cols=3):

    FIG_WIDTH = 12
    FIG_HEIGHT = 4

    nodes_num = len(gene_nodes)

    # 获取共同路径, 起始节点, 节点大小
    common_nodes = set()
    start_nodes = []
    node_sizes = {}
    
    for gene_node in gene_nodes:
        start_nodes.append(gene_node.get_node(0).name)
        for i in range(gene_node.get_node_size()):
            node = gene_node.get_node(i)
            if node.name in node_sizes:
                common_nodes.add(node.name)
            node_sizes[node.name] = node.size


    # 创建主图和子图的布局
    rows = math.ceil(nodes_num / cols) if subplot else FIG_HEIGHT
    cols = cols if subplot else 1
    fig = plt.figure(figsize=(FIG_WIDTH, max(2, rows * 2)+1))
    gs = gridspec.GridSpec(rows, cols + 1, width_ratios=[3] + [1]*cols)
    main_ax = fig.add_subplot(gs[:, 0])    
    G = nx.DiGraph()


    # 计算并重用节点位置
    G.add_nodes_from(node_sizes.keys())
    for gene_node in gene_nodes:
        nodes_list  = gene_node.get_nodes()
        if len(nodes_list ) > 1:
            edges = [(nodes_list[i].name, nodes_list[i + 1].name) for i in range(len(nodes_list) - 1)]
            G.add_edges_from(edges)

    pos = nx.spring_layout(G)  # 仅计算一次布局

    # 绘制主图
    nx.draw(G, pos, ax=main_ax, with_labels=True, node_size=500, node_color='skyblue')
    
    #nx.draw_networkx_labels(G, pos, labels=node_sizes, ax=main_ax, font_size=8)


    if highlight_common:
        nx.draw_networkx_nodes(G, pos, ax=main_ax, nodelist=common_nodes, node_size=600, node_color='yellow')
    if highlight_start:
        nx.draw_networkx_nodes(G, pos, ax=main_ax, nodelist=start_nodes, node_size=600, node_color='green')
    
    main_ax.text(0.5, -0.1, f"Gene ID: {gene_id}", transform=main_ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='center')

    # 绘制子图
    if subplot:
        for i, gene_node in enumerate(gene_nodes):
            row, col = divmod(i, cols)
            ax = fig.add_subplot(gs[row, col + 1], aspect='equal')
            G_seq = nx.DiGraph()
            nodes_list = gene_node.get_nodes()
            G_seq.add_nodes_from([node.name for node in nodes_list])
            if len(nodes_list) > 1:
                edges = [(nodes_list[i].name, nodes_list[i + 1].name) for i in range(len(nodes_list) - 1)]
                G_seq.add_edges_from(edges)

            nx.draw(G_seq, pos, ax=ax, with_labels=True, node_size=300, node_color='lightgreen')
            ax.set_title(f"{'->'.join(map(str, nodes_list))}", loc="left")

            ax.set_xlim(main_ax.get_xlim())
            ax.set_ylim(main_ax.get_ylim())
            ax.axis('off')
    
    # 如果不是子图模式，添加路径信息框
    if not subplot:
        info_ax = fig.add_subplot(gs[0, -1])
        info_ax.axis('off')
        #paths_str = ['->'.join(map(str, gene_node.get_nodes())) for gene_node in gene_nodes]
        paths_str = ['->'.join(str(node.name) for node in gene_node.get_nodes()) for gene_node in gene_nodes ]
        paths_text = "\n".join(paths_str)
        info_ax.text(0.5, 0.5, paths_text, fontsize=9, verticalalignment='center', horizontalalignment='left')
        
        # 添加节点大小信息框
        info_ax = fig.add_subplot(gs[1, -1])
        info_ax.axis('off')
        node_sizes_str = [f"{name}:{size}" for name,size in node_sizes.items() ]
        node_sizes_text = "\n".join(node_sizes_str)
        info_ax.text(0.5, 0.5, node_sizes_text, fontsize=9, verticalalignment='top' , horizontalalignment='left')

    #plt.tight_layout()
    plt.show()

class Gene:

    def __init__(self, name, isoforms:List[Isoform]=None ) -> None:
        self.name = name # gene name
        self.isoforms = isoforms # isoform list
        self.size = len(isoforms) # isoform number

        # best isoform of gene
        self.best_isoform = None

    
    def add(self, isoform: Isoform):
        self.isoforms.append(isoform)
        self.size += 1

    def get_isoform(self, index:int):
        if index >= self.size:
            raise IndexError(f"index out of range, size: {self.size}, index: {index}")
        return self.isoforms[index]
    
    def set_best_isoform(self, isoform:Isoform):
        self.best_isoform = isoform

    def __repr__(self) -> str:
        return f"{self.name}: size={self.size}"

def read_fasta(fasta_file):
    fasta_dict = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[seq_record.id] = seq_record
    return fasta_dict

def read_salmon_sf(salmon_output_sf):
    # read salmon output and merged them 
    expr_df = pd.DataFrame()

    for sf in salmon_output_sf:
        sample_name = sf.split("/")[-1].split("_")[0]
        tmp = pd.read_csv(sf, sep="\t")
        # add sample name column
        tmp["sample"] = sample_name
        expr_df = pd.concat([expr_df, tmp])


    pivot_df = expr_df.pivot_table(index='Name', columns='sample', values='TPM', aggfunc='max')
    # 使用idxmax获取TPM值最大的sample
    pivot_df['max_tpm'] = pivot_df.max(axis=1)
    tpm_dict = pivot_df["max_tpm"].to_dict()
    
    return tpm_dict

def filter_orf(tx_orf_dict, tpm_dict, min_tpm=1.0, min_orf_length=100, min_orf_coverage=0.5):
    # 筛选标准
    # 1. isoform with complete orf(even if tpm < 1.0)
    # 2. the tpm of isoform > 1.0 (even no complete orf)) 

    remove_id = set()
    min_tpm = 1.0

    for tx, predictions in tx_orf_dict.items():
        
        if tpm_dict[tx] > min_tpm:
            continue
        
        flag = True
        for prediction in predictions:
            if "complete" in prediction.orf_type:
                flag = False
                break
        if flag:
            remove_id.add(tx)

    # delete the transcript
    for tx in remove_id:
        if tx in tx_orf_dict:
            del tx_orf_dict[tx]
    
    return tx_orf_dict


def rank_gene(genes: List[Gene]):

    # 对于多个multi_isoform_gene的情况
    best_gene_list = []
    uncomplete_gene_list = []

    RATIO = 1.5
    for gene in genes:
        best_isoform = gene.best_isoform

        threshold = best_isoform.score / RATIO

        # 如果最佳转录本是完整的, 只保留 > score/1.5的转录本
        # 并添加到best_gene_list,

        if best_isoform.rank == IsoformRank.single_complete:
            best_gene_list.append(gene)
            # 重新构建isoform list
            complete_isoform = []
            for isoform in gene.isoforms:
                if isoform.rank == IsoformRank.single_complete and isoform.score > threshold:
                    complete_isoform.append(isoform)
            gene.isoforms = complete_isoform
        else:
            complete_isoform = []
            for isoform in gene.isoforms:
                if isoform.rank == IsoformRank.single_complete and isoform.score > threshold:
                    complete_isoform.append(isoform)

            # 如果存在完整的转录本, 重新构建isoform list   
            if len(complete_isoform) > 0:
                gene.isoforms = complete_isoform
                best_gene_list.append(gene)
            else:
                uncomplete_gene_list.append(gene)
    return gene, best_gene_list, uncomplete_gene_list


def main(fata_file, salmon_output_sf, bed_file ):

    # read trinity fasta
    fasta_dict = read_fasta(fasta_file)

    # filter base on tpm and orf prediction
    tpm_dict = read_salmon_sf(salmon_output_sf)

    # tx : [orf1, orf2, orf3]
    tx_orf_dict = {}

    with open(bed_file, 'r') as file:
        # Skip the header
        next(file)
        
        # Parse each line
        for line in file:
            orf = ORF.parse_bed_line(line)
            transcript_name = orf.transcript_name
            if transcript_name not in tx_orf_dict:
                tx_orf_dict[transcript_name] = []
            tx_orf_dict[transcript_name].append(orf)

    
    tx_orf_dict = filter_orf(tx_orf_dict, tpm_dict)


    # 解析orf_dict, 得到isoform_dict
    # gene: [tx1, tx2, tx3 ....]
    gene_tx_dict = {}

    for tx, orfs in tx_orf_dict.items():
        gene_name = "_".join(tx.split("_")[:-1])
        isoform = Isoform(tx, orfs)
        if gene_name not in gene_tx_dict:
            gene_tx_dict[gene_name] = []
        gene_tx_dict[gene_name].append(isoform)


    # 解析isoform_dict, 得到gene_dict
    gene_dict = {}
    for gene_name, isoforms in gene_tx_dict.items():
        
        gene_dict[gene_name] = Gene(gene_name, isoforms)
    
    # 将基因按照isoform数目分成两类
    single_isoform_gene = [ v for v in gene_dict.values() if v.size == 1]
    multi_isoform_gene = [ v for v in gene_dict.values() if v.size > 1]

    print(f"single isoform gene: {len(single_isoform_gene)}")
    print(f"multi isoform gene: {len(multi_isoform_gene)}")

    # 对于多个multi_isoform_gene的情况
    good_gene_list = []
    bad_gene_list = []
    for gene in multi_isoform_gene:
        isoforms = gene.isoforms
        # 获取每个转录本的预测结果, 从中选取得分最高的
        highest_isoform = ""
        # 无穷小
        highest_score = float("-inf") 

        for isoform in isoforms:
            if isoform.score > highest_score:
                highest_score = isoform.score 
                highest_isoform = isoform
        
        if highest_isoform.rank == IsoformRank.single_complete:
            good_gene_list.append(gene)
        else:
            bad_gene_list.append(gene)



if __name__ == "__main__":
    fasta_file = "/data5/xzg_data/arabidopsis_single_cell_atalas/rna-seq/metrics/trinity_ss_flnc.fasta"
    # read salmon output 
    salmon_output_sf = [
        "/home/xzg/project2/arabidopsis_single_cell_atalas/rna-seq/metrics/salmon/SRR25073470_quant.sf",
        "/home/xzg/project2/arabidopsis_single_cell_atalas/rna-seq/metrics/salmon/SRR25073471_quant.sf",
        "/home/xzg/project2/arabidopsis_single_cell_atalas/rna-seq/metrics/salmon/SRR25073472_quant.sf"
    ]
    bed_file = "/data5/xzg_data/arabidopsis_single_cell_atalas/rna-seq/metrics/transdecoer.standard.code/trinity_ss_flnc.fasta.transdecoder.bed"
    main(fasta_file, salmon_output_sf, bed_file)