
import io
from logging import warning
import os

import concurrent.futures
import copy

from subprocess import run as _run
from operator import methodcaller


from typing import List, Dict
from pathlib import Path
from collections import defaultdict

from Bio import (
    SeqIO,
    SeqRecord,
    Phylo
)

import argparse  


"""脚本说明
为MCMCTree整理输入文件 和 CODEML分析

输入文件如下

- 物种树, 如OrthoFinder的输出结果, Species_Tree/SpeciesTree_rooted.txt (OrthoFinder输出结果)
- 校准点信息, 来源于文献, 可以在timetree.org中查询
- 单拷贝基因编号文件, 第一行记录物种名, 后续的每一行是对应一个单拷贝基因, 制表符分割
- cds和pep文件


依赖软件

- mafft
- pal2nal
- PAML
"""



def update_tree(tree_file, calibrations = None, verbose=False):
    """更新树文件, 返回物种名
    :param tree_file newick 格式
    :param calibration 校准点信息, 时间单位 100Myr, 格式's1,s2:>0.01;s3,s4:<0.02;s5,s6:>.01 <0.2'
    """

    tree = Phylo.read(tree_file, 'newick')
    # show tree structure
    if verbose:
        Phylo.draw_ascii(tree)
    # remove the branch length
    species_name = []
    for clade in tree.find_clades():
        if clade.name is not None:
            species_name.append(clade.name)
        
        clade.branch_length = None

    # add fossial calibration
    ## format: 'x,y: >.06 <.08'
    if calibrations is not None:
        for calibration in calibrations.split(';'):
            species, time_pt = calibration.split(':')
            s1, s2 = species.split(',')
            clade = tree.common_ancestor(s1, s2)
            clade.name = time_pt

    tree_buf = io.StringIO()
    Phylo.write(tree, tree_buf, 'newick')
    tree_buf.seek(0)

    species_num = len(species_name)

    tree_data =  ['{}\t1\n'.format(species_num) ] + tree_buf.readlines()
    tree_data[1] = tree_data[1].replace(':0.00000','').replace('1.00','')

    return species_name,tree_data

def prepare_seq_db(files: List[Path]):
    """以字典形式记录每个物种的序列, 如cds, pep
    :param files 一组cds, pep的文件路径
    """

    seq_dict = defaultdict(dict)

    for file in files:
        if file is not Path:
            file = Path(file)
        species_name = file.stem
        seq_dict[species_name] = {record.id : record.seq for record in SeqIO.parse(file, 'fasta') }

    return seq_dict    

class Pal2Nal:

    def __init__(self, og_id, pep_file, cds_file):
        """
        :param pep_file aa fasta
        :param cds_file nt fasta
        """
        # input
        self.og_id = og_id
        self.pep_file = pep_file
        self.cds_file = cds_file
        
        # 是否需要将输出结果记录在这里呢？
        self.nal_dict = None

    def run(self):
        """运行mafft and pal2nal.pl
        """
        aln_file = f'tmp/{self.og_id}.aln'
        
        if not os.path.exists(aln_file):
            mafft_cmd = ['mafft', '--localpair', '--maxiterate', '1000', '--anysymbol', self.pep_file  ]
            output = _run(mafft_cmd, capture_output=True)
            
            
            with open(aln_file, 'w') as handler:
                handler.writelines(output.stdout.decode('ascii'))
                
        # -nogap 
        # -nomismatch

        cmd = ['pal2nal.pl', aln_file, self.cds_file, '-output', 'fasta', '-nogap']

        output = _run(cmd, capture_output=True)
        string_io = io.StringIO(output.stdout.decode('ascii'))

        nal_dict = {}
        for record in SeqIO.parse(string_io, 'fasta'):
            species_name = record.id.split('|')[1]
            nal_dict[species_name] = record

        if len(nal_dict) == 0:
            # TO DO: add some warnning message
            pass

        self.nal_dict =  nal_dict
    
    def cleanup(self, delete=True):
        """临时文件
        """
        # delete file
        if delete:
            os.unlink(self.pep_file)
            os.unlink(self.cds_file)
       

    @classmethod
    def from_seq(cls,og_id, pal_seq_list, cds_seq_list):

        if not os.path.exists('tmp'):
            os.mkdir('tmp')
        
        # cds and coresponding protein
        cds_file = os.path.join('tmp', og_id + '.cds')
        pep_file = os.path.join('tmp', og_id + '.pep')

        with open(pep_file, 'w') as handler:
            SeqIO.write(pal_seq_list, handler, 'fasta')

        with open(cds_file, 'w') as handler:
            SeqIO.write(cds_seq_list, handler, 'fasta')

        return  cls(og_id, pep_file, cds_file )
    
def prepare_pal2nal(og_id, gene_dict, pep_dict, cds_dict):
    """构建单个pal2nal
    :param og_id 单拷贝基因家族的编号
    :param gene_dict 单拷贝基因家族对应的基因的编号
    :param pep_dict aa序列
    :param cds_dict cds序列
    """

    pep_seq_list = []
    cds_seq_list = []

    for species,gene_id in gene_dict.items():
        #print(species,gene_id)
        pep_seq = pep_dict[species][gene_id]
        cds_seq = cds_dict[species][gene_id]

        pep_seq_list.append(
            SeqRecord.SeqRecord(id=gene_id + "|" + species, seq = pep_seq, name="", description="")
        )

        cds_seq_list.append(
            SeqRecord.SeqRecord(id=gene_id + "|" + species, seq = cds_seq, name="", description="")
        )
    
    return Pal2Nal.from_seq(og_id, pep_seq_list, cds_seq_list)

def prepare_pal2nal_list(scg_groups, pep_dict, cds_dict):
    """基于编号生成任务列表
    :param scg_groups 一组单拷贝基因的字典, 里面的元素对应各个og
    :param pep_dict aa序列
    :param cds_dict cds序列
    """
    pal2nal_list = []
    for og_id,gene_dict in scg_groups.items():
        # og_id: OG0012728
        pal2nal_list.append(
            prepare_pal2nal(og_id, gene_dict, pep_dict, cds_dict)
        )

    return pal2nal_list

def run_all_pal2nal(scg_groups, pep_dict, cds_dict, n_job=1):
    """运行所有的pal2nal, 返回列表, 
    :param scg_group each group contain dict of single copy gene 
    :param pepe_dict protein dictionary
    :param cds_dict cds dictionary
    :param n_job number of njobs
    """

    all_nal = []
    

    # 任务并发
    if n_job == 1:
        for scg_id,gene_dict in scg_groups.items():
            pal2nal = prepare_pal2nal(scg_id, gene_dict, pep_dict, cds_dict)
            pal2nal.run()
            all_nal.append(  
                copy.deepcopy(pal2nal.nal_dict)
            )
            pal2nal.cleanup(True)

    else:
        tot_jobs = len(scg_groups)
        scg_ids = list(scg_groups.keys())
        for i in range(0, tot_jobs, n_job):
            idx_start = i
            idx_end = i + n_job if (i + n_job) < tot_jobs else tot_jobs

            cur_scg_ids = scg_ids[idx_start:idx_end]
            
            sub_scg_groups = { og_id : scg_groups[og_id] for og_id in cur_scg_ids }

            pal2nal_list = prepare_pal2nal_list(sub_scg_groups,pep_dict,cds_dict)
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_job) as executor:
                executor.map(methodcaller('run'), pal2nal_list)

            # 汇总结果
            for pal2nal in pal2nal_list:
                all_nal.append(copy.deepcopy(pal2nal.nal_dict))
            
            # 清理临时文件
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_job) as executor:
                executor.map(methodcaller('cleanup'), pal2nal_list)

    return all_nal
  
class McmcTreeCtl:

    def __init__(self, seqfile, treefile,**kwargs):

        # initialize
        self.seed = -1
        self.seqfile = seqfile  # seqfile in phylip format
        self.treefile = treefile # treefile with fossial calibration
        self.outfile = "mcmc.txt"  # 

        self.ndata = 3        # partition number
        self.seqtype = 0      # 0: nucleotides; 1:codons; 2:AAs
        self.usedata = 1      # 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        self.clock = 1        # 1: global clock; 2: independent rates; 3: correlated rates
        self.RootAge = '<1.0' # safe constraint on root age, used if no fossil for root.

        self.model = 0        # 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        self.alpha = 0        # alpha for gamma rates at sites
        self.ncatG = 5        # No. categories in discrete gamma

        self.cleandata = 0    # remove sites with ambiguity data (1:yes, 0:no)?

        self.BDparas = '1 1 0'        # birth, death, sampling
        self.kappa_gamma = '6 2'      # gamma prior for kappa
        self.alpha_gamma = '1 1'      # gamma prior for alpha

        self.rgene_gamma = '2 2'        # gamma prior for overall rates for genes
        self.sigma2_gamma = '1 10'      # gamma prior for sigma^2     (for clock=2 or 3)

        self.finetune = '1: 0.1  0.1  0.1  0.01 .5'  # auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

        self.print = 1
        self.burnin = 8000
        self.sampfreq = 2
        self.nsample = 200000

        # Resetting parameters
        for k,v in kwargs.items():
            if k in self.__dict__:
                self.__dict__[k] = v
            # TO DO: warnning

    def __repr__(self) -> str:
        lines = ""
        for k,v in self.__dict__.items():
            lines += f'{k} = {v} \n'
        
        return lines
    
    def write_ctl(self, file_name):
        
        with open(file_name, 'w') as handler:
            for k,v in self.__dict__.items():
                line = f'{k} = {v} \n'
                handler.write(line)

def export_all_nal(scg_ids, all_nal, outdir="SCG" ,format='fasta'):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pair_iter = zip(scg_ids, all_nal)

    for og,nal_dict in pair_iter:
        nal_list = []
        for _,v in nal_dict.items():
            gene_id,species_name = v.name.split('|')
            v.name = species_name + '_' +  gene_id
            v.id = species_name + '_' +  gene_id
            nal_list.append(v)

        # some dictionary may empty due to inconsistent  
        if len(nal_list) == 0:
            continue
        file_name = os.path.join(outdir, og + ".fas" )
        with open(file_name, 'w') as handler:
            SeqIO.write(nal_list, handler, format)

def export_phylip(nal_dict, file_name):
    """保存为phylip格式
    :param nal_dict 记录Seq的字典
    """
    species_num = len(nal_dict)
    # tricky: using iter to obtain iterator  of dict_values([...])
    nuc_num = len(next(iter(nal_dict.values()))) // 3
    header = '   {} {}\n'.format(species_num, nuc_num)
    for i in range(3):
        with open(file_name, 'a') as f:
            f.write(header)
            
        with open(file_name, 'a') as f:

            for key, value in nal_dict.items():
                f.write('{}\n{}\n'.format(key, value[i::3]  ) )


def main(args):
    """workflow:
    1. load pep and cds dict
    2. load single copy gene orthogroups
    3. mafft align the protein and translate the pep to cds
    """
    
    # parse argument
    cds_dir = args.cds_dir
    cds_suffix = args.cds_suffix
    pep_dir = args.pep_dir
    pep_suffix = args.pep_suffix
    scg_file = args.scg_file
    tree_file =args.tree_file
    calibrations = args.calibrations


    # load all cds
    import glob
    import pickle

    if not os.path.exists('all_nal_bak.pickle'):
        # load the pep and cds dict
        cds_files = glob.glob( os.path.join(cds_dir , "*." + cds_suffix ) )
        cds_dict = prepare_seq_db( cds_files )
        
        pep_files = glob.glob( os.path.join(pep_dir , "*." + pep_suffix ) )
        pep_dict = prepare_seq_db( pep_files )

        # load the single copy gene family
        scg_file_handler = open(scg_file, 'r')
        species_list = next(scg_file_handler).strip().split('\t')[1:]

        scg_groups = {}
        for line in scg_file_handler:
            fields = line.strip().split('\t')
            og_id = fields[0]
            gene_list = fields[1:]
            gene_dict = dict(zip(species_list, gene_list))

            scg_groups[og_id] = gene_dict

        # scg_ids all_nal 是配对

        with open('scg.txt', 'w') as f:
            for scg in scg_groups.keys():
                f.write(scg + '\n')

        # 输出记录各个species codon的字典
        # concatenate to a single object
        # have efficiency improvement potential
        all_nal = run_all_pal2nal(scg_groups, pep_dict, cds_dict, args.threads)
        #print(all_nal)
    
        with open('all_nal_bak.pickle', 'wb') as f:
            pickle.dump(all_nal, f)
    else:
        with open('all_nal_bak.pickle', 'rb') as f:
            all_nal = pickle.load(f)

    # export codon alignment alone
    scg_ids = [x.strip() for x in open('scg.txt')]

    export_all_nal(scg_ids, all_nal)

    # export to phylip format       
    nal_dict = {}
    for nal in all_nal:
        for k,v in nal.items():
            if k in nal_dict:
                nal_dict[k] += v.seq
            else:
                nal_dict[k] = v.seq

    export_phylip(nal_dict, 'mcmctree.phylip')

    # 整合treefile和校准点信息 
    species_name,tree_data = update_tree(tree_file, calibrations)
    print(species_name)

    with open('mcmctree.txt', 'w') as f:
        f.writelines(tree_data)



def test_main():
    import glob
    cds_dir = "."
    cds_suffix = "cds"
    pep_dir = "."
    pep_suffix = "faa"
    threads = 20

    cds_files = glob.glob( os.path.join(cds_dir , "*." + cds_suffix ) )
    cds_dict = prepare_seq_db( cds_files )
    
    pep_files = glob.glob( os.path.join(pep_dir , "*." + pep_suffix ) )
    pep_dict = prepare_seq_db( pep_files )

    scg_file_handler = open("SingleCopyOrthologues.txt", 'r')
    species_list = next(scg_file_handler).strip().split('\t')[1:]

    scg_groups = {}
    for line in scg_file_handler:
        fields = line.strip().split('\t')
        og_id = fields[0]
        gene_list = fields[1:]
        gene_dict = dict(zip(species_list, gene_list))

        scg_groups[og_id] = gene_dict

    all_nal = run_all_pal2nal(scg_groups, pep_dict, cds_dict, threads)


def arg_parser():

    parser = argparse.ArgumentParser(description="Orthofinder result to MCMC")
    
    parser.add_argument('--threads', default=20, type=int)
    parser.add_argument('--cds_dir', required=True) # directory for cds file
    parser.add_argument('--cds_suffix', default = 'cds') # suffix for pattern matching

    parser.add_argument('--pep_dir', required=True) # directory for pep file
    parser.add_argument('--pep_suffix', default = 'faa') # suffix for pattern matching

    parser.add_argument('--calibrations')
    parser.add_argument('scg_file', help = 'eg. Orthogroups/Orthogroups_SingleCopyOrthologues.txt') # 
    parser.add_argument('tree_file', help = 'eg. Species_Tree/SpeciesTree_rooted.txt')


    return parser


if __name__ == "__main__":
    args = arg_parser().parse_args()

    main(args)

