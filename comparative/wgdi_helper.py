
import argparse
import os
import sys
import subprocess

#用于根据一组输出, 生成WGDI的配置文件 


class BaseConf:

    def __init__(self, blast_file, gff1, gff2, lens1, lens2):
        self.blast = blast_file
        self.gff1 = gff1
        self.gff2 = gff2
        self.lens1 = lens1
        self.lens2 = lens2

    def write_ctl(self, file_name):
        
        with open(file_name, 'a') as handler:
            for k,v in self.__dict__.items():
                if k == 'header':
                    line = f'[{v}]\n'
                else:
                    line = f'{k} = {v} \n'
                handler.write(line)
    
    def __repr__(self) -> str:
        lines = ""
        for k,v in self.__dict__.items():
            lines += f'{k} = {v} \n'
        
        return lines

class DotPlot(BaseConf):

    def __init__(self, blast_file, gff1, gff2, lens1, lens2, genome1, genome2, prefix, suffix = "svg", **kwargs):

        self.header ='dotplot'

        BaseConf.__init__(self,blast_file, gff1, gff2, lens1, lens2)

        self.genome1_name = genome1
        self.genome2_name = genome2

        self.multiple  = 1
        self.score = 100
        self.evalue = "1e-5"
        self.repeat_number = 10
        self.position = 'order'
        self.blast_reverse = 'false'
        self.ancestor_left =  'none'
        self.ancestor_top =  'none'
        self.markersize = 0.5
        self.figsize = "10,10"

        self.savefig = prefix + "_dotplot." + suffix

        for k,v in kwargs.items():
            if k in self.__dict__:
                self.__dict__[k] = v
            # TO DO: warnning

class IclConf(BaseConf):
    def __init__(self, blast_file, gff1, gff2, lens1, lens2, prefix, **kwargs):

        self.header = 'collinearity'
        BaseConf.__init__(self,blast_file, gff1, gff2, lens1, lens2)

        self.blast_reverse = "false"
        self.multiple  = 1
        self.process = 8
        self.evalue = "1e-5"
        self.score = 100
        self.grading = "50,40,25"
        self.mg = "25,25"
        self.pvalue = 1
        self.repeat_number = 10
        self.positon = "order"
        self.savefile = prefix + "_collinearity.txt"

        for k,v in kwargs.items():
            if k in self.__dict__:
                self.__dict__[k] = v
            # TO DO: warnning


class KsConf(BaseConf):

    def __init__(self, cds_file, pep_file, icl_file, prefix) -> None:
        self.header = 'ks'
        self.cds_file =  cds_file
        self.pep_file =  pep_file
        self.align_software = 'mafft'
        self.pairs_file = icl_file
        self.ks_file = prefix + ".ks"


class BlockinfoConf(BaseConf):

    def __init__(self, blast_file, gff1, gff2, lens1, lens2, icl_file, ks_file, prefix,**kwargs  ):
        self.header = 'blockinfo'
        BaseConf.__init__(self,blast_file, gff1, gff2, lens1, lens2)
        
        self.collinearity = icl_file
        self.score = 100
        self.evalue = "1e-5"
        self.repeat_number = 20
        self.position = "order"
        self.ks = ks_file
        self.ks_col = "ks_NG86"
        self.savefile = prefix + "_blockinfo.csv"   

class KsPeakConf(BaseConf):

    def __init__(self, blockinfo, block_length, prefix, format="png",**kwargs):

        self.header = "kspeaks"
        self.blockinfo = blockinfo
        self.pvalue = 0.2
        self.tandem = "true"
        self.block_length = block_length
        self.ks_area = "0,10"
        self.multiple  = 1
        self.homo = "0,1"
        self.fontsize = 9
        self.area = "0,3"
        self.figsize = "10,6.18"
        self.savefig = prefix + "_ks_peak." + format
        self.savefile = prefix + "_ks_peak.txt"

class BlockKs(BaseConf):

    def __init__(self, lens1, lens2, genome1, genome2, blockinfo, block_length, prefix, format = "png", **kwargs):

        self.header = 'blockks'
        self.lens1 = lens1
        self.lens2 = lens2
        self.genome1_name =  genome1
        self.genome2_name =  genome2
        self.blockinfo = blockinfo
        self.pvalue = 0.2
        self.tandem = "true"
        self.tandem_length = 200
        self.marksersize = 1
        self.area = "0,2"
        self.block_length = block_length
        self.figsize = "8,8"
        self.savefig = prefix + "_blockks." + format

class PeaksFit(BaseConf):

    def __init__(self, blockinfo, prefix, format="png"):

        self.header = "peaksfit"
        self.blockinfo = blockinfo
        self.mode = "median"
        self.bins_number = 200
        self.ks_area = "0,10"
        self.fontsize = 9
        self.area = "0,3"
        self.figsize = "10,6.18"
        self.shadow = "true"
        self.savefig = prefix + "peaks_fit." + format


class KsFigure(BaseConf):

    def __init__(self, ksfit, prefix, format="png"):

        self.header = 'ksfigure'
        self.ksfit = ksfit
        self.labelfontsize = 15
        self.legendfontsize = 15
        self.xlabel = "none"
        self.ylabel = "none"
        self.title = "none"
        self.area = "0,2"
        self.figsize = "10,6.18"
        self.shadow = "true" #(true/false)
        self.savefig =  prefix + "_ks_figure." + format

def prepare_conf(args):

    gff1,gff2,lens1,lens2 = args.gff1, args.gff2, args.lens1,args.lens2
    genome1,genome2 = args.genome_name1, args.genome_name2

    if len(genome1) == 0:
        genome1 = gff1.split('.')[0]
    if len(genome2) == 0:
        genome2 = gff2.split('.')[0]
    
    
    cds_file,pep_file = args.cds_file, args.pep_file
    block_length = args.block_length
    prefix = args.prefix
    blast_file, conf_file = args.blast_file, args.conf_file

    icl_file = prefix + "_collinearity.txt"
    ks_file  = prefix + ".ks"
    blockinfo = prefix + "_blockinfo.csv"
    ksfit = prefix + '_ks_peak.txt'

    if os.path.exists(conf_file):
        os.unlink(conf_file)

    DotPlot(blast_file, gff1, gff2, lens1, lens2, genome1, genome2, prefix).write_ctl(conf_file)
    IclConf(blast_file, gff1, gff2, lens1, lens2, prefix).write_ctl(conf_file)

    KsConf(cds_file, pep_file, icl_file, prefix).write_ctl(conf_file)
    BlockinfoConf(blast_file, gff1, gff2, lens1, lens2, icl_file, ks_file, prefix).write_ctl(conf_file)

    BlockKs(lens1, lens2, genome1, genome2, blockinfo, block_length, prefix).write_ctl(conf_file)

    KsPeakConf(blockinfo, block_length, prefix).write_ctl(conf_file)
    PeaksFit(blockinfo, prefix).write_ctl(conf_file)
    KsFigure(ksfit, prefix).write_ctl(conf_file)

def run(args):

    conf = args.conf_file

    print('Running dotplot', file = sys.stderr)
    subprocess.run(['wgdi', '-d', conf])
    print('Running collinearity', file = sys.stderr)
    subprocess.run(['wgdi', '-icl', conf])
    print('Running ks ', file = sys.stderr)
    subprocess.run(['wgdi', '-ks', conf])
    print('Running blokcinfo', file = sys.stderr)
    subprocess.run(['wgdi', '-bi', conf])
    print('Running ks dotplot', file = sys.stderr)
    subprocess.run(['wgdi', '-bk', conf])
    

def arg_parser():
    parser = argparse.ArgumentParser(description="WGDI helper")

    subparser = parser.add_subparsers(help="sub command help")

    conf_parser = subparser.add_parser('conf', help = 'generate conf')
    
    conf_parser.add_argument('--gff1', required=True)
    conf_parser.add_argument('--gff2', required=True)
    conf_parser.add_argument('--lens1', required=True)
    conf_parser.add_argument('--lens2', required=True)
    conf_parser.add_argument('--genome_name1', required=False, default="")
    conf_parser.add_argument('--genome_name2', required=False, default="")
    conf_parser.add_argument('--cds_file')
    conf_parser.add_argument('--pep_file')
    conf_parser.add_argument('--block_length', default='5')
    conf_parser.add_argument('--prefix', default="output")
    conf_parser.add_argument('blast_file')
    conf_parser.add_argument('conf_file')

    conf_parser.set_defaults(func = prepare_conf)

    run_parser = subparser.add_parser('run', help = 'run all WGDI')
    run_parser.add_argument('conf_file')
    run_parser.set_defaults(func = run)

    return parser


if __name__ == "__main__":
    
    args = arg_parser().parse_args()
    args.func(args)    
    # try:
    #     args.func(args)
    # except:
    #     arg_parser().print_help()
    
    
    






