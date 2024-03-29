from flask import Flask
from flask import request
import tempfile
import subprocess
from os.path import basename

app = Flask(__name__)


"""
注意:
    1. 修改ref_dict中的ref_name1和ref_name2fa的位置
	2. 修改query_seq_db中的 ref1 和 ref2, 改成你的物种里的特有关键词, 例如A.thaliana 我用的是thaliana

"""


ref_dict = {
    "ref_name1":  "/path/to/your/reference1.fa",
    "ref_name2":  "/path/to/your/reference2.fa"
}

def query_seq_db(db):
    #if "\\" in db:
    #    dbname = db.split("\\")[-1]
    #else:
    #    dbname = basename(db)
    #dbname = dbname.split(".")[0]
    #ref = ref_dict[dbname]
    #return ref
    db = db.lower()
    if "ref1" in db:
        return ref_dict["ref_name1"]
    elif "ref2" in db:
        return ref_dict["ref_name2"]
    else:
        return ref_dict["ref_name1"]
   
    

def blat(seq,db,type):
    # 需要改成你实际的blat路径
    blat = "/opt/biosoft/ucsctools/blat/blat"
    # 需要改成你实际的参考位置
    ref = query_seq_db(db)

    fd1, fa_file = tempfile.mkstemp(suffix=".fa")
    fd2, psl_file = tempfile.mkstemp(suffix=".psl")

    with open(fa_file,"w") as f:
        f.write(f">query\n{seq}\n")
    process = subprocess.Popen([blat, ref, fa_file, psl_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    
    stdout,stderr = process.communicate()

    psl_ret = [line for line in open(psl_file, "r") ]
    psl_rec = [line.strip().split('\t') for line in psl_ret[5:] ]

    return psl_rec

@app.route("/myBlat", methods=['GET'])
def myBlat():

    seq = request.args.get('userSeq','')
    db = request.args.get('db','')
    output = request.args.get('output','json')
    type = request.args.get('type', 'DNA')

    psl_rec = blat(seq, db, type)
    data =  {
        "track": "blat",
        "genome": "ath",
        "fields": ["matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", \
                "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", \
                "tSize", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"],
        "blat" : psl_rec
        }
    return data

    
