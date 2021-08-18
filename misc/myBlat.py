from flask import Flask
from flask import request
import tempfile
import subprocess
from os.path import basename

app = Flask(__name__)

ref_dict = {
    "Athaliana":  "/data/reference/blat/Athaliana.2bit",
    "Alyrata": "/data/reference/blat/Alyrata_384_v1.2bit"
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
    if "thaliana" in db:
        return ref_dict["Athaliana"]
    elif "lyrata" in db:
        return ref_dict["Alyrata"]
    else:
        return ref_dict["Athaliana"]
   
    

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

    
