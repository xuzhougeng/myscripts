

FASTA=$1
DB=$2
target=Viridiplantae
THREADS=100

blastn -query $FASTA -db $DB -num_alignments 10 -outfmt 6 -num_threads 100 > round1.blastn.txt

# find best match
python -m jcvi.formats.blast best round1.blastn.txt

# find the best-matched accession
cut -f 2 round1.blastn.txt.best | sort -u > accesions.txt

# 根据accession分析它的lineage
blastdbcmd -db $DB -dbtype nucl -entry_batch accesions.txt -outfmt '%a %T' > accesion2taxid.txt

# 根据accesion分析lineage
cut -d ' ' -f 2 accesion2taxid.txt | sort -u | taxonkit lineage > round1_lineage.txt

# 挑选绿色植物门
grep $target round1_lineage.txt | cut -f 1  > filter_lineage.txt

# 根据taxonomy挑选accession
grep -wf filter_lineage.txt accesion2taxid.txt | cut -d ' ' -f 1 > filter_accesion.txt

# 根据accesion挑选对应的query
grep -wf filter_accesion.txt round1.blastn.txt.best  | cut -f 1 > filter_seq_id.txt
