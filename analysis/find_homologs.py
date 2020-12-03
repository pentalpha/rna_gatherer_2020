import sys
import os
import subprocess
import pandas as pd
import numpy as np

query = sys.argv[1]
db = sys.argv[2]
outdir = sys.argv[3]
names_dump_path = sys.argv[4]
gff_annotation = sys.argv[5]
threads = int(sys.argv[6])

def read_name_dump(path):
    taxid_to_name = {}
    with open(path, 'r') as stream:
        for line in stream:
            cells = [cell.rstrip().lstrip() for cell in line.rstrip("\n").split("|")]
            taxid = cells[0]
            name = cells[1]
            taxid_to_name[cells[0]] = cells[1]
    return taxid_to_name

def get_gff_attributes(attrs):
    try:
        parts = attrs.split(";")
    except:
        print("Could not split the following attributes:")
        print(str(attrs))
        raise ValueError("Attribute spliting error")
    last_named_part = 0
    for i in range(1,len(parts)):
        if "=" in parts[i]:
            last_named_part = i
        else:
            parts[last_named_part] += ";"+parts[i]
            parts[i] = ""
    values = {}
    for part in parts:
        if len(part) > 0:
            ab = part.split("=")
            if len(ab) == 2:
                name = ab[0]
                val = ab[1].lstrip("'").rstrip("'")
                values[name] = val
    return values

def read_annotation(gff_path):
    details = {}
    with open(gff_path,'r') as stream:
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            source = cells[1]
            attrs = get_gff_attributes(cells[-1])
            attrs['source'] = source
            details[attrs['ID']] = attrs
    return details

def runCommand(cmd, print_cmd=True):
    if print_cmd:
        print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def align(query, db, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = outdir+"/alignments.paf"
    if not os.path.exists(outfile):
        cmd = " ".join(["minimap2",
            "-x splice:hq -uf -t", str(threads),
            db, query,
            ">", outfile])
        code = runCommand(cmd)
        if code != 0:
            print("Error: Alignment unsucessful.")
            os.remove(outfile)
            return None
    return outfile

def get_best_mapping(df):
    sorted = df.sort_values(["identity","qcovs","matchs"],
            ascending=[False,False,False])
    '''if len(sorted) > 1:
        print("Alignments: ")
        print(sorted)
        print("Selecting:", sorted.iloc[0])'''
    return sorted.iloc[0]

paf_file = align(query, db, outdir)
if paf_file == None:
    quit()

print("Reading taxon names")
#taxid = read_name_dump(names_dump_path)

from ete3 import NCBITaxa
ncbi = NCBITaxa()

print("Reading gff annotation for ncRNA details")
gene_details = read_annotation(gff_annotation)

print('Reading paf file')
minimap_df = pd.read_csv(paf_file, sep='\t', header=None, index_col=False,
            names=["qseqid","qseq_len","qstart","qend","strand",
                "sseqid","sseq_len","sstart","send","matchs",
                "block_len","quality","13th","14th","15th","16th","17th","18th"])
minimap_df = minimap_df.astype({"qstart": 'int32', "qend": 'int32', "qseq_len": "int32",
            "sstart": 'int32', "send": 'int32', "sseq_len": "int32",
            "quality": 'int32', "block_len": "int32", "matchs": "int32"})

print("Filtering...")

minimap_df["qcovs"] = minimap_df.apply(
    lambda row: (row["qend"]-row["qstart"]) / row["qseq_len"], axis=1)
minimap_df["identity"] = minimap_df.apply(
    lambda row: row["matchs"] / row["block_len"], axis=1)

print(str(minimap_df.head()))
print(str(len(minimap_df)) + " alignments")

low_th = 0.80
medium_th = 0.85
high_th = 0.90
def classify(value):
    if value >= high_th:
        return "HIGH"
    elif value >= medium_th:
        return "MEDIUM"
    elif value >= low_th:
        return "LOW"
    else:
        return "INVALID"

minimap_df["result"] = minimap_df.apply(
    lambda row: classify(min(row['identity'], row['qcovs'])), axis=1)
minimap_df = minimap_df[minimap_df['result'] != "INVALID"]
'''minimap_df = minimap_df[minimap_df["identity"] >= min_id]
print(str(len(minimap_df)) + " alignments after filtering by identity")
minimap_df = minimap_df[minimap_df["qcovs"] >= min_cov]
print(str(len(minimap_df)) + " alignments after filtering by coverage")'''

def get_tax_name(id):
    ids = ncbi.translate_to_names([int(id)])
    return str(ids[0])

print("Making new rows")
minimap_df["taxid"] = minimap_df.apply(
    lambda row: get_tax_name(str(row['sseqid'].split("_")[-1])), axis=1)

minimap_df = minimap_df[minimap_df['taxid'] != "Arapaima gigas"]

def get_rna_type_str(rna_name):
    tp = gene_details[rna_name]['type'].split(";")
    if len(tp) == 1 or tp[0] == 'Cis-reg':
        return tp[0]
    else:
        return tp[1]

minimap_df["type"] = minimap_df.apply(
    lambda row:  get_rna_type_str(row['qseqid']), axis=1)

minimap_df["source"] = minimap_df.apply(
    lambda row:  gene_details[row['qseqid']]['source'], axis=1)

i = 0
def get_next_id(i):
    i += 1
    return i

minimap_df["id"] = np.arange(len(minimap_df))

print(minimap_df["id"])

best_hits = set()
for name, hits in minimap_df.groupby(["qseqid"]):
    hit = get_best_mapping(hits)
    best_hits.add(hit['id'])

minimap_df["best_hit"] = minimap_df.apply(
    lambda row: row['id'] in best_hits, axis=1)

minimap_df = minimap_df[minimap_df['best_hit'] == True]

print(str(len(minimap_df)) + " total ncRNA with homologs.")

print("Calculating group stats for sumarry")
def group_stats(col_name, df):  
    data = []
    others = [0,0,0,0]
    for name, homologs in df.groupby([col_name]):
        high = len(homologs[homologs['result'] == "HIGH"])
        medium = len(homologs[homologs['result'] == "MEDIUM"])
        low = len(homologs[homologs['result'] == "LOW"])
        total = len(homologs)
        if total < 10:
            others[0] += low
            others[1] += medium
            others[2] += high
            others[3] += total
        else:
            data.append([name, low, medium, high, total])
    data.sort(key=lambda x: x[-1], reverse=True)
    if others[-1] > 0:
        data.append(['Others'] + others)
    return data

total_data = [len(minimap_df[minimap_df['result'] == val]) for val in ['LOW','MEDIUM','HIGH']] + [len(minimap_df)]
species_data = group_stats('taxid', minimap_df)
source_data = group_stats('source', minimap_df)
type_data = group_stats('type', minimap_df)

species_str = "\n".join(["\t".join([str(cell) for cell in row]) for row in species_data])+"\n"
source_str = "\n".join(["\t".join([str(cell) for cell in row]) for row in source_data])+"\n"
type_str = "\n".join(["\t".join([str(cell) for cell in row]) for row in type_data])+"\n"

with open(outdir+"/sumarry.tsv",'w') as stream:
    stream.write("\t".join(['','Low Similarity','Medium Similarity','High Similarity','Total'])+"\n")
    stream.write("All Homologs\t"+'\t'.join([str(cell) for cell in total_data])+"\n")
    stream.write("Homologs by Species:\n"+species_str)
    stream.write("Homologs by Annotation Strategy:\n"+source_str)
    stream.write("Homologs by ncRNA type:\n"+type_str)

with open(outdir+"/homologs.tsv", 'w') as stream:
    print(minimap_df)
    cols = ['id','qseqid','sseqid','result','taxid','identity','qcovs']
    stream.write("\t".join(cols)+"\n")
    for name, hit in minimap_df.iterrows():
        '''try:'''
        to_write = "\t".join([str(hit[cell_name]) for cell_name in cols])+"\n"
        stream.write(to_write)
        '''except Exception:
            print(hit)
            print(Exception)'''
        
