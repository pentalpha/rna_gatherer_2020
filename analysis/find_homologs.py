import sys
import os
import subprocess
import pandas as pd
import numpy as np
from ete3 import NCBITaxa

query = sys.argv[1]
db = sys.argv[2]
gff_annotation = sys.argv[3]
outdir = sys.argv[4]
threads = int(sys.argv[5])

ncbi = NCBITaxa()

def evol_sim(taxid1, taxid2):
    l1 = set(ncbi.get_lineage(taxid1))
    l2 = set(ncbi.get_lineage(taxid2))
    common_taxid = l1.intersection(l2)
    return len(common_taxid)

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
    sorted = df.sort_values(["identity","qcovs","common_taxid"],
            ascending=[False,False,False])
    return sorted.iloc[0]

paf_file = align(query, db, outdir)
if paf_file == None:
    quit()

print("Reading taxon names")
#taxid = read_name_dump(names_dump_path)

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

def get_tax_name(id):
    ids = ncbi.translate_to_names([id])
    return str(ids[0])

print("Making new columns")
minimap_df["taxid"] = minimap_df.apply(
    lambda row: int(row['sseqid'].split("_")[-1]), axis=1)

minimap_df["species"] = minimap_df.apply(
    lambda row: get_tax_name(row['taxid']), axis=1)

minimap_df = minimap_df[minimap_df['species'] != "Arapaima gigas"]

def get_rna_type_str(rna_name):
    tp = gene_details[rna_name]['type'].split(";")
    if len(tp) == 1 or tp[0] == 'Cis-reg':
        return tp[0]
    else:
        return tp[1]

minimap_df["type"] = minimap_df.apply(
    lambda row:  get_rna_type_str(row['qseqid']), axis=1)

types_annotated = set(get_rna_type_str(ID) for ID in gene_details.keys())
def count_type(tp_str):
    count = 0
    for ID, details in gene_details.items():
        if get_rna_type_str(ID) == tp_str:
            count += 1
    return count
type_size = {tp_str: count_type(tp_str) for tp_str in types_annotated}

minimap_df["source"] = minimap_df.apply(
    lambda row:  gene_details[row['qseqid']]['source'], axis=1)

minimap_df["id"] = np.arange(len(minimap_df))

print("Calculating taxonomic closeness")
minimap_df["common_taxid"] = minimap_df.apply(lambda row: evol_sim(row['taxid'], 113544), axis=1)

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
minimap_filtered = minimap_df[minimap_df['result'] != "INVALID"]
'''minimap_df = minimap_df[minimap_df["identity"] >= min_id]
print(str(len(minimap_df)) + " alignments after filtering by identity")
minimap_df = minimap_df[minimap_df["qcovs"] >= min_cov]
print(str(len(minimap_df)) + " alignments after filtering by coverage")'''

print("Finding best hits")
best_hits = set()
for name, hits in minimap_filtered.groupby(["qseqid"]):
    hit = get_best_mapping(hits)
    best_hits.add(hit['id'])

minimap_filtered["best_hit"] = minimap_filtered.apply(
    lambda row: row['id'] in best_hits, axis=1)

minimap_bests = minimap_filtered[minimap_filtered['best_hit'] == True]

print(str(len(minimap_bests)) + " total ncRNA with homologs.")

print("Calculating group stats for sumarry")
def group_stats(col_name, df):  
    data = []
    #others = [0,0,0,0]
    for name, homologs in df.groupby([col_name]):
        high = len(homologs[homologs['result'] == "HIGH"])
        medium = len(homologs[homologs['result'] == "MEDIUM"])
        low = len(homologs[homologs['result'] == "LOW"])
        total = len(homologs)
        '''if total < 10:
            others[0] += low
            others[1] += medium
            others[2] += high
            others[3] += total
        else:'''
        data.append([name, low, medium, high, total])
    data.sort(key=lambda x: x[-1], reverse=True)
    '''if others[-1] > 0:
        data.append(['Others'] + others)'''
    return data

total_data = [len(minimap_bests[minimap_bests['result'] == val]) for val in ['LOW','MEDIUM','HIGH']] + [len(minimap_bests)]
species_data = group_stats('taxid', minimap_bests)
#source_data = group_stats('source', minimap_df)
type_data = group_stats('type', minimap_bests)

types_conserved = []
for i in range(len(type_data)):
    tp, low, mid, high, total = type_data[i]
    conserved = round((total / type_size[tp])*100, 2)
    type_data[i] += [str(conserved)]
    types_conserved.append(tp)

for tp, count in type_size.items():
    if not tp in types_conserved:
        type_data.append([tp,0,0,0,0,0.0])

min_genes_in_taxa = 13
#arapaima_lineage = set(ncbi.get_lineage(113544)[:-1])
min_height = 16
def group_species(species_data):
    common_names = ncbi.get_common_names([taxid for taxid, low, med, high, total in species_data])
    common_names[41665] = "bony fishes"
    main_data = [row for row in species_data if row[-1] >= min_genes_in_taxa]
    species_names = ncbi.translate_to_names([taxid for taxid, low, med, high, total in main_data])
    assert len(main_data) == len(species_names)
    for i in range(len(main_data)):
        taxid = main_data[i][0]
        other_name = "NCBI:taxid"+str(taxid)
        if taxid in common_names:
            other_name = common_names[taxid]
        main_data[i][0] = species_names[i] + " (" + other_name + ")\t"

    others = [row for row in species_data if row[-1] < min_genes_in_taxa]
    
    others_sumarry = {}
    others_total = [0,0,0,0]
    for taxid, low, med, high, total in others:
        l = ncbi.get_lineage(taxid)
        parent_taxid = l[-2]
        if len(l) > min_height:
            parent_taxid = l[min_height]
        if 2 in l:
            #group bacteria
            parent_taxid = 2
        elif 10239 in l:
            #group viruses
            parent_taxid = 10239
        elif 2787823 in l or 2787854 in l:
            #group unclassified
            parent_taxid = 2787823
        elif 6656 in l:
            #group arthropods
            parent_taxid = 6656

        if not parent_taxid in others_sumarry:
            others_sumarry[parent_taxid] = [0,0,0,0]
            common_name_attempt = ncbi.get_common_names([parent_taxid])
            if parent_taxid in common_name_attempt:
                common_names[parent_taxid] = common_name_attempt[parent_taxid]
        others_sumarry[parent_taxid][0] += low
        others_total[0] += low
        others_sumarry[parent_taxid][1] += med
        others_total[1] += med
        others_sumarry[parent_taxid][2] += high
        others_total[2] += high
        others_sumarry[parent_taxid][3] += total
        others_total[3] += total
    
    others_list = []
    for taxid, data in others_sumarry.items():
        others_list.append([taxid]+data)
    others_list.sort(key=lambda x: x[-1], reverse=True)
    species_names_others = ncbi.translate_to_names([taxid for taxid, low, med, high, total in others_list])
    assert len(others_list) == len(species_names_others)
    for i in range(len(others_list)):
        taxid = others_list[i][0]
        other_name = "NCBI:txid"+str(taxid)
        if taxid in common_names:
            other_name = common_names[taxid]
        others_list[i][0] = "\t" + species_names_others[i] + " (" + other_name + ")"
    main_data.append(["Others:\t"] + others_total)
    main_data += others_list
    return main_data

print("Grouping less represented species")
species_data = group_species(species_data)
#species_data.sort(key=lambda x: x[-1], reverse=True)
print("Outputing")

main_taxa = [row for row in species_data if row[-1] >= min_genes_in_taxa-2]
less_frequent_taxa = ["\tLess Frequent Taxa",0,0,0,0]
for taxid, low, med, high, total in [row for row in species_data if row[-1] < min_genes_in_taxa-2]:
    less_frequent_taxa[1] += low
    less_frequent_taxa[2] += med
    less_frequent_taxa[3] += high
    less_frequent_taxa[4] += total
main_taxa.append(less_frequent_taxa)

def cells_to_str(cell_list):
    cell_list_str = [str(x) for x in cell_list]
    return cell_list_str

species_str = "\n".join(["\t".join([str(cell) for cell in row]) for row in species_data])+"\n"
species_str_short = "\n".join(["\t".join([str(cell) for cell in row]) for row in main_taxa])+"\n"
#source_str = "\n".join(["\t".join(cells_to_str(row)) for row in source_data])+"\n"
type_str = "\n".join(["\t".join(cells_to_str(row)) for row in type_data])+"\n"

with open(outdir+"/by_type.tsv",'w') as stream:
    stream.write("\t".join(['','Low Similarity','Medium Similarity','High Similarity','Total',
                    'Predictions With An Homolog (%)'])+"\n")
    #stream.write("All Homologs\t"+'\t'.join([str(cell) for cell in total_data])+"\n")
    stream.write(type_str)

with open(outdir+"/by_taxon.tsv",'w') as stream:
    stream.write("\t".join(['','','Low Similarity','Medium Similarity','High Similarity','Total'])+"\n")
    stream.write("All Homologs\t"+'\t'.join([str(cell) for cell in total_data])+"\n")
    stream.write(species_str_short)

with open(outdir+"/by_taxon_full_list.tsv",'w') as stream:
    stream.write("\t".join(['','','Low Similarity','Medium Similarity','High Similarity','Total'])+"\n")
    stream.write("All Homologs\t"+'\t'.join([str(cell) for cell in total_data])+"\n")
    stream.write(species_str)

with open(outdir+"/homologs.tsv", 'w') as stream:
    print(minimap_bests)
    cols = ['id','qseqid','sseqid','result','taxid','species','type','identity','qcovs']
    stream.write("\t".join(cols)+"\n")
    for name, hit in minimap_bests.iterrows():
        '''try:'''
        to_write = "\t".join([str(hit[cell_name]) for cell_name in cols])+"\n"
        stream.write(to_write)
        '''except Exception:
            print(hit)
            print(Exception)'''

with open(outdir+"/alignments.tsv", 'w') as stream:
    print(minimap_df)
    cols = ['id','qseqid','sseqid','result','taxid','species','type','identity','qcovs']
    stream.write("\t".join(cols)+"\n")
    for name, hit in minimap_df.iterrows():
        '''try:'''
        to_write = "\t".join([str(hit[cell_name]) for cell_name in cols])+"\n"
        stream.write(to_write)
        
# %%
print("Getting bacterial best hits")
bacterial_types = {}
is_bacterial = {x: 2 in ncbi.get_lineage(x) for x in set(minimap_df['taxid'].tolist())}
bacterial_taxids = set([x for x,y in is_bacterial.items() if y == True])
print(len(bacterial_taxids), " bacterial species found")
bacterial_df = minimap_bests[minimap_bests['taxid'].isin(bacterial_taxids)]
with_bact_homolog = bacterial_df['qseqid'].tolist()
possible_hgt_df = minimap_filtered[minimap_filtered['qseqid'].isin(with_bact_homolog)]
for query, data in possible_hgt_df.groupby(['qseqid']):
    bact_homos = set(data[data['taxid'].isin(bacterial_taxids)]['species'].tolist())
    euc_homos = set(data[~data['taxid'].isin(bacterial_taxids)]['species'].tolist())
    if len(euc_homos) > 0:
        print(query, " is homolog to these bacteria:")
        print("\t"+"\n\t".join(bact_homos))
        print("and these eucariots:")
        print("\t"+"\n\t".join(euc_homos))
        print("\n")
'''eucariot_df = possible_hgt_df[~possible_hgt_df['taxid'].isin(bacterial_taxids)]
print(len(eucariot_df))
print(eucariot_df)'''

#artificial sequences (vectors and etc)
artificial_id = 29278