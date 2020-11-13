#%%
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import subprocess
import sys
import obonet
import networkx as nx
import os

gene_list_dir = sys.argv[1]
#gene_list_dir = "/home/pitagoras/main/dev/on-going/rna_nexus/results/gigas_tissue_specific_lncrna/"
obo_path = sys.argv[2]
predictions_file = sys.argv[3]
#predictions_file = "gigas_predictions/lnc_rna_prediction_normal.tsv"
population_file = sys.argv[4]
#population_file = "/home/pitagoras/main/dev/on-going/rna_nexus/test_data/lnc_list/gigas_lnc.txt"
max_pval = 0.01
if len(sys.argv) > 5:
    max_pval = float(sys.argv[5])
#outdir = "goatools_results"
outdir = gene_list_dir + "/enrichment_analysis"
#obo_path = "/home/pitagoras/main/data/go/go.obo"

def get_first_words(p):
    return set([l.rstrip("\n").split()[0] for l in open(p,'r').readlines()])

def runCommand(cmd, print_cmd=True):
        if print_cmd:
                print("\t> " + cmd)
        process = subprocess.call(cmd, shell=True)
        return process

def make_id2gos(outdir, id2go_file):
    id2gos = {}
    with open(id2go_file, 'r') as stream:
        print("Reading " + id2go_file)
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            transcript_name = cells[0]
            go_name = cells[1]
            if not transcript_name in id2gos:
                id2gos[transcript_name] = set()
            id2gos[transcript_name].add(go_name)
    out_file = outdir + "/associations.tsv"
    with open(out_file, 'w') as stream:
        print("\tWriting in id2gos format")
        for name, gos in id2gos.items():
            stream.write(name+"\t"+";".join(gos)+"\n")
    return id2gos, out_file

runCommand("mkdir " + outdir)
tissue_names = ['Brain', 'Gonad', 'Heart', 
                'Kidney', 'Liver', 'Lung',
                'Muscle', 'Skin']

graph = obonet.read_obo(obo_path)
roots = ['GO:0005575', 'GO:0003674', 'GO:0008150']
print("Reading graph")
obo_nodes = graph.nodes(data=True)
go_ids = [id_ for id_, data in obo_nodes]
correct_id = {}
descriptions={}
print("Solving redundant GO ids")
for ID in tqdm(go_ids):
    if "alt_id" in obo_nodes[ID]:
        for alt_id in obo_nodes[ID]["alt_id"]:
            correct_id[alt_id] = ID
    correct_id[ID] = ID
    descriptions[ID] = obo_nodes[ID]['name']
print("Solved " + str(len(correct_id.keys())))
#%%
#prefix = "gigas_tpm."
print("Parsing gene lists")
lists_to_enrich = []
for tissue in tissue_names:
    lists_to_enrich.append((tissue,gene_list_dir+"/"+tissue+".tissue.txt"))

for sex in ["female", "male"]:
    lists_to_enrich.append(("Skin-"+sex, gene_list_dir+"/Skin."+sex+"_expressed.txt"))
    lists_to_enrich.append((sex+"-diff", gene_list_dir+"/"+sex+"_diff.txt"))

lists_to_enrich.append(("housekeeping",gene_list_dir+"/housekeeping.txt"))
lists_to_enrich.append(("growth",gene_list_dir+"/involved_in_growth.txt"))
lists_to_enrich.append(("maturation",gene_list_dir+"/involved_in_maturation.txt"))

associations, associations_file = make_id2gos(outdir, predictions_file)
# = make_population_from_associations(outdir, associations)
enrichments_dir = outdir + "/enrichments"
runCommand("mkdir " + enrichments_dir)

# %%
for name, list_file in lists_to_enrich:
    if not os.path.exists(list_file):
        print("Could not find " + list_file)

cmds = {name: " ".join(["find_enrichment.py",
                          "--pval="+str(max_pval)+" --indent",
                          "--obo", obo_path,
                          "--outfile", 
                          enrichments_dir+"/"+name+".tsv",
                          list_file, population_file,
                          associations_file])
        for name, list_file in lists_to_enrich}

for tissue_name, cmd in tqdm(cmds.items()):
    runCommand(cmd)