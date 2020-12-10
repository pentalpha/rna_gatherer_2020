#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 09:04:37 2020

@author: pitagoras
"""
import pandas as pd
import numpy as np
import os

gatherer_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))+"/../")
print("RNA Gatherer dir = ", gatherer_dir)
df_path = gatherer_dir + "/result_data/gigas_expression_analysis/tissue_analysis.tsv"

associations_name = gatherer_dir + "/result_data/gigas_geneset_enrichment_analysis/associations.tsv"
correlations_name = gatherer_dir + "/result_data/gigas_annotation/gigas-lnc_coding-correlations-spearman.tsv"
growth_genes_path = gatherer_dir + "/data/function_sets/growth_functions.txt"

coding_ass_path = gatherer_dir + "/test_data/annotation/identified_annotated_coding_mRNA.id2go.tsv"
coding_realnames_path = gatherer_dir + "/test_data/martins_annotation/gene_names.tsv"

analysis_dir = gatherer_dir + '/result_data/gigas_expression_analysis'

df = pd.read_csv(df_path,sep="\t")
df = df[df['Classification'] == "Housekeeping"]
df = df[df['Involved_in_Growth'] == True]
gene_names = df['Name'].tolist()

print(len(gene_names), " housekeeping lncRNA")

growth_funcs = {}

with open(growth_genes_path, 'r') as stream:
    for line in stream:
        name, description = line.rstrip('\n').split("\t")
        growth_funcs[name] = description

print(len(growth_funcs.keys()), 'growth terms')

growth_coding_genes = {}
with open(coding_ass_path,'r') as stream:
    for line in stream:
        go, name, aspect = line.rstrip('\n').split("\t")
        if go in growth_funcs:
            if not name in growth_coding_genes:
                growth_coding_genes[name] = []
            growth_coding_genes[name].append(go)

print(len(growth_coding_genes.keys()), 'growth coding genes')

associations = {}
with open(associations_name,'r') as stream:
    for line in stream:
        name, gos = line.rstrip('\n').split("\t")
        if name in gene_names:
            gos = gos.split(";")
            new_funcs = []
            for go in gos:
                if go in growth_funcs:
                    new_funcs.append((go, growth_funcs[go]))
            if len(new_funcs) > 0:
                associations[name] = new_funcs

print(len(associations.keys()), 'growth lncrna')            
            
correlations = {name: [] for name in associations.keys()}
with open(correlations_name, 'r') as stream:
    for line in stream.readlines():
        cod_name, lnc_name, corr = line.rstrip('\n').split("\t")
        if lnc_name in correlations:
            if cod_name in growth_coding_genes:
                correlations[lnc_name].append((cod_name, float(corr)))

for name in correlations.keys():
    correlations[name].sort(key=lambda x: x[1], reverse=True)

real_names = {}
coding_descriptions = {}
with open(coding_realnames_path, 'r') as stream:
    stream.readline()
    #print(stream.readline().split('\t'))
    for line in stream.readlines():
        cols = line.rstrip('\n').split("\t")
        name, real_name, descript = (cols[1], cols[3], cols[4])
        if name in growth_coding_genes:
            real_names[name] = real_name
            coding_descriptions[name] = descript

#%%

corr_output = "\t".join(['lncRNA transcript', 'Correlated Growth Proteins', 'Unique Protein Annotations', 'lncRNA Growth functions', 'Protein Descriptions'])+"\n"
for lnc_name in associations:
    protein_transcripts = len(correlations[lnc_name])
    correlated = set([real_names[x]+": "+coding_descriptions[x].replace(' [Scleropages formosus]', '')
                      for x,y in correlations[lnc_name]])
    correlated_str = "; ".join(correlated)
    protein_genes = len(correlated)
    funcs = [go_id+": "+descript
             for go_id, descript in associations[lnc_name]]
    funcs_str = "; ".join(funcs)
    corr_output += "\t".join([lnc_name,str(protein_transcripts),str(protein_genes),
                              funcs_str, correlated_str])+"\n"
print(corr_output)
open(analysis_dir + "/growth_housekeeping_analysis.tsv",'w').write(corr_output)