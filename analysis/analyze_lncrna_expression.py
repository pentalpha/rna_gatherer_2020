# -*- coding: utf-8 -*-

import pandas as pd
from tqdm import tqdm
import sys
import numpy as np
import math
import os
import obonet
import networkx as nx
'''
usage:
    analyze_lncrna_expression.py <counts>.tsv <annotation>.tsv \
        <growth_names>.txt <maturation_names>.txt \
        <output_path>
'''

df_path = sys.argv[1]
associations_path = sys.argv[2]
growth_names_path = sys.argv[3]
maturation_names_path = sys.argv[4]
peridot_results = sys.argv[5]
obo_path = sys.argv[6]
output_path = sys.argv[7]

if not os.path.exists(output_path):
    os.mkdir(output_path)

def get_go_list(p):
    return set([l.rstrip("\n").split()[0] for l in open(p,'r').readlines()])

def read_id2gos(id2gos_file):
    id2gos = {}
    with open(id2gos_file, 'r') as stream:
        print("Reading " + id2gos_file)
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            transcript_name = cells[0]
            go_names = cells[1].split(";")
            if not transcript_name in id2gos:
                id2gos[transcript_name] = set()
            for go_name in go_names:
                id2gos[transcript_name].add(go_name)
    return id2gos

#df_path = "../test_data/counts/gigas-tpm.tsv"
#output_path = "../results/gigas_tissue_specific_lncrna/gigas_tpm"
print("Reading counts")
df = pd.read_csv(df_path, sep='\t', header=0)

def not_coding(name):
    return name.startswith('jg') and len(name.split(".")) == 2

df['is_coding'] = df.apply(lambda row: not_coding(row['Name']), axis=1)
df = df[df['is_coding'] == False]
del df['is_coding']

associations = read_id2gos(associations_path)
print("Reading GO lists")
growth_names = get_go_list(growth_names_path)
maturation_names = get_go_list(maturation_names_path)

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

groups = {'Heart': ['Heart_Male'], 
          'Brain': ['Brain_Female', 'Brain_Male'], 
          'Liver': ['Liver_Female', 'Liver_Male'], 
          'Skin': ['Skin_Female', 'Skin_Male'],
          'Muscle':['Muscle_Female'],
          'Gonad': ['Gonad_Female', 'Gonad_Male'], 
          'Lung': ['Lung_Female', 'Lung_Male'],
          'Kidney': ['Kidney_Female', 'Kidney_Male']}

sample_names = ['Heart_Male', 'Muscle_Female',
    'Brain_Female', 'Brain_Male',
    'Liver_Female', 'Liver_Male',
    'Skin_Female', 'Skin_Male',
    'Gonad_Female', 'Gonad_Male',
    'Lung_Female', 'Lung_Male',
    'Kidney_Female', 'Kidney_Male']

def genes_annotations_in_set(gene_id, associations, interest_list):
    if gene_id in associations:
        ann = associations[gene_id]
        interest_ann = ann.intersection(interest_list)
        return [go + ": " + descriptions[correct_id[go]] for go in interest_ann]
        #return list(interest_ann)
    else:
        return []

def is_specific(sample_count, other_samples, counts):
    for other_sample in other_samples:
        if counts[other_sample] > 0:
            if sample_count / counts[other_sample] < 5:
                return False
    return True

def calc_tissue_specificity(counts):
    higher_tpm = sample_names[0]
    lower_tpm = sample_names[0]
    for sample in sample_names:
        if counts[sample] > counts[higher_tpm]:
            higher_tpm = sample
        if counts[sample] < counts[lower_tpm]:
            lower_tpm = sample
    lower_tpm = counts[lower_tpm]

    male_mean = 0.0
    males = 0
    female_mean = 0.0
    females = 0

    relevant = False
    tissues_expressed = []
    sample_counts = {}
    for sample in sample_names:
        if counts[sample] >= 1.0:
            relevant = True
            tissues_expressed.append(sample)
            sample_counts[sample] = counts[sample]
        if "Male" in sample:
            male_mean += counts[sample]
            males += 1
        elif "Female" in sample:
            female_mean += counts[sample]
            females += 1

    sex_specific = np.nan
    if male_mean > 0 and female_mean == 0.0:
        sex_specific = "Male"
    elif female_mean > 0 and male_mean == 0.0:
        sex_specific = "Female"

    male_mean = male_mean / males if males > 0 else 0
    female_mean = female_mean / females if females > 0 else 0
    raw_fc = np.inf
    if male_mean > 0 and female_mean > 0:
        raw_fc = female_mean/male_mean
    elif (not male_mean > 0) and (not female_mean > 0):
        raw_fc = 1.0
    fold_change = abs(math.log(raw_fc,2))

    if not relevant:
        return 'Not Expressed', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific

    other_samples = sample_names[:]
    
    #other_samples.remove(higher_tpm)
    #specific = is_specific(counts[higher_tpm], other_samples, counts)
    group_name = higher_tpm.split('_')[0]
    for s in groups[group_name]:
        other_samples.remove(s)
    
    specifics_bool = [is_specific(counts[s], other_samples, counts)
                       for s in groups[group_name]]
    n_specifics = sum(specifics_bool)
    tissue_specific = group_name if n_specifics > 0 else None
    
    if tissue_specific != None:
        return 'Tissue Specific', group_name, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific
    else:
        expressed_in_all = True
        for sample in sample_names:
            if counts[sample] < 1.0:
                expressed_in_all = False
                break
        if expressed_in_all:
            return 'Housekeeping', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific
        else:
            return 'Mixed Expression', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific

def get_max_expression(name):
    return max([val 
                for val in df[df['Name'] == name].values.tolist()[0] 
                if isinstance(val, float)])

def write_names_to_file(names, f):
    with open(f,'w') as stream:
        for name in names:
            stream.write(name+"\n")
    return names

df2_path = output_path+"/tissue_analysis.tsv"

print("Analyzing expression of each lncrna")
results = {}
bar = tqdm(total=len(df))
for index, row in tqdm(df.iterrows()):
    results[row['Name']] = calc_tissue_specificity(row)
    bar.update(1)
bar.close()

diff_sex_lnc = {}
#read differential genes
consensus_file = peridot_results + "/VennDiagram.PostAnalysisModule/1-Intersect.tsv"
for line in open(consensus_file, 'r'):
    cells = line.rstrip("\n").split("\t")
    name = cells[0].lstrip('"').rstrip('"')
    diff_sex_lnc[name] = len(cells[1].split(","))

# %%
print("Creating new columns with results")
df['Classification'] = df.apply(lambda row: results[row['Name']][0], axis=1)
df['Specific_Tissue'] = df.apply(lambda row: results[row['Name']][1], axis=1)
df['Lowest_Expression_Count'] = df.apply(lambda row: results[row['Name']][2], axis=1)
df['Male_Mean_Expression'] = df.apply(lambda row: results[row['Name']][3], axis=1)
df['Female_Mean_Expression'] = df.apply(lambda row: results[row['Name']][4], axis=1)
df['Log_FC_Expression'] = df.apply(lambda row: results[row['Name']][5], axis=1)
df['Samples_With_Expression'] = df.apply(lambda row: results[row['Name']][6], axis=1)
df['Expressed_Only_In_This_Sex'] = df.apply(lambda row: results[row['Name']][7], axis=1)

print("\tDE Columns")
df['Number_Of_DE_Packages'] = df.apply(lambda row: diff_sex_lnc[row['Name']] if row['Name'] in diff_sex_lnc else 0, axis=1)
df['Diff_Sex'] = df.apply(lambda row: row['Number_Of_DE_Packages'] > 1, axis=1)

print("\tGrowth columns")
df['Growth_Functions'] = df.apply(
    lambda row: len(genes_annotations_in_set(row['Name'],associations, growth_names)), axis=1)
df['Growth_Functions_Percent'] = df.apply(
    lambda row: (row['Growth_Functions'] / len(associations[row['Name']])) if row['Name'] in associations else 0.0,
    axis=1)
df['Involved_in_Growth'] = df.apply(
    lambda row: row['Growth_Functions'] > 0, axis=1)

print("\tMaturation Columns")
df['Maturation_Functions'] = df.apply(
    lambda row: len(genes_annotations_in_set(row['Name'],associations, maturation_names)), axis=1)
df['Maturation_Functions_Percent'] = df.apply(
    lambda row: (row['Maturation_Functions'] / len(associations[row['Name']])) if row['Name'] in associations else 0.0,
    axis=1)
df['Involved_in_Maturation'] = df.apply(
    lambda row: row['Maturation_Functions'] > 0, axis=1)

df2 = df.copy(deep=True)
for sample_name in sample_names:
    del df2[sample_name]

df2.to_csv(df2_path, sep='\t', index=False, header=True)

print("Writing gene lists")
housekeeping_df = df2[df2['Classification'] == 'Housekeeping'][['Name', 'Lowest_Expression_Count', 'Samples_With_Expression']]
housekeeping_df.sort_values(by=['Samples_With_Expression', 'Lowest_Expression_Count'], 
                            ascending = [False, False], inplace=True)
housekeeping_df.to_csv(output_path +"/housekeeping.tsv", sep="\t")
write_names_to_file(df2[df2['Classification'] == 'Housekeeping']['Name'].tolist(), 
                    output_path +"/housekeeping.txt")

sex_diff_df = df2[df2['Diff_Sex'] == True][['Name', 'Classification', 'Specific_Tissue', 'Male_Mean_Expression', 'Female_Mean_Expression', 'Log_FC_Expression', 'Number_Of_DE_Packages']]
sex_diff_df.sort_values(by=['Number_Of_DE_Packages', 'Log_FC_Expression'], 
                            ascending = [False, False], inplace=True)
sex_diff_df.to_csv(output_path +"/sex_diff.tsv", sep="\t")
write_names_to_file(df2[df2['Diff_Sex'] == True]['Name'].tolist(), 
                    output_path +"/sex_diff.txt")

growth_df = df2[df2['Involved_in_Growth'] == True][['Name', 'Classification', 'Specific_Tissue', 'Growth_Functions', 'Growth_Functions_Percent']]
growth_df['Functions'] = growth_df.apply(lambda row: str(genes_annotations_in_set(row['Name'], 
                                                        associations, growth_names)),
                                        axis=1) 
growth_df.sort_values(by=['Growth_Functions', 'Growth_Functions_Percent'], 
                            ascending = [False, False], inplace=True)
growth_df.to_csv(output_path +"/involved_in_growth.tsv", sep="\t")
write_names_to_file(df2[df2['Involved_in_Growth'] == True]['Name'].tolist(), 
                    output_path +"/involved_in_growth.txt")

growth_hk_df = growth_df[growth_df['Classification'] == 'Housekeeping']
growth_hk_df.to_csv(output_path +"/involved_in_growth-housekeeping.tsv", sep="\t")

maturation_df = df2[df2['Involved_in_Maturation'] == True][['Name', 'Classification', 'Specific_Tissue', 'Maturation_Functions', 'Maturation_Functions_Percent', 'Diff_Sex', 'Male_Mean_Expression', 'Female_Mean_Expression']]
maturation_df['Functions'] = maturation_df.apply(lambda row: str(genes_annotations_in_set(row['Name'], 
                                                        associations, maturation_names)),
                                        axis=1)
maturation_df.sort_values(by=['Maturation_Functions', 'Maturation_Functions_Percent'], 
                            ascending = [False, False], inplace=True)
maturation_df.to_csv(output_path +"/involved_in_maturation.tsv", sep="\t")
write_names_to_file(df2[df2['Involved_in_Maturation'] == True]['Name'].tolist(), 
                    output_path +"/involved_in_maturation.txt")
#%%

print("Reading results to print summary")
df2 = pd.read_csv(df2_path, sep="\t", header=0)
tissue_data = []
summary_path = output_path+"/tissue_sumarry.tsv"
name_lists = {}
for tissue_name, tissue_df in tqdm(
        list(df2[df2['Classification'] == 'Tissue Specific'].groupby('Specific_Tissue'))):
    print(tissue_name)
    names = set(tissue_df['Name'].tolist())
    sex_diff = set(tissue_df[tissue_df['Diff_Sex'] == True]['Name'].tolist())
    growth = set(tissue_df[tissue_df['Involved_in_Growth'] == True]['Name'].tolist())
    maturation = set(tissue_df[tissue_df['Involved_in_Maturation'] == True]['Name'].tolist())
    tissue_data.append(
        {'Tissue Name': tissue_name, 
         'Tissue Specific': len(tissue_df),
         'Differentially Expressed by Sex': len(sex_diff),
         'Growth Genes': len(growth),
         'Maturation Genes': len(maturation)
         })
    
    name_list_prefix = output_path +"/"+tissue_name+"."
    print("\tWriting genes")
    names = write_names_to_file(names, name_list_prefix+"tissue.txt")

    if "Skin" in tissue_name:
        male_df = tissue_df[tissue_df['Expressed_Only_In_This_Sex'] == 'Male']
        female_df = tissue_df[tissue_df['Expressed_Only_In_This_Sex'] == 'Female']
        male_specific = write_names_to_file(male_df['Name'].tolist(), name_list_prefix+"male_expressed.txt")
        female_specific = write_names_to_file(female_df['Name'].tolist(), name_list_prefix+"female_expressed.txt")
    print("\tDone")
    #name_lists[tissue_name] = names
print("Calculating for Mixed Expression")
print("\Mixed Expression")
tissue_data.append({'Tissue Name': "Mixed Expression", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Mixed Expression']),
         'Differentially Expressed by Sex': 
            len(df2[df2['Classification'] == 'Mixed Expression'][df2['Diff_Sex'] == True]['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Mixed Expression'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Mixed Expression'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
print("\tHouse keeping")
tissue_data.append({'Tissue Name': "Housekeeping", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Housekeeping']),
         'Differentially Expressed by Sex': 
            len(df2[df2['Classification'] == 'Housekeeping'][df2['Diff_Sex'] == True]['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Housekeeping'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Housekeeping'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
print("\tNot expressed")
tissue_data.append({'Tissue Name': "Not Expressed", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Not Expressed']),
         'Differentially Expressed by Sex': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Diff_Sex'] == True]['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
print("Writing dataframe")
summary_df = pd.DataFrame(data=tissue_data,columns=tissue_data[0].keys())
summary_df.to_csv(summary_path, sep='\t', index=False, header=True)
print(summary_df)