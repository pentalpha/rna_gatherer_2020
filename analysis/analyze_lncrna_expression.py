# -*- coding: utf-8 -*-

import pandas as pd
from tqdm import tqdm
import sys
import numpy as np
import math
import os
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
output_path = sys.argv[5]

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
        return list(interest_ann)
    else:
        return []

def is_specific(sample_count, other_samples, counts):
    for other_sample in other_samples:
        if counts[other_sample] > 0:
            if sample_count / counts[other_sample] < 5:
                return False
    return True

def calc_tissue_specificity(counts):
    relevant = False
    tissues_expressed = []
    sample_counts = {}
    for sample in sample_names:
        if counts[sample] >= 1.0:
            relevant = True
            tissues_expressed.append(sample)
            sample_counts[sample] = counts[sample]
    if not relevant:
        return 'Not Expressed', np.nan, np.nan, np.nan, np.nan, np.nan
    
    male_expressions = [x for x in tissues_expressed if "Male" in x]
    female_expressions = [x for x in tissues_expressed if "Female" in x]

    diff_sex = np.nan
    log2fc = np.nan
    min_expressed_tissues = 4
    sex_with_enough_expressions = None
    if (len(male_expressions) >= min_expressed_tissues 
        or len(female_expressions) >= min_expressed_tissues):
            male_counts = [c 
                for sample_name, c in sample_counts.items()
                if "Male" in sample_name]
            female_counts = [c 
                for sample_name, c in sample_counts.items()
                if "Female" in sample_name]
            avg_male = sum(male_counts) / len(male_counts) if len(male_counts) > 0 else 0.0
            avg_female = sum(female_counts) / len(female_counts) if len(female_counts) > 0 else 0.0
            if avg_male == 0:
                diff_sex = "Female"
                log2fc = float('inf')
            elif avg_female == 0:
                diff_sex = "Male"
                log2fc = float('inf')
            else:
                log2fc = abs(math.log((avg_male/avg_female),2))
                if log2fc >= 2.5:
                    diff_sex = "Male" if avg_male > avg_female else "Female"

    higher_tpm = sample_names[0]
    for sample in sample_names:
        if counts[sample] > counts[higher_tpm]:
            higher_tpm = sample
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
    specific_sex = np.nan
    if len(groups[group_name]) == 2:
        male_counts = counts[group_name+"_Male"]
        female_counts = counts[group_name+"_Female"]
        if n_specifics == 2:
            if male_counts * 10 < female_counts:
                specific_sex = "Female"
            elif female_counts * 10 < male_counts:
                specific_sex = "Male"
        elif n_specifics == 1:
            if male_counts > female_counts:
                specific_sex = "Male"
            else:
                specific_sex = "Female"
    
    if tissue_specific != None:
        return 'Tissue Specific', group_name, specific_sex, diff_sex, log2fc
    else:
        expressed_in_all = True
        for sample in sample_names:
            if counts[sample] < 1.0:
                expressed_in_all = False
                break
        if expressed_in_all:
            return 'Expressed in All', np.nan, specific_sex, diff_sex, log2fc
        else:
            return 'Others', np.nan, specific_sex, diff_sex, log2fc

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

# %%
print("Creating new columns with results")
df['Classification'] = df.apply(lambda row: results[row['Name']][0], axis=1)
df['Specific_Tissue'] = df.apply(lambda row: results[row['Name']][1], axis=1)
df['Sex_Expression'] = df.apply(lambda row: results[row['Name']][2], axis=1)
df['Diff_Sex'] = df.apply(lambda row: results[row['Name']][3], axis=1)
df['log2_FC'] = df.apply(lambda row: results[row['Name']][4], axis=1)

df['Involved_in_Maturation'] = df.apply(
    lambda row: len(genes_annotations_in_set(row['Name'],associations, maturation_names)) > 0, axis=1)
df['Involved_in_Growth'] = df.apply(
    lambda row: len(genes_annotations_in_set(row['Name'],associations, growth_names)) > 0, axis=1)

df2 = df.copy(deep=True)
for sample_name in sample_names:
    del df2[sample_name]

df2.to_csv(df2_path, sep='\t', index=False, header=True)
print("Writing more gene lists")
write_names_to_file(df2[df2['Classification'] == 'Expressed in All']['Name'].tolist(), 
                    output_path +"/housekeeping.txt")
write_names_to_file(df2[df2['Diff_Sex'] == 'Male']['Name'].tolist(), 
                    output_path +"/male_diff.txt")
write_names_to_file(df2[df2['Diff_Sex'] == 'Female']['Name'].tolist(), 
                    output_path +"/female_diff.txt")
write_names_to_file(df2[df2['Involved_in_Growth'] == True]['Name'].tolist(), 
                    output_path +"/involved_in_growth.txt")
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
    male_specific = set(tissue_df[tissue_df['Sex_Expression'] == 'Male']['Name'].tolist())
    female_specific = set(tissue_df[tissue_df['Sex_Expression'] == 'Female']['Name'].tolist())
    
    male_diff = set(tissue_df[tissue_df['Diff_Sex'] == 'Male']['Name'].tolist())
    female_diff = set(tissue_df[tissue_df['Diff_Sex'] == 'Female']['Name'].tolist())
    growth = set(tissue_df[tissue_df['Involved_in_Growth'] == True]['Name'].tolist())
    maturation = set(tissue_df[tissue_df['Involved_in_Maturation'] == True]['Name'].tolist())
    tissue_data.append(
        {'Tissue Name': tissue_name, 
         'Tissue Specific': len(tissue_df),
         'Male Expressed': len(male_specific),
         'Female Expressed': len(female_specific),
         'Male Differential': len(male_diff),
         'Female Differential': len(female_diff),
         'Growth Genes': len(growth),
         'Maturation Genes': len(maturation)
         })
    
    name_list_prefix = output_path +"/"+tissue_name+"."
    print("\tWriting genes")
    names = write_names_to_file(names, name_list_prefix+"tissue.txt")

    if "Skin" in tissue_name:
        print("\tWriting male expressed genes")
        male_specific = write_names_to_file(male_specific,
            name_list_prefix+"male_expressed.txt")
        print("\tWriting female expressed genes")
        female_specific = write_names_to_file(female_specific,
            name_list_prefix+"female_expressed.txt")
    print("\tDone")
    #name_lists[tissue_name] = names
print("Calculating for others")
print("\tOthers")
tissue_data.append({'Tissue Name': "Others", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Others']),
         'Male Expressed:': np.nan,
         'Female Expressed': np.nan,
         'Male Differential': 
            len(df2[df2['Classification'] == 'Others'][df2['Diff_Sex'] == 'Male']['Name'].tolist()),
         'Female Differential': 
            len(df2[df2['Classification'] == 'Others'][df2['Diff_Sex'] == 'Female']['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Others'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Others'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
print("\tHouse keeping")
tissue_data.append({'Tissue Name': "House-keeping", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Expressed in All']),
         'Male Expressed': np.nan,
         'Female Expressed': np.nan,
         'Male Differential': 
            len(df2[df2['Classification'] == 'Expressed in All'][df2['Diff_Sex'] == 'Male']['Name'].tolist()),
         'Female Differential': 
            len(df2[df2['Classification'] == 'Expressed in All'][df2['Diff_Sex'] == 'Female']['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Expressed in All'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Expressed in All'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
print("\tNot expressed")
tissue_data.append({'Tissue Name': "Not Expressed", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Not Expressed']),
         'Male Expressed': np.nan,
         'Female Expressed': np.nan,
         'Male Differential': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Diff_Sex'] == 'Male']['Name'].tolist()),
         'Female Differential': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Diff_Sex'] == 'Female']['Name'].tolist()),
         'Growth Genes': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
         'Maturation Genes': 
            len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})


print("Writing dataframe")
summary_df = pd.DataFrame(data=tissue_data,columns=tissue_data[0].keys())
summary_df.to_csv(summary_path, sep='\t', index=False, header=True)

print(summary_df)
#%%
'''from math import log
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from random import shuffle
def get_counts_list(name):
    row = df[df['Name'] == name].iloc[0]
    return [row[sample_name] for sample_name in sample_names]

def convert_to_log(vec, base=2):
    return [log(x+1, base) for x in vec]

sking_others, skin_male, skin_female = name_lists['Skin']
#shuffle(sking_others)
#shuffle(skin_male)
#shuffle(skin_female)
counts = [convert_to_log(get_counts_list(name)) 
          for name in skin_male+skin_female]

name_axis_labels = skin_male+skin_female
sample_name_axis = sample_names

fig, ax = plt.subplots(figsize=(8,14))
im = ax.imshow(counts, interpolation='nearest', aspect='auto')
#ax.set_title("Heatmap of expression counts for sex-specific skin lncRNA")

ax.set_xticks(np.arange(len(sample_name_axis)))
ax.set_yticks(np.arange(len(name_axis_labels)))
# ... and label them with the respective list entries
ax.set_xticklabels(sample_name_axis)
ax.set_yticklabels(name_axis_labels, fontsize=9)

cbarlabel_0 = "log2(TPM+1)"
cbar = ax.figure.colorbar(im, ax=ax, 
                          cmap="YlGn")
cbar.ax.set_ylabel(cbarlabel_0, rotation=-90, va="bottom")

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

fig.tight_layout()
plt.show()
fig.savefig(output_path+"_sking_sex-heatmap.png")'''