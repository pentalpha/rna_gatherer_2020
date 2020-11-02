# -*- coding: utf-8 -*-

import pandas as pd
from tqdm import tqdm
import sys
import numpy as np
#df_path = sys.argv[1]
#output_path = sys.argv[2]
df_path = "../test_data/counts/gigas-tpm.tsv"
output_path = "../results/gigas_tissue_specific_lncrna/gigas_tpm"
df = pd.read_csv(df_path, sep='\t', header=0)

def not_coding(name):
    return name.startswith('jg') and len(name.split(".")) == 2

df['is_coding'] = df.apply(lambda row: not_coding(row['Name']), axis=1)
df = df[df['is_coding'] == False]
del df['is_coding']

#%%

groups = {'Heart': ['Heart_Male'], 
          'Brain': ['Brain_Female', 'Brain_Male'], 
          'Liver': ['Liver_Female', 'Liver_Male'], 
          'Skin': ['Skin_Female', 'Skin_Male'],
          'Muscle':['Muscle_Female'],
          'Gonad': ['Gonad_Female', 'Gonad_Male'], 
          'Lung': ['Lung_Female', 'Lung_Male'],
          'Kidney': ['Kidney_Female', 'Kidney_Male']}

sample_names = df.columns.tolist()[1:]

def is_specific(sample_count, other_samples, counts):
    for other_sample in other_samples:
        if counts[other_sample] > 0:
            if sample_count / counts[other_sample] < 5:
                return False
    return True

def calc_tissue_specificity(counts):
    relevant = False
    for sample in sample_names:
        if counts[sample] >= 1.0:
            relevant = True
            break
    if not relevant:
        return 'Not Expressed', np.nan, np.nan, np.nan
    
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
    
    n_specifics = sum([is_specific(counts[s], other_samples, counts)
                       for s in groups[group_name]])
    tissue_specific = group_name if n_specifics > 0 else None
    specific_sex = np.nan
    if n_specifics == 2:
        male_counts = counts[group_name+"_Male"]
        female_counts = counts[group_name+"_Female"]
        if male_counts * 5 < female_counts:
            specific_sex = "Female"
        elif female_counts * 5 < male_counts:
            specific_sex = "Male"
    
    if tissue_specific != None:
        return 'Tissue Specific', group_name, specific_sex
    else:
        expressed_in_all = True
        for sample in sample_names:
            if counts[sample] < 1.0:
                expressed_in_all = False
                break
        if expressed_in_all:
            return 'Expressed in All', np.nan, specific_sex
        else:
            return 'Others', np.nan, specific_sex
        
results = {}
bar = tqdm(total=len(df))
for index, row in df.iterrows():
    results[row['Name']] = calc_tissue_specificity(row)
    bar.update(1)
bar.close()

# %%
df['Classification'] = df.apply(lambda row: results[row['Name']][0], axis=1)
df['Specific_Tissue'] = df.apply(lambda row: results[row['Name']][1], axis=1)
df['Specific_Sex'] = df.apply(lambda row: results[row['Name']][2], axis=1)

df2 = df.copy(deep=True)
for sample_name in sample_names:
    del df2[sample_name]
df2_path = output_path+".tissue_analysis.tsv"
df2.to_csv(df2_path, sep='\t', index=False, header=True)

#%%
def get_max_expression(name):
    return max([val 
                for val in df[df['Name'] == name].values.tolist()[0] 
                if isinstance(val, float)])

def write_names_to_file(names, f):
    name_and_expression = [(name, get_max_expression(name)) for name in names]
    name_and_expression.sort(key=lambda x: x[1], reverse=True)
    with open(f,'w') as stream:
        for name, max_expression in name_and_expression:
            stream.write(name+"\t"+str(max_expression)+"\n")

df2 = pd.read_csv(df2_path, sep="\t", header=0)
tissue_data = []
summary_path = output_path+".tissue_sumarry.tsv"
for tissue_name, tissue_df in tqdm(
        list(df2[df2['Classification'] == 'Tissue Specific'].groupby('Specific_Tissue'))):
    print(tissue_name)
    names = set(tissue_df['Name'].tolist())
    male_specific = set(tissue_df[tissue_df['Specific_Sex'] == 'Male']['Name'].tolist())
    female_specific = set(tissue_df[tissue_df['Specific_Sex'] == 'Female']['Name'].tolist())
    #sex_indiferent = tissue_df[tissue_df['Name'].isin(tissue_and_sample)]
    sex_indiferent = tissue_df[tissue_df['Specific_Sex'] != 'Female']
    sex_indiferent = set(sex_indiferent[sex_indiferent['Specific_Sex'] != 'Male']['Name'].tolist())
    tissue_data.append(
        {'Tissue Name': tissue_name, 
         'Tissue Specific': len(tissue_df),
         'Male Specific': len(male_specific),
         'Female Specific': len(female_specific),
         'Not Sex Specific': len(sex_indiferent)})
    
    name_lists_prefix = output_path +"."+tissue_name+"-"
    print("\tWriting genes")
    write_names_to_file(names, name_lists_prefix+"all.txt")
    print("\tWriting male genes")
    write_names_to_file(male_specific, name_lists_prefix+"male.txt")
    print("\tWriting female genes")
    write_names_to_file(female_specific, name_lists_prefix+"female.txt")

tissue_data.append({'Tissue Name': "Others", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Others']),
         'Male Specific:': np.nan,
         'Female Specific': np.nan,
         'Not Sex Specific': np.nan})
tissue_data.append({'Tissue Name': "House-keeping", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Expressed in All']),
         'Male Specific': np.nan,
         'Female Specific': np.nan,
         'Not Sex Specific': np.nan})
tissue_data.append({'Tissue Name': "Not Expressed", 
         'Tissue Specific': len(df2[df2['Classification'] == 'Not Expressed']),
         'Male Specific': np.nan,
         'Female Specific': np.nan,
         'Not Sex Specific': np.nan})
write_names_to_file(df2[df2['Classification'] == 'Expressed in All']['Name'].tolist(), 
                    output_path +".housekeeping.txt")
summary_df = pd.DataFrame(data=tissue_data,columns=tissue_data[0].keys())
summary_df.to_csv(summary_path, sep='\t', index=False, header=True)
