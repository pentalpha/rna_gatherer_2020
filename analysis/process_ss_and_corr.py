import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import dask
import dask.dataframe as dd

ontologies = ["molecular_function",
    "biological_process",
    "cellular_component"]

def associations_to_dict(associations):
    d = {x: set() for x in [a for a,b in associations]}
    for a,b in associations:
        d[a].add(b)
    return d

def read_gaf(gaf_path):
    names_dict = {}
    names = set()
    all_associations = []
    mf = []
    bp = []
    cc = []
    with open(gaf_path, 'r') as stream:
        i = 0
        for line in stream.readlines():
            if (not line.startswith("!")) and not(line.startswith("#")):
                cells = line.rstrip("\n").split("\t")
                db_obj_id = cells[1]
                names.add(db_obj_id)
                go_id = cells[4]
                aspect = cells[8]
                if aspect == "F":
                    mf.append((db_obj_id, go_id))
                elif aspect == "P":
                    bp.append((db_obj_id, go_id))
                elif aspect == "C":
                    cc.append((db_obj_id, go_id))
                all_associations.append((db_obj_id, go_id))
    
    all_anno = associations_to_dict(all_associations)
    mf_anno  = associations_to_dict(mf)
    bp_anno  = associations_to_dict(bp)
    cc_anno  = associations_to_dict(cc)
    names_list = list(names)
    for i in range(len(names_list)):
        names_dict[names_list[i]] = i
    
    return all_anno, mf_anno, bp_anno, cc_anno, names_dict, names_list

def make_histogram(data, title, axe):
    axe.hist(x=data, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
    axe.set_title(title)

def name_to_id(gene_name, n_dict):
    if gene_name in n_dict:
        return n_dict[gene_name]
    else:
        new_id = len(n_dict)+1
        n_dict[gene_name] = new_id
        return new_id

if __name__ == "__main__":
    ss_path = sys.argv[1]
    corr_file_paths = sys.argv[2].split(",")
    gaf_path = sys.argv[3]
    output_dir = sys.argv[4]

    if len(corr_file_paths) == 1:
        sys.argv[2] == "None"
        corr_file_paths = []

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    print("Reading annotation")
    all_anno, mf_anno, bp_anno, cc_anno, names_dict, names = read_gaf(gaf_path)
    max_id = max(names_dict.values())
    print(str(max_id))
    next_new_id = max_id + 1
    print("Reading SS dataframe")
    dask.config.set(scheduler='single-threaded')
    ss_df = dd.read_csv(ss_path, sep="\t", header=None, names=["gene_a","gene_b",
                        "mf", "bp", "cc"])

    print("Creating new ID columns")
    ss_df["gene_id_a"] = ss_df.apply(lambda row: name_to_id(row["gene_a"], 
                                    names_dict),
                                    axis=1, meta=int)
    ss_df["gene_id_b"] = ss_df.apply(lambda row: name_to_id(row["gene_b"],
                                    names_dict),
                                    axis=1, meta=int)
    print("Removing old ones")
    del ss_df["gene_a"]
    del ss_df["gene_b"]
    print("Writing it to disk")
    ss_df.compute().to_csv(output_dir+"/ss.new_indexes.tsv", sep="\t")
    ss_df.compute().to_hdf(output_dir+"/ss.new_indexes.hdf", '/data')
    
    del ss_df

    print("Reading correlation dataframes")
    for corr_file_path in tqdm(corr_file_paths):
        corr_df = dd.read_csv(corr_file_path, sep="\t", header=None, 
                    names=["gene_a","gene_b","corr"])
        print("Creating new ID columns for " + corr_file_path)
        corr_df["gene_id_a"] = corr_df.apply(lambda row: name_to_id(row["gene_a"], 
                                        names_dict),
                                        axis=1, meta=int)
        corr_df["gene_id_b"] = corr_df.apply(lambda row: name_to_id(row["gene_b"],
                                        names_dict),
                                        axis=1, meta=int)
        print("Removing old ones")
        del corr_df["gene_a"]
        del corr_df["gene_b"]

        print("Writing it to disk")
        corr_df.compute().to_hdf(output_dir+"/"+(corr_file_path.split("/")[-1])+".hdf", 
                                '/data')

    print("Writing dictonary to disk")
    with open(output_dir+"/gene_id.tsv", 'w') as stream:
        for gene_name, gene_id in names_dict.items():
            stream.write(gene_name + "\t" + str(gene_id) + "\n")

'''print("Allocating space for data")
max_i = max(names_dict.values())
all_data = {}

ss_name = ss_path.split("/")[-1].split(".")[0]
out_basename = output_dir+"/"+ss_name+"_"

print("Reading " + ss_path)
names_not_found = set()
with open(ss_path,'r') as in_stream:
    for line in in_stream:
        cells = line.rstrip("\n").split("\t")
        if cells[0] in names_dict and cells[1] in names_dict:
            pair_key = (cells[0], cells[1])
            values = [None if x == "None" else float(x) 
                        for x in [cells[2],cells[3],cells[4]]]
            all_data[pair_key] = values
#print(len(names_not_found) + " semantic similarity names not found in annotation")

print("Joining values for distribution of semantic similarities")
all_mf_values = [mf for mf,bp,cc in all_data.values() if mf]
all_cc_values = [cc for mf,bp,cc in all_data.values() if cc]
all_bp_values = [bp for mf,bp,cc in all_data.values() if bp]

print("Making distribution of semantic similarities graphs")
fig, axes = plt.subplots(nrows=4, figsize=(8, 12))
plot_all, plot_mf, plot_cc, plot_bp = axes.flatten()
#all_ss_values = np.concatenate((all_mf_values,all_cc_values,all_bp_values))
make_histogram(np.concatenate([all_mf_values, all_cc_values, all_bp_values]),
    "Similarity Values For All Ontologies", plot_all)
make_histogram(all_mf_values, 
    "Molecular Function Similarities", plot_mf)
make_histogram(all_cc_values,
    "Celullar Component Similarities", plot_cc)
make_histogram(all_bp_values, 
    "Biological Process Similarities", plot_bp)
fig.tight_layout()
fig.savefig(out_basename+"ss_values.png", bbox_inches='tight')

del fig
del axes
del all_mf_values
del all_bp_values
del all_cc_values

print("Grouping values by annotation size difference")
values_by_annotation_diff = {}
for id_pair, ss_values in all_data.items():
    
    a, b = id_pair
    annos = [mf_anno, bp_anno, cc_anno]
    for i in range(len(annos)):
        onto_anno = annos[i]
        if a in onto_anno and b in onto_anno:
            diff = abs(len(onto_anno[a])-len(onto_anno[b]))
            value = ss_values[i]

            if not diff in values_by_annotation_diff:
                values_by_annotation_diff[diff] = list()
            values_by_annotation_diff[diff].append(value)

a = []
b = []
c = []
for diff, values in values_by_annotation_diff.items():
    a.append(diff)
    not_none = [x for x in values if x]
    n_nones = len(values)-len(not_none)
    b.append((n_nones/len(values))*100.0)
    c.append(np.mean(not_none))

print("Plotting it")
fig, axes = plt.subplots(nrows=2, figsize=(6, 6))
plot_b, plot_c = axes.flatten()

plot_b.plot(a, b, 'o')
plot_b.set_title("% of None values by difference in annotation lengths")
plot_c.plot(a, c, 'o')
plot_c.set_title("Average semantic similarity by difference in annotation lengths")

fig.tight_layout()
fig.savefig(out_basename+"ss_values_by_diff.png", bbox_inches='tight')'''