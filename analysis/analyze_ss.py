import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os

ss_path = sys.argv[1]
corr_file_paths = sys.argv[2].split(",")
gaf_path = sys.argv[3]
output_dir = sys.argv[4]

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

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

print("Reading annotation")
all_anno, mf_anno, bp_anno, cc_anno, names_dict, names = read_gaf(gaf_path)
print("Allocating space for data")
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
fig.savefig(out_basename+"ss_values_by_diff.png", bbox_inches='tight')