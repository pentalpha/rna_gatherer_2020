import sys
import obonet
import networkx
from gatherer.functional_prediction import get_ancestors, get_descendants
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import statistics as stc

reference_path = sys.argv[1]
go_obo = sys.argv[2]
output = sys.argv[3]
annotation_path_dir = [sys.argv[4]]
if len(sys.argv) >= 6:
    dirs = sys.argv[5:]
else:
    dirs = []

def get_annotations(dir, name_part):
    files = []
    # r=root, d=directories, f = files
    for f in os.listdir(dir):
        if name_part in f:
                files.append(os.path.join(dir, f))
    return files

def getFilesWith(directory, name_part):
    files = []
    # r=root, d=directories, f = files
    for f in os.listdir(directory):
        if name_part in f:
                files.append(os.path.join(directory, f))
    return files

annotation_paths = []
for path in annotation_path_dir:
    annotation_paths += get_annotations(path, ".tsv")
#print("Annotations to compare:" + str(annotation_paths))

print("Loading GO network.")
graph = obonet.read_obo(go_obo)

association_sets = {}
gos = {}
antecesors = {}
descendants = {}

id_2_ontology = {}
def plot_scores(data_points):
    print("Sorting")
    data_points.sort(key=lambda x: stc.mean(x[1]))
    print("Extracting data")
    values = [vals for name,vals in data_points]
    names = [name for name,vals in data_points]
    sns.set(style="whitegrid")
    sns.plotting_context("talk", font_scale=1.5)
    sns.set_context("talk")
    print("Making plot")
    f, ax = plt.subplots(figsize=(19, 40))
    #ax.set_xscale("log")
    ax = sns.boxplot(data=values, orient='h', palette="vlag")
    #ax = seaborn.swarmplot(data=values, color=".2")
    ax.set(yticklabels=names)
    axes = ax.axes
    axes.set_xlim(0.0,1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
    #ax = (ax.set_axis_labels("Z-Scores")).set(xlim=(-1,1))
    plt.title("Similarities")
    #plt.show(ax)
    print("Saving")
    plt.savefig(output+"similarities.png",bbox_inches='tight')
    #plt.savefig("similarities.png")


def read_data(path):
    data = []
    with open(path,'r') as stream:
        for line in stream.readlines():
            line = line.rstrip("\n")
            if not ("NA" in line):
                data.append(float(line))
            #else:
            #    data.append(0.0)
    return data

def file_name(path):
    return ".".join(os.path.basename(path).split(".")[:-1])

def read_id2gos(filepath,ontology_type="All"):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            if not rfam_id in gos_dict:
                gos_dict[rfam_id] = set()
            go_ids = gos_str.split(";")
            for go in go_ids:
                gos_dict[rfam_id].add(go)
                if ontology_type != "All":
                    id_2_ontology[go] = ontology_type
                else:
                    if len(cols) > 2:
                        ontology = cols[2]
                        id_2_ontology[go] = ontology
    return gos_dict

reference = read_id2gos(reference_path)
reference_name = reference_path.split("/")[-1]

def get_lower_terms(go_id):
    if not go_id in descendants:
        if go_id in graph:
            descendants[go_id] = get_ancestors(graph, go_id)
        else:
            descendants[go_id] = set()
    return descendants[go_id]

def get_upper_terms(go_id):
    if not go_id in antecesors:
        if go_id in graph:
            antecesors[go_id] = get_descendants(graph, go_id)
        else:
            antecesors[go_id] = set()
    return antecesors[go_id]

def is_part_of(go_id, ref_gos):
    for ref_go in ref_gos:
        if go_id in get_lower_terms(ref_go):
            return True
    return False

def is_upper_to(go_id, ref_gos):
    for ref_go in ref_gos:
        if go_id in get_upper_terms(ref_go):
            return True
    return False

def overview(id2gos, name):
    association_sets[name] = set()
    associations = 0
    genes = len(id2gos.keys())
    gos[name] = set()

    len_dict = {}

    for key in id2gos.keys():
        terms = id2gos[key]
        for term in terms:
            associations += 1
            association_sets[name].add((key,term))
            gos[name].add(term)
        
        value = id2gos[key]
        l = len(value)
        if not l in len_dict:
            len_dict[l] = 1
        else:
            len_dict[l] += 1

    '''print(name+":")
    print("\tGene->GO associations: " + str(associations))
    print("\tGenes annotated:" + str(genes))
    print("\tGOs annotated:" + str(len(gos[name])))'''
    keys_list = []
    for key in len_dict.keys():
        keys_list.append(key)
    #print(str(len_dict))
    #print(str(keys_list))
    keys_list.sort()
    strs = [str(key)+": "+str(len_dict[key]) for key in keys_list]
    #print("\tNumber of annotations per gene:\n"
    #    +", ".join(strs))
    return len_dict

def analyze_annotation(p):
    if "CC" in p or "cellular_component" in p:
        ont = "cellular_component"
    elif "MF" in p or "molecular_function" in p:
        ont = "molecular_function"
    else:
        ont = "biological_process"

    id2gos = read_id2gos(p,ontology_type=ont)
    #print("This one is a " + ont + " annotation")
    file_name_parts = p.split("/")
    file_name = file_name_parts[-1].replace(".LNC","").replace(".tsv", "")
    if len(file_name_parts) > 1:
        file_name = file_name_parts[-2] + "/" + file_name
    
    lens_dict = overview(id2gos, file_name)
    
    #print("Calculating similarity to reference...")
    
    total = len(association_sets[file_name])
    in_reference = 0
    part_of_reference = 0
    upper_to_reference = 0
    unrelated = 0
    
    for association in association_sets[file_name]:
        name = association[0]
        term = association[1]
        terms_in_reference = set()
        if name in reference:
            terms_in_reference = reference[name]

        if (name,term) in association_sets[reference_name]:
            in_reference += 1
        elif is_part_of(term,terms_in_reference):
            part_of_reference += 1
        elif is_upper_to(term,terms_in_reference):
            upper_to_reference += 1
        else:
            unrelated += 1

    reference_total = 0
    reference_found = 0
    partially_found = 0
    completelly_found = 0
    for name in reference.keys():
        if name in id2gos:
            reference_set = reference[name]
            reference_set_ont = set()
            contains_right_ontology = False
            for go_id in reference_set:
                if id_2_ontology[go_id] == ont:
                    reference_set_ont.add(go_id)
            if len(reference_set_ont) > 0:
                reference_total += 1
                new_set = id2gos[name]
                n_found = len(reference_set_ont.intersection(new_set))
                if n_found == reference_set_ont:
                    completelly_found += 1
                elif n_found > 0:
                    reference_found += 1
                else:
                    for go_id in reference_set_ont:
                        lower = set(get_lower_terms(go_id))
                        if len(lower.intersection(new_set)) > 1:
                            partially_found += 1
                            break
                            

    in_reference_perc       = (float(in_reference) / total) * 100
    part_of_reference_perc  = (float(part_of_reference) / total) * 100
    upper_to_reference_perc = (float(upper_to_reference) / total) * 100
    unrelated_perc          = (float(unrelated) / total) * 100
    completelly_found_perc    = (float(completelly_found) / reference_total) * 100
    reference_found_perc    = (float(reference_found) / reference_total) * 100
    partially_found_perc    = (float(partially_found) / reference_total) * 100
    
    '''print("Total associations: " + str(total)
            +"\n\tIn ref: "+str(in_reference_perc)
            +"\n\tPart of: "+str(part_of_reference_perc)
            +"\n\tUpper to: "+str(upper_to_reference_perc)
            +"\n\tNot in: "+str(unrelated_perc)
            +"\n\tComplete found: "+str(completelly_found_perc)
            +"\n\tReference found: "+str(reference_found_perc)
            +"\n\tReference found (partially): "+str(partially_found_perc))'''
    
    return (file_name,[in_reference,part_of_reference,upper_to_reference,unrelated],lens_dict,
             [completelly_found_perc,reference_found_perc,partially_found_perc],ont)

overview(reference, reference_path.split("/")[-1])

labels = []
ontology_types = []
counts = []
lens_dicts = []
ref_founds = []
part_labels = ["In reference", "Part of reference", "Reference is part of", "Unknown"]
part_labels_b = ["% of reference genes entirelly annotated", 
                "% of reference genes annotated", "% of reference genes partially annotated"]

for annotation_path in tqdm(annotation_paths):
    #print(annotation_path)
    name, parts, lens, ref_found, ont_type = analyze_annotation(annotation_path)
    ontology_types.append(ont_type)
    labels.append(name)
    counts.append(parts)
    lens_dicts.append(lens)
    ref_founds.append(ref_found)

in_reference = [x[0] for x in counts]
part_of_reference = [x[1] for x in counts]
upper_to_reference = [x[2] for x in counts]
unrelated = [x[3] for x in counts]

'''print("Plotting Similarities")
print(str(counts))
print(str(labels))'''
N = len(labels)
ind = np.arange(N)
width = 0.7

def sum(x,y):
    return [x[i]+y[i] for i in range(len(x))]

on_reference = sum(part_of_reference, in_reference)

f, ax = plt.subplots(figsize=(16, 1.2*N + 1))

bars1 = ax.barh(ind, in_reference, width,
                color="#0400ff")
left = in_reference
bars2 = ax.barh(ind, part_of_reference, width, 
                left=left,color="#8052ff")
left = sum(left,part_of_reference)
bars3 = ax.barh(ind, upper_to_reference, width, 
                left=left,color="#ff61f7")
left = sum(left,upper_to_reference)
bars4 = ax.barh(ind, unrelated, width, 
                left=left,color="#ff2b59")

#print("Formating")
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
plt.xlabel('Number of Associations')
plt.title('Similarity to Reference Annotation')
plt.yticks(ind, labels)
#plt.yticks(np.arange(0, 100, 10))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend((bars1,bars2,bars3,bars4), part_labels,loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".1.png",bbox_inches='tight')

print("Plotting reference found")
#complete_found = [x[0] for x in ref_founds]
reference_found = [x[1] for x in ref_founds]
partially_found = [x[2] for x in ref_founds]

f, ax = plt.subplots(figsize=(16, 1.2*N + 1))

'''barsb1 = ax.barh(ind, complete_found, width,
                color="#0000a0")
left = reference_found'''
ax.set_axisbelow(True)
plt.grid(axis='x')
barsb1 = ax.barh(ind, reference_found, width,
                color="#0400ff")
left = reference_found
barsb2 = ax.barh(ind, partially_found, width, 
                left=left,color="#8052ff")
plt.xlabel('%')
plt.title('Reference Annotations Recovered')
plt.yticks(ind, labels)
ax.xaxis.set_major_locator(ticker.MultipleLocator(3))
#plt.xticks(np.arange(0, 100, 10))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend((barsb1,barsb2), part_labels_b[1:],loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".2.png",bbox_inches='tight')

def short_method_name(name):
    name = name.replace("molecular_function","MF")
    name = name.replace("biological_process","BP")
    name = name.replace("cellular_component","CC")
    name = name.replace("-pval=0.05-fdr=0.05", "-pval=fdr=0.05")
    name = name.replace("-pval=0.01-fdr=0.01", "-pval=fdr=0.01")
    return name

##########
methods = {}
for i in range(len(labels)):
    label = labels[i]
    ont_type = ontology_types[i]
    label_parts = label.split("/")[-1].split(".")
    method = ont_type+"-"+label_parts[0]+"-pval=0."+label_parts[4]+"-fdr=0."+label_parts[6]
    method = short_method_name(method)
    if not method in methods:
        methods[method] = []
    threshold_str = "0." + label.split(".th0.")[-1].split(".")[0]
    threshold_val = float(threshold_str)
    ref = reference_found[i]
    methods[method].append((threshold_val,ref,on_reference[i],unrelated[i]))

vecs = []
for method in methods:
    pairs = methods[method]
    pairs.sort(key=lambda x: x[0])
    ts = [t for t,r,o,u in pairs]
    rs = [r for t,r,o,u in pairs]
    os = [o for t,r,o,u in pairs]
    us = [u for t,r,o,u in pairs]
    vecs.append([method,ts,rs,os,us])

def get_method_color(method):
    colors = {"PRS":'r',"MIC":'b',"DC":'g',"SPR":'k'}
    for method_name, color in colors.items():
        if method_name in method:
            return color

#print(vecs)
size = (8,4)
plt.rcParams.update({'font.size': 13})
f, ax = plt.subplots(figsize=size)
plt.style.use('seaborn-whitegrid')
plt.title('Reference Annotations Recovered 2')
for method, thresholds, reference_perc, on_reference, unknown in vecs:
    ax.plot(thresholds, reference_perc,label=method,color=get_method_color(method))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('% of the reference genes annotated')
plt.xlabel('Correlation coefficient threshold')
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".3.png",bbox_inches='tight')



f, ax = plt.subplots(figsize=size)
plt.style.use('seaborn-whitegrid')
plt.title('Prediction size over different thresholds')
methods_to_print = list(set([m for m,t,r,o,u in vecs]))
colors = ["#ebffa2","#c99800","#8d0000","#400040","008080"]
method_color = {}
for i in range(len(methods_to_print)):
    method_color[methods_to_print[i]] = colors[i]
for method, thresholds, reference_perc, on_reference, unknown in vecs:
    ax.plot(thresholds, on_reference,color=get_method_color(method),label=method+", Part of Reference")
    ax.plot(thresholds, unknown,color=get_method_color(method),label=method+", Unknown",linestyle='dashed')
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Number of associations')
plt.xlabel('Correlation coefficient threshold')
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".4.png",bbox_inches='tight')


if len(dirs) > 0:
    print("Plotting semantic similarities")
    print("Reading files")
    data = []
    for d in dirs:
        for f in getFilesWith(d, ".tsv"):
            d = read_data(f)
            name = f.split("/")[-1].replace("-Wang-BMA", "-") + " ("+str(len(d))+")"
            data.append([name,d])

    plot_scores(data)