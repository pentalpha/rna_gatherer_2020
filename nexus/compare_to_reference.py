import sys
import obonet
import networkx
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import statistics as stc

reference_path = sys.argv[1]
go_obo = sys.argv[2]
aspect = sys.argv[3].split(",")
output = sys.argv[4]
annotation_path_dir = [sys.argv[5]]
if len(sys.argv) >= 7:
    dirs = sys.argv[6:]
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
    f, ax = plt.subplots(figsize=(14, (len(values)*3) + 2))
    #ax.set_xscale("log")
    ax = sns.boxplot(data=values, orient='h', palette="vlag")
    #ax = seaborn.swarmplot(data=values, color=".2")
    ax.set(yticklabels=names)
    axes = ax.axes
    axes.set_xlim(0.0,1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    #ax = (ax.set_axis_labels("Z-Scores")).set(xlim=(-1,1))
    plt.title("Semantic Similarity to Reference Annotation")
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

def read_id2gos(filepath,ontology_type=["molecular_function",
                                        "biological_process",
                                        "cellular_component"]):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            aspect = cols[2]
            if aspect in ontology_type:
                if not rfam_id in gos_dict:
                    gos_dict[rfam_id] = set()
                go_ids = gos_str.split(";")
                for go in go_ids:
                    gos_dict[rfam_id].add(go)
                    id_2_ontology[go] = aspect
    return gos_dict

reference = read_id2gos(reference_path, ontology_type = aspect)
reference_name = reference_path.split("/")[-1]


def get_ancestors(graph, parent_id):
    ancestors = networkx.ancestors(graph, parent_id)
    return ancestors

def get_descendants(graph, parent_id):
    descendants = networkx.descendants(graph, parent_id)
    return descendants

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
    id2gos = read_id2gos(p, ontology_type=aspect)
    #print("This one is a " + ont + " annotation")
    file_name_parts = p.split("/")
    file_name = file_name_parts[-1].replace(".LNC","").replace(".tsv", "")
    if len(file_name_parts) > 1:
        file_name = file_name_parts[-2] + "/" + file_name
    
    lens_dict = overview(id2gos, file_name)
    
    #print("Calculating similarity to reference...")
    
    association_found = 0
    association_children_found = 0
    association_parent_found = 0
    association_missing = len(association_sets[reference_name])
    for association in tqdm(association_sets[reference_name]):
        if association in association_sets[file_name]:
            association_found += 1
            association_missing -= 1
        else:
            name = association[0]
            term = association[1]
            matching_predictions = set()
            if name in id2gos:
                matching_predictions.update(id2gos[name])
            '''for pred_name,pred_term in association_sets[file_name]:
                if pred_name == name:
                    matching_predictions.add(pred_term)'''
            if len(matching_predictions) > 0:
                if is_upper_to(term,matching_predictions):
                    association_children_found += 1
                    association_missing -= 1
                elif is_part_of(term,matching_predictions):
                    association_parent_found += 1
                    association_missing -= 1
        

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
            reference_total += 1
            new_set = id2gos[name]
            n_found = len(reference_set.intersection(new_set))
            if n_found == reference_set:
                completelly_found += 1
            elif n_found > 0:
                reference_found += 1
            else:
                for go_id in reference_set:
                    lower = set(get_lower_terms(go_id))
                    if len(lower.intersection(new_set)) > 1:
                        partially_found += 1
                        break
                            
    if total > 0:
        in_reference_perc       = (float(in_reference) / total) * 100
        part_of_reference_perc  = (float(part_of_reference) / total) * 100
        upper_to_reference_perc = (float(upper_to_reference) / total) * 100
        unrelated_perc          = (float(unrelated) / total) * 100
        completelly_found_perc    = (float(completelly_found) / reference_total) * 100
        reference_found_perc    = (float(reference_found) / reference_total) * 100
        partially_found_perc    = (float(partially_found) / reference_total) * 100
        print("Total associations: " + str(total)
                +"\n\tIn ref: "+str(in_reference_perc)
                +"\n\tPart of: "+str(part_of_reference_perc)
                +"\n\tUpper to: "+str(upper_to_reference_perc)
                +"\n\tNot in: "+str(unrelated_perc)
                +"\n\tComplete found: "+str(completelly_found_perc)
                +"\n\tReference found: "+str(reference_found_perc)
                +"\n\tReference found (partially): "+str(partially_found_perc)
                +"\n\tAssociation found: "+str(association_found)
                +"\n\tAssociation children found: "+str(association_children_found)
                +"\n\tAssociation parent found: "+str(association_parent_found)
        )
        
        return (file_name.split("/")[-1],[in_reference,part_of_reference,upper_to_reference,unrelated],lens_dict,
                [completelly_found_perc,reference_found_perc,partially_found_perc],
                [association_found,association_children_found,association_parent_found,association_missing])
    else:
        return (file_name.split("/")[-1],[0,0,0,0],lens_dict,
                [0,0,0],
                [association_found,association_children_found,association_parent_found,association_missing])

overview(reference, reference_name)

labels = []
counts = []
lens_dicts = []
ref_founds = []
ref_association_search = []
#part_labels = ["Na referência", "Parte da referência", "Referência é parte de", "Desconhecidos"]
#part_labels_b = ["% of reference genes entirelly annotated", 
#                "% of reference genes annotated", "% of reference genes partially annotated"]
part_labels_c = ["Recuperada", "Termo descendente recuperado", "Termo antepassado recuperado", "Não recuperada"]
#part_labels_c = ["a", "b", "c", "d"]

results = []
for annotation_path in tqdm(annotation_paths):
    results.append((analyze_annotation(annotation_path), annotation_path))
results.sort(key=lambda result: result[0][-1][-1])

for result, annotation_path in results:
    name, parts, lens, ref_found, ref_associations = result
    labels.append(name)
    counts.append(parts)
    lens_dicts.append(lens)
    ref_founds.append(ref_found)
    ref_association_search.append(ref_associations)

'''print("Plotting Similarities")
print(str(counts))
print(str(labels))'''
N = len(labels)
ind = np.arange(N)
width = 0.7

def sum(x,y):
    return [x[i]+y[i] for i in range(len(x))]

association_found = [x[0] for x in ref_association_search]
association_children_found = [x[1] for x in ref_association_search]
association_parent_found = [x[2] for x in ref_association_search]
association_missing = [x[3] for x in ref_association_search]
x = 8
y = 0.8*N + 1
f, ax = plt.subplots(figsize=(x, y))

bars1 = ax.barh(ind, association_found, width,
                color="#0400ff",label=part_labels_c[0])
left = association_found
bars2 = ax.barh(ind, association_children_found, width, 
                left=left,color="#8052ff",label=part_labels_c[1])
left = sum(left,association_children_found)
bars3 = ax.barh(ind, association_parent_found, width, 
                left=left,color="#ff61f7",label=part_labels_c[2])
left = sum(left,association_parent_found)
bars4 = ax.barh(ind, association_missing, width, 
                left=left,color="#ff2b59",label=part_labels_c[3])

#print("Formating")
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
plt.xlabel('Número de Associações')
plt.title('Associações de Referência Recuperadas' + " ("+",".join(aspect)+")")
plt.yticks(ind, labels)
#plt.yticks(np.arange(0, 100, 10))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".0.png",bbox_inches='tight')
'''
in_reference = [x[0] for x in counts]
part_of_reference = [x[1] for x in counts]
upper_to_reference = [x[2] for x in counts]
unrelated = [x[3] for x in counts]

on_reference = sum(part_of_reference, in_reference)

f, ax = plt.subplots(figsize=(x, y))

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
plt.title('Similarity to Reference Annotation' + "("+",".join(aspect)+")")
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

f, ax = plt.subplots(figsize=(x, y))'''

'''barsb1 = ax.barh(ind, complete_found, width,
                color="#0000a0")
left = reference_found''''''
ax.set_axisbelow(True)
plt.grid(axis='x')
barsb1 = ax.barh(ind, reference_found, width,
                color="#0400ff")
left = reference_found
barsb2 = ax.barh(ind, partially_found, width, 
                left=left,color="#8052ff")
plt.xlabel('% of reference genes')
plt.title('Reference Annotations Recovered' + "("+",".join(aspect)+")")
plt.yticks(ind, labels)
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
#plt.xticks(np.arange(0, 100, 10))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend((barsb1,barsb2), part_labels_b[1:],loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".2.png",bbox_inches='tight')
'''
def short_method_name(name):
    name = name.replace("molecular_function","MF")
    name = name.replace("biological_process","BP")
    name = name.replace("cellular_component","CC")
    name = name.replace("-pval=0.05-fdr=0.05", "-pval=fdr=0.05")
    name = name.replace("-pval=0.01-fdr=0.01", "-pval=fdr=0.01")
    return name

##########
'''methods = {}
for i in range(len(labels)):
    label = labels[i]
    ont_type = ontology_types[i]
    label_parts = label.split("/")[-1].split(".")
    method = ont_type+"-"+label_parts[0]
    if len(label_parts) >= 5:
        method += "-pval=0."+label_parts[4]
        if len(label_parts) >= 7:
            method += "-fdr=0."+label_parts[6]
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
'''