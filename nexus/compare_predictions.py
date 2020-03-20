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
confidence_file = sys.argv[3]
aspect = sys.argv[4].split(",")
output = sys.argv[5]
annotation_path_dir = [sys.argv[6]]

def load_confidences(confidence_file):
    conf_by_method = {}
    with open(confidence_file, 'r') as stream:
        lines = [line.rstrip("\n").split(",") for line in stream.readlines()[2:]]
        for line in lines:
            m = line[0]
            ths = [float(x) for x in line[1:]]
            conf_by_method[m] = ths
    return conf_by_method

confidence_by_method = load_confidences(confidence_file)

#print(str(confidence_by_method))
def get_annotations(dir, name_part):
    files = []
    # r=root, d=directories, f = files
    for f in os.listdir(dir):
        if name_part in f:
                files.append(os.path.join(dir, f))
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

path_cache = {}
def path_length(s, t):
    if not (s,t) in path_cache:
        try:
            path_len = len(networkx.shortest_path(graph,
                                          source=s,
                                          target=t))
            path_cache[(s,t)] = path_len
        except networkx.NetworkXNoPath:
            path_cache[(s,t)] = None
    return path_cache[(s,t)]

def path_value(s,t):
    if s == t:
        return 1.0
    else:
        path_len = path_length(s,t)
        if path_len != None:
            return 1.0/pow(float(path_len),3)
        else:
            return 0.0

def term_value(ref_terms, new_term):
    term_values = [path_value(new_term, ref_term) 
                   for ref_term in ref_terms]
    return max(term_values)

def read_id2gos(filepath,ontology_type=["ALL"]):
    if ontology_type[0] == "ALL":
        ontology_type = ["molecular_function",
                        "biological_process",
                        "cellular_component"]
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

def confidence_closer_to_threshold(threshold, confidences):
    diffs = []
    for conf in confidences:
        diffs.append(abs(threshold-conf))
    min_i = 0
    for i in range(len(diffs)):
        if diffs[i] < diffs[min_i]:
            min_i = i
    return min_i

def get_label_and_confidence(filename):
    methods_str = filename.split(".")[0]
    confs = confidence_by_method[methods_str]
    if ".th=0." in filename:
        th_str = "0." + filename.split(".th=0.")[1].split(".")[0]
        return methods_str,confidence_closer_to_threshold(float(th_str),
            confs)
    elif ".th0." in filename:
        th_str = "0." + filename.split(".th0.")[1].split(".")[0]
        return methods_str,confidence_closer_to_threshold(float(th_str),
            confs)
    elif ".c" in filename:
        c_str = filename.split(".c")[1].split(".")[0]
        return methods_str,int(c_str)
    else:
        print("No threshold for " + filename)
        methods_str, None

def analyze_annotation(p):
    id2gos = read_id2gos(p, ontology_type=aspect)
    #print(str(len(id2gos)) + " annotated genes in " + p)
    file_name = p.split("/")[-1]
    lens_dict = overview(id2gos, file_name)
    
    association_found = 0
    association_children_found = 0
    association_parent_found = 0
    association_missing = len(association_sets[reference_name])
    for association in association_sets[reference_name]:
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
    
    not_found = 0.0
    association_values = []
    for association in association_sets[file_name]:
        could_not_find = True
        name = association[0]
        term = association[1]
        terms_in_reference = set()
        if name in reference:
            terms_in_reference = reference[name]
            try:
                association_values.append(term_value(terms_in_reference,term))
                could_not_find = False
            except networkx.exception.NodeNotFound:
                could_not_find = True
                association_values.append(0)
        else:
            association_values.append(0)
        if could_not_find:
            not_found += 1
        '''if (name,term) in association_sets[reference_name]:
            in_reference += 1
        elif is_part_of(term,terms_in_reference):
            part_of_reference += 1
        elif is_upper_to(term,terms_in_reference):
            upper_to_reference += 1
        else:
            unrelated += 1'''
    not_found = not_found / float(total)
    #print("For " + file_name + ": " + str(not_found*100) + "% go_pairs not found.")
                            
    if total > 0:
        in_reference_perc       = (float(in_reference) / total) * 100
        part_of_reference_perc  = (float(part_of_reference) / total) * 100
        upper_to_reference_perc = (float(upper_to_reference) / total) * 100
        unrelated_perc          = (float(unrelated) / total) * 100
        '''print("Total associations: " + str(total)
                +"\n\tIn ref: "+str(in_reference_perc)
                +"\n\tPart of: "+str(part_of_reference_perc)
                +"\n\tUpper to: "+str(upper_to_reference_perc)
                +"\n\tNot in: "+str(unrelated_perc)
                +"\n\tAssociation found: "+str(association_found)
                +"\n\tAssociation children found: "+str(association_children_found)
                +"\n\tAssociation parent found: "+str(association_parent_found)
        )'''
        reference_total = len(association_sets[reference_name])
        association_missing_w = 1-(association_missing/reference_total)
        '''association_children_w = 1-(association_children_found/reference_total)
        association_parent_w = 1-(association_parent_found/reference_total)
        completion = ((association_found 
                    + association_children_found*association_children_w 
                    + association_parent_found*pow(association_parent_w,2)) 
                        / reference_total)
        
        unrelated_w = 1-(unrelated/total)
        prediction_power = ((in_reference + part_of_reference*0.5) 
                        / reference_total)'''
        prediction_quality = np.sum(association_values)
        completion = ((association_found)/ reference_total)
        label, conf_level = get_label_and_confidence(file_name.split("/")[-1])
        if conf_level != None:
            return (label, conf_level, prediction_quality, completion)
    else:
        print(file_name + " is empty")
    return None

overview(reference, reference_name)

labels = []
conf_levels = []
prediction_powers = []
completions = []

results = []
for annotation_path in tqdm(annotation_paths):
    vals = analyze_annotation(annotation_path)
    if vals != None:
        results.append(vals)
results.sort(key=lambda result: result[1])

for result in results:
    label, conf_level, prediction_power, completion = result
    labels.append(label)
    conf_levels.append(conf_level)
    prediction_powers.append(prediction_power)
    completions.append(completion)

methods = {}
for i in range(len(labels)):
    methods[labels[i]] = []
for i in range(len(labels)):
    label = labels[i]
    conf_level = conf_levels[i]
    prediction_power = prediction_powers[i]
    completion = completions[i]
    methods[label].append((conf_level,prediction_power,completion))

vecs = []
for method in methods:
    pairs = methods[method]
    pairs.sort(key=lambda x: x[0])
    ts = [t for t,r,o in pairs]
    rs = [r for t,r,o in pairs]
    os = [o for t,r,o in pairs]
    vecs.append([method,ts,rs,os])

def get_method_color(method):
    colors = {"PRS":'#ff0000',"MIC":'#dd0074',"DC":'#006600',
            "SPR":'#000099',"SOB":'#000000',"FSH":'#00ff00'}
    for method_name, color in colors.items():
        if method_name in method:
            return color

print("Prediction quality")
major_locator = 1
minor_locator = 1
size = (8,4)
plt.rcParams.update({'font.size': 13})
f, ax = plt.subplots(figsize=size)
plt.style.use('seaborn-whitegrid')
plt.title('Poder de Predição das Métricas')
for method, conf_levels, prediction_powers, completions in vecs:
    ax.plot(conf_levels, prediction_powers, label=method,
            color=get_method_color(method),linewidth=3)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_locator))
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_locator))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Poder de Predição')
plt.xlabel('Nível de confiança')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+"prediction_power.png",bbox_inches='tight')

print("Plotting completion measure")

f, ax = plt.subplots(figsize=size)
plt.style.use('seaborn-whitegrid')
plt.title('Completude das Predições Funcionais')
for method, conf_levels, prediction_powers, completions in vecs:
    ax.plot(conf_levels, completions, label=method,
            color=get_method_color(method),linewidth=3)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_locator))
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_locator))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Completude')
plt.xlabel('Nível de confiança')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+"completude.png",bbox_inches='tight')