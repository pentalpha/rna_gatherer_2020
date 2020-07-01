import sys
import obonet
import networkx
from tqdm import tqdm
import numpy as np
import os
import statistics as stc
import pandas as pd

'''
usage:
python predictions_dotplot.py [reference_path] [go_obo] [confidence_file] [aspect] [output] [annotation_path_dir]
'''

reference_path = sys.argv[1]
go_obo = sys.argv[2]
confidence_file = sys.argv[3]
aspect = sys.argv[4].split(",")
output = sys.argv[5]
annotation_path_dir = [sys.argv[6]]
external_anno = None
if len(sys.argv) >= 8:
    external_anno = sys.argv[7]

def load_confidences(intervals_file):
    with open(intervals_file, 'r') as stream:
        lines = [line.rstrip("\n").split(",") for line in stream.readlines()]
        th_lines = lines[2:]
        confidences = [{} for i in range(len(th_lines[0])-1)]
        for cells in th_lines:
            metric = cells[0]
            for i in range(1,len(cells)):
                if cells[i] == "None":
                    confidences[i-1][metric] = None
                else:
                    confidences[i-1][metric] = float(cells[i])
        return confidences
    return []

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
            if len(cols) >= 3:
                gene_id = cols[0]
                go_str = cols[1]
                aspect = cols[2]
                if aspect in ontology_type:
                    if not gene_id in gos_dict:
                        gos_dict[gene_id] = set()
                    gos_dict[gene_id].add(go_str)
                    id_2_ontology[go_str] = aspect
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
    #confs = confidence_by_method[methods_str]
    if ".c" in filename:
        c_str = filename.split(".c")[1].split(".")[0]
        return methods_str,int(c_str)
    else:
        print("No threshold for " + filename)
        return methods_str, -1

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

    for association in association_sets[file_name]:
        name = association[0]
        term = association[1]
        terms_in_reference = set()
        if name in reference:
            terms_in_reference = reference[name]
            if term in terms_in_reference:
                in_reference += 1
            elif is_part_of(term,terms_in_reference):
                part_of_reference += 1
            elif is_upper_to(term,terms_in_reference):
                upper_to_reference += 1
            else:
                unrelated += 1
        else:
            unrelated += 1
    
    reference_total = len(association_sets[reference_name])
    del association_sets[file_name]
    del gos[file_name]
    if total > 0:
        unrelated_perc = (float(unrelated) / total) * 100
        related_perc = 100.0 - unrelated_perc
        size = float(total)/float(reference_total)
        
        missing_coef = ((association_missing)/ reference_total)
        completion = ((association_found)/ reference_total)
        
        label, conf_level = get_label_and_confidence(file_name.split("/")[-1])
        number_of_metrics = len(label.split("-"))
        if conf_level != None:
            return (label, conf_level, completion, missing_coef, number_of_metrics, related_perc, size)
    else:
        print(file_name + " is empty")
    return None

overview(reference, reference_name)

results = []
if external_anno != None:
    annotation_paths = [external_anno] + annotation_paths
for annotation_path in tqdm(annotation_paths):
    vals = analyze_annotation(annotation_path)
    if vals != None:
        results.append(vals)
results.sort(key=lambda result: result[2])

id_series = pd.Series(name='index')
conf_series = pd.Series(name='confidenceLevel')
metric_comb_series = pd.Series(name='metrics')
completion_series = pd.Series(name='Completude')
missing_series = pd.Series(name='Incompletude')
has_path_series = pd.Series(name="hasPath")
size_series = pd.Series(name="size")
distance_series = pd.Series(name="distance")
quality_series = pd.Series(name="quality")
number_of_metrics_series = pd.Series(name='Número de Métricas')

def euclids(a, b):
    return np.sqrt(np.sum((b-a)**2))
top_point = np.array([0.0,1.0])

print("Calculating DF")
i = 0
prediction_labels = []
qualities_by_conf = {}
for label, conf_level, completion, missing_coef, number_of_metrics, has_path_perc, size in results:
    prediction_label = label+"_conf="+str(conf_level)
    if conf_level == -1:
        prediction_label = label
        metric_comb_series[prediction_label] = "PRS,SPR,SOB,FSH"
    else:
        metric_comb_series[prediction_label] = label
    id_series[prediction_label] = i
    conf_series[prediction_label] = conf_level
    
    
    completion_series[prediction_label] = completion*100.0
    missing_series[prediction_label] = missing_coef*100.0
    number_of_metrics_series[prediction_label] = int(number_of_metrics)
    has_path_series[prediction_label] = has_path_perc
    size_series[prediction_label] = size

    dist = euclids(np.array([completion, missing_coef]), top_point)*100.0
    distance_series[prediction_label] = dist
    quality = dist * (has_path_perc/100.0)
    quality_series[prediction_label] = quality
    if conf_level >= 0:
        if not conf_level in qualities_by_conf:
            qualities_by_conf[conf_level] = []
        qualities_by_conf[conf_level].append((prediction_label,quality))

    prediction_labels.append(prediction_label)
    i+=1

def custom_max(tuples):
    max_i = 0
    for i in range(len(tuples)):
        if tuples[i][1] >= tuples[max_i][1]:
            max_i = i
    return tuples[i][0]

bests = []
for conf_level, qualities in qualities_by_conf.items():
    bests.append(custom_max(qualities))

best_series = pd.Series(name='best')
for label in prediction_labels:
    if label in bests:
        best_series[label] = "BEST"
    else:
        best_series[label] = "NOT_BEST"

df = pd.concat([conf_series, metric_comb_series, completion_series,
                missing_series, has_path_series, size_series,
                number_of_metrics_series,
                distance_series, quality_series, best_series],
            axis=1)

df.to_csv(output + "-predictions_stats.tsv",sep="\t",header=True,index=True)
best_df = df[df['best']=='BEST']
best_df.to_csv(output + "-best_predictions.tsv",sep="\t",header=True,index=True)