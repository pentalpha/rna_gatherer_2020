import subprocess
import datetime
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
import shutil

now = datetime.datetime.now()
#mode = conf_mode or th_mode
go_obo = sys.argv[1]
mode = "conf_mode"
if len(sys.argv) > 2:
    mode = sys.argv[2]

def runCommand(cmd, print_cmd=True):
        if print_cmd:
                print("\t> " + cmd)
        process = subprocess.call(cmd, shell=True)
        return process

def get_annotations(dir, name_part):
    files = []
    # r=root, d=directories, f = files
    for f in os.listdir(dir):
        if name_part in f:
                files.append(os.path.join(dir, f))
    return files

def short_method_name(name):
    name = name.replace("molecular_function","MF")
    name = name.replace("biological_process","BP")
    name = name.replace("cellular_component","CC")
    name = name.replace("-pval=0.05-fdr=0.05", "-pval=fdr=0.05")
    name = name.replace("-pval=0.01-fdr=0.01", "-pval=fdr=0.01")
    return name

predictor = "../predict.py"
counts_file = "../test_data/counts/mus_musculus_tpm.tsv"
annotation_file = "../test_data/annotation/mgi_genes_annotation.tsv"
reference_file = "../test_data/reference/mgi_lncRNA_annotation.tsv"
#reg_file = "../test_data/lnc_list/mus_musculus_lncRNA.txt"
reg_file = "../test_data/lnc_list/mus_musculus_lncRNA.txt"

onto_types = ["molecular_function"]
mic_thresholds = [0.7435552598651071,
                0.7435552598651071,0.7511610977726769,
                0.7511610977726769,0.7511610977726769,
                0.9980008838722987,0.9980008838722987]
mic_levels = [0,1,2,3,4,5,6]

test_name = (mode+"-"+now.strftime("%y")
            +"_"+now.strftime("%j")
            +"_"+now.strftime("%H")
            +"."+now.strftime("%M")
            +"."+now.strftime("%S"))

if os.path.exists(test_name):
    shutil.rmtree(test_name)

#test_name = "conf_mode-20_089_06.50.34"
ths_str = "-conf " + ",".join([str(x) for x in mic_levels])
if mode == "th_mode":
    ths_str = "-th " + ",".join([str(x) for x in mic_thresholds])

cmd = " ".join(["python",predictor,"-cr",counts_file,
                "-reg",reg_file,"-ann",annotation_file,
                "-o",test_name,ths_str,"-ont",",".join(onto_types),
                "-met","MIC","-p","6"])

code = runCommand(cmd)
#code = 0
if code != 0:
    print("The test " + test_name + " failed with exit code "+str(code))
    quit()

annotation_paths = []
for onto_type in onto_types:
    annotation_paths += get_annotations(
                        test_name + "/" + short_method_name(onto_type) + "/",
                        ".tsv")

print("Loading GO network.")
graph = obonet.read_obo(go_obo)

association_sets = {}
gos = {}
antecesors = {}
descendants = {}
id_2_ontology = {}

def file_name(path):
    return ".".join(os.path.basename(path).split(".")[:-1])

def read_id2gos(filepath, name, ontology_type=["molecular_function",
                                        "biological_process",
                                        "cellular_component"]):
    association_sets[name] = set()
    gos[name] = set()
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            gene_id = cols[0]
            go_str = cols[1]
            ont_type = cols[2]
            if ont_type in ontology_type:
                if not gene_id in gos_dict:
                    gos_dict[gene_id] = set()
                gos_dict[gene_id].add(go_str)
                id_2_ontology[go_str] = ont_type

    for key in gos_dict.keys():
        terms = gos_dict[key]
        for term in terms:
            association_sets[name].add((key,term))
            gos[name].add(term)

    return gos_dict

reference_name = reference_file.split("/")[-1]
reference = read_id2gos(reference_file, reference_name,
                        ontology_type = onto_types)


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

def analyze_annotation(p):
    
    #print("This one is a " + ont + " annotation")
    file_name_parts = p.split("/")
    file_name = file_name_parts[-1].replace(".LNC","").replace(".tsv", "")
    if len(file_name_parts) > 1:
        file_name = file_name_parts[-2] + "/" + file_name
    id2gos = read_id2gos(p, file_name, ontology_type=onto_types)
    
    #print("Calculating similarity to reference...")
    
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
        
        return (file_name.split("/")[-1],[in_reference,part_of_reference,upper_to_reference,unrelated],
                [completelly_found_perc,reference_found_perc,partially_found_perc],
                [association_found,association_children_found,association_parent_found,association_missing])
    else:
        return (file_name.split("/")[-1],[0,0,0,0],
                [0,0,0],
                [association_found,association_children_found,association_parent_found,association_missing])

labels = []
counts = []
ref_founds = []
ref_association_search = []
part_labels_c = ["Recuperada", "Termo descendente recuperado", 
                "Termo antepassado recuperado", "Não recuperada"]

results = []
for annotation_path in tqdm(annotation_paths):
    results.append((analyze_annotation(annotation_path), annotation_path))
results.sort(key=lambda result: result[0][-1][-1])

for result, annotation_path in results:
    name, parts, ref_found, ref_associations = result
    labels.append(name)
    counts.append(parts)
    ref_founds.append(ref_found)
    ref_association_search.append(ref_associations)

table_str = "\t".join(["prediction","association_found","association_children_found",
                        "association_parent_found","association_missing"]) + "\n"
for i in range(len(labels)):
    table_str += "\t".join([labels[i],str(ref_association_search[i][0]),
                            str(ref_association_search[i][1]),
                            str(ref_association_search[i][2]),
                            str(ref_association_search[i][3])]) + "\n"

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
plt.title('Associações de Referência Recuperadas' + " ("+",".join(onto_types)+")")
plt.yticks(ind, labels)
#plt.yticks(np.arange(0, 100, 10))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(test_name+"/associations.png",bbox_inches='tight')
with open(test_name+"/associations.tsv",'w') as stream:
    stream.write(table_str)