import sys
import obonet
import networkx
from nexus.functional_prediction import get_ancestors, get_descendants
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

reference_path = sys.argv[1]
go_obo = sys.argv[2]
output = sys.argv[3]
annotation_paths = sys.argv[4:]
print("Annotations to compare:" + str(annotation_paths))

print("Loading GO network.")
graph = obonet.read_obo(go_obo)

association_sets = {}
gos = {}
antecesors = {}
descendants = {}

def read_id2gos(filepath):
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
    return gos_dict

reference = read_id2gos(reference_path)
reference_name = reference_path.split("/")[-1]

def get_lower_terms(go_id):
    if not go_id in descendants:
        descendants[go_id] = get_ancestors(graph, go_id)
    return descendants[go_id]

def get_upper_terms(go_id):
    if not go_id in antecesors:
        antecesors[go_id] = get_descendants(graph, go_id)
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

    print(name+":")
    print("\tGene->GO associations: " + str(associations))
    print("\tGenes annotated:" + str(genes))
    print("\tGOs annotated:" + str(len(gos[name])))
    keys_list = []
    for key in len_dict.keys():
        keys_list.append(key)
    #print(str(len_dict))
    #print(str(keys_list))
    keys_list.sort()
    strs = [str(key)+": "+str(len_dict[key]) for key in keys_list]
    print("\tNumber of annotations per gene:\n"
        +", ".join(strs))
    return 

def analyze_annotation(p):
    id2gos = read_id2gos(p)
    file_name = p.split("/")[-1]
    
    overview(id2gos, file_name)
    
    print("Calculating similarity to reference...")
    
    total = len(association_sets[file_name])
    in_reference = 0
    part_of_reference = 0
    upper_to_reference = 0
    unrelated = 0
    
    for association in tqdm(association_sets[file_name]):
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

    in_reference_perc       = (float(in_reference) / total) * 100
    part_of_reference_perc  = (float(part_of_reference) / total) * 100
    upper_to_reference_perc = (float(upper_to_reference) / total) * 100
    unrelated_perc          = (float(unrelated) / total) * 100
    
    print("Total associations: " + str(total)
            +"\n\tIn ref: "+str(in_reference_perc)
            +"\n\tPart of: "+str(part_of_reference_perc)
            +"\n\tUpper to: "+str(upper_to_reference_perc)
            +"\n\tNot in: "+str(unrelated_perc))
    
    return (file_name,[in_reference,part_of_reference,upper_to_reference,unrelated])

overview(reference, reference_path.split("/")[-1])

labels = []
counts = []
part_labels = ["In reference", "Part of reference", "Reference is part of", "Unknown"]

for annotation_path in annotation_paths:
    print(annotation_path)
    name, parts = analyze_annotation(annotation_path)
    labels.append(name)
    counts.append(parts)

in_reference = [x[0] for x in counts]
part_of_reference = [x[1] for x in counts]
upper_to_reference = [x[2] for x in counts]
unrelated = [x[3] for x in counts]

print("Plotting")
print(str(counts))
print(str(labels))
N = len(labels)
ind = np.arange(N)
width = 0.7

def sum(x,y):
    return [x[i]+y[i] for i in range(len(x))]

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

print("Formating")
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
plt.xlabel('Number of Associations')
plt.title('Similarity to Reference Annotation')
plt.yticks(ind, labels)
#plt.yticks(np.arange(0, 100, 10))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend((bars1,bars2,bars3,bars4), part_labels,loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output,bbox_inches='tight')