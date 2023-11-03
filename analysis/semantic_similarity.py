import networkx as nx
from tqdm import tqdm
import math
import numpy as np

computational_evi = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']

def load_full_gaf(filepath: str, go_tocorrect):
    print('Reading', filepath)
    stream = open(filepath, 'r')
    ann = {'MF': {}, 'BP': {}, 'CC': {}}
    
    for rawline in stream:
        if rawline[0] != '!':
            cells = rawline.rstrip('\n').split('\t')
            db = cells[0]
            name = cells[1]

            goid = cells[4]
            if goid in go_tocorrect:
                goid = go_tocorrect[goid]
            ont = cells[8].replace('F', 'MF').replace('C', 'CC').replace('P', 'BP')
            if not goid in ann[ont]:
                ann[ont][goid] = set()
            ann[ont][goid].add(name)
    
    return ann

def count_freq(goid: list, ann: dict, desc_lists: dict, total_annots: dict):
    if goid in desc_lists:
        descendant_goids = desc_lists[goid]
        annots = set()
        annots.update(ann[goid])

        for descendant_goid in descendant_goids:
            if not descendant_goid in total_annots:
                count_freq(descendant_goid, ann, desc_lists, total_annots)
            annots.update(total_annots[descendant_goid])
        total_annots[goid] = annots
    else:
        #go id not in annotation, zero frequency
        total_annots[goid] = set()

def calc_IC(ann, graph: nx.DiGraph):
    print('Calculating IC of', len(ann.keys()), 'GO terms')
    total_ann = sum(len(genes) for goid, genes in ann.items())

    print('Getting descendants lists')
    desc_lists = {goid: [anc for anc in nx.ancestors(graph, goid)] 
                  for goid, genes in ann.items()}
    #GOIDs with least descendants first
    all_goids = list(desc_lists.keys())
    all_goids.sort(key = lambda goid: len(desc_lists[goid]))

    print('Calculating expanded frequencies')
    total_annots = {}
    for goid in tqdm(all_goids):
        count_freq(goid, ann, desc_lists, total_annots)
    print("total_annots: ")
    all_freq = [len(x) for x in total_annots.values()]
    print('Min:', min(all_freq), 'Max:', max(all_freq), 'Mean:', np.mean(all_freq))

    
    relative_freqs = {goid: len(x)/total_ann for goid, x in total_annots.items()}
    print("relative_freqs: ")
    all_freq = [x for x in relative_freqs.values()]
    print('Min:', min(all_freq), 'Max:', max(all_freq), 'Mean:', np.mean(all_freq))

    max_log = math.log2(total_ann)
    ics = {goid: math.log2(1/freq)/max_log for goid, freq in relative_freqs.items() if freq > 0}
    print("ics: ")
    all_freq = [x for x in ics.values()]
    print('Min:', min(all_freq), 'Max:', max(all_freq), 'Mean:', np.mean(all_freq))

    return ics

def calc_IC_indexes(filepath: str, go_tocorrect, graph):
    ontologies = load_full_gaf(filepath, go_tocorrect)
    ic_indexes = {}
    for ontoname, ann in ontologies.items():
        new_ics = calc_IC(ann, graph)
        ic_indexes[ontoname] = {goid: ic for goid, ic in new_ics.items()}
    #quit()
    return ic_indexes

def simgic(goids_a, goids_b, ontology_ic: dict):
    a = set([x for x in goids_a if x in ontology_ic])
    b = set([x for x in goids_b if x in ontology_ic])

    go_union = a.union(b)
    go_intersect = a.intersection(b)

    union_sum = sum([ontology_ic[goid] for goid in go_union])
    intersect_sum = sum([ontology_ic[goid] for goid in go_intersect])

    return intersect_sum / union_sum

