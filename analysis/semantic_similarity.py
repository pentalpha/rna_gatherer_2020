import networkx as nx
from tqdm import tqdm
import math
import numpy as np
from os import path

computational_evi = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']

def load_full_gaf(filepath: str, go_tocorrect, confirmed_lnc):

    rna_central_lncrnas_path = path.join(path.dirname(filepath), 'rna_central_lncrna_list.txt')
    rna_central_lncrnas = open(rna_central_lncrnas_path,'r').read().split('\n')
    rna_central_lncrnas = set(rna_central_lncrnas)

    print('Reading', filepath)
    stream = open(filepath, 'r')
    ann = {'MF': {}, 'BP': {}, 'CC': {}}
    
    all_gos = set()
    for rawline in stream:
        if rawline[0] != '!':
            cells = rawline.rstrip('\n').split('\t')
            db = cells[0]
            name = cells[1]

            if db in ['MGI', 'RNAcentral']:
                is_nc_rna = False
                if db == 'RNAcentral':
                    #islncRNA = 'lnc' in cells[2] or 'long' in cells[2]
                    is_nc_rna = name in rna_central_lncrnas
                elif db == 'MGI':
                    is_nc_rna = 'lnc' in cells[-6]

                if not is_nc_rna and name in confirmed_lnc:
                    is_nc_rna = True
                
                if is_nc_rna:
                    goid = cells[4]
                    if goid in go_tocorrect:
                        goid = go_tocorrect[goid]
                    ont = cells[8].replace('F', 'MF').replace('C', 'CC').replace('P', 'BP')
                    if not goid in ann[ont]:
                        ann[ont][goid] = set()
                        all_gos.add(goid)
                    ann[ont][goid].add(name)

    return ann, all_gos

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

def calc_IC_indexes(filepath: str, go_tocorrect, graph, confirmed_lnc):
    ontologies, all_gos = load_full_gaf(filepath, go_tocorrect, confirmed_lnc)
    sub_graph = graph.subgraph(all_gos).copy()
    ic_indexes = {}
    for ontoname, ann in ontologies.items():
        new_ics = calc_IC(ann, sub_graph)
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

def resnik_pair(go_a, go_b, ontology_ic: dict, graph: nx.DiGraph):
    #https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html#resnik-method
    ancestors_a = set([d for d in nx.descendants(graph, go_a)])
    ancestors_b = set([d for d in nx.descendants(graph, go_b)])
    common_ancestors = ancestors_a.intersection(ancestors_b)
    if len(common_ancestors) > 0:
        common_ancestors_ics = [ontology_ic[goid] for goid in common_ancestors if goid in ontology_ic]
        if len(common_ancestors_ics) > 0:
            return max(common_ancestors_ics)
        else:
            return 0.0
    else:
        return 0.0

def resnik_pair_ancestor(go_a, go_b, ontology_ic: dict, graph: nx.DiGraph):
    #https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html#resnik-method
    ancestors_a = set([d for d in nx.descendants(graph, go_a)])
    ancestors_b = set([d for d in nx.descendants(graph, go_b)])
    common_ancestors = ancestors_a.intersection(ancestors_b)
    if len(common_ancestors) > 0:
        common_ancestors_ics = [(goid, ontology_ic[goid]) for goid in common_ancestors if goid in ontology_ic]
        common_ancestors_ics.sort(key=lambda tp: tp[1])
        if len(common_ancestors_ics) > 0:
            return common_ancestors_ics[-1][0]
        else:
            return ''
    else:
        return ''
    
def resnik_combine(a, b, ontology_ic: dict, graph: nx.DiGraph):
    a = sorted(a)
    b = sorted(b)
    #Calculate Resnik for two lists of GO IDs, using the BMA strategy
    #https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html#bma
    resnik_matrix = [[resnik_pair(go_a, go_b, ontology_ic, graph)
            for go_b in b] for go_a in a]
    resnik_matrix_b = [[resnik_pair_ancestor(go_a, go_b, ontology_ic, graph)
            for go_b in b] for go_a in a]
    a_best_matchs = []
    for i in range(len(resnik_matrix)):
        resniks = resnik_matrix[i]
        a_best_matchs.append(max(resniks))
    
    b_best_matchs = []
    for j in range(len(resnik_matrix[0])):
        resniks = [resnik_matrix[i][j]
            for i in range(len(resnik_matrix))]
        best_match = max(resniks)
        b_best_matchs.append(best_match)

    
    resnik_bma = (np.mean(a_best_matchs) + np.mean(b_best_matchs)) / 2
    '''if len(b) >= 4 and len(a) <= 6 and len(a) >= 5:
        print(a)
        print(b)
        print('a', a_best_matchs)
        print('b', b_best_matchs)
        print(resnik_bma)
        print('\t'+'\t'.join(b))
        for i in range(len(resnik_matrix)):
            row  = resnik_matrix_b[i]
            #print('\t'.join([a[i]] + [str(round(x, 3)) for x in row]))
            print('\t'.join([a[i]] + row))
        quit()'''
    return resnik_bma
