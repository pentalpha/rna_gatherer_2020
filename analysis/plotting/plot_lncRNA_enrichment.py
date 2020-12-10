# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import subprocess
import obonet
import networkx as nx
import math
from scipy.spatial import distance
import sys

gene_list_dir = sys.argv[1]
#gene_list_dir = "/home/pitagoras/main/dev/on-going/rna_gatherer/results/gigas_tissue_specific_lncrna/"
obo_path = sys.argv[2]
population_file = sys.argv[3]
#population_file = "/home/pitagoras/main/dev/on-going/rna_gatherer/test_data/lnc_list/gigas_lnc.txt"
#outdir = "goatools_results"
outdir = gene_list_dir + "/../gigas_geneset_enrichment_analysis"
graphs_dir = outdir+"/graphs"
enrichment_dir = outdir+"/enrichments"
tissue_names = ['Brain', 'Gonad', 'Heart', 
                'Kidney', 'Liver', 'Lung',
                'Muscle', 'Skin']

roots = ['GO:0005575', 'GO:0003674', 'GO:0008150']

min_depth = 2
min_genes_involved = 3
min_fold_change = 1.5

cols = ['go', 'onto', 'in_study', 'in_population', 
        'depth', 'log2fc', 'fdr', 'n_genes', 
        'description', 'genes']

field_index = {cols[i]: i for i in range(len(cols))}
print("Columns:",field_index)

description_minimizer_dict = {'biological_process': "BP",
                            'molecular_function': "MF",
                            'cellular_component': "CC",
                            'regulation': "reg.",
                            'function': "funct.",
                            'cellular': "cell.",
                            'response': "resp.",
                            'positive': "pos.",
                            'negative': "neg."}

def runCommand(cmd, print_cmd=True):
        if print_cmd:
                print("\t> " + cmd)
        process = subprocess.call(cmd, shell=True)
        return process

def get_lines(p):
    return [l.rstrip("\n") for l in open(p,'r').readlines()]

def read_id2gos(id2gos_file):
    id2gos = {}
    with open(id2gos_file, 'r') as stream:
        print("Reading " + id2gos_file)
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            transcript_name = cells[0]
            go_names = cells[1].split(";")
            if not transcript_name in id2gos:
                id2gos[transcript_name] = set()
            for go_name in go_names:
                id2gos[transcript_name].add(go_name)
    return id2gos

#gene_list_dir = "/home/pitagoras/main/dev/on-going/rna_gatherer/results/gigas_tissue_specific_lncrna/"
#population_file = "/home/pitagoras/main/dev/on-going/rna_gatherer/test_data/lnc_list/gigas_lnc.txt"
#predictions_file = "gigas_predictions/lnc_rna_prediction_normal.tsv"

runCommand("mkdir " + graphs_dir)
runCommand("mkdir " + enrichment_dir)

graph = obonet.read_obo(obo_path)
print("Reading graph")
obo_nodes = graph.nodes(data=True)
go_ids = [id_ for id_, data in obo_nodes]
correct_id = {}
descriptions={}
print("Solving redundant GO ids")
for ID in tqdm(go_ids):
    if "alt_id" in obo_nodes[ID]:
        for alt_id in obo_nodes[ID]["alt_id"]:
            correct_id[alt_id] = ID
    correct_id[ID] = ID
    text = obo_nodes[ID]['name']
    words = text.split(" ")
    for i in range(len(words)):
        if words[i] in description_minimizer_dict:
            words[i] = description_minimizer_dict[words[i]]
    text2 = " ".join(words)
    '''if text != text2:
        print("Minimized:", text, "|", text2)'''
    descriptions[ID] = text2
print("Solved " + str(len(correct_id.keys())))

associations = read_id2gos(outdir+"/associations.tsv")

#%%
def filter_enrichment(lines, min_depth, min_genes, min_fc, small_enrichment=False):
    if small_enrichment:
        min_genes = 2
    try:
        return [line for line in lines 
                if line[field_index['depth']] >= min_depth 
                and line[field_index['in_study']] >= min_genes
                and line[field_index['log2fc']] >= min_fc]
    except IndexError:
        print("Invalid line. Example:", lines[0])
        quit()


def highlight_gos(lines, words_to_look_for=[],
                  min_names=6,
                  sortby=5):
    lower_values = sortby == 6
    highlighted = set()
    #for onto, lines in cells_by_ontology.items():
    if len(words_to_look_for)>0:
        for cells in lines:
            for word in words_to_look_for:
                if word.lower() in cells[-2].lower():
                    highlighted.add(cells[0])
    else:
        if not lower_values:
            highlighted.update(
                [cells[0] for cells in lines[-min_names:]])
        else:
            highlighted.update(
                [cells[0] for cells in lines[:min_names]])
    return list(highlighted)
        
#descriptions = {}
def values_from_line(raw_line):
    cells = raw_line.rstrip("\n").split("\t")
    go = cells[0].replace(".", "")
    depth = cells[0].count(".")
    if "GO:" in go:
        onto = cells[1]
        description = cells[3]
        in_study = int(cells[4].split("/")[0])
        if in_study >= 1:
            in_population = int(cells[5].split("/")[0])
            study_total = int(cells[4].split("/")[1])
            pop_total = int(cells[5].split("/")[1])
            study_fraction = 0.0 if in_study == 0 else in_study/study_total
            pop_fraction = 0.0 if in_population == 0 else in_population/pop_total
            if study_fraction > 0.0:
                log2fc = abs(math.log(study_fraction/pop_fraction,2))
                fdr = float(cells[-2])
                n_genes = len(cells[-1].split(","))
                genes = cells[-1].replace(" ", "").split(",")
                #descriptions[go] = description
                return (correct_id[go], onto, in_study, in_population,
                        depth, log2fc, fdr, n_genes, description,
                        genes)

def sep_enrich_by_onto(lines):
    ontos = {"BP": [], "MF": [], "CC": []}
    for line in lines:
        ontos[line[1]].append(line)
    return ontos

def get_edges(indented_nodes):
    edges = []
    #print("Getting edges")
    all_gos = set()
    has_parent = set()
    for i in range(len(indented_nodes)):
        go_id, depth = indented_nodes[i]
        all_gos.add(go_id)
        upper_nodes = graph.neighbors(go_id)
        for node in upper_nodes:
            has_parent.add(go_id)
            edges.append((node, go_id, 1))
            all_gos.add(node)
    found_parents = True
    while found_parents:
        no_parent = all_gos-has_parent
        for go_id in no_parent:
            found_parents = False
            upper_nodes = list(graph.neighbors(go_id))
            if len(upper_nodes) > 0:
                for node in upper_nodes:
                    has_parent.add(go_id)
                    edges.append((node, go_id, 1))
                    all_gos.add(node)
                has_parent.add(go_id)
                found_parents = True
    edges.append(("Root","GO:0008150", 1))
    edges.append(("Root","GO:0003674", 1))
    edges.append(("Root","GO:0005575", 1))
    return list(set(edges))

def enriched_list_from_outfile(outdir, outfile):
    non_filtered_lines = []
    None
    total_lines = 0
    #print("Parsing " + outfile)
    with open(outfile, 'r') as stream:
        previous_line = values_from_line(stream.readline())
        while previous_line == None:
            previous_line = values_from_line(stream.readline())
        total_lines += 1
        for line in stream:
            new_values = values_from_line(line)
            if new_values != None:
                '''new_is_children = 
                if not (new_values[4] > previous_line[4]):
                    #new line is children of last line
                    lines.append(previous_line + [False])'''
                if new_values[4] <= previous_line[4]:
                    #last line is a leaf
                    non_filtered_lines.append(previous_line)
                else:
                    non_filtered_lines.append(previous_line)
                previous_line = new_values
                total_lines += 1
    non_filtered_lines.append(previous_line)
    #lines = filter_enrichment(non_filtered_lines, 0, 0)
    #print("\tUsing", len(lines),"out of", total_lines, "lines.")
    return sep_enrich_by_onto(non_filtered_lines)

def write_parsed_lines(enrich_lines, outfile):
    parsed_enrichment = open(outfile,'w')
    parsed_enrichment.write("\t".join(cols)+"\n")
    for line in enrich_lines:
        new_l = "\t".join([str(x) for x in line])+"\n"
        parsed_enrichment.write(new_l)
    '''parsed_enrichment.write("\n".join(
        ["\t".join([str(cell) for cell in cells]) 
         for cells in enrich_lines]))'''
    parsed_enrichment.close()

def write_parsed_funcs(lines_by_onto, outfile):
    enrich_lines = []
    enrich_lines += lines_by_onto["BP"]
    enrich_lines += lines_by_onto["MF"]
    enrich_lines += lines_by_onto["CC"]
    write_parsed_lines(enrich_lines, outfile)

print("Reading enrichment files")
enrichment_names = (tissue_names + ['Skin-male','Skin-female', 'housekeeping',
                                    'sex_diff','growth_housekeeping'])

raw_enrichment_lists = {name: enriched_list_from_outfile(outdir, 
                                                enrichment_dir+"/"+name+".tsv")
                    for name in tqdm(enrichment_names)}

for name, enrichments_by_ont in tqdm(raw_enrichment_lists.items()):
    for ont, enrichments in enrichments_by_ont.items():
        if len(enrichments) == 0:
            print("No enrichments read for " + name + " " + ont 
                    +  " in " + enrichment_dir+"/"+name+".tsv")

print("Filtering enrichments")
enrichment_lists = {name: {ont: filter_enrichment(enrichments, 
                                            min_depth, min_genes_involved,
                                            min_fold_change, small_enrichment=(name in ['growth_housekeeping']))
                        for ont, enrichments in enrichments_by_ont.items()}  
                    for name, enrichments_by_ont in tqdm(raw_enrichment_lists.items())}

for name, enrichments_by_ont in tqdm(raw_enrichment_lists.items()):
    for ont, enrichments in enrichments_by_ont.items():
        if len(enrichments) == 0:
            print("No enrichments remaining after filtering for " + name + " " + ont 
                    +  " in " + enrichment_dir+"/"+name+".filtered.tsv")

print("Creating edges")
edge_lists = {name: get_edges([(cells[0],cells[4]) for cells in enrichments_by_ont["MF"]]
                            + [(cells[0],cells[4]) for cells in enrichments_by_ont["BP"]]
                            + [(cells[0],cells[4]) for cells in enrichments_by_ont["CC"]])
            for name, enrichments_by_ont in tqdm(enrichment_lists.items())}

sumarry_lines = ["\t".join(["Set of Genes", "MF", "BP", "CC"])]
print("Writing parsed enrichments")
for name, lines_by_onto in enrichment_lists.items():
    lines_by_onto["BP"].sort(key=lambda x: x[-5], reverse=True)
    lines_by_onto["MF"].sort(key=lambda x: x[-5], reverse=True)
    lines_by_onto["CC"].sort(key=lambda x: x[-5], reverse=True)
    
    write_parsed_funcs(lines_by_onto, 
                       enrichment_dir+"/"+name+".filtered.tsv")
    sumarry_lines += ["\t".join([name, str(len(lines_by_onto["MF"])),
                                str(len(lines_by_onto["BP"])),
                                str(len(lines_by_onto["CC"]))])]
open(outdir+"/enrichment_summary.tsv",'w').write("\n".join(sumarry_lines))
#%%
def simplify_graph(g, main_nodes):
    simplified = True
    while simplified:
        nodes = list(g.nodes)
        simplified = False
        for node in nodes:
            if not (node in main_nodes):
                next_nodes = list(g.neighbors(node))
                parent_nodes = list(g.predecessors(node))
                if len(parent_nodes) == 1 and len(next_nodes) >= 1:
                    if parent_nodes[0] != "Root":
                        #print("Removing "+node)
                        g.remove_node(node)
                        for next_node in next_nodes:
                            g.add_edge(parent_nodes[0], 
                                       next_node, 
                                       weight=1)
                        simplified = True
                        break

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def reduce_arrow_len(p1, p2, reducer):
    d1, d2 = (p2[0]-p1[0], p2[1]-p1[1])
    arrow_len, angle = cart2pol(d1, d2)
    new_len = arrow_len - reducer
    new_d1, new_d2 = pol2cart(new_len, angle)
    return new_d1, new_d2

def get_graph(board, tissue_to_analyze, favorite_functions, n_best,
              labels_for_all=True, sorting_col='n_genes',
              simplify=True, ontologies_to_plot=['BP','MF','CC']):
    enrich_lines = []
    for ont_to_plot in ontologies_to_plot:
        enrich_lines += enrichment_lists[tissue_to_analyze][ont_to_plot]
    G = nx.DiGraph()
    indexed_lines = {cells[0]: cells for cells in enrich_lines}
    tissue_edges = edge_lists[tissue_to_analyze]
    for go1, go2, dist in tissue_edges:
        G.add_edge(go1, go2, weight=1)
    highlight_goids = highlight_gos(enrich_lines, 
                                words_to_look_for=favorite_functions,
                                min_names=n_best, 
                                sortby=field_index[sorting_col])
    main_nodes = [node for node in highlight_goids if node in G]
    
    all_nodes = list(G.nodes)
    nodes_to_erase = []
    for node in all_nodes:
        if node in main_nodes:
            continue
        has_path = False
        for main_node in main_nodes:
            if main_node in G:
                if nx.has_path(G, node, main_node):
                    has_path = True
                    break
        if not has_path:
            nodes_to_erase.append(node)
    #print("Erasing",len(nodes_to_erase),"nodes of",len(all_nodes))
    for to_erase in nodes_to_erase:
        G.remove_node(to_erase)
    for root_go in roots:
        if root_go in G:
            if len(list(G.successors(root_go))) > 0:
                main_nodes.append(root_go)
                #print("added",descriptions[root_go])
    n_roots = 0
    for root_go in roots:
        if root_go in G:
            n_roots += 1
    if n_roots == 1:
        G.remove_node("Root")
    if simplify:
        simplify_graph(G, main_nodes)
    #print("Getting positions")
    pos = nx.nx_agraph.graphviz_layout(G)
    
    names = []
    x = []
    y = []
    
    for name, xy in pos.items():
        if name in main_nodes and not(name in roots):
            names.append(name)
            x.append(xy[0])
            y.append(xy[1])

    first_secondary = len(names)
    for name, xy in pos.items():
        if not (name in main_nodes) or (name in roots):
            names.append(name)
            x.append(xy[0])
            y.append(xy[1])
    
    labels = []
    label_sizes = []
    max_label_len = 32
    for name in names:
        if labels_for_all:
            cond = name in descriptions
        else:
            cond = name in main_nodes
        if cond:
            if name in descriptions:
                raw = descriptions[name]
                if len(raw) > max_label_len:
                    descrip_words = raw.split()
                    total = 0
                    raw = descrip_words[0]
                    for word in descrip_words[1:]:
                        if (len(raw)+len(word)+1) <= max_label_len - 5:
                            raw += " " + word
                    raw += "[...]"
                labels.append(raw)
            else:
                labels.append(name)
        else:
            labels.append("")
        if name in main_nodes:
            label_sizes.append(9)
        else:
            label_sizes.append(7)
    min_size = 50
    max_size = 160
    '''n_genes = [indexed_lines[name][-3] if name in indexed_lines else min_size 
               for name in names]
    sizes_unorm = [genes_n+min_size for genes_n in n_genes]
    max_n = max(sizes_unorm)
    sizes = [(max_size*x)/max_n for x in sizes_unorm]'''
    fdr = [indexed_lines[name][-4] if name in indexed_lines else 0.01 
           for name in names[:first_secondary]]
    log2fc = [indexed_lines[name][field_index['log2fc']] 
              if name in indexed_lines else 1.0 
              for name in names[:first_secondary]]
    if len(log2fc) == 0:
        print("Error plotting " + tissue_to_analyze + " not enough items in 'names'")
        return None, highlight_goids
    sizes = [min_size*fc for fc in log2fc]
    indexed_sizes = {names[i]: sizes[i] for i in range(len(sizes))}
    norm_fdr = [-math.log(f,10) 
                if f > 0 
                else -math.log(math.pow(10,-320),10)
                for f in fdr]
    '''for i in range(first_secondary, len(names)):
        sizes[i] = min_size/2
        norm_fdr[i] = '''
        
    lines = []
    for n1, n2 in G.edges:
        if (n1 in pos) and (n2 in pos):
            a = pos[n1]
            b = pos[n2]
            c = reduce_arrow_len(a, b, 15.0)
            lines.append((a,c))
            '''d1, d2 = (n2_x-n1_x, n2_y-n1_y)
            angle = math.atan(d2/d1)
            arrow_len = math.sqrt(d1**2 + d2**2)
            new_len = arrow_len - 0.0
            new_d1 = new_len * math.cos(angle)
            new_d2 = new_len * math.sin(angle)'''
            #arrow_len = distance.euclidean(pos[n1], pos[n2])
            
    
    
        '''#ax1.arrow(start[0], start[1], 
                  end[0], end[1],
                  ,
                  color='blue',#head_width=0.1,
                  head_length=0.075,
                  length_includes_head=True)'''
    handle_fcs = [min(log2fc), 
                  (max(log2fc)+min(log2fc))/2, max(log2fc)]
                          #(max(log2fc)-min(log2fc))/5)
    handle_sizes = [min_size*0.75]+[fc*min_size for fc in handle_fcs]
    handle_labels = ["Not Enriched"] + [str(round(fc,2)) for fc in handle_fcs]
    patchs = []
    for i in range(len(handle_sizes)):
        if i == 0:
            p = plt.scatter([1], [1], min_size*0.75, 
                            'white', edgecolors='black',
                            linewidths=1,
                            alpha=0.6)
            patchs.append(p)
        else:
            p = plt.scatter([1], [1], handle_sizes[i], 
                            'blue',
                            alpha=0.9)
            patchs.append(p)
    
    scatter = board.scatter(x[:first_secondary], y[:first_secondary],
                            sizes[:first_secondary], 
                            norm_fdr[:first_secondary], 
                          cmap="winter", label="Function", alpha=1.0)
    scatter = board.scatter(x[first_secondary:], y[first_secondary:],
                            min_size*0.75, 
                            color="white",
                            linewidths=1,
                            edgecolors="black",
                            label="Function", alpha=1.0)
    for start, end in lines:
        board.arrow(start[0], start[1], end[0], end[1],
                alpha=0.3, width=6, color='blue',zorder=0,
                length_includes_head=True, head_length=10.0)
    '''print(labels)
    print(log2fc)
    print(sizes)'''
    bbox_props = dict(boxstyle="round", fc="w", 
                      ec="0.5", alpha=0.6, edgecolor="white")
    label_down = [True for i in range(len(labels))]
    def label_pos(x, y, is_down):
        return (x, y - 11) if is_down else (x, y + 11)
    dist_th_y = 10
    dist_th_x = 80
    def change_label_positions():
        for i in range(len(labels)):
            positions = [label_pos(x[j],y[j],label_down[j])
                               for j in range(len(labels))]
            current_pos = positions[i]
            dists_y = [distance.euclidean(current_pos[1], other_pos[1])
                     for other_pos in positions]
            dists_x = [distance.euclidean(current_pos[0], other_pos[0])
                     for other_pos in positions]
            smaller_dist = -1
            for j in range(len(dists_y)):
                if j != i:
                    if dists_x[j] <= dist_th_x:
                        #print("small x distance:", labels[i], labels[j])
                        if dists_y[j] <= dists_y[smaller_dist]:
                            smaller_dist = j
            if smaller_dist >= 0:
                if i != smaller_dist:
                    if dists_y[smaller_dist] <= dist_th_y:
                        #print(labels[i])
                        #print("\tsmaller distance is to", labels[smaller_dist])
        
                        #print("\t",labels[i],"and",labels[smaller_dist],
                        #      "are too close")
                        #print(i,smaller_dist)
                        alternative_pos = label_pos(x[i],y[i],
                                                    not label_down[i])
                        if distance.euclidean(alternative_pos[1], 
                                    positions[smaller_dist][1]) > dists_y[smaller_dist]:
                            label_down[i] = not label_down[i]
                            #print("\tSwitching position for " + labels[i])
    change_label_positions()
    change_label_positions()    
    def paint_label(index):
        label_x, label_y = label_pos(x[index],y[index],
                                     label_down[index])
        board.text(label_x, label_y,
                   labels[index], 
                   ha="center", va="center", 
                   size=label_sizes[index],
                   bbox=bbox_props)
    for i in range(len(labels)):
        if not names[i] in main_nodes:
            paint_label(i)
    for i in range(len(labels)):
        if names[i] in main_nodes:
            paint_label(i)
    board.axis("off")
    '''handles, size_labels = scatter.legend_elements(prop="sizes",
                                                   alpha=0.6,
                                                  num=4)'''
    #from matplotlib.patches import Circle
    '''#from matplotlib.patches import Patch
    from matplotlib.patches import Ellipse'''
    
    plt.scatter([1], [1], handle_sizes[-1]+100,
                "white", edgecolors='white',
                alpha=1.0)
    '''size_labels = [label.split('{')[0] 
                   + "{" + str(
                       round(
                           float(label.split('{')[-1].split("}")[0])/min_size
                       ,2))
                   + "}" + label.split('}')[-1] 
                   for label in size_labels]'''
    '''legend2 = ax1.legend(handles, size_labels, loc="lower center", 
                         title="Genes")'''
    #patchs = [ for s in handle_sizes]
    legend2 = board.legend(
        patchs, handle_labels,
        bbox_to_anchor=(0., -0.15, 1., .102), 
        loc='upper center', title="abs(log2(foldChange))",
        ncol=len(patchs))
    
    cbar = plt.colorbar(scatter, ax=board, orientation='horizontal')
    cbar.set_label('FDR')
    #print(list(cbar.get_ticks()))
    if tissue_to_analyze == "Lungug":
        cbar.ax.set_xticklabels(["{:.1e}".format(math.pow(10,-tick)) 
                             for tick in [50,100,150,200,250,300,320]])
    else:
        cbar.ax.set_xticklabels(["{:.1e}".format(math.pow(10,-tick)) 
                                 for tick in cbar.get_ticks()])
    plt.setp(cbar.ax.get_xticklabels(), rotation=25, ha="right",
         rotation_mode="anchor")
    cbar.ax.tick_params(axis='both', which='major', labelsize=8)
    return cbar, highlight_goids

fig, ax1 = plt.subplots(1, figsize=(7,7))
cbar, highlight_goids = get_graph(ax1, "sex_diff", ['pigment', 'melanosome', 'melanin', 'melanocyte'], 6,
          labels_for_all=True,
          sorting_col='log2fc',
          simplify = False)
print(highlight_goids)
ax1.set_title("Sex Differential lncRNA, Enriched Pigmentation Functions")
fig.tight_layout()
#fig.show()
fig.savefig(graphs_dir+"/pigmentation_graph.png")

transcripts = set()

for ont, enrichments in enrichment_lists['sex_diff'].items():
    for enrich_line in enrichments:
        if enrich_line[0] in highlight_goids:
            transcripts.update(enrich_line[-1])
print(len(transcripts), "pigmentation genes.")
transcripts_out = "\n".join(transcripts)
open(outdir+"/pigmentation.txt",'w').write(transcripts_out)
#%%
'''maturation_base_functions = ['GO:0061458', 'GO:0007548']
def get_preds(data, go_id, depth = 0):
    pres = list(data.predecessors(go_id))
    if len(pres) == 0 or depth == 2:
        return [go_id]
    else:
        results = [go_id]
        for pred in pres:
            results += get_preds(data, pred, depth=depth+1)
        return results

maturation_funcs = set()
for func in maturation_base_functions:
    maturation_funcs.update(get_preds(graph, func))
maturation_funcs = list(maturation_funcs)
mat_descriptions = [descriptions[func] for func in maturation_funcs]
with open(outdir+"/maturation_funcs.txt", 'w') as stream:
    for i in range(len(maturation_funcs)):
        stream.write(maturation_funcs[i] + " " 
                     + mat_descriptions[i] + "\n")

for genename, genefuncs in associations.items():
    found = 0
    for func in maturation_funcs:
        if func in genefuncs:
            found += 1
    if found > 0:
        print(found," maturation functions annotated in", genename)'''
#%%

#%%
to_plot=['Brain', 'Gonad', 'Heart', 'Kidney', 
            'Liver', 'Lung', 'Muscle', 'Skin']

for field in ['n_genes', 'fdr', 'log2fc']:
    fig, axes = plt.subplots(4,2, figsize=(14,28))
    i = 0
    ax_list = []
    for ax_row in axes:
        for ax in ax_row:
            ax_list.append((to_plot[i],ax))
            i+=1
    figure = {}

    for tissue_name, ax in tqdm(ax_list):
        cb, hg_go = get_graph(ax, tissue_name, [], 12, 
                    labels_for_all=False, 
                    sorting_col=field, simplify = False)
        ax.set_title(tissue_name)
        figure[tissue_name] = (ax,cb)
    fig.tight_layout()
    #fig.show()
    fig.savefig(graphs_dir+"/tissues_graph_"+field+".png")

#%%
other_analysis = ['housekeeping', 'sex_diff', 'growth_housekeeping']
cool_name = ['Housekeeping Genes', 'Sex Differential Genes', 'Housekeeping Growth Genes']
onts_to_plt = ['MF', 'BP', 'CC']
for field in ['fdr']:
    for i in range(len(other_analysis)):
        other_key = other_analysis[i]
        c_name = cool_name[i]
        for ont in onts_to_plt:
            fig, ax = plt.subplots(1, figsize=(7,7))
            cb, hg_go = get_graph(ax, other_key, [], 10, 
                    labels_for_all=False, 
                    sorting_col=field, simplify = False,
                    ontologies_to_plot=[ont])
            ax.set_title(c_name + " - " + ont)
            fig.tight_layout()
            #fig.show()
            fig.savefig(graphs_dir+"/"+other_key+"."+ont+".png")