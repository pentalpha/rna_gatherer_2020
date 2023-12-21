import obonet
import networkx
import sys
from os import path
from tqdm import tqdm
import numpy as np
import pandas as pd

from semantic_similarity import calc_IC_indexes, resnik_combine, simgic

current_dir = path.dirname(path.abspath(__file__))
proj_dir = path.dirname(current_dir)
test_dir = path.join(proj_dir, 'test')

go_roots = ['GO:0005575', 'GO:0003674', 'GO:0008150']

def read_rna_central_trans(filepath):
    print('Reading', filepath)
    stream = open(filepath, 'r')
    rnacentral2mgi = {}
    mgi2rnacentral = {}

    for rawline in stream:
        cells = rawline.rstrip('\n').split('\t')
        rnacentral = cells[0] + '_' + cells[3]
        mgi_name = cells[2]
        rnacentral2mgi[rnacentral] = mgi_name
        mgi2rnacentral[mgi_name] = rnacentral
    return rnacentral2mgi, mgi2rnacentral

def read_id2go_annotation(filepath, go_tocorrect, translator=None):
    print('Reading', filepath)
    stream = open(filepath, 'r')
    ann = {}
    not_translated = set()
    for rawline in stream:
        cells = rawline.rstrip('\n').split('\t')
        name = cells[0]
        if translator:
            if name in translator:
                name = translator[name]
            '''else:
                #print(name, 'not in translator')
                not_translated.add(name)'''
        if not name in ann:
            ann[name] = []
        pval = float(cells[3]) if len(cells) >= 4 else 0.0
        fdr = float(cells[4]) if len(cells) >= 5 else 0.0
        goid = cells[1]
        if goid in go_tocorrect:
            goid = go_tocorrect[goid]
        ann[name].append((goid, pval, fdr))
    
    return ann



def read_gaf_annotation(filepath: str, go_tocorrect, translator = None, 
                        surely_are_lnc = set()):
    rna_central_lncrnas_path = path.join(path.dirname(filepath), 'rna_central_lncrna_list.txt')
    rna_central_lncrnas = open(rna_central_lncrnas_path,'r').read().split('\n')
    rna_central_lncrnas = set(rna_central_lncrnas)
    print('Reading', filepath)
    stream = open(filepath, 'r')
    ann = {'MF': {}, 'BP': {}, 'CC': {}}
    
    evi_not_allowed = []
    computational_evi = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']
    enough_evi = 0
    not_enough_evi = 0
    not_lncrna = 0
    is_a_lncrna = 0
    
    from_rnacentral = 0
    not_lnc_from_rnacentral = 0
    
    for rawline in stream:
        if rawline[0] != '!':
            cells = rawline.rstrip('\n').split('\t')
            db = cells[0]
            name = cells[1]
            if translator:
                if name in translator:
                    name = translator[name]

            #if db in ['RNAcentral', 'MGI']:
            if db in ['MGI', 'RNAcentral']:
                is_nc_rna = False
                if db == 'RNAcentral':
                    #islncRNA = 'lnc' in cells[2] or 'long' in cells[2]
                    is_nc_rna = name in rna_central_lncrnas
                    if is_nc_rna:
                        from_rnacentral += 1
                    else:
                        not_lnc_from_rnacentral += 1
                elif db == 'MGI':
                    is_nc_rna = 'lnc' in cells[-6]
                
                if not is_nc_rna and name in surely_are_lnc:
                    is_nc_rna = True
                
                if is_nc_rna:
                    is_a_lncrna += 1
                    
                    evi_code = cells[6]
                    has_evi = True
                    if len(evi_code) <= 3:
                        if evi_code in evi_not_allowed:
                            has_evi = False
                            not_enough_evi += 1
                    else:
                        print('Long evidence code:', evi_code)
                    
                    if has_evi:
                        enough_evi += 1
                        goid = cells[4]
                        if goid in go_tocorrect:
                            goid = go_tocorrect[goid]
                        ont = cells[8].replace('F', 'MF').replace('C', 'CC').replace('P', 'BP')
                        if not name in ann[ont]:
                            ann[ont][name] = set()
                        ann[ont][name].add((goid, 0, 0))
                else:
                    not_lncrna += 1
    print('enough_evi', enough_evi)
    print('not_enough_evi', not_enough_evi)
    print('is_a_lncrna', is_a_lncrna)
    print('not_lncrna', not_lncrna)
    
    print('from_rnacentral', from_rnacentral)
    print('not_lnc_from_rnacentral', not_lnc_from_rnacentral)
    
    return ann['MF'], ann['BP'], ann['CC']

'''ref_mf_ann, ref_bp_ann, ref_cc_ann = read_gaf_annotation(ref_ann_path, 
    go_tocorrect, translator=mgi2rnacentral)'''

def filter_ann(ann_dict, to_erase_path):
    to_erase = set()

    for rawline in open(to_erase_path, 'r'):
        if rawline.startswith('id: GO'):
            go = rawline.lstrip('id: ').rstrip('\n')
            to_erase.add(go)

    to_erase.update(go_roots)

    for genename in ann_dict.keys():
        go_list = ann_dict[genename]
        ann_dict[genename] = [x for x in go_list if not x in to_erase]

    empty_anns = [genename for genename in ann_dict.keys() if len(ann_dict[genename]) == 0]
    for genename in empty_anns:
        del ann_dict[genename]
    print('Erased annotations to', len(empty_anns), 'genes')


def expand_ann(go_graph, ann_dict, ont_name, term_to_namespace):
    ont = (ont_name.replace('MF', 'molecular_function')
        .replace('BP', 'biological_process')
        .replace('CC', 'cellular_component'))
    print('Expanding', ont)
    expand_list = []
    for genename, go_list in ann_dict.items():
        all_superterms = set()
        goids = [go_id for go_id, pval, fdr in go_list]
        for go_id in goids:
            #if go_id in go_graph:
            superterms = [superterm for superterm in networkx.descendants(go_graph, go_id)
                          if not superterm in go_list]
            same_ont_terms = [term for term in superterms if term in term_to_namespace]
            same_ont_terms = [term for term in same_ont_terms if term_to_namespace[term] == ont]
            if len(same_ont_terms) > 4:
                #print(genename, go_id, 'could be expanded with:')
                #print(same_ont_terms)
                all_superterms.update(same_ont_terms)
        if len(all_superterms) > 0:
            new_goids = [(term, 0.0, 0.0) for term in all_superterms]
            expand_list.append([genename, new_goids])
    
    print(len(expand_list), 'expansions')
    print('Potential new associations:', 
          sum([len(expansion) for gene, expansion in expand_list]))
    
    new_associations = 0
    for gene, expansion in expand_list:
        before = len(ann_dict[gene])
        if type(ann_dict[gene]) == set:
            ann_dict[gene].update(expansion)
        else:
            existing_gos = set([go for go, fdr, pval in ann_dict[gene]])
            sub_expansion = [go for go in expansion if not go in existing_gos]
            for go in sub_expansion:
                ann_dict[gene].append((go, -1, -1))
        after = len(ann_dict[gene])
        new_associations += after - before
        
    print('New associations:', new_associations)

    for genename in ann_dict.keys():
        go_list = ann_dict[genename]
        ann_dict[genename] = [x[0] if len(x) == 3 else x for x,y,z in go_list]

def load_obo(obo_path):
    print('Loading go graph')
    go_graph = obonet.read_obo(obo_path, ignore_obsolete=False)
    obo_nodes = go_graph.nodes(data=True)
    go_ids = [id_ for id_, data in obo_nodes]
    alt_ids = {}
    term_to_namespace = {}
    print("Solving redundant GO ids")
    for ID in tqdm(go_ids):
        if "alt_id" in obo_nodes[ID]:
            for alt_id in obo_nodes[ID]["alt_id"]:
                alt_ids[alt_id] = ID
        if "namespace" in obo_nodes[ID]:
            term_to_namespace[ID] = obo_nodes[ID]["namespace"]
    print("Solved " + str(len(alt_ids.keys())))
    return go_graph, alt_ids, term_to_namespace

def make_filtered_obo(obo_path):
    new_obo_path = obo_path.replace('.obo', '.filtered.obo')
    input_obo = open(obo_path, 'r')
    output = open(new_obo_path, 'w')
    
    erased = {}
    kept = 0
    for rawline in input_obo:
        if rawline.startswith('relationship:'):
            parts = rawline.split()
            if parts[1] == 'part_of':
                output.write(rawline)
                kept += 1
            else:
                if not parts[1] in erased:
                    erased[parts[1]] = 0
                erased[parts[1]] += 1
        else:
            output.write(rawline)
        
    output.close()
    print('Relationships erased:')
    print(erased)
    print('Relationships kept (part_of):')
    print(kept)
    return new_obo_path
    

def compare_to_reference(new_ann_original, ref_ann_original, ic, go_graph):
    new_ann = {}
    for key, values in new_ann_original.items():
        #key_tuples = [go if type(go) == str else go[0] for go, fdr, pval in values]
        #print(key_tuples)
        new_ann[key] = set(values)
    ref_ann = {}
    for key, values in ref_ann_original.items():
        #key_tuples = [go if type(go) == str else go[0] for go, fdr, pval in values]
        #print(key_tuples)
        ref_ann[key] = set(values)
    new_genes = set(new_ann.keys())
    ref_genes = set(ref_ann.keys())
    #print('\tOnly new annotation:', len(new_genes-ref_genes))
    #print('\tNot predicted:', len(ref_genes-new_genes))
    predicted_and_annotated = ref_genes.intersection(new_genes)
    predicted_and_annotated_perc = len(predicted_and_annotated) / len(ref_genes)
    #print('\t% of Annotated covered:', predicted_and_annotated_perc*100)
    new_ass_total = sum([len(vals) for key, vals in new_ann.items()])
    ref_ass_total = sum([len(vals) for key, vals in ref_ann.items()])
    size = (new_ass_total/ref_ass_total)*100
    #print('\tSize compared to reference:', size)
    
    tp = 0
    fp = 0
    for gene, goids in new_ann.items():
        for goid in goids:
            if gene in ref_ann:
                if goid in ref_ann[gene]:
                    tp += 1
                else:
                    fp += 1
            else:
                fp += 1
    #print('\tTP', tp)
    #print('\tFP', fp)

    simgics = []
    resniks = []
    for genename in predicted_and_annotated:
        new_goids_set = new_ann[genename]
        '''for a in new_goids_set:
            assert type(a) == str, 'a is' + str(a)'''
        ref_goids_set = ref_ann[genename]
        simgics.append(simgic(new_goids_set, ref_goids_set, ic))
        resniks.append(resnik_combine(new_goids_set, ref_goids_set, ic, go_graph))
    simgics.sort()
    resniks.sort()
    
    #print('\tSIMGIC:', np.quantile(simgics, 0.25), np.mean(simgics), np.quantile(simgics, 0.75))
    #print('\tResniks:', np.quantile(resniks, 0.25), np.mean(resniks), np.quantile(resniks, 0.75))
    
    return [tp, round(tp/(tp+fp), 4), round(np.mean(simgics), 4), round(np.mean(resniks), 4), new_ass_total, ref_ass_total]
    

#%%
if __name__ == '__main__':
    "python new_article_benchmark.py <analysis_dir>"
    
    analysis_dir = "/home/pitagoras/main/experiments/mgi_simgic_tpm-exon-backup/"
    #analysis_dir = sys.argv[1]
    obo_path = path.join(analysis_dir, 'go.obo')
    to_erase_path = path.join(analysis_dir, 'gocheck_do_not_annotate.obo')
    
    new_obo_path = make_filtered_obo(obo_path)
    go_graph, go_tocorrect, term_to_namespace = load_obo(new_obo_path)

    #%%

    rnacentral_translation_path = path.join(analysis_dir, 'mgi.tsv')
    rnacentral2mgi, mgi2rnacentral = read_rna_central_trans(rnacentral_translation_path)

    new_ann_names = [['MF', 'MF.Spearman.min0.929.tsv'], 
                     ['BP', 'BP.Spearman.min0.936.tsv'],
                     ['CC', 'CC.Spearman.min0.91.tsv']]
    lncrna2goa_names = [['MF', 'LNCRNA2GOA_MF.tsv'],
                        ['BP', 'LNCRNA2GOA_BP.tsv']]
    ref_ann_path = path.join(analysis_dir, 'mus_musculus.gaf')

    new_ann = []
    untranslated = set()
    for ont_name, ann_name in new_ann_names:
        print(ont_name)
        ann_path = path.join(analysis_dir, ann_name)
        annotation = read_id2go_annotation(ann_path, go_tocorrect)
        new_ann.append([ont_name, annotation])
        print(len(annotation.keys()), 'genes with', ont_name)
    
    lncrna2goa_ann = []
    for ont_name, ann_name in lncrna2goa_names:
        print(ont_name)
        ann_path = path.join(test_dir, ann_name)
        annotation = read_id2go_annotation(ann_path, go_tocorrect, translator=mgi2rnacentral)
        lncrna2goa_ann.append([ont_name, annotation])
        print(len(annotation.keys()), 'genes with', ont_name)
        
    all_genes = list(new_ann[0][1].keys()) + list(new_ann[1][1].keys()) + list(new_ann[2][1].keys())
    all_genes = set(all_genes)
    print(len(all_genes), 'ncRNA genes annotated in all ontologies')

    information_contents = calc_IC_indexes(ref_ann_path, go_tocorrect, go_graph, all_genes)
    
    ref_mf_ann, ref_bp_ann, ref_cc_ann = read_gaf_annotation(ref_ann_path, 
        go_tocorrect, translator=mgi2rnacentral, surely_are_lnc=all_genes)
    
    all_ref_genes = list(ref_mf_ann.keys()) + list(ref_bp_ann.keys()) + list(ref_cc_ann.keys())
    all_ref_genes = set(all_ref_genes)
    print(len(all_ref_genes), 'ncRNA genes with reference annotations')
    no_ref = all_genes-all_ref_genes
    print('No ref annotations:', len(no_ref))
    for ont_name, ann in new_ann + lncrna2goa_ann:
        for genename in no_ref:
            if genename in ann:
                del ann[genename]
    reference_annotations = [['MF',ref_mf_ann], ['BP', ref_bp_ann], ['CC', ref_cc_ann]]
    #%%
    '''print('Expanding references')
    for ont_name, ann in reference_annotations:
        #expand_ann(go_graph, ann, ont_name, term_to_namespace)
        filter_ann(ann, to_erase_path)
    
    print('Expanding new annotations')
    for ont_name, ann in new_ann:
        #expand_ann(go_graph, ann, ont_name, term_to_namespace)
        filter_ann(ann, to_erase_path)'''
    
    '''print('Expanding lncrna2goa')
    for ont_name, ann in lncrna2goa_ann:
        expand_ann(go_graph, ann, ont_name, term_to_namespace)
        filter_ann(ann, to_erase_path)'''
        
    #%%
    '''print('Annotation', 'True Positives', 'Precision', 'SIMGIC', 'Resnik', 
          'Annotations Predicted', 'Ref. Annotations', sep='\t')
    #print('RNA Gatherer')
    for i in range(len(new_ann)):
        ont_name, ann = new_ann[i]
        #print(ont_name)
        _, ref_ann = reference_annotations[i]
        
        result = compare_to_reference(ann, ref_ann, information_contents[ont_name], go_graph)
        result_str = '\t'.join([str(x) for x in result])
        print('RNA Gatherer '+ont_name + '\t' + result_str)
    #print('LNCRNA2GOA')
    for i in range(len(lncrna2goa_ann)):
        ont_name, ann = lncrna2goa_ann[i]
        #print(ont_name)
        _, ref_ann = reference_annotations[i]
        
        result = compare_to_reference(ann, ref_ann, information_contents[ont_name], go_graph)
        result_str = '\t'.join([str(x) for x in result])
        print('LNCRNA2GOA '+ont_name + '\t' + result_str)'''
    
    #%%
    enrichments_dir = ("/home/pitagoras/main/dev/rna_gatherer/result_data/gigas_geneset_enrichment_analysis/"
                     +"enrichments/")
    sex_diff_path = enrichments_dir+"sex_diff.filtered.tsv"
    growth_housekeeping_path = enrichments_dir+"growth_housekeeping.filtered.tsv"
    pigmentation_gos = "/home/pitagoras/main/dev/rna_gatherer/dissertation/go_sets/all_pigmentation.tsv"
    pigmentation_gos = [rawline.split()[0] for rawline in open(pigmentation_gos, 'r').readlines()]
    
    sex_diff_df = pd.read_csv(sex_diff_path, sep='\t')
    sex_diff_df = sex_diff_df[sex_diff_df['go'].isin(pigmentation_gos)]
    print(sex_diff_df)
    allowed_gos = set()
    
    for index, row in sex_diff_df.iterrows():
        goid = row['go']
        allowed_gos.update([superterm for superterm in networkx.descendants(go_graph, goid)])
        allowed_gos.add(goid)

        go_graph.nodes[goid]['enriched'] = True
        go_graph.nodes[goid]['fdr'] = row['fdr']
    
    pigmentation_graph = go_graph.subgraph(list(allowed_gos)).copy()

    networkx.readwrite.gml.write_gml(pigmentation_graph, 
        ("/home/pitagoras/main/dev/rna_gatherer/result_data/gigas_geneset_enrichment_analysis/"+
         "enrichments/pigmentation_graph.gml"))
    
    growth_df = pd.read_csv(growth_housekeeping_path, sep='\t')
    growth_df = growth_df[growth_df['onto'] == 'BP']
    growth_df = growth_df[growth_df['fdr'] <= 0.05]
    print(growth_df)
    allowed_gos2 = set()
    
    for index, row in growth_df.iterrows():
        print(row['go'], row['onto'], row['fdr'], row['description'])
        goid = row['go']
        allowed_gos2.update([superterm for superterm in networkx.descendants(go_graph, goid)])
        allowed_gos2.add(goid)

        go_graph.nodes[goid]['enriched'] = True
        go_graph.nodes[goid]['fdr'] = row['fdr']
        go_graph.nodes[goid]['enriched_name'] = row['description']
    
    growth_graph = go_graph.subgraph(list(allowed_gos2)).copy()

    networkx.readwrite.gml.write_gml(growth_graph, 
        ("/home/pitagoras/main/dev/rna_gatherer/result_data/gigas_geneset_enrichment_analysis/"+
         "enrichments/growth_housekeeping_graph.gml"))