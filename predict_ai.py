"""Function predictor for lncRNA that uses gene coexpression.

This script is used to apply statistics in order to figure out
possible Gene Ontology terms for lncRNAs, based on gene 
expression counts and an annotation for coding genes. The gene
expression is read from a counts table (see test_data/counts)
and the expression of lncRNA and coding genes is compared
by calculating correlation coefficients through several
differente metrics. The pairs of (lncRNA,gene) with 
coefficients passing the required minimum value to be 
considered inside the confidence level are taken as
coexpressed.
The associations between coding genes and GO terms in the 
functional annotation (see test_data/annotation) of the
coding genes are passed as possible associations for their
lncRNA coexpressed pairs. Statistical tests are used to filter
which of these possible associations are statistically
relevant, calculating P-Values and FDRs.
The associations passing the filtering are grouped together in
a output file, the functional prediction.

Author: Pitágoras Alves (github.com/pentalpha)
"""

import sys
from nexus.functional_prediction import *
from nexus.util import *
from nexus.confidence_levels import *
from nexus.bioinfo import load_metrics, short_ontology_name
from nexus.interaction_predictor import *
import obonet
import networkx
import numpy as np
import argparse
import multiprocessing

'''
Loading configurations and command line arguments
'''
from config import configs, require_files
mandatory_files = ["go_obo"]
require_files(mandatory_files)

#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

available_species = get_available_species()

display_cache = get_cache()

all_methods = ["MIC","DC","PRS","SPR","SOB","FSH"]
all_ontologies = ["molecular_function","cellular_component","biological_process"]
default_methods = load_metrics(confs["metrics_table"])
highest_confidence = str(max([int(conf) for conf in default_methods.keys()]))

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-cr", "--count-reads", required=True,
        help=("Count reads table path. TSV file with where the first column must be the gene names and the " +
            " following columns are the FPKM normalized counts of each sample."))
    ap.add_argument("-reg", "--regulators", required=True,
        help=("A list of regulators, where each line contains the name of one gene."))
    ap.add_argument("-ann", "--annotation", required=True,
        help=("Functional annotation file of the genes."))
    ap.add_argument("-o", "--output-dir", required=True, help=("Output directory."))
    ap.add_argument("-ont", "--ontology-type", required=False,
        default="molecular_function", help=("One of the following: molecular_function (default),"
        +" cellular_component or biological_process."))
    ap.add_argument("-pre", "--precision-level", help = ("Used to increase the interaction "
        + "preditor precision. Possible values: " + ",".join(prob_diff_values.keys())+". Default=3."), 
        default="3")

    default_threads = max(2, multiprocessing.cpu_count()-1)
    ap.add_argument("-p", "--processes", required=False,
        default=default_threads, help=("CPUs to use. Default: " + str(default_threads)+"."))
    ap.add_argument("-K", "--k-min-coexpressions", required=False,
        default=1, help=("The minimum number of ncRNAs a Coding Gene must be coexpressed with."
                        +" Increasing the value improves accuracy of functional assignments, but"
                        +" may restrict the results. Default: 1."))
    ap.add_argument("-pv", "--pvalue", required=False,
        default=0.05, help=("Maximum pvalue. Default: 0.05"))
    ap.add_argument("-fdr", "--fdr", required=False,
        default=0.05, help=("Maximum FDR. Default: 0.05"))
    ap.add_argument("-m", "--min-m", required=False,
        default=1, help=("Minimum m value. Default: 1"))
    ap.add_argument("-M", "--min-M", required=False,
        default=1, help=("Minimum M value. Default: 1"))
    ap.add_argument("-n", "--min-n", required=False,
        default=1, help=("Minimum n value. Default: 1"))

    ap.add_argument("-chu", "--cache-usage", required=False,
        default=0.6, help=("Portion of the cache memory to use for storing the counts table."))
    ap.add_argument("-ch", "--cache-size", required=False,
        default=display_cache, help=("Sets the size of the cache memory. Default: auto-detection of CPU cache size."))
    return vars(ap.parse_args())

cmdArgs = getArgs()
count_reads_path = cmdArgs["count_reads"]
regulators_path = cmdArgs["regulators"]
coding_gene_ontology_path = cmdArgs["annotation"]
tempDir = cmdArgs["output_dir"]
go_path = confs["go_obo"]
interaction_model_path = confs["interaction_model"]
precision_level = cmdArgs["precision_level"]
prob_diff = 0.0
if precision_level in prob_diff_values.keys():
    prob_diff = prob_diff_values[precision_level]
else:
    print("Invalid precision level.")
    quit()
print("Loading interaction predictor")
interaction_model = InteractionPredictor(interaction_model_path, diff=prob_diff)

ontology_types_arg = cmdArgs["ontology_type"].split(",")
if ontology_types_arg[0] == "ALL":
    ontology_types_arg = all_ontologies

threads = int(cmdArgs["processes"])

cache_usage = float(cmdArgs["cache_usage"])
available_cache = get_cache(usage=cache_usage)
if "cache_size" in cmdArgs:
    available_cache = int(int(cmdArgs["cache_size"]) * cache_usage)
print("Available cache memory: " + str(int(available_cache/1024)) + "KB")

pval = float(cmdArgs["pvalue"])
fdr = float(cmdArgs["fdr"])
min_n = int(cmdArgs["min_n"])
min_M = int(cmdArgs["min_M"])
min_m = int(cmdArgs["min_m"])
K = int(cmdArgs["k_min_coexpressions"])

regulators_max_portion = 0.4

'''
Creating output directory
'''

if not os.path.exists(tempDir):
    os.mkdir(tempDir)

def find_correlated(reads, regulators, threads, predictor, output_stream):
    """Find coexpressed pairs using a set of metrics."""
    if len(reads) < threads*2:
        threads = len(reads)/2
    coding_noncoding_pairs = []
    func = ai_find_coexpression_process
    genes_per_process = int(len(reads) / threads)
    
    limit = len(reads)-1
    end = 0
    last = -1
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    processes = []
    last_pid = 0
    for i in range(threads-1):
        start = last+1
        end = start + genes_per_process
        if end >= limit:
            end = limit
        parcial_df = reads.iloc[start:end]
        p = multiprocessing.Process(target=func, 
            args=(i, parcial_df, regulators, predictor, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(start) + " to " + str(end))
        p.start()
        last = end

        last_pid = i
        if end == limit:
            break
    if end < limit:
        parcial_df = reads.iloc[end:limit]
        p = multiprocessing.Process(target=func, 
            args=(last_pid+1, parcial_df, regulators, predictor, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(end) + " to " + str(limit))
        p.start()

    for p in processes:
        p.join()

    #print("Merging results")
    for value in return_dict.values():
        coding_noncoding_pairs += value

    #print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    
    for coding_name, noncoding_name in coding_noncoding_pairs:
        output_stream.write("\t".join([coding_name,noncoding_name]) + "\n")
    manager._process.terminate()
    manager.shutdown()
    del manager

'''
Pre-processing of the count-reads table
'''
print("Ontology type is " + str(ontology_types_arg))
reads = pd.read_csv(count_reads_path, sep='\t')
print(str(len(reads)))
reads["constant"] = reads.drop([reads.columns[0]], axis=1).apply(
        lambda row: is_constant(np.array(row.values,dtype=np.float32)),axis=1
    )
mask = reads["constant"] == False
reads = reads[mask]
del reads["constant"]
print(reads.head())
print(str(len(reads)))

print("Reading regulators")
regulators = []
with open(regulators_path,'r') as stream:
    for line in stream.readlines():
        regulators.append(line.rstrip("\n").lstrip(">"))

interaction_list = tempDir+"/interacting_pairs.tsv"
if not os.path.exists(interaction_list):
    out_stream = open(interaction_list, 'a+')
    print("Separating regulators from regulated.")
    mask = reads[reads.columns[0]].isin(regulators)
    regulators_reads = reads.loc[mask]
    non_regulators_reads = reads.loc[~mask]
    print(str(len(non_regulators_reads)) + " regulated.")
    print(str(len(regulators_reads)) + " regulators.")
    #print("Generating correlations.")

    available_size = available_cache

    '''
    Split the table into cache-sized smaller parts.
    '''
    max_for_regulators = available_size*regulators_max_portion
    regs_size = getsizeof(regulators_reads)
    regulator_dfs = [regulators_reads]
    if regs_size > max_for_regulators:
        regulator_dfs = split_df_to_max_mem(regulators_reads, max_for_regulators)
        available_size -= max_for_regulators
    else:
        available_size -= regs_size

    dfs = split_df_to_max_mem(non_regulators_reads, available_size)

    '''print("Chunks for regulated: " + str(len(dfs)) 
    + "\nChunks for regulators: " + str(len(regulator_dfs)))'''

    df_pairs = []
    for df in dfs:
        for regulator_df in regulator_dfs:
            df_pairs.append((df,regulator_df))

    '''
    Calculate the correlations for each part of the table
    '''
    i = 1
    for df,regulator_df in tqdm(df_pairs):
        find_correlated(df, regulators_reads, int(threads), interaction_model, out_stream)
        i += 1

    out_stream.close()

coding_genes = {}
genes_coexpressed_with_ncRNA = {}
total_interactions = 0

'''
Load predicted interactions
'''
with open(interaction_list, 'r') as stream:
    lines = 0
    invalid_lines = 0
    
    for raw_line in stream.readlines():
        cells = raw_line.rstrip("\n").split("\t")
        if len(cells) == 2:
            gene = cells[0]
            rna = cells[1]
            if not gene in coding_genes:
                coding_genes[gene] = 0
            coding_genes[gene] += 1
            if not rna in genes_coexpressed_with_ncRNA:
                genes_coexpressed_with_ncRNA[rna] = set()
            genes_coexpressed_with_ncRNA[rna].add(gene)
            total_interactions += 1
        else:
            invalid_lines += 1
        lines += 1
    if lines == 0:
        print("Fatal error, no correlations could be loaded from "
                + interaction_list + "\n(The file may be "
                + "corrupted or just empty)")
        quit()
    else: 
        print(str(float(invalid_lines)/lines)
            + " lines without proper number of columns (2 columns)")
    print("Regulators interacting = "+str(len(genes_coexpressed_with_ncRNA.keys())))
    print("Coding genes interacting = "+str(len(coding_genes.keys())))
    print("Predicted interactions = "+str(total_interactions))

N = len(reads)

print("Loading GO network.")
graph = obonet.read_obo(go_path)

onto_id2gos = {"biological_process":{},"molecular_function":{},"cellular_component":{}}
onto_genes_annotated_with_term = {"biological_process":{},"molecular_function":{},"cellular_component":{}}

print("Reading annotation of coding genes from " + coding_gene_ontology_path)
with open(coding_gene_ontology_path,'r') as stream:
    lines = 0
    invalid_lines = 0
    associations = 0
    for raw_line in stream.readlines():
        cells = raw_line.rstrip("\n").split("\t")
        if len(cells) == 3 or len(cells) == 4:
            gene = cells[0].upper()
            go = cells[1]
            onto = cells[2]
            if onto in onto_id2gos.keys():
                id2gos = onto_id2gos[onto]
                genes_annotated_with_term = onto_genes_annotated_with_term[onto]
                #coding_genes.add(gene)
                if gene in coding_genes:
                    if not gene in id2gos:
                        id2gos[gene] = set()
                    id2gos[gene].add(go)
                    if not go in genes_annotated_with_term:
                        genes_annotated_with_term[go] = set()
                    genes_annotated_with_term[go].add(gene)
                    associations += 1
            else:
                invalid_lines += 1
        else:
            invalid_lines += 1
        lines += 1
    print(str(float(invalid_lines)/lines)
            + " lines without proper number of columns (3 or 4 columns)")
    print(str(associations) + " valid associations loaded.")

def predict(tempDir, ontology_type="molecular_function",
            k_min_coexpressions=1,
            pval_threshold=0.05, fdr_threshold=0.05,
            min_n=5, min_M=5, min_m=1):
    """Predict functions based on the loaded correlation coefficients."""
    
    K = k_min_coexpressions
    ontology_type_mini = short_ontology_name(ontology_type)
    
    print("Ontology type = " + ontology_type
            + ", pvalue = " + str(pval_threshold) + ", fdr = " + str(fdr_threshold))
    
    out_dir = tempDir+"/"+ontology_type_mini
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    longer = []
    if K != 1 or min_n != 1 or min_M != 1 or min_m != 1:
        longer = ["K"+str(K),"n"+str(min_n),"M"+str(min_M),"m"+str(min_m)]
    run_mode = "LNC"
    output_file = (out_dir + "/" +".".join(["pval"+str(pval_threshold),
            "fdr"+str(fdr_threshold)]
            + longer + ["tsv"]))

    if os.path.exists(output_file):
        print(output_file + " already created, skiping.")
        return output_file

    print("Selecting significant genes for processing")
    valid_coding_genes = {}
    valid_genes_coexpressed_with_ncRNA = {}
    total = 0
    
    for rna in genes_coexpressed_with_ncRNA.keys():
        for gene in genes_coexpressed_with_ncRNA[rna]:
            if not gene in valid_coding_genes:
                valid_coding_genes[gene] = 0
            valid_coding_genes[gene] += 1
            if not rna in valid_genes_coexpressed_with_ncRNA:
                valid_genes_coexpressed_with_ncRNA[rna] = set()
            valid_genes_coexpressed_with_ncRNA[rna].add(gene)
            total += 1
    
    #print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
    #print("Discarding coding genes with too little correlations with regulators.")
    genes_to_discard = set()
    for coding_gene in valid_coding_genes.keys():
        if valid_coding_genes[coding_gene] < K:
            genes_to_discard.add(coding_gene)

    for gene in genes_to_discard:
        del valid_coding_genes[gene]
        for rna in valid_genes_coexpressed_with_ncRNA.keys():
            if gene in valid_genes_coexpressed_with_ncRNA[rna]:
                valid_genes_coexpressed_with_ncRNA[rna].remove(gene)
    #print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))

    valid_id2gos = {}
    valid_genes_annotated_with_term = {}

    print("Reading annotation of coding genes from " + coding_gene_ontology_path)
    id2gos = onto_id2gos[ontology_type]
    for _id in id2gos.keys():
        if _id in valid_coding_genes.keys():
            valid_id2gos[_id] = id2gos[_id]
    
    genes_annotated_with_term = onto_genes_annotated_with_term[ontology_type]
    for term in genes_annotated_with_term.keys():
        for _id in genes_annotated_with_term[term]:
            if _id in valid_coding_genes.keys():
                if not term in valid_genes_annotated_with_term.keys():
                    valid_genes_annotated_with_term[term] = set()
                valid_genes_annotated_with_term[term].add(_id)

    genes_annotated_with_term2 = {}

    #print("Extending associations of terms to genes by including children")
    found = 0
    for go in valid_genes_annotated_with_term.keys():
        genes = set()
        genes.update(valid_genes_annotated_with_term[go])
        before = len(genes)
        if go in graph:
            found += 1
            childrens = get_ancestors(graph, go)
            for children_go in childrens:
                if children_go in valid_genes_annotated_with_term:
                    genes.update(valid_genes_annotated_with_term[children_go])
        genes_annotated_with_term2[go] = genes
    #print(str((float(found)/len(valid_genes_annotated_with_term.keys()))*100) + "% of the GO terms found in network.")
    valid_genes_annotated_with_term = genes_annotated_with_term2

    print("Listing possible associations between rnas and GOs")
    possible_gene_term = []
    for rna in valid_genes_coexpressed_with_ncRNA.keys():
        for go in valid_genes_annotated_with_term.keys():
            possible_gene_term.append((rna, go))
    if len(possible_gene_term) == 0:
        print("No possible association to make, under current parameters and data."
            + " Suggestions: Try a different correlation threshold or a different method.")
        return
    #print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
    #print("len(valid_genes_annotated_with_term)= " + str(len(valid_genes_annotated_with_term)))
    #print("Possible gene,term = " + str(len(possible_gene_term)))
    valid_gene_term, n_lens, M_lens, m_lens = get_valid_associations(valid_genes_coexpressed_with_ncRNA,
                                            valid_genes_annotated_with_term,
                                            possible_gene_term,
                                            min_n=min_n, min_M=min_M, min_m=min_m)

    print("Calculating p-values")
    gene_term_pvalue = parallel_pvalues(N, possible_gene_term, 
                                        valid_gene_term, n_lens, M_lens, m_lens, 
                                        threads, available_cache)
    print("Calculating corrected p-value (FDR)")
    pvalues = [pval for gene, term, pval in gene_term_pvalue]
    reject, fdrs, alphacSidak, alphacBonf = multitest.multipletests(pvalues, alpha=0.05, method='fdr_by')
    #print("Finished calculating pvalues, saving now")
    with open(tempDir + "/association_pvalue.tsv", 'w') as stream:
        for i in range(len(gene_term_pvalue)):
            if valid_gene_term[i]:
                rna, term, pvalue = gene_term_pvalue[i]
                stream.write(rna+"\t"+term+"\t"+str(pvalue)+"\t" + str(fdrs[i]) + "\n")

    print("Selecting relevant pvalues and fdr")
    relevant_pvals = []
    rna_id2gos = {}
    pval_passed = 0
    fdr_passed = 0
    for i in tqdm(range(len(gene_term_pvalue))):
        rna, term, pvalue = gene_term_pvalue[i]
        fdr = fdrs[i]
        if pvalue <= pval_threshold:
            pval_passed += 1
            if fdr <= fdr_threshold:
                fdr_passed += 1
                relevant_pvals.append((rna, term, pvalue, fdr))
                if not rna in rna_id2gos:
                    rna_id2gos[rna] = set()
                rna_id2gos[rna].add(term)
    print(str(pval_passed) + " rna->go associations passed p-value threshold ("
            + str((pval_passed/len(gene_term_pvalue))*100) + "%)")
    print(str(fdr_passed) + " rna->go associations passed fdr threshold ("
            + str((fdr_passed/len(gene_term_pvalue))*100) + "%)")
    
    print("Writing results")
    print("Output annotation is " + output_file)
    with open(output_file, 'w') as stream:
        for rna, term, pvalue, fdr in relevant_pvals:
            stream.write("\t".join([rna,term,ontology_type,str(pvalue),str(fdr)])+"\n")
    return output_file

output_files = []
for onto in ontology_types_arg:
    out_file = predict(tempDir, ontology_type=onto,
            conf_arg=conf, k_min_coexpressions=K,
            pval_threshold=pval, fdr_threshold=fdr,
            min_n=min_n, min_M=min_M, min_m=min_m)
    output_files.append((out_file,onto))

print("Writing annotation file with all ontologies")
if len(output_files) > 1:
    lines = []
    ontos = set()
    for output_file,onto_value in output_files:
        with open(output_file,'r') as stream:
            new_lines = [line for line in stream.readlines()]
            lines += new_lines
        ontos.add(onto_value)
    ontos = list(ontos)
    ontos.sort()
    ontos_str = "_".join([short_ontology_name(str(onto)) 
                        for onto in ontos])
    if len(ontos) == 3:
        ontos_str = "ALL"
    onto_dir = tempDir + "/" + ontos_str
    if not os.path.exists(onto_dir):
        os.mkdir(onto_dir)
    output_file = (onto_dir + "/" + ".".join(["pval"+str(pval),"fdr"+str(fdr),"tsv"]))
    with open(output_file,'w') as stream:
        for line in lines:
            stream.write(line)
