import sys
from nexus.functional_prediction import *
from nexus.util import *
#from bioinfo import *
import obonet
import networkx
import numpy as np
import argparse
import multiprocessing
from config import configs
#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

print("Available cache memory: " + str(int(get_cache()/1024)) + "KB")

all_methods = ["MIC","DC","PRS","SPR"]
all_ontologies = ["molecular_function","cellular_component","biological_process"]

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

    ap.add_argument("-met", "--method", required=False,
        default="MIC", help=("Correlation coefficient calculation method: MIC (default), "+
        "DC (Distance Correlation), SPR (Spearman Coefficient) or PRS (Pearson Coefficient)."))
    ap.add_argument("-bh", "--benchmarking", required=False,
        default=False, help=("Enables the (much slower) leave one out strategy: regulators can be regulated too."
        +"Default: False."))
    default_threads = max(2, multiprocessing.cpu_count()-1)
    ap.add_argument("-p", "--processes", required=False,
        default=default_threads, help=("CPUs to use. Default: " + str(default_threads)+"."))

    ap.add_argument("-th", "--threshold", required=False,
        default=0.85, help=("Correlation coefficient threshold: 0.5 <= x <= 0.99. Default=0.85."))
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
    ap.add_argument("-H", "--high-precision", required=False,
        default=False, help=("Set of arguments to find only high quality annotations. "
        +" Number of annotations is greatly reduced. Default: False.\n\t"
        +"Overrides other arguments: threshold=0.85, K=2, pval=0.05, fdr=0.05, n=3, M=8, m=1"))

    ap.add_argument("-ch", "--cache-usage", required=False,
        default=0.6, help=("Portion of the cache memory to use for storing the counts table."))
    ap.add_argument("-st", "--threshold-step", required=False,
        default=None, help=("Run the prediction for many thresold values, up to a maximum. Example: "+
                            "step=0.01,maximum=0.95: 0.01,0.95"))

    return vars(ap.parse_args())

cmdArgs = getArgs()
count_reads_path = cmdArgs["count_reads"]
regulators_path = cmdArgs["regulators"]
coding_gene_ontology_path = cmdArgs["annotation"]
tempDir = cmdArgs["output_dir"]
go_path = confs["go_obo"]

ontology_type = cmdArgs["ontology_type"].split(",")
if ontology_type[0] == "ALL":
    ontology_type = all_ontologies

method = cmdArgs["method"].split(",")
if method[0] == "ALL":
    method = all_methods
benchmarking = cmdArgs["benchmarking"]
threads = int(cmdArgs["processes"])

thresholds = [float(th_str) for th_str in cmdArgs["threshold"].split(",")]
min_threshold = thresholds[0]
for th in thresholds:
    if th < min_threshold:
        min_threshold = th

threshold_step = cmdArgs["threshold_step"]
if threshold_step != None:
    parts = threshold_step.split(",")
    step = float(parts[0])
    maximum = float(parts[1])
    thresholds = [float(x) for x in list(np.arange(min_threshold,maximum,step).round(3))]

print("Correlation coefficient thresholds: " + str(thresholds))

pval = float(cmdArgs["pvalue"])
fdr = float(cmdArgs["fdr"])
min_n = int(cmdArgs["min_n"])
min_M = int(cmdArgs["min_M"])
min_m = int(cmdArgs["min_m"])
K = int(cmdArgs["k_min_coexpressions"])

cache_usage = float(cmdArgs["cache_usage"])
regulators_max_portion = 0.4

high_precision = {"threshold": [0.94], "K": 2, "pval":0.05,"fdr":0.05,"n":3,"M":8,"m":1}
is_high = cmdArgs["high_precision"]
if is_high:
    thresholds = high_precision["threshold"]
    min_threshold = thresholds[0]
    pval = high_precision["pval"]
    fdr = high_precision["fdr"]
    min_n = high_precision["n"]
    min_M = high_precision["M"]
    min_m = high_precision["m"]
    K = high_precision["K"]

def find_correlated(reads, regulators, tempDir, methods, separate_regulators = False):
    coding_noncoding_pairs = []
    func = try_find_coexpression_process
    if separate_regulators:
        print("Running 'leave one out' (benchmarking) mode.")
        func = leave_one_out
    genes_per_process = int(len(reads) / threads)
    limit = len(reads)-1
    last = -1
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    processes = []
    showing = True
    last_pid = 0
    for i in range(threads-1):
        start = last+1
        end = start + genes_per_process
        if end >= limit:
            end = limit
        parcial_df = reads.iloc[start:end]
        p = multiprocessing.Process(target=func, 
            args=(i, parcial_df, regulators, showing, methods, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(start) + " to " + str(end))
        p.start()
        if showing:
            showing = False
        last = end

        last_pid = i
        if end == limit:
            break
    if end < limit:
        parcial_df = reads.iloc[end:limit]
        p = multiprocessing.Process(target=func, 
            args=(last_pid+1, parcial_df, regulators, showing, methods, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(end) + " to " + str(limit))
        p.start()

    for p in processes:
        p.join()

    print("Merging results")
    for value in return_dict.values():
        coding_noncoding_pairs += value

    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    output = tempDir+"/correlated.tsv"
    with open(output,'a+') as stream:
        for coding_name, noncoding_name, corr, method in coding_noncoding_pairs:
            stream.write("\t".join([coding_name,noncoding_name,str(corr),method]) + "\n")
    manager.shutdown()
    del manager

print("Ontology type is " + str(ontology_type))
print("Method is " + str(method))

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

if not os.path.exists(tempDir):
    os.mkdir(tempDir)

print("Reading regulators")
regulators = []
with open(regulators_path,'r') as stream:
    for line in stream.readlines():
        regulators.append(line.rstrip("\n"))

correlations_file_path = tempDir + "/correlated.tsv"
print("Loading correlations from " + correlations_file_path + ".")

if not os.path.exists(correlations_file_path):
    print("Separating regulators from regulated.")
    mask = reads[reads.columns[0]].isin(regulators)
    regulators_reads = reads.loc[mask]

    non_regulators_reads = reads
    if not benchmarking:
        non_regulators_reads = reads.loc[~mask]
    print(str(len(non_regulators_reads)) + " regulated.")
    print(str(len(regulators_reads)) + " regulators.")
    print("Generating correlations.")
    
    available_size = get_cache(usage=cache_usage)

    max_for_regulators = available_size*regulators_max_portion
    regs_size = getsizeof(regulators_reads)
    regulator_dfs = [regulators_reads]
    if regs_size > max_for_regulators:
        regulator_dfs = split_df_to_max_mem(regulators_reads, max_for_regulators)
        available_size -= max_for_regulators
    else:
        available_size -= regs_size

    dfs = split_df_to_max_mem(non_regulators_reads, available_size)
    
    print("Chunks for regulated: " + str(len(dfs)) 
    + "\nChunks for regulators: " + str(len(regulator_dfs)))

    df_pairs = []
    for df in dfs:
        for regulator_df in regulator_dfs:
            df_pairs.append((df,regulator_df))

    i = 1
    for df,regulator_df in tqdm(df_pairs):
        print("Processing chunk " + str(i) + "/" + str(len(df_pairs)))
        find_correlated(df, regulators_reads, tempDir, method, separate_regulators = benchmarking)
        i += 1

coding_genes = {}
    #keys = rnas, values = sets of genes
genes_coexpressed_with_ncRNA = {}
#keys = (rna,gene), values = (correlation,method)
correlation_values = {}
#correlation_lines = []
with open(correlations_file_path,'r') as stream:
        lines = 0
        invalid_lines = 0
        for raw_line in stream.readlines():
            cells = raw_line.rstrip("\n").split("\t")
            if len(cells) == 4:
                #correlation_lines.append([cells[0].upper(),cells[1].upper(),float(cells[2]),cells[3]])
                gene = cells[0]
                rna = cells[1]
                corr = cells[2]
                corr_val = float(corr)
                if corr_val >= min_threshold or corr_val <= -min_threshold:
                    m = cells[3]
                    if m in method:
                        if not gene in coding_genes:
                            coding_genes[gene] = 0
                        coding_genes[gene] += 1
                        
                        if not rna in genes_coexpressed_with_ncRNA:
                            genes_coexpressed_with_ncRNA[rna] = set()
                        genes_coexpressed_with_ncRNA[rna].add(gene)
                        correlation_values[(rna,gene,m)] = corr
            else:
                invalid_lines += 1
            lines += 1
        print(str(float(invalid_lines)/lines)
            + " lines without proper number of columns (4 columns)")
print("correlation_values = "+str(len(correlation_values.keys())))
print("genes_coexpressed_with_ncRNA = "+str(len(genes_coexpressed_with_ncRNA.keys())))
print("coding_genes = "+str(len(coding_genes.keys())))
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

def predict(tempDir,ontology_type="molecular_function",current_method="MIC",
            threshold=0.9,benchmarking=False,k_min_coexpressions=1,
            pval_threshold=0.05,fdr_threshold=0.05,
            min_n=5, min_M=5, min_m=1):
    print("Current method = " + current_method + ", threshold = " + str(threshold)
            + ", pvalue = " + str(pval_threshold) + ", fdr = " + str(fdr_threshold))
    K = k_min_coexpressions
    
    print("Selecting significant genes for processing")
    valid_coding_genes = {}
    valid_genes_coexpressed_with_ncRNA = {}
    total = 0
    correct_method = 0
    valid_corr = 0
    for rna in genes_coexpressed_with_ncRNA.keys():
        for gene in genes_coexpressed_with_ncRNA[rna]:
            key = (rna,gene,current_method)
            if key in correlation_values.keys():
                corr = float(correlation_values[key])
                correct_method += 1
                if corr >= threshold or corr <= -threshold:
                    valid_corr += 1
                    if not gene in valid_coding_genes:
                        valid_coding_genes[gene] = 0
                    valid_coding_genes[gene] += 1
                    if not rna in valid_genes_coexpressed_with_ncRNA:
                        valid_genes_coexpressed_with_ncRNA[rna] = set()
                    valid_genes_coexpressed_with_ncRNA[rna].add(gene)
            total += 1
    '''print(str((float(correct_method)/total)*100) + "% of coexpressions with correct method ("+method+").")
    print(str((float(valid_corr)/total)*100) + "% of coexpressions with relevant coefficient.")
    print("valid_coding_genes = "+str(len(valid_coding_genes.keys())))
    print("valid_genes_coexpressed_with_ncRNA = "+str(len(valid_genes_coexpressed_with_ncRNA.keys())))'''
    print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
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
    '''print("valid_coding_genes = "+str(len(valid_coding_genes.keys())))
    print("valid_genes_coexpressed_with_ncRNA = "+str(len(valid_genes_coexpressed_with_ncRNA.keys())))'''
    print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))

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
    print("len(valid_genes_annotated_with_term)= " + str(len(valid_genes_annotated_with_term)))
    '''print("id2gos = "+str(len(id2gos.keys())))
    print("valid_id2gos = "+str(len(valid_id2gos.keys())))
    print("genes_annotated_with_term = "+str(len(genes_annotated_with_term.keys())))
    print("valid_genes_annotated_with_term = "+str(len(valid_genes_annotated_with_term.keys())))'''

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
        '''if len(genes) > before:
            print(go + " genes increased by " + str(len(genes)-before))'''
    #print(str((float(found)/len(valid_genes_annotated_with_term.keys()))*100) + "% of the GO terms found in network.")
    valid_genes_annotated_with_term = genes_annotated_with_term2

    print("Listing possible associations between rnas and GOs")
    possible_gene_term = []
    for rna in valid_genes_coexpressed_with_ncRNA.keys():
        for go in valid_genes_annotated_with_term.keys():
            possible_gene_term.append((rna, go))
    print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
    print("len(valid_genes_annotated_with_term)= " + str(len(valid_genes_annotated_with_term)))
    print("Possible gene,term = " + str(len(possible_gene_term)))
    valid_gene_term, n_lens, M_lens, m_lens = get_valid_associations(valid_genes_coexpressed_with_ncRNA,
                                            valid_genes_annotated_with_term,
                                            possible_gene_term,
                                            min_n=min_n, min_M=min_M, min_m=min_m)

    print("Calculating p-values")
    #N = len(coding_genes) #population size
    gene_term_pvalue = parallel_pvalues(N, possible_gene_term, 
                                        valid_gene_term, n_lens, M_lens, m_lens, 
                                        threads, get_cache(usage=cache_usage))
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
    ontology_type_mini = ontology_type.replace("biological_process","BP")
    ontology_type_mini = ontology_type_mini.replace("molecular_function","MF")
    ontology_type_mini = ontology_type_mini.replace("cellular_component","CC")
    run_mode = "LNC"
    if benchmarking:
        run_mode = "BENCH"
    out_dir = tempDir+"/"+ontology_type_mini
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    longer = []
    if K != 1 or min_n != 1 or min_M != 1 or min_m != 1:
        longer = ["K"+str(K),"n"+str(min_n),"M"+str(min_M),"m"+str(min_m)]
    output_file = (out_dir + "/" +".".join([current_method,"th"+str(threshold),
            run_mode,"pval"+str(pval_threshold),"fdr"+str(fdr_threshold)]
            + longer + ["tsv"])).replace(".LNC.",".")
    print("Output annotation is " + output_file)
    with open(output_file, 'w') as stream:
        for rna, term, pvalue, fdr in relevant_pvals:
            stream.write("\t".join([rna,term,ontology_type,str(pvalue),str(fdr)])+"\n")

for onto in ontology_type:
    for meth in method:
        for th in thresholds:
            predict(tempDir,ontology_type=onto,current_method=meth,
                    threshold=th,benchmarking=benchmarking,k_min_coexpressions=K,
                    pval_threshold=pval,fdr_threshold=fdr,
                    min_n=min_n, min_M=min_M, min_m=min_m)
