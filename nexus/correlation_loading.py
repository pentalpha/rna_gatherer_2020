import os
from tqdm import tqdm
from nexus.functional_prediction import *
import multiprocessing
import pandas as pd

def readlines(stream,n=200000):
    line = stream.readline()
    lines = []
    while n > 0 and line:
        lines.append(line.rstrip("\n").split("\t"))
        line = stream.readline()
        n -= 1
    return lines

def filepath_to_metric(filepath):
    return input_file.split("/")[-1].split(".")[0]

def discover_right_files(output_dir, correlations):
    metrics = {}
    for gene_a,gene_b,coef,metric in correlations:
        if not metric in metrics:
            metrics[metric] = []
        metrics[metric].append((gene_a,gene_b,coef))
    return [(output_dir+"/correlations/"+metric+".tsv", coefs) 
            for metric,coefs in metrics.items()]

def write_unmixed_correlations(output_dir, correlations):
    unmixed_correlations = discover_right_files(output_dir, correlations)
    for filepath, coefs in unmixed_correlations:
        with open(filepath,'a+') as stream:
            stream.write("\n".join(["\t".join([a,b,c]) for a,b,c in coefs])+"\n")

def unmix(output_dir):
    legacy_correlations_file = output_dir + "/correlated.tsv"
    with open(legacy_correlations_file,'r') as stream:
        total_lines = 0
        lines = readlines(stream)
        while len(lines) > 0:
            write_unmixed_correlations(output_dir,lines)
            total_lines += len(lines)
            print("Unmixed " + str(total_lines))
            lines = readlines(stream)
    #os.remove(legacy_correlations_file)

def read_correlations(input_file, metric, threshold,
                    correlation_values, coding_genes, genes_coexpressed_with_ncRNA):
    #print(input_file + "'s metric is " + metric)
    with open(input_file,'r') as stream:
        lines = 0
        invalid_lines = 0
        for raw_line in stream.readlines():
            cells = raw_line.rstrip("\n").split("\t")
            if len(cells) == 3:
                #correlation_lines.append([cells[0].upper(),cells[1].upper(),float(cells[2]),cells[3]])
                gene = cells[0]
                rna = cells[1]
                corr = cells[2]
                corr_val = float(corr)
                if (corr_val >= threshold or corr_val <= -threshold):
                    if not gene in coding_genes:
                        coding_genes[gene] = 0
                    coding_genes[gene] += 1
                    
                    if not rna in genes_coexpressed_with_ncRNA:
                        genes_coexpressed_with_ncRNA[rna] = set()
                    genes_coexpressed_with_ncRNA[rna].add(gene)
                    correlation_values[(rna,gene,metric)] = corr
            else:
                invalid_lines += 1
            lines += 1
        if lines == 0:
            print("Fatal error, no correlations could be loaded from "
                    + input_file + "\n(The file may be "
                    + "corrupted or just empty)")
            os.remove(input_file)
            quit()

def load_correlations(output_dir, metrics, thresholds, 
    correlation_values, coding_genes, genes_coexpressed_with_ncRNA):
    correlations_dir = output_dir + "/correlations"
    if not os.path.exists(correlations_dir):
        os.mkdir(correlations_dir)
    legacy_correlations_file = output_dir + "/correlated.tsv"
    if os.path.exists(legacy_correlations_file):
        print("Unmixing legacy correlations file.")
        unmix(output_dir)
    
    missing_metrics = []
    for metric in tqdm(metrics):
        input_file = correlations_dir + "/" + metric + ".tsv"
        if os.path.exists(input_file):
            read_correlations(input_file, metric, thresholds[metric],
                correlation_values, coding_genes, genes_coexpressed_with_ncRNA)  
        else:
            missing_metrics.append(metric)
            print("Cannot find the following correlations file: " + input_file)
    return missing_metrics

def save_correlations(output_dir, corrs_by_metric):
    correlations_dir = output_dir + "/correlations"
    if not os.path.exists(correlations_dir):
        os.mkdir(correlations_dir)
    for metric, corrs in corrs_by_metric.items():
        output_file = correlations_dir + "/" + metric + ".tsv"
        with open(output_file,'a+') as stream:
            for gene_a,gene_b,corr in corrs:
                stream.write(gene_a+"\t"+gene_b+"\t"+str(corr)+"\n")

def find_correlated(reads, regulators, tempDir, methods, minimum_coef, 
                    threads, separate_regulators = False):
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
            args=(i, parcial_df, regulators, showing, 
            methods, minimum_coef, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(start) + " to " + str(end))
        p.start()
        if showing:
            showing = False
        last = end

        last_pid = i
        if end == limit:
            break
    if threads > 1:
        if end < limit:
            parcial_df = reads.iloc[end:limit]
            p = multiprocessing.Process(target=func, 
                args=(last_pid+1, parcial_df, regulators, showing, 
                methods, minimum_coef, return_dict, ))
            processes.append(p)
            #print("Spawned process from gene " + str(end) + " to " + str(limit))
            p.start()
    else:
        parcial_df = reads
        p = multiprocessing.Process(target=func, 
            args=(last_pid+1, parcial_df, regulators, showing, 
            methods, minimum_coef, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(end) + " to " + str(limit))
        p.start()

    for p in processes:
        p.join()
    method_corrs = {}
    for metric_name in methods:
        method_corrs[metric_name] = []
    print("Merging results")
    corrs_found = 0
    for value in return_dict.values():
        for coding_name, noncoding_name, corr, method in value:
            method_corrs[method].append((coding_name, noncoding_name, corr))
        corrs_found += len(value)

    print(str(corrs_found) + " correlation pairs found.")
    save_correlations(tempDir, method_corrs)
    '''output = tempDir+"/correlated.tsv"
    with open(output,'a+') as stream:
        for coding_name, noncoding_name, corr, method in coding_noncoding_pairs:
            stream.write("\t".join([coding_name,noncoding_name,str(corr),method]) + "\n")'''
    manager.shutdown()
    del manager

def calc_correlation_in_chunks(count_reads_path, regulators_path, benchmarking, available_size,
    method, tempDir, minimum_coef, threads):

    regulators_max_portion = 0.4

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
            regulators.append(line.rstrip("\n"))

    print("Separating regulators from regulated.")
    mask = reads[reads.columns[0]].isin(regulators)
    regulators_reads = reads.loc[mask]

    non_regulators_reads = reads
    if not benchmarking:
        non_regulators_reads = reads.loc[~mask]
    print(str(len(non_regulators_reads)) + " regulated.")
    print(str(len(regulators_reads)) + " regulators.")
    print("Generating correlations.")

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
        find_correlated(df, regulators_reads, tempDir, method, minimum_coef, 
            threads, separate_regulators = benchmarking)
        i += 1