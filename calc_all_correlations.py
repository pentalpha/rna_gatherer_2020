import sys
from nexus.functional_prediction import *
from nexus.util import getFilesWith
import argparse
import multiprocessing
from config import configs
from tqdm import tqdm
#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

display_cache = get_cache()

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-cr", "--count-reads", required=True,
        help=("Count reads table path. TSV file with where the first column must be the gene names and the " +
            " following columns are the FPKM normalized counts of each sample."))
    ap.add_argument("-o", "--output-dir", required=True, help=("Output directory."))
    default_threads = max(2, multiprocessing.cpu_count()-1)
    ap.add_argument("-p", "--processes", required=False,
        default=default_threads, help=("CPUs to use. Default: " + str(default_threads)+"."))
    ap.add_argument("-m", "--metrics", required=False, default=",".join(all_methods),
        help="Metrics to calculate. Options are FSH, SOB, MIC, DC, PRS, SPR. Default: All.")   
    ap.add_argument("-chu", "--cache-usage", required=False,
        default=0.6, help=("Portion of the cache memory to use for storing the counts table."))
    ap.add_argument("-ch", "--cache-size", required=False,
        default=display_cache, help=("Sets the size of the cache memory."
        +" Default: auto-detection of CPU cache size."))

    return vars(ap.parse_args())

cmdArgs = getArgs()
count_reads_path = cmdArgs["count_reads"]
tempDir = cmdArgs["output_dir"]
threads = int(cmdArgs["processes"])
cache_usage = float(cmdArgs["cache_usage"])
available_cache = get_cache(usage=cache_usage)
if "cache_size" in cmdArgs:
    available_cache = int(int(cmdArgs["cache_size"]) * cache_usage)
print("Available cache memory: " + str(int(available_cache)))
available_cache = int(available_cache/threads)
metrics_used = cmdArgs["metrics"].split(",")

regulators_max_portion = 0.5

if not os.path.exists(tempDir):
    os.mkdir(tempDir)

correlations_dir = tempDir + "/correlations"
if not os.path.exists(correlations_dir):
    os.mkdir(correlations_dir)

def get_temp_metric_file(base_path, metric_name):
    return open(base_path + "." + metric_name, 'a+')

def get_metric_file(metric_name):
    return open(correlations_dir + "/" + metric_name + ".tsv", 'w')

def find_correlated(reads, regulators, tempDir, method_streams):
    coding_noncoding_pairs = []
    func = calc_all
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
            args=(i, parcial_df, regulators, return_dict, ))
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
            args=(last_pid+1, parcial_df, regulators, return_dict, ))
        processes.append(p)
        #print("Spawned process from gene " + str(end) + " to " + str(limit))
        p.start()

    for p in processes:
        p.join()

    #print("Merging results")
    for value in return_dict.values():
        coding_noncoding_pairs += value

    #print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    output = tempDir+"/correlated.tsv"
    #with open(output,'a+') as stream:
    for coding_name, noncoding_name, corr, method_name in coding_noncoding_pairs:
        method_streams[method_name].write("\t".join([coding_name,noncoding_name,str(corr)]) + "\n")
    manager._process.terminate()
    manager.shutdown()
    del manager

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

regulator_names = reads[reads.columns[0]].tolist()

correlation_files = {method_name:correlations_dir+"/"+method_name+".tsv" for method_name in metrics_used}

for m,f in correlation_files.items():
    delete_if_empty(f,min_cells=3,sep="\t")

available_size = available_cache

'''max_for_regulators = available_size*regulators_max_portion
regs_size = getsizeof(regulators_reads)
regulator_dfs = [regulators_reads]
if regs_size > max_for_regulators:
    regulator_dfs = split_df_to_max_mem(regulators_reads, max_for_regulators)
    available_size -= max_for_regulators
else:
    available_size -= regs_size'''


dfs = split_df_to_max_mem(reads, int(available_size/2))
print("Splitting up dataframe into smaller DFs: " + str(len(dfs)))
for i in range(len(dfs)):
    path_to_write = tempDir + "/" + str(i) + "-counts_part.tsv"
    if not os.path.exists(path_to_write):
        dfs[i].to_csv(path_to_write, sep="\t", header=True, index=False)
del dfs

print("Reading smaller DFs")
def get_corrs_paths(df_path):
    corrs_paths = [(metric, df_path+"."+metric) for metric in metrics_used]
    for metric,p in corrs_paths:
        delete_if_empty(p)
    final_corrs_paths = []
    for metric, path in corrs_paths:
        if os.path.exists(path):
            final_corrs_paths.append((metric, path))
    if len(final_corrs_paths) != len(corrs_paths):
        for metric, path in final_corrs_paths:
            os.remove(path)
        return []
    else:
        return corrs_paths

small_df_paths = getFilesWith(tempDir, "-counts_part.tsv")
dfs = {p: pd.read_csv(p, sep="\t") for p in small_df_paths}
mini_corrs_paths = {p: get_corrs_paths(p) 
                    for p in small_df_paths}
should_run = {p: True if len(corr_paths) == len(metrics_used) else False 
                    for p, corr_paths in mini_corrs_paths.items()}

print("Processing counts")

pbar = tqdm(total=len(small_df_paths)*len(small_df_paths))
for i in range(len(small_df_paths)):
    current_df_path = small_df_paths[i]
    current_df = dfs[current_df_path]
    if not should_run[current_df_path]:
        method_streams = {metric_name:get_temp_metric_file(current_df_path,metric_name) 
                    for metric_name in metrics_used}
        for j in range(len(small_df_paths)):
            find_correlated(current_df, dfs[small_df_paths[j]], 
                tempDir, method_streams)
            pbar.update(1)
        for method_name, stream in method_streams.items():
            stream.close()
    else:
        print("Skiping " + current_df_path)
pbar.close()
        
print("Joining results")
mini_corrs_paths = {p: get_corrs_paths(p) 
                    for p in small_df_paths}
final_method_streams = {metric_name:get_metric_file(metric_name) 
                    for metric_name in metrics_used}

for p, corrs_paths in tqdm(mini_corrs_paths.items()):
    if len(corrs_paths) != len(metrics_used):
        print("Not all correlations calculated for " + p)
        quit()
    for metric, corr_path in corrs_paths:
        with open(corr_path, 'r') as corr_stream:
            lines_read = corr_stream.read()
            if lines_read[-1] != "\n":
                lines_read += "\n"
            final_method_streams[metric].write(lines_read)