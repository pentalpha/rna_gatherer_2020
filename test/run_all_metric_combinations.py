from itertools import combinations
import subprocess
import multiprocessing
import sys
import numpy as np

'''
usage:  
    python run_all_metric_combinations.py <prediction_dir> <counts_file> <annotation_file> <lnc_list>
'''

prediction_dir = sys.argv[1]
counts_file = sys.argv[2]
annotation_file = sys.argv[3]
lnc_list_file = sys.argv[4]
start, end = [int(x) for x in sys.argv[5].split(",")]
model_name = sys.argv[6]
threads = max(2, int(multiprocessing.cpu_count()*0.8))

def runCommand(cmd, print_cmd=True):
    if print_cmd:
            print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

metrics = ['DC', 'SPR', 'PRS', 'FSH','SOB']
ns = [5,1,4,2,3]
#ns = [1,2,3,4,5,6]
def combs_n(n):
    combs = set()
    for vals in combinations(metrics,n):
        new_set = list(set(vals))
        new_set.sort()
        combs.add(",".join(new_set))
    return list(combs)

comb_n = [combs_n(n) for n in ns]
n_predictions = 0
for combs in comb_n:
    n_predictions += len(combs)
conf_levels = list(range(start,end+1))
n_predictions = n_predictions*len(conf_levels)
confs_str = ",".join([str(x) for x in conf_levels])

print(str(n_predictions))
base_command = ("cd ../ && "
            +"python predict.py -cr " + counts_file
            +" -reg " + lnc_list_file
            +" -ann " + annotation_file
            +" -o " + prediction_dir + " -conf " + confs_str
            +" -ont molecular_function,biological_process,cellular_component"
            +" -met [METRICS] -p " + str(threads) + " -md " + model_name)

for combs_list in comb_n:
    print("combinations of " + str(len(combs_list[0].split(','))) 
          + " elements: (" + str(len(combs_list)) + ")")
    for comb in combs_list:
        ms = comb.split(",")
        ms.sort()
        #print("\t"+"-".join(ms))
        metrics_str = ",".join(ms)
        runCommand(base_command.replace("[METRICS]", metrics_str))
