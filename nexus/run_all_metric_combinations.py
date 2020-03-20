from itertools import combinations
import subprocess

def runCommand(cmd, print_cmd=True):
    if print_cmd:
            print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

metrics = ['MIC', 'DC', 'SPR', 'PRS', 'FSH','SOB']
ns = [6,1,5,2,4,3]
#ns = [1,2,3,4,5,6]
def combs_n(n):
    combs = set()
    for vals in combinations(metrics,n):
        new_set = list(set(vals))
        new_set.sort()
        combs.add(",".join(new_set))
    return list(combs)

comb_n = [combs_n(n) for n in ns]

base_command = ("cd ../ && "
            +"python predict.py -cr test_data/counts/mus_musculus_tpm.tsv"
            +" -reg test_data/lnc_list/mus_musculus_lncRNA.txt"
            +" -ann test_data/annotation/mgi_genes_annotation.tsv"
            +" -o output/mgi_tpm_combs -conf 3,4,5,6"
            +" -ont molecular_function,biological_process"
            +" -met [METRICS] -p 5")

for combs_list in comb_n:
    print("combinations of " + str(len(combs_list[0].split(','))) 
          + " elements: (" + str(len(combs_list)) + ")")
    for comb in combs_list:
        ms = comb.split(",")
        ms.sort()
        #print("\t"+"-".join(ms))
        runCommand(base_command.replace("[METRICS]",",".join(ms)))