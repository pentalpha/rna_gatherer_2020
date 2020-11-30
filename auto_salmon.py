import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
#from gatherer.util import *
import math
from scipy.stats.stats import pearsonr
import multiprocessing
import sys

usage = ("python find_tissue_specifics.py <transcriptome_path> "
        +"<sample_list_path> <output_dir> <threads>\n"
        +"Lines in the sample list file should be in the following form\n:"
        +"tissue_name\tsample_1_p1,sample_1_p2,sample_1_singles;sample_2_singles;sample_n_p1,sample_n_p2\n"
        +"\tSamples are separated by ';' and files are separated by ','"
)

if len(sys.argv) != 5:
    print("Example usage:")
    print(usage)

fasta_path = sys.argv[1]
sample_list_path = sys.argv[2]
output_dir = sys.argv[3]
threads = int(sys.argv[4])

def run_salmon(salmon_path, index_path, tissues, threads, output_dir):
    if threads > 30:
        threads = 30
        print("Limiting salmon processes to 30 cores.")
    tmpDir = output_dir
    threads_per_process = max(3, math.ceil(float(threads) / float(len(tissues))))
    semaphor_slots = max(1, threads / threads_per_process)
    print("Threads available: " + str(threads)
        + "; Threads p/ process: " + str(threads_per_process)
        + "; Semaphor slots: " + str(semaphor_slots))

    semaphore = threading.Semaphore(semaphor_slots)
    thread_codes = {}
    thread_cmds = {}
    thread_objs = []

    def run_with_semaphore(cmd):
        with semaphore:
            return runCommand(cmd)

    def thread_function(tissue):
        name, sample_list = tissue
        pair1 = ""
        pair2 = ""
        single_only = ""

        for sample in sample_list:
            if len(sample) == 1:
                single_only += " " + sample[0]
            elif len(sample) >= 2:
                pair1 += " " + sample[0]
                pair2 += " " + sample[1]
                if len(sample) == 3:
                    single_only += " " + sample[2]

        sample_str = ""
        single_sample_str = ""
        if len(pair1) > 0:
            sample_str = "-1" + pair1 + " -2" + pair2 + " "
        if len(single_only) > 0:
            single_sample_str += "-r"+single_only

        stdout = tmpDir + "/salmon_"+name+".stdout"
        stderr = tmpDir + "/salmon_"+name+".stderr"
        out_dir = tmpDir+"/salmon_"+name
        if os.path.exists(out_dir):
            runCommand("rm -Rf " + out_dir)
        cmd = " ".join([salmon_path,"quant","-l A",sample_str,"-o", out_dir,
                    "-p",str(int(threads_per_process)),"--validateMappings -i",index_path,
                    ">", stdout, "2>", stderr])
        cmd_final = cmd
        '''if len(single_sample_str) > 0:
            stdout_single = tmpDir + "/salmon_"+name+".singles.stdout"
            stderr_single = tmpDir + "/salmon_"+name+".singles.stderr"
            out_dir_single = tmpDir+"/salmon_"+name+"singles."
            if os.path.exists(out_dir_single):
                runCommand("rm -Rf " + out_dir_single)
            cmd2 = " ".join([salmon_path,"quant","-l A",single_sample_str,"-o", out_dir_single,
                        "-p",str(int(threads_per_process)),"--validateMappings -i",index_path,
                        ">", stdout_single, "2>", stderr_single])
            cmd_final = cmd+" && "+cmd2'''

        thread_cmds[name] = cmd_final
        thread_codes[name] = run_with_semaphore(cmd_final)
        if thread_codes[name] != 0:
            print(out_dir + " failed, removing it.")
            runCommand("rm -Rf " + out_dir)
            '''if len(single_sample_str) > 0:
                runCommand("rm -Rf " + out_dir_single)'''

    for i in range(len(tissues)):
        tissue_samples = tissues[i]
        print("Creating thread for " + tissue_samples[0])
        thread_objs.append(threading.Thread(target=thread_function, args=(tissue_samples, )))

    for th in thread_objs:
        th.start()

    for th in thread_objs:
        th.join()
        print("Thread finished")

    some_failed = False
    for key in thread_codes.keys():
        if thread_codes[key] != 0:
            print(key + " failed with code " + str(thread_codes[key]))
            print(key + " trying again " + str(thread_cmds[key]))
            #some_failed = True
            try_cmd = thread_cmds[key]
            code = runCommand(try_cmd)
            if code != 0:
                print("Could definetly not run quantification for " + key)
                some_failed = True
    if some_failed:
        return False
    else:
        return True

def merge_quants(salmon_path, input_dir, output_dir):
    salmon_dirs = get_subdirs(input_dir)
    if len(salmon_dirs) == 0:
        print("No salmon dirs in " + input_dir + ", skiping...")
        return True
    quants = []
    names = []
    for dir in salmon_dirs:
        name = dir.split("salmon_")[-1]
        quants.append(dir)
        names.append(name)
    quants_arg = "--quants " + " ".join(quants) + " --names " + " ".join(names)
    out_reads = output_dir + "/raw_reads.tsv"
    out_tpm = output_dir + "/tpm.tsv"
    cmd = " ".join([salmon_path, 'quantmerge', quants_arg, "--column numreads", "-o", out_reads])
    cmd2 = " ".join([salmon_path, 'quantmerge', quants_arg, "--column tpm", "-o", out_tpm])
    code1 = runCommand(cmd)
    code2 = runCommand(cmd2)
    return code1 == 0 and code2 == 0

def quantify_coding(fasta_path, sample_list, output_dir):
    runCommand("mkdir " + output_dir)
    indexing_dir = output_dir + "/indexing/"
    runCommand("mkdir " + indexing_dir)
    quantifying_dir = output_dir + "/quantifying/"
    runCommand("mkdir " + quantifying_dir)

    path = fasta_path
    index = indexing_dir + "/index"
    cmd = " ".join(['salmon',"index","-t",path,"-i",index,"-k",str(31)])
    code = runCommand(cmd)
    #code = 0
    if code != 0:
        return False

    success1 = run_salmon('salmon', index, sample_list,
                        threads, quantifying_dir)
    success1 = True
    if not success1:
        print("Some of the quantifications have to have failed")

    success2 = merge_quants('salmon', quantifying_dir, output_dir)
    return success1 and success2

def read_tissues(sample_list_path):
    tissue_list = []
    with open(sample_list_path, 'r') as stream:
        for raw_line in stream.readlines():
            cells = raw_line.rstrip("\n").split("\t")
            if len(cells) != 2:
                print("All lines must have two collumns: ")
                print(raw_line)
                break
            name = cells[0]
            sample_list = []
            samples = cells[1].split(";")
            print("Adding to " + name + ":")
            for sample in samples:
                sample_files = sample.split(",")
                print("\t" + " ".join(sample_files))
                sample_list.append(sample_files)
            tissue_list.append((name, sample_list))
    return tissue_list

tissue_list = read_tissues(sample_list_path)
for name, samples in tissue_list:
    print(name+"={"+str(samples)+"}")
success = quantify_coding(fasta_path, tissue_list, output_dir)
