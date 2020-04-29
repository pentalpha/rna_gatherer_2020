import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral, header_to_id
from nexus.util import *
import math
from scipy.stats.stats import pearsonr
from scipy.special import comb
import multiprocessing
import networkx
import obonet
import statsmodels.stats.multitest as multitest

def read_ids2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            go_ids = gos_str.split(";")
            gos_dict[rfam_id] = set(go_ids)
    return gos_dict

def read_rfam2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split()
            rfam_id = cols[0].split(":")[-1]
            if not rfam_id in gos_dict:
                gos_dict[rfam_id] = set()
            go_str = cols[-1]
            gos_dict[rfam_id].add(go_str)
    return gos_dict

def write_id2go(filepath, gos_dict):
    with open(filepath, 'w') as stream:
        for key, gos in gos_dict.items():
            for go in gos:
                stream.write(key + "\t" + go + "\n")

def write_transcriptome(args, confs, tmpDir, stepDir):
    print("Loading annotation")
    annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
                names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Loading genome: " + args['genome_link'])
    genome_dict = seqListToDict(readSeqsFromFasta(args['genome_link']), header_to_name = header_to_id)
    transcriptome = []
    print("Creating transcriptome file")
    for index, row in annotation.iterrows():
        #print(fasta_header)
        s = genome_dict[str(row["seqname"])] #cant find key PGUA01000001.1 #TODO
        new_header = get_gff_attributes(row["attribute"])["ID"]
        from_seq = int(row["start"])
        to_seq = int(row["end"])
        begin = min(from_seq,to_seq)-1
        up_to = max(from_seq,to_seq)
        new_seq = s[begin:up_to]
        transcriptome.append((new_header, new_seq))
    print("Writing transcriptome")
    writeFastaSeqs(transcriptome, tmpDir + "/transcriptome.fasta")
    return True

def make_id2go(args, confs, tmpDir, stepDir):
    id2go_path = confs["rfam2go"]
    if os.path.exists(id2go_path):
        print("Loading ids2go associations")
        global_ids2go = read_rfam2go(id2go_path)
        print("Loading annotation")
        annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
            names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        local_ids2go = {}
        print("Associating IDs to GO terms")
        ids = []
        for index, row in annotation.iterrows():
            attrs = get_gff_attributes(row["attribute"])
            ID = attrs["ID"]
            ids.append(ID)
            if "rfam" in attrs:
                rfam_id = attrs["rfam"]
                if rfam_id in global_ids2go:
                    go_list = global_ids2go[rfam_id]
                    local_ids2go[ID] = go_list
        write_id2go(tmpDir + "/id2go.tsv", local_ids2go)
        print("Writing population: " + str(len(ids)) + " ids")
        with open(tmpDir + "/ids.txt", 'w') as stream:
            for ID in ids:
                stream.write(ID + "\n")
        return True
    else:
        print(id2go_path + " does not exist.")
        return False

def index_transcriptome(args, confs, tmpDir, stepDir):
    path = stepDir['write_transcriptome'] + "/transcriptome.fasta"
    output = tmpDir + "/transcriptome_index"
    cmd = " ".join([confs['salmon'],"index","-t",path,"-i",output,"-k",str(31)])
    code = runCommand(cmd)
    return code == 0

def separatePairedAndSingles(files):
    print("Separating:\n\t"+str(files))
    paired_1 = list()
    paired_2 = list()
    singles = list()

    while len(files) > 0:
        print("Deciding on " + files[0])
        f = files[0]

        to_repl = ''
        repl = ''
        if '1' in f:
            possible_2 = replace_last(f,'1','2')
            if possible_2 in files:
                find_third = replace_last(f, '1', '3')
                if find_third in files or find_third in singles:
                    singles.append(f)
                    print(f + " is a single.")
                else:
                    paired_1.append(f)
                    paired_2.append(possible_2)
                    print(f + " is part of a pair")
                    files.remove(possible_2)
                #files.remove(f)
        elif '2' in f:
            possible_1 = replace_last(f,'2','1')
            if possible_1 in files:
                find_third = replace_last(f, '2', '3')
                if find_third in files or find_third in singles:
                    singles.append(f)
                    print(f + " is a single.")
                else:
                    paired_1.append(possible_1)
                    paired_2.append(f)
                    print(f + " is part of a pair")
                    files.remove(possible_1)
                #files.remove(f)
        else:
            singles.append(f)
            print(f + " is a single.")
        files.remove(f)
    return paired_1, paired_2, singles

def getFastqs(dir):
    files = list(set(getFilesWith(dir, ".fq") + getFilesWith(dir, ".fastq")
            + getFilesWith(dir, ".fq.gz") + getFilesWith(dir, ".fastq.gz")))
    return separatePairedAndSingles(files)

def run_salmon(salmon_path, index_path, pairs1, pairs2, singles, threads, output_dir):
    if threads > 30:
        threads = 30
        print("Limiting salmon processes to 30 cores.")
    tmpDir = output_dir
    threads_per_process = max(3, math.ceil(float(threads) / float(len(pairs1)+len(singles))))
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

    def thread_function(fq1, fq2):
        fq_name = file_name(fq1)
        if fq2 != None:
            fq_name = file_name(fq1).replace("1", "")
        stdout = tmpDir + "/salmon"+fq_name+".stdout"
        stderr = tmpDir + "/salmon"+fq_name+".stderr"
        output_dir = tmpDir+"/salmon_"+fq_name
        cmd = " ".join([salmon_path,"quant","-l A -r", fq1,"-o", output_dir,
                    "-p",str(int(threads_per_process)),"-i",index_path,
                    ">", stdout, "2>", stderr])
        thread_cmds[fq_name] = cmd
        if fq2 != None:
            cmd = " ".join([salmon_path,"quant","-l A -1", fq1, "-2", fq2,
                    "-o",tmpDir+"/salmon_"+fq_name,"-p",str(int(threads_per_process))
                    ,"-i",index_path, ">", stdout, "2>", stderr])
            thread_cmds[fq_name] = cmd
        thread_codes[fq_name] = run_with_semaphore(cmd)
        if thread_codes[fq_name] != 0:
            runCommand("rm -Rf " + output_dir)

    for i in range(len(pairs1)):
        fq1 = pairs1[i]
        fq2 = pairs2[i]
        print("Creating thread for " + fq1)
        thread_objs.append(threading.Thread(target=thread_function, args=(fq1, fq2, )))

    for i in range(len(singles)):
        fq1 = singles[i]
        print("Creating thread for " + fq1)
        thread_objs.append(threading.Thread(target=thread_function, args=(fq1, None, )))

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

def quantify(args, confs, tmpDir, stepDir):
    index_path = stepDir['index_transcriptome']+ "/transcriptome_index"
    if "fastq_directory" in args:
        fastq_directory = args['fastq_directory']
        paired1, paired2, singles = getFastqs(fastq_directory)
        print("Pairs1:\n\t" + str(paired1))
        print("Pairs2:\n\t" + str(paired2))
        print("Singles:\n\t" + str(singles))

        success = run_salmon(confs['salmon'], index_path, paired1, paired2, singles,
                            confs['threads'], tmpDir)
        return success
    return True

def quantmerge(args, confs, tmpDir, stepDir):
    success = merge_quants(confs['salmon'], stepDir['quantify'], tmpDir)
    return success

def quantify_coding(args, confs, tmpDir, stepDir):
    if "coding_transcriptome" in args:
        indexing_dir = tmpDir + "/indexing/"
        runCommand("mkdir " + indexing_dir)
        quantifying_dir = tmpDir + "/quantifying/"
        runCommand("mkdir " + quantifying_dir)

        path = args['coding_transcriptome']
        index = indexing_dir + "/index"
        cmd = " ".join([confs['salmon'],"index","-t",path,"-i",index,"-k",str(31)])
        code = runCommand(cmd)
        if code != 0:
            return False

        if "fastq_directory" in args:
            fastq_directory = args['fastq_directory']
            paired1, paired2, singles = getFastqs(fastq_directory)
            print("Pairs1:\n\t" + str(paired1))
            print("Pairs2:\n\t" + str(paired2))
            print("Singles:\n\t" + str(singles))

            success1 = run_salmon(confs['salmon'], index, paired1, paired2, singles,
                                confs['threads'], quantifying_dir)
            if not success1:
                print("Some of the quantifications have to have failed")

            success2 = merge_quants(confs['salmon'], quantifying_dir, tmpDir)
            return success1 and success2
    else:
        print("No coding transcriptome, skiping...")
    return True

def complete_id2go(args, confs, tmpDir, stepDir):
    new_id2gos = read_ids2go(stepDir['assign_terms'] + "/new_terms.id2gos")
    old_id2gos = read_ids2go(stepDir['make_id2go'] + "/ids2go.tsv")
    for ID in new_id2gos.keys():
        if ID in old_id2gos:
            old_id2gos[ID].update(new_id2gos[ID])
        else:
            old_id2gos[ID] = new_id2gos[ID]
    write_ids2go(tmpDir+"/ncRNA.id2gos", old_id2gos)
    return True