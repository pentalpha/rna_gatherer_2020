import sys
import subprocess
import os
from tqdm import tqdm

def runCommand(cmd):
    print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def get_fastq_zipped_data(fqc_file):
    runCommand("unzip " + fqc_file)
    data_file = fqc_file.replace(".zip", "/fastqc_data.txt")
    lines = [line.rstrip("\n") for line in open(data_file,'r').readlines()]
    runCommand("rm -Rf " + fqc_file.replace(".zip", "/"))
    if len(lines) == 0:
        print(data_file, " is empty")
        quit()
    return lines

def get_gc(lines):
    for line in lines:
        if "%GC" in line:
            gc = line.split()[-1]
            return float(gc)
    return None

def get_n_sequences(lines):
    for line in lines:
        if "Total Sequences" in line:
            n = line.split()[-1]
            return int(n)
    return None

def get_module(lines, module_key):
    mod_lines = []
    inside = False
    for line in lines:
        if not inside:
            if module_key in line:
                inside = True
        else:
            if "END_MODULE" in line:
                break
            elif not line.startswith('#'):
                mod_lines.append(line)
    return mod_lines

def get_lens(lines):
    lens_module = get_module(lines, "Sequence Length Distribution")
    lens = {}
    base_sum = 0
    total_seqs = 0
    for line in lens_module:
        l, n = line.split("\t")
        sub_lens = l.split('-')
        n = int(float(n))/len(sub_lens)
        for sub_len in sub_lens:
            lens[int(sub_len)] = n
            total_seqs += n
            base_sum += int(sub_len)*n
    return total_seqs, lens, base_sum

def get_n(lines, total, lens):
    lens_module = get_module(lines, "Per base N content")
    n_percs = [0.0]
    for line in lens_module:
        cells = line.split("\t")
        sub_ns = cells[0].split('-')
        for x in sub_ns:
            n_percs.append(float(cells[-1])/100)
    n_totals = [0]
    for i in range(len(n_percs)):
        if i > 0:
            n_seqs = sum([total for seq_len, total in lens.items() if seq_len >= i])
            n_totals.append(n_percs[i]*n_seqs)
            print('len=',i,'n_seqs=',n_seqs,'n_percs[i]=',n_percs[i],'n_totals[-1]=',n_totals[-1])
    total_n = sum(n_totals)
    n_perc = (total_n/total)*100
    print('total_n=',total_n,'n_perc=',n_perc)
    return n_perc

def get_Q25(lines, total):
    Q25_module = get_module(lines, "Per sequence quality scores")
    Q25_total = 0
    for line in Q25_module:
        cells = line.split("\t")
        quality = int(cells[0])
        if quality >= 25:
            pass_count = int(float(cells[1]))
            Q25_total += pass_count
    return (Q25_total / total)*100

def get_fastqc_stats(fqc_file):
    print(fqc_file)
    lines = get_fastq_zipped_data(fqc_file)
    n_seqs = get_n_sequences(lines)
    if n_seqs == None:
        print("Error, could not read total sequences")
        quit()
    gc = get_gc(lines)
    if gc == None:
        print("Error, could not gc content")
        quit()
    total_seqs, lens, base_sum = get_lens(lines)
    if total_seqs != n_seqs:
        print(total_seqs, "!=", n_reads)
        quit()
    min_len = min(lens.keys())
    max_len = max(lens.keys())
    n = get_n(lines, int(base_sum), lens)
    Q25 = get_Q25(lines, n_seqs)

    return n_seqs, int(base_sum), gc, n, Q25

def sum_len_dicts(d1, d2):
    for key in d2.keys():
        if key in d1:
            d1[key] += d2[key]
        else:
            d1[key]
    return d1

def mean_stats(s1, s2):
    return s1[0], s1[1] + s2[1], (s1[2]+s2[2])/2,  (s1[3]+s2[3])/2,  (s1[4]+s2[4])/2

sample_names_file = sys.argv[1]
fastqc_dir = sys.argv[2]
outfile = sys.argv[3]

with open(outfile,'w') as out_stream:
    out_stream.write("\t".join(['Sample Name','Total Reads','Base Sum','GC (%)','N (%)','Q25(%)']) + "\n")
    with open(sample_names_file, 'r') as in_stream:
        for line in tqdm(in_stream.readlines()):
            cells = line.rstrip("\n").split()
            if len(cells) == 2:
                out_stream.write(cells[0]+" "+cells[1]+":\n")
            elif len(cells) == 3:
                sample_name = cells[0].upper().replace("_"," ")
                sample_id = cells[1]
                sample_extension = cells[2]
                left_sample = fastqc_dir+"/"+sample_id + "1" + sample_extension + "_fastqc.zip"
                right_sample = fastqc_dir+"/"+sample_id + "2" + sample_extension +"_fastqc.zip"
                if not os.path.exists(left_sample):
                    print("Could not find " + left_sample)
                else:
                    l_stats, r_stats = (get_fastqc_stats(left_sample),get_fastqc_stats(right_sample))
                    n_reads, base_sum, gc, n, Q25 = mean_stats(l_stats, r_stats)
                    line = "\t".join([sample_name, "{:.4E}".format(n_reads), "{:.4E}".format(base_sum), 
                                    str(round(gc, 3)), str(round(n, 3)), str(round(Q25, 3))])
                    out_stream.write(line+"\n")

