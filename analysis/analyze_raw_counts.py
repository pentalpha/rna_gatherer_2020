from tqdm import tqdm
import sys
import numpy as np
import math
import os

gatherer_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))+"/../")
print("RNA Gatherer dir = ", gatherer_dir)
counts_file = gatherer_dir + "/test_data/counts/gigas-raw-lnc.tsv"
outdir = gatherer_dir + "/result_data/gigas_expression_analysis/"
if not os.path.exists(outdir):
    os.mkdir(outdir)
outfile = gatherer_dir + "/result_data/gigas_expression_analysis/lnc_reads_stats.txt"

reads_aligned = []
reads_by_sample = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0]
low_reads = 0
with open(counts_file,'r') as stream:
    n_samples = len(stream.readline().split("\t"))-1
    for raw_line in tqdm(stream.readlines()):
        cells = raw_line.rstrip("\t").split('\t')
        name = cells[0]
        for i in range(len(reads_by_sample)):
            reads_by_sample[i] += float(cells[i+1])
        reads = [float(x) for x in cells[1:]]
        reads_sum = sum(reads)
        max_count = max(reads)
        if max_count < 5.0:
            low_reads += 1
        reads_aligned.append((name,reads_sum))

reads_aligned.sort(key=lambda x: x[1], reverse=True)
reads_vec = np.array([y for x, y in reads_aligned])
total_reads_aligned = sum(reads_vec)
mean = np.median(reads_vec)
p25 = np.percentile(reads_vec,25)
p75 = np.percentile(reads_vec,75)
avg_reads_per_sample = np.mean(reads_by_sample)
min_reads_per_sample = min(reads_by_sample)
max_reads_per_sample = max(reads_by_sample)

with open(outfile,'w') as stream:
    stream.write("lncRNAs with estimated expression: "+str(len(reads_aligned))+"\n")
    stream.write("Total Reads Aligned: "+str(total_reads_aligned)+"\n")
    stream.write("Median Reads Aligned: "+str(mean)+"\n")
    stream.write("\t25% Percentile: "+str(p25)+"\n")
    stream.write("\t75% Percentile: "+str(p75)+"\n")
    stream.write(str((len(reads_aligned)-low_reads)/len(reads_aligned)*100)+ "% of reads with expression >= 5\n")
    stream.write("Avg Reads Per Sample: "+str(avg_reads_per_sample)+"\n")
    stream.write("\tMin Reads Per Sample: "+str(min_reads_per_sample)+"\n")
    stream.write("\tMax Reads Per Sample: "+str(max_reads_per_sample)+"\n")

    stream.write("\nTop Expressed lncRNA:\n")
    for lnc_name, reads in reads_aligned[0:10]:
        percent = round((reads/total_reads_aligned)*100,4)
        stream.write(lnc_name+"\t"+str(reads)+"\t"+str(percent)+"%\n")

    stream.write("\n")
    for i in range(len(reads_aligned)):
        lnc_name, reads = reads_aligned[i]
        pos = i + 1
        if not "Transcript" in lnc_name:
            stream.write("In the position "+str(pos)+" is: " + lnc_name + "\t"+str(reads)+"\n")
