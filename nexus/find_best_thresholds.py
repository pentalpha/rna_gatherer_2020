import sys
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import statistics as stc

correlations_path = sys.argv[1]
sims_path = sys.argv[2]

print("Reading similarities")
sims = {}
with open(sims_path,'r') as stream:
    #lines = [line.rstrip("\n") for line in stream.readlines()]
    #lines_splited = [line.split("\t") for line in lines]
    #sim_tuples = [(gene1,gene2,float(sim)) for gene1,gene2,sim in lines_splited]
    for line in stream:
        cells = line.rstrip("\n").split("\t")
        if len(cells) == 3:
            gene1,gene2,sim = cells
            sim = float(sim)
            if sim > 0.0:
                if not gene1 in sims:
                    sims[gene1] = {}
                if not gene2 in sims:
                    sims[gene2] = {}
                #sims[gene1][gene2] = sim
                sims[gene2][gene1] = sim
        elif len(cells) > 0:
            print("Line with too many values:\n\t"+line)

print("Reading correlations")
sim_correlations = {}
corrs_with_sim = 0
total_corrs = 0
with open(correlations_path, 'r') as stream:
    for line in stream:
        cells = line.rstrip("\n").split("\t")
        method = cells[-1]
        if cells[0] in sims:
            gene_sims = sims[cells[0]]
            if cells[1] in gene_sims:
                if float(cells[2]) >= 0.7:
                    if not method in sim_correlations:
                        sim_correlations[method] = []
                    sim_correlations[method].append((cells[0],cells[1],float(cells[2]),gene_sims[cells[1]]))
                    
print("Sorting " + str(len(sim_correlations)) + " items")
for method in list(sim_correlations.keys()):
    sim_correlations[method].sort(key=lambda x: x[2])
    #sim_correlations[method] = sim_correlations[method][::-1]
print("Sorted")

def filter_by_threshold(min_corr,corrs_vec):
    i = 0
    for corr,sim in corrs_vec:
        if corr >= min_corr:
            break
        i += 1
    return corrs_vec[i:]

def avg_similarity(min_corr,corrs_vec):
    filtered_vec = filter_by_threshold(min_corr,corrs_vec)
    if len(filtered_vec) > 0:
        sims = np.array([y for x,y in filtered_vec])
        return np.mean(sims)
    else:
        return None

def median_similarity(min_corr,corrs_vec):
    filtered_vec = filter_by_threshold(min_corr,corrs_vec)
    if len(filtered_vec) > 0:
        sims = np.array([y for x,y in filtered_vec])
        return np.quantile(sims,0.50)
    else:
        return None

def high_similarity(min_corr,corrs_vec):
    filtered_vec = filter_by_threshold(min_corr,corrs_vec)
    if len(filtered_vec) > 0:
        sims = np.array([y for x,y in filtered_vec])
        return np.quantile(sims,0.25)
    else:
        return None

def find_threshold(vec, sim_calculator, 
                initial_threshold=0.7, initial_step=0.01,
                good_similarity=0.97,execs=4):
    current_threshold = initial_threshold
    step = initial_step
    sim = sim_calculator(current_threshold,vec)
    if sim == None:
        return None
    while sim < good_similarity:
        current_threshold += step
        sim = sim_calculator(current_threshold,vec)
        if sim == None:
            return None

    if execs > 0:
        return find_threshold(vec, sim_calculator,
                            initial_threshold=current_threshold-step,
                            initial_step=step/10,
                            good_similarity=good_similarity,
                            execs=execs-1)
    else:
        return current_threshold

def count_corrs(min_value, corrs):
    if min_value == None:
        return []
    new_vec = []
    for gene_a,gene_b,corr,sim in corrs:
        if corr >= min_value:
            new_vec.append((gene_a,gene_b))
    return new_vec

result_sets = {}

for method, vec in sim_correlations.items():
    corrs_vec = [(c,d) for a,b,c,d in vec]
    print("Finding median high similarity threshold for " + method)
    avg_threshold = find_threshold(corrs_vec, avg_similarity)
    print("\tAverage threshold: " + str(avg_threshold))
    avg_set = set(count_corrs(avg_threshold, vec))
    print("\t\tNumber of correlations: " + str(len(avg_set)))
    if len(avg_set) > 0:
        result_sets[method] = avg_set
    '''median_threshold = find_threshold(corrs_vec, median_similarity)
    print("\tMedian threshold: " + str(median_threshold))
    median_set = set(count_corrs(median_threshold, vec))
    print("\t\tNumber of correlations: " + str(len(median_set)))
    high_threshold = find_threshold(corrs_vec, high_similarity,
                        initial_threshold=median_threshold)
    print("\tHigh threshold: " + str(high_threshold))
    high_set = set(count_corrs(high_threshold, vec))
    print("\t\tNumber of correlations: " + str(len(high_set)))'''
    

print("\t"+"\t".join([method for method in result_sets.keys()]))
for method_a in result_sets.keys():
    corrs_a = result_sets[method_a]
    sims = method_a + ":\t"
    for method_b in result_sets.keys():
        corrs_b = result_sets[method_b]
        intersect_len = len(corrs_a.intersection(corrs_b))
        percent = float(intersect_len) / float(len(corrs_a))
        sims += (str(percent)[0:4]+"\t")
    print(sims)
