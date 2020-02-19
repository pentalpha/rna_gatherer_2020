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
min_threshold,step,max_threshold = [float(x) for x in sys.argv[3].split(",")]
thresholds = np.arange(min_threshold,max_threshold,step).round(3)
methods = [x for x in sys.argv[4].split(",")]
output = sys.argv[5]

print("Reading similarities")
sims = {}
with open(sims_path,'r') as stream:
    #lines = [line.rstrip("\n") for line in stream.readlines()]
    #lines_splited = [line.split("\t") for line in lines]
    #sim_tuples = [(gene1,gene2,float(sim)) for gene1,gene2,sim in lines_splited]
    for line in stream:
        gene1,gene2,sim = line.rstrip("\n").split("\t")
        sim = float(sim)
        if sim > 0.0:
            if not gene1 in sims:
                sims[gene1] = {}
            if not gene2 in sims:
                sims[gene2] = {}
            #sims[gene1][gene2] = sim
            sims[gene2][gene1] = sim

print("Reading correlations")
sim_correlations = []
with open(correlations_path, 'r') as stream:
    #line = stream.readline()
    for line in stream:
        cells = line.rstrip("\n").split("\t")
        if cells[0] in sims:
            gene_sims = sims[cells[0]]
            if cells[1] in gene_sims:
                if float(cells[2]) >= thresholds[0]:
                    sim_correlations.append((cells[-1],float(cells[2]),gene_sims[cells[1]]))
        else:
            sim_correlations.append((cells[-1],float(cells[2]),0.0))

print("Sorting " + str(len(sim_correlations)) + " items")
sim_correlations.sort(key=lambda x: x[1])
current_threshold = 0
locations = []
for i in tqdm(range(len(sim_correlations))):
    method, corr, sim = sim_correlations[i]
    if corr >= thresholds[current_threshold]:
        locations.append(i)
        current_threshold += 1
        if current_threshold >= len(thresholds):
            break

print(str(locations))

results = []

for i in tqdm(range(len(locations))):
    loc = locations[i]
    threshold = thresholds[i]
    print("Threshold " + str(threshold))
    sub_vec = sim_correlations[loc:]
    for method in methods:
        print("Method " + method)
        total = 0.0
        n = 0.0
        for m,corr,sim in sub_vec:
            if m == method:
                total += sim
                n += 1
        avg = total/n
        results.append((method,threshold,n,avg))

print(str(results))

data_points = []

for method in methods:
    x = []
    y = []
    for m,threshold,n,avg in results:
        if m == method:
            x.append(threshold)
            y.append(avg)
    data_points.append([method,x,y])

print(str(data_points))

f, ax = plt.subplots(figsize=(20, 6))
plt.style.use('seaborn-whitegrid')
plt.title('Relation between correlation coefficient and GO similarity')
for method, thresholds, avgs in data_points:
    ax.plot(thresholds, avgs,label=method)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Average similarity between genes with correlation above threshold')
plt.xlabel('Correlation coefficient threshold')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output+".png",bbox_inches='tight')