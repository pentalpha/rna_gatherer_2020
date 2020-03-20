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
#min_threshold,step,max_threshold = [float(x) for x in sys.argv[3].split(",")]
#thresholds = np.arange(min_threshold,max_threshold,step).round(3)
max_coefs = int(sys.argv[3])
coefs_step = int(sys.argv[4])
methods = [x for x in sys.argv[5].split(",")]
output = sys.argv[6]

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
for metric in methods:
    sim_correlations[metric] = []
corrs_with_sim = 0
total_corrs = 0

with open(correlations_path, 'r') as stream:
    #line = stream.readline()
    for line in stream:
        cells = line.rstrip("\n").split("\t")
        if cells[0] in sims:
            gene_sims = sims[cells[0]]
            if cells[1] in gene_sims:
                #if float(cells[2]) >= 0.6:
                sim_correlations[cells[-1]].append((float(cells[2]),gene_sims[cells[1]]))
                corrs_with_sim += 1
        total_corrs += 1
        #else:
        #    sim_correlations.append((cells[-1],float(cells[2]),0.0))
for metric in methods:
    if len(sim_correlations[metric]) == 0:
        del sim_correlations[metric]

print(str(float((corrs_with_sim)/float(total_corrs))*100.0) 
        + "% of the correlations had a corresponding similarity")
for metric, vec in sim_correlations.items():
    print("Sorting " + str(len(vec)) + " items")
    vec.sort(key=lambda x: x[0])
'''current_threshold = 0
locations = []
for i in tqdm(range(len(sim_correlations))):
    method, corr, sim = sim_correlations[i]
    if corr >= thresholds[current_threshold]:
        locations.append(i)
        current_threshold += 1
        if current_threshold >= len(thresholds):
            break'''

#print(str(locations))

results = []
total = 0
while total < max_coefs:
    total += coefs_step
    for metric, vec in sim_correlations.items():
        if len(vec) >= total:
            tops = vec[-total::]
            top_sims = np.array([sim for corr,sim in tops])
            th = tops[0][0]
            avg = top_sims.mean()
            results.append((metric,total,avg,th))
    print(str(total) + "/" + str(max_coefs))
'''    loc = locations[i]
    threshold = thresholds[i]
    #print("Threshold " + str(threshold))
    sub_vec = sim_correlations[loc:]
    method_vecs = {"PRS":[],"SPR":[],"MIC":[],"DC":[],"SOB":[],"FSH":[]}
    for m,corr,sim in sub_vec:
        method_vecs[m].append(sim)
    for method in methods:
        method_vec = np.array(method_vecs[method])
        n = len(method_vec)
        if n > 0:
            if method != "SOB" or threshold >= 0.9:
                lower_quantile = np.quantile(method_vec,0.25)
                middle_quantile = np.quantile(method_vec,0.50)
                upper_quantile = np.quantile(method_vec,0.75)
                avg = np.mean(method_vec)
                results.append((method,threshold,n,lower_quantile,middle_quantile,upper_quantile,avg))'''

#print(str(results))

data_points = []

for method in methods:
    #x = []
    #lower_quantiles = []
    #middle_quantiles = []
    #upper_quantiles = []
    #avgs = []
    dots = []
    for m,n,avg,th in results:
        if m == method:
            #lower_quantiles.append(lower_quantile)
            #middle_quantiles.append(middle_quantile)
            #upper_quantiles.append(upper_quantile)
            dots.append((n,avg,th))
    dots.sort(key=lambda x: x[0],reverse=False) 
    data_points.append([method,[a for a,b,c in dots],[b for a,b,c in dots],
                                [c for a,b,c in dots]])
data_points.sort(key=lambda vals: vals[0])
#printprint(str(data_points))


f, ax = plt.subplots(figsize=(12, 7))
#plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 15})
plt.title('Média de similaridade semântica nas melhores correlações')

def get_method_color(method):
    colors = {"PRS":'#ff0000',"MIC":'#dd0074',"DC":'#006600',
            "SPR":'#000099',"SOB":'#000000',"FSH":'#00ff00'}
    for method_name, color in colors.items():
        if method_name in method:
            return color

#colors = ['r','b','g','k','m']
#color_i = 0
for method, ns, avgs, ths in data_points:
    print("Plotting " + method)
    method_color = get_method_color(method)
    '''ax.plot(thresholds, middle_quantiles, 
        label=method+", median",color=method_color)
    ax.plot(thresholds, lower_quantiles, 
        label=method+", 0.25 quantile",marker='v',color=method_color)
    ax.plot(thresholds, upper_quantiles, 
        label=method+", 0.75 quantile",marker='^',color=method_color)'''
    ax.plot(ns, avgs, 
        label=method+", average",color=method_color)
    #color_i += 1
    #break
ax.set_xlim(max_coefs, 0)
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))

ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Similaridade semântica média')
plt.xlabel('N melhores correlações')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output,bbox_inches='tight')

print("Searching confidence intervals")
intervals = [0.7,0.75,0.8,0.85,0.90,0.95,1.0]
interval_start = [ [] for i in range(len(intervals)) ]
interval_size = [ [] for i in range(len(intervals)) ]
for i in tqdm(range(len(intervals))):
    min_sim = intervals[i]
    for method, ns, avgs, ths in data_points:
        min_th = None
        n = None
        for j in range(len(avgs))[::-1]:
            if avgs[j] < min_sim:
                if j < len(avgs)-1:
                    min_th = ths[j+1]
                    n = ns[j+1]
        interval_start[i].append((method,min_th))
        interval_size[i].append((method,n))

'''table = ("Average Semantic Similarity\tInterval Name\t" 
    + "\t".join([m for m,min_th in interval_start[0]])+"\n")
for i in range(len(intervals)):
    table += str(intervals[i])+"\t"
    table += str(i)+"\t"
    table += "\t".join([str(min_coef) for metric, min_coef in interval_start[i]])
    table += "\n"'''
sep = ","
table = "Average Semantic Similarity"+sep+sep.join([str(x) for x in intervals])+"\n"
table += "Interval Name"+sep+sep.join([str(x) for x in range(len(intervals))])+"\n"
method_ths = {}
for method, ns, avgs, ths in data_points:
    method_ths[method] = []
for inter in interval_start:
    for m, min_th in inter:
        method_ths[m].append(min_th)
for method, min_ths in method_ths.items():
    table += method + sep + sep.join([str(min_th) for min_th in min_ths]) + "\n"
print(table)

table = ("Average Semantic Similarity\tInterval Name\t" 
    + "\t".join([m for m,min_th in interval_start[0]])+"\n")
for i in range(len(intervals)):
    table += str(intervals[i])+"\t"
    table += str(i)+"\t"
    table += "\t".join([str(n) for metric, n in interval_size[i]])
    table += "\n"
print(table)