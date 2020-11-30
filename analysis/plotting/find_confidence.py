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
#methods = [x for x in sys.argv[5].split(",")]
output = sys.argv[5]

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
            if sim >= 0.0:
                if not gene1 in sims:
                    sims[gene1] = {}
                if not gene2 in sims:
                    sims[gene2] = {}
                sims[gene1][gene2] = sim
                sims[gene2][gene1] = sim
        elif len(cells) > 0:
            print("Line with too many values:\n\t"+line)

print("Reading correlations")
sim_correlations = {}
corrs_with_sim = 0
total_corrs = 0

corr_files = [correlations_path + "/" + f for f in os.listdir(correlations_path)]
methods = []
for corr_file in corr_files:
    metric_name = corr_file.split("/")[-1].split(".")[0]
    methods.append(metric_name)
    sim_correlations[metric_name] = []
    with open(corr_file, 'r') as stream:
        #line = stream.readline()
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            if cells[0] in sims:
                gene_sims = sims[cells[0]]
                if cells[1] in gene_sims:
                    #if float(cells[2]) >= 0.6:
                    sim_correlations[metric_name].append((float(cells[2]),gene_sims[cells[1]]))
                    corrs_with_sim += 1
            total_corrs += 1
            #else:
            #    sim_correlations.append((cells[-1],float(cells[2]),0.0))

for metric in methods:
    l = len(sim_correlations[metric])
    if l == 0:
        del sim_correlations[metric]
    print("For " + metric + ": " + str(l) + " correlations with similarities.")

print(str(float((corrs_with_sim)/float(total_corrs))*100.0) 
        + "% of the correlations had a corresponding similarity")
for metric, vec in sim_correlations.items():
    print("Sorting " + str(len(vec)) + " items")
    vec.sort(key=lambda x: x[0])

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
    #print(str(total) + "/" + str(max_coefs))
#print(str(results))

data_points = []

for method in methods:
    dots = []
    for m,n,avg,th in results:
        if m == method:
            dots.append((n,avg,th))
    dots.sort(key=lambda x: x[0],reverse=True) 
    data_points.append([method,[a for a,b,c in dots],[b for a,b,c in dots],
                                [c for a,b,c in dots]])
data_points.sort(key=lambda vals: vals[0])
#printprint(str(data_points))

print("Searching confidence intervals")
intervals = [0.5,0.6,0.7,0.8,0.9,1.0]
interval_start = [ [] for i in range(len(intervals)) ]
interval_size = [ [] for i in range(len(intervals)) ]
threshold_points = {}
for i in tqdm(range(len(intervals))):
    min_sim = intervals[i]
    for method, ns, avgs, ths in data_points:
        found = False
        if not method in threshold_points:
            threshold_points[method] = [[],[]]
        for j in range(len(avgs)):
            if avgs[j] >= min_sim:
                interval_start[i].append((method,ths[j]))
                interval_size[i].append((method,ns[j]))
                found = True
                threshold_points[method][0].append(ns[j])
                threshold_points[method][1].append(avgs[j])
                break
        if not found:
            interval_start[i].append((method,None))
            interval_size[i].append((method,None))

f, ax = plt.subplots(figsize=(8, 5))
#plt.style.use('seaborn-whitegrid')
#plt.rcParams.update({'font.size': 15})
plt.title('Média de similaridade semântica nas melhores correlações')

def get_method_color(method):
    colors = {"PRS":'#ff0000',"MIC":'#dd0074',"DC":'#006600',
            "SPR":'#000099',"SOB":'#000000',"FSH":'#00ff00'}
    for method_name, color in colors.items():
        if method_name in method:
            return color

for method, ns, avgs, ths in data_points:
    ax.plot(ns, avgs, 
        label=method, color=get_method_color(method), linewidth=2.0)

ax.axhline(y=intervals[0], color="#555555", linestyle='--', linewidth=1.0, 
        label="Nível de Confiança")
for interval in intervals[1:]:
    ax.axhline(y=interval, color="#555555", linestyle='--', linewidth=1.0)

ax.set_xlim(max_coefs, 0)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))

ax.set_axisbelow(True)
plt.grid(axis='x',which='both')

for method, xy in threshold_points.items():
    m_color = get_method_color(method)
    if m_color == get_method_color("SOB"):
        ax.scatter(xy[0], xy[1], s=100, c=get_method_color(method), marker='>',
            edgecolors= ["#000000"]*len(xy[0]),label="Ponto de Confiança Equivalente")
    else:
        ax.scatter(xy[0], xy[1], s=100, c=get_method_color(method), marker='>', 
            edgecolors= ["#000000"]*len(xy[0]))

plt.ylabel('Similaridade semântica média')
plt.xlabel('N melhores correlações')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output,bbox_inches='tight')

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