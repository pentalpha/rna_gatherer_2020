import sys
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import statistics as stc

correlations_path = sys.argv[1]
min_threshold,step,max_threshold = [float(x) for x in sys.argv[2].split(",")]
thresholds = np.arange(min_threshold,max_threshold,step).round(3)
print("Thresholds: " + str(thresholds))
methods = [x for x in sys.argv[3].split(",")]
output = sys.argv[4]

print("Reading correlations")
correlations = []
with open(correlations_path, 'r') as stream:
    #line = stream.readline()
    for line in stream:
        cells = line.rstrip("\n").split("\t")
        value = float(cells[2])
        value = abs(value)
        if float(cells[2]) >= thresholds[0]:
            correlations.append((cells[-1],value))

print("Sorting " + str(len(correlations)) + " items")
correlations.sort(key=lambda x: x[1])
current_threshold = 0
locations = []
for i in tqdm(range(len(correlations))):
    method, corr = correlations[i]
    if corr >= thresholds[current_threshold]:
        locations.append(i)
        current_threshold += 1
        if current_threshold >= len(thresholds):
            break

#print(str(locations))

corr_results = []
for i in tqdm(range(len(locations))):
    loc = locations[i]
    threshold = thresholds[i]
    #print("Threshold " + str(threshold))
    sub_vec = correlations[loc:]
    #if i < len(locations)-1:
    #    sub_vec = correlations[loc:locations[i+1]]
    #print("from " + str(locations[i]) + " to " + str(locations[i]) + ": " + str(len(sub_vec)))
    method_cors = {"PRS":[],"SPR":[],"MIC":[],"DC":[],"SOB":[],"FSH":[]}
    for m,corr in sub_vec:
        method_cors[m].append(corr)
    #for method,vec in method_cors.items():
    #    print(method + " has " + str(len(vec)) + " elements")
    for method in methods:
        method_corrs = len(method_cors[method])
        corr_results.append((method,threshold,method_corrs))

#print(str(results))

corr_points = []
for method in methods:
    x = []
    y = []
    for m,threshold,method_corrs in corr_results:
        if m == method:
            x.append(threshold)
            y.append(method_corrs)
            print(m + " has " + str(method_corrs) + " correlations at " + str(threshold))
    corr_points.append([method,x,y])

#printprint(str(data_points))



#plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 15})


def get_method_color(method):
    colors = {"PRS":'#ff0000',"MIC":'#ff0074',"DC":'#008800',
            "SPR":'#0000ff',"SOB":'#000000',"FSH":'#00ff00'}
    for method_name, color in colors.items():
        if method_name in method:
            return color

f, ax = plt.subplots(figsize=(15, 6))
for method, x, y in corr_points:
    print("Plotting " + method)
    method_color = get_method_color(method)
    ax.plot(x, y, label=method,color=method_color)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax.set_axisbelow(True)
plt.grid(axis='x',which='both')
plt.ylabel('Number of correlations above a certain threshold')
plt.yscale('symlog')
plt.xlabel('Correlation coefficient threshold')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(output,bbox_inches='tight')