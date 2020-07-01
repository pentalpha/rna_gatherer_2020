import os
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statistics as stc

dirs = sys.argv[1:]

def getFilesWith(directory, name_part):
    files = []
    # r=root, d=directories, f = files
    for f in os.listdir(directory):
        if name_part in f:
                files.append(os.path.join(directory, f))
    return files

def plot_scores(data_points):
    print("Sorting")
    data_points.sort(key=lambda x: stc.mean(x[1]))
    print("Extracting data")
    values = [vals for name,vals in data_points]
    names = [name for name,vals in data_points]
    sns.set(style="whitegrid")
    sns.plotting_context("talk", font_scale=1.5)
    sns.set_context("talk")
    print("Making plot")
    f, ax = plt.subplots(figsize=(19, 40))
    #ax.set_xscale("log")
    ax = sns.boxplot(data=values, orient='h', palette="vlag")
    #ax = seaborn.swarmplot(data=values, color=".2")
    ax.set(yticklabels=names)
    axes = ax.axes
    axes.set_xlim(0.0,1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
    #ax = (ax.set_axis_labels("Z-Scores")).set(xlim=(-1,1))
    plt.title("Similarities")
    #plt.show(ax)
    print("Saving")
    plt.savefig("similarities.png",bbox_inches='tight')
    #plt.savefig("similarities.png")


def read_data(path):
    data = []
    with open(path,'r') as stream:
        for line in stream.readlines():
            line = line.rstrip("\n")
            if not ("NA" in line):
                data.append(float(line))
            #else:
            #    data.append(0.0)
    return data

def file_name(path):
    return ".".join(os.path.basename(path).split(".")[:-1])

print("Reading files")
data = []
for d in dirs:
    for f in getFilesWith(d, ".tsv"):
        d = read_data(f)
        name = f.split("/")[-1].replace("-Wang-BMA", "-") + " ("+str(len(d))+")"
        data.append([name,d])

plot_scores(data)