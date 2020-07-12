import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import pandas as pd
import math
import matplotlib.colors as colors

translate_dict_en = {"cc": "Cellular Component", "CC": "Cellular Component",
                "mf": "Molecular Function", "MF": "Molecular Function",
                "bp": "Biological Process", "BP": "Biological Process",
                "fsh": "Fisher", "FSH": "Fisher",
                "sob": "Sobolev", "SOB": "Sobolev",
                "mic": "Maximal Information\nCoefficient", 
                "MIC": "Maximal Information\nCoefficient",
                "prs": "Pearson", "PRS": "Pearson",
                "spr": "Spearman", "SPR": "Spearman",
                "dc": "Distance\nCorrelation", "DC": "Distance\nCorrelation"}

translate_dict = {"cc": "Componente Celular", "CC": "Componente Celular",
                "mf": "Função Molecular", "MF": "Função Molecular",
                "bp": "Processo Biologico", "BP": "Processo Biologico",
                "max": "Similaridade Semântica Máxima", 
                "avg": "Similaridade Semântiac Média",
                "fsh": "Fisher", "FSH": "Fisher",
                "sob": "Sobolev", "SOB": "Sobolev",
                "mic": "Maximal Information\nCoefficient", 
                "MIC": "Maximal Information\nCoefficient",
                "prs": "Pearson", "PRS": "Pearson",
                "spr": "Spearman", "SPR": "Spearman",
                "dc": "Distance\nCorrelation", "DC": "Distance\nCorrelation"}

full_pallete = ['#ff0000', '#dd0074', '#006600',
                '#000099', '#000000', '#00ff00']

def get_filename(full_path):
    last = full_path.split("/")[-1]
    #name = ".".join(last.split(".")[:-1])
    return last

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
            yield lst[i:i + n]

def make_pallete(names):
    pal_dict = {}
    for i in range(len(names)):
        pal_dict[names[i]] = full_pallete[i]
    return pal_dict

def translator(name):
    name.replace("_"," ")
    words = name.split()
    new_words = [(translate_dict[word] if word in translate_dict else word)
                    for word in words]
    new_name = " ".join(new_words)
    return new_name

def get_min_max(df_path):
    current_line = 1
    invalid_values = 0
    with open(df_path, 'r') as stream:
        header_line = stream.readline().rstrip("\n").split("\t")
        current_line += 1
        mins = []
        maxes = []
        for i in range(len(header_line)):
            mins.append(100.)
            maxes.append(0.0)
        raw_line = stream.readline()
        while raw_line:
            cells = raw_line.rstrip("\n").split("\t")
            current_line += 1
            for i in range(len(cells)):
                if cells[i] != "nan" and cells[i] != "NaN":
                    try:
                        value = float(cells[i])
                    except Exception:
                        invalid_values += 1
                    if value > maxes[i]:
                        maxes[i] = value
                    if value < mins[i]:
                        mins[i] = value
            raw_line = stream.readline()
    print(str(invalid_values/current_line) + " invalid lines.")
    return mins, maxes, header_line

def get_lines(stream, limit = 100000, min_len=1):
    lines = []
    raw_line = stream.readline()
    while raw_line:
        cells = raw_line.rstrip("\n").split("\t")
        if len(cells) >= min_len:
            lines.append(cells)
            if len(lines) >= limit:
                return lines, True
        raw_line = stream.readline()
    return lines, False

def get_frequencies(file_ss, index, converter, my_bins):
    frequencies = np.array([])
    line_limit = 10000
    failed_lines = 0
    with open(file_ss, 'r') as stream:
        header_line = stream.readline().rstrip("\n").split("\t")
        line_chunk, more_lines = get_lines(stream, min_len=index+1)
        while more_lines:
            values = []
            for value_str in [cells[index] for cells in line_chunk]:
                try:
                    value = converter(value_str)
                    values.append(value)
                except ValueError:
                    failed_lines += 1
            new_frequencies, bins = np.histogram(values, bins=my_bins)
            if len(frequencies) == 0:
                frequencies = new_frequencies
            else:
                frequencies = frequencies + new_frequencies
            line_chunk, more_lines = get_lines(stream, min_len=index+1)
    return frequencies

def legend_without_duplicate_labels(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    fig.legend(*zip(*unique), loc='center left', bbox_to_anchor=(1, 0.5))

if __name__ == "__main__":
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    ss_col_names = ["mf", "bp", "cc", "avg", "max"]
    mins, maxes, header = get_min_max(file_ss)
    min_max_path = output_dir+"/"+name+"min_max.tsv"
    with open(min_max_path, 'w') as stream:
        stream.write("\t".join(["Column Name","Min","Max"])+"\n")
        for i in range(len(mins)):
            stream.write("\t".join([str(x) for x in [header[i], mins[i],maxes[i]]]) + "\n")
    
    print(open(min_max_path, 'r').read())
    mins = [math.floor(x) for x in mins]
    maxes = [math.ceil(x) for x in maxes]
    bins = [np.linspace(mins[i], maxes[i], 100) for i in range(len(mins))]

    print("Plot semantic similarities")
    fig, ax = plt.subplots(1, 
                        figsize=(5, 7))
    print("\tCalculating frequencies")
    freqs = [get_frequencies(file_ss, i, float, bins[i])
            for i in tqdm([1,2,3,4,5])]
    x_values = [((bin_edges[:-1] + bin_edges[1:]) / 2) 
                for bin_edges in [bins[1],bins[2],bins[3],bins[4],bins[5]]]
    print("Plotting lines")
    pallete = make_pallete(ss_col_names)
    for i in tqdm(range(len(ss_col_names))):
        col = ss_col_names[i]
        freq = freqs[i]
        x = x_values[i]
        ax.plot(x, freq, label=translator(col), color=pallete[col])
    ax.set_yscale('symlog')
    ax.set_title("Distribuição dos Valores de Similaridade Semântica")
    legend_without_duplicate_labels(fig, ax)
    fig.tight_layout()
    fig.savefig(output_dir+"/semantic_similarity_hist.png", bbox_inches='tight')

    bins = [np.linspace(mins[i], maxes[i], 300) for i in range(len(mins))]

    print("Plotting correlations")
    fig, ax = plt.subplots(figsize=(6, 5))
    print("\tCalculating frequencies")
    cols = header[6:]
    print("\t\tGetting frequencies for columns:")
    print("\t\t\t"+str(cols))
    col_indexes = list(range(4,6+len(cols)))
    print("\t\t\t"+str(col_indexes))
    freqs = [get_frequencies(file_ss, i, float, bins[i])
            for i in tqdm(col_indexes)]
    x_values = [((bin_edges[:-1] + bin_edges[1:]) / 2) 
                for bin_edges in [bins[i] for i in col_indexes]]
    
    print("\tPlotting")
    pallete = make_pallete(cols)
    for i in tqdm(range(len(cols))):
        col = cols[i]
        freq = freqs[i]
        x = x_values[i]
        ax.plot(x, freq, label=translator(col), color=pallete[col])
    legend_without_duplicate_labels(fig, ax)
    ax.set_title("Distribuição dos Valores de Coeficientes de Correlação")
    fig.tight_layout()
    fig.savefig(output_dir+"/correlations_hist.png", bbox_inches='tight')

    


