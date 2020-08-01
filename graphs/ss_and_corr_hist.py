import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import pandas as pd
import math
import matplotlib.colors as colors
import dcor
from scipy.stats.stats import pearsonr, spearmanr

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
    #print("Index " + str(index) + " of " + file_ss)
    frequencies = np.array([])
    failed_lines = 0
    with open(file_ss, 'r') as stream:
        header_line = stream.readline().rstrip("\n").split("\t")
        more_lines = True
        while more_lines:
            line_chunk, more_lines = get_lines(stream, min_len=index+1)
            #print("Processing chunk of " + str(len(line_chunk)) + " lines.")
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
    return frequencies

def legend_without_duplicate_labels(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    fig.legend(*zip(*unique), loc='center left', bbox_to_anchor=(1, 0.5))

def convert_to_ranges(ranges, y_points):
    new_x = []
    new_y = []
    for i in range(len(ranges)):
        a, b = ranges[i]
        y = y_points[i]
        new_x.append(a)
        new_y.append(y)
        new_x.append(b)
        new_y.append(y)
    return new_x, new_y

def get_percentiles_by_bins(ss_column, corr_column, invalid_lines, n_bins=500):
    bin_vals = [[] for i in range(n_bins)]
    #bin_ids  = [[] for i in range(n_bins)]
    for i in range(len(ss_column)):
        if not i in invalid_lines:
            ss_val   = ss_column[i]
            corr_val = corr_column[i]
            corr_bin = int(round(corr_val*n_bins,0))
            if corr_bin >= n_bins:
                corr_bin = n_bins - 1
            #bin_ids[corr_bin].append(i)
            bin_vals[corr_bin].append(ss_val)

    x               = []
    medians         = []
    top_quantile    = []
    bottom_quantile = []
    for i in range(n_bins):
        if len(bin_vals[i]) > 0:
            x1 = i/n_bins
            x2 = (i+1)/n_bins
            x.append((x1,x2))

            bin_vals[i].sort()
            top_quantile.append(np.percentile(bin_vals[i], 75))
            medians.append(np.percentile(bin_vals[i], 50))
            bottom_quantile.append(np.percentile(bin_vals[i], 25))
        '''else:
            top_quantile.append(0.0)
            medians.append(0.0)
            bottom_quantile.append(0.0)'''
    corr = dcor.distance_correlation(np.array([a for a,b in x]), np.array(medians))

    return x, top_quantile, medians, bottom_quantile, corr

def load_df_values(file_ss):
    min_cells = 1+5+6
    columns = []
    for n_cell in range(min_cells):
        columns.append([])
    failed_lines = 0
    aprox_lines_to_read = int((25*1024*1024*1024)/195)
    stream = open(file_ss, 'r')
    progress_bar = tqdm(total=aprox_lines_to_read)
    header_line = stream.readline().rstrip("\n").split("\t")
    raw_line = stream.readline()
    progress_bar.update(2)
    while raw_line:
        try:
            np_values = np.fromstring(raw_line.rstrip("\n"),
                                    sep='\t',dtype=float)
            if len(np_values) >= min_cells:
                for i in range(len(np_values)):
                    columns[i].append(np_values[i])
            else:
                failed_lines += 1
        except ValueError:
            failed_lines += 1
        raw_line = stream.readline()
        progress_bar.update(1)
    print("Failed lines: " + str(failed_lines))
    stream.close()
    progress_bar.close()
    return columns[1:], header_line[1:], [int(x) for x in columns[0]]

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (not np.isfinite(x[i]) or np.isnan(x[i])) or x[i] < min_cor:
            f.add(i)
    return f

if __name__ == "__main__":
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]
    
    if not os.path.exists(output_dir):
        print("creating output dir")
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    ss_col_names = ["mf", "bp", "cc", "avg", "max"]
    print("Reading minimum and maximum values for columns")
    unorm_mins, unorm_maxes, header = get_min_max(file_ss)
    
    min_max_path = output_dir+"/"+name+"min_max.tsv"
    with open(min_max_path, 'w') as stream:
        stream.write("\t".join(["Column Name","Min","Max"])+"\n")
        for i in range(len(unorm_mins)):
            stream.write("\t".join([str(x) for x in [header[i], unorm_mins[i], unorm_maxes[i]]]) + "\n")
    
    print(open(min_max_path, 'r').read())
    mins  = [math.floor(x) for x in unorm_mins]
    maxes = [math.ceil(x) for x in unorm_maxes]
    bins = [np.linspace(mins[i], maxes[i], 100) for i in range(len(mins))]

    '''print("Plot semantic similarities")
    fig, ax = plt.subplots(1, 
                        figsize=(5, 7))
    print("\tCalculating frequencies")
    freqs = [get_frequencies(file_ss, i, float, bins[i])
            for i in tqdm([1,2,3])]
    x_values = [((bin_edges[:-1] + bin_edges[1:]) / 2) 
                for bin_edges in [bins[1],bins[2],bins[3]]]
    #print(str([len(freq) for freq in freqs]))
    #print(str(freqs[0]))
    #print(str([len(x) for x in x_values]))
    print("Plotting lines")
    pallete = make_pallete(ss_col_names)
    for i in tqdm(range(len(ss_col_names[:3]))):
        col = ss_col_names[i]
        freq = freqs[i]
        x = x_values[i]
        ax.plot(x, freq, label=translator(col), color=pallete[col])
    ax.set_yscale('symlog')
    ax.set_title("Distribuição dos Valores de Similaridade Semântica")
    legend_without_duplicate_labels(fig, ax)
    fig.tight_layout()
    fig.savefig(output_dir+"/semantic_similarity_hist.png", bbox_inches='tight')

    range_len = [min(int(abs(unorm_mins[i]-unorm_maxes[i])*25),2000)
            for i in range(len(unorm_mins))]
    bins = [np.linspace(unorm_mins[i], unorm_maxes[i], range_len[i]) 
            for i in range(len(unorm_mins))]
    print("Plotting correlations")
    fig, ax = plt.subplots(figsize=(6, 5))
    print("\tCalculating frequencies")
    cols = header[6:]
    print("\t\tGetting frequencies for columns:")
    print("\t\t\t"+str(cols))
    col_indexes = [6,7,8,9,10,11]
    print("\t\t\t"+str(col_indexes))
    print("\t\t\t"+str(range_len))
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
    fig.savefig(output_dir+"/correlations_hist.png", bbox_inches='tight')'''

    mins  = [round(x,2) for x in unorm_mins]
    maxes = [round(x,2) for x in unorm_maxes]

    mins = mins[1:]
    maxes = maxes[1:]
    column_values, header, ids = load_df_values(file_ss)
    corr_col_names = header[5:]
    #Normalizing negatives and high values to 0.0 -> 1.0
    norm_negatives = [mins[i] < 0 for i in range(len(mins))]
    norm_by_max = [maxes[i] > 1.0 for i in range(len(maxes))]
    for i in range(len(column_values)):
        if norm_negatives[i]:
            print("Normalizing negative values of " + header[i])
            column_values[i] = [abs(x) for x in column_values[i]]
        if norm_by_max[i]:
            print("Normalizing values to range 0.0-1.0, for " + header[i])
            column_values[i] = [x/maxes[i] for x in column_values[i]]
        if header[i] in ["SOB","FSH"]:
            column_values[i] = [max(1.0-x, 0) for x in column_values[i]]
    print("Creating filters...")
    all_corrs_filter = set.union(*[filter_column(x, 0.0) 
                                          for x in tqdm(column_values[5:])])
    ss_filters = [filter_column(x, 0.0) for x in tqdm(column_values[:5])]
    ss_filters = [ss_filter.union(all_corrs_filter) for ss_filter in ss_filters]

    print("Plot average SS values")
    x_plots = len(corr_col_names)
    y_plots = len(ss_col_names[:-2])
    fig, axes = plt.subplots(y_plots, x_plots,
                            figsize=(x_plots*8, y_plots*3))
    print(str([len(axes),len(axes[0])]))
    bar = tqdm(total=x_plots*y_plots)
    axes_j = 0
    higher_y = 0.0
    for j in range(len(ss_col_names[:-2])):
        corr_plots = axes[axes_j]
        ss_filter = ss_filters[j]
        i = 0
        #bizarres = set()
        for corr_i in [5,6,7,8,9,10]:
            sub_plot = corr_plots[i]
            x, top_quantile, medians, bottom_quantile, corr = get_percentiles_by_bins(
                                        column_values[j],
                                        column_values[corr_i],
                                        ss_filter)
            
            x_median, y_median = convert_to_ranges(x, medians)
            x_top, y_top = convert_to_ranges(x, top_quantile)
            x_bottom, y_bottom = convert_to_ranges(x, bottom_quantile)

            sub_plot.fill_between(x_top, y_top, y_bottom, color='red')
            #sub_plot.xaxis.set_label_position('top')
            if j == 2:
                sub_plot.set_xlabel(translator(header[corr_i]))
            if i == 0:
                sub_plot.text(-0.2, 0.5, translator(header[j]),
                    rotation=90,
                    horizontalalignment='center',
                    verticalalignment='center',
                    multialignment='center', 
                    transform=sub_plot.transAxes)
            sub_plot.set_ylabel("Percentis de S.S. por segmento")
            sub_plot.set_title("Correlação das Medianas = " + str(round(corr,4)))

            sub_plot.plot(x_median, y_median, color="blue")
            sub_plot.set_ylim(0.0, 0.6)

            '''avg_ax = sub_plot.twinx()
            avg_ax.plot(x, y, color="blue")
            avg_ax.set_ylabel("Média por segmento", color="blue")
            avg_ax.tick_params(axis='y', colors='blue')'''
            bar.update(1)
            i += 1
        '''bizarre_rows = [[ids[index]] + [column_values[j][index]] + [col[index] for col in column_values[5:]] 
                        for index in bizarre_indexes]
        bizarre_lines= ["\t".join([str(val) for val in row]) for row in bizarre_rows]
        with open(output_dir+"/bizarre_ss-"+ss_col_names[j]+".tsv", 'w') as stream:
            stream.write("\t".join(["pair_id"]+[ss_col_names[j]]+header[5:])+"\n")
            for line in bizarre_lines:
                stream.write(line+"\n")'''
        axes_j += 1
    bar.close()
    fig.tight_layout()
    fig.savefig(output_dir+"/correlations_ranges.png", bbox_inches='tight')


