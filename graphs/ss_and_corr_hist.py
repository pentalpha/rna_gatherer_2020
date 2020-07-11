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

def get_plotcoordinates(number_of_files, max_cols):
    print("For " + str((number_of_files, max_cols)) + ": ")
    coords = []
    x, y = (0, 0)
    for i in range(number_of_files*max_cols):
        coords.append((x,y))
        if x == (max_cols-1):
            x = 0
            y += 1
        else:
            x += 1
    print(str(coords))
    return coords, coords[-1][1]+1

def add_col(plot_ax, df, col_name, label_name, my_bins, color):
    df_col = df[col_name].dropna()
    counts, bin_edges = np.histogram(df_col, bins=my_bins)
    center = (bin_edges[:-1] + bin_edges[1:]) / 2
    plot_ax.plot(center, counts, label=label_name, color=color)

def plot_df_ss_nan(plot_ax, df, name, n_termos, col_names= ['mf','bp','cc']):
    pallete = make_pallete(col_names)
    #bottom_val = min([df[np.isfinite(df[col])][col].min() for col in col_names])
    #top_val = max([df[np.isfinite(df[col])][col].max() for col in col_names])
    onto_results = {'mf': {"y": [], "x": [], "total": 0},
        'bp': {"y": [], "x": [], "total": 0},
        'cc': {"y": [], "x": [], "total": 0}}
    
    total_rows = len(df)
    print("\tTotal rows:" + str(total_rows))
    print("\tGrouping by the difference in number of terms")
    #df = df[df['mf'].isna() or df['cc'].isna() or df['bp'].isna()]
    def annotation_len(gene, termos):
        gene_int = int(gene)
        if not (gene_int in n_terms):
            termos[gene_int] = 0
        return termos[gene_int]

    df["terms_diff"] = df.apply(
        lambda row: abs(annotation_len(row["gene_a"], n_termos) 
                        - annotation_len(row["gene_b"], n_termos)),
        axis=1, meta=int)

    def stats(df_group, col_names= ['mf','bp','cc']):
        l = len(df_group)
        stat = {}
        for col in col_names:
            stat[col] = {}
            y = len(df[df[col].isna()])
            stat[col]["y_abs"] = y
            stat[col]["y"] = (y/l)*100.0
        return stat

    diff_res = df.groupby("terms_diff").apply(stats, meta=object).compute()
    print(str(diff_res))

    print("\tComputing stats for groups")
    for diff, res in diff_res.iteritems():
        for col, stats in res.items():
            onto_results[col]['x'].append(diff)
            onto_results[col]['y'].append(stats['y'])
            onto_results[col]['total'] += stats['y_abs']
    
    print(str(onto_results))
    '''for diff, entries in tqdm(diff_groups):
        l = len(entries)
        for col in col_names:
            y = len(nan_df = entries[entries[col].isna()])
            onto_results[col]["total"] += y
            onto_results[col]["y"].append((y/l)*100.0)
            onto_results[col]["x"].append(l)'''

    print("\tPlotting")
    for onto_name, results in tqdm(onto_results.items()):
        print("plotting " + str(len(results['x'])))
        plot_ax.plot(results["x"], results["y"], 
                    label=translator(onto_name),
                    color=pallete[onto_name])
        print("\tNaN rate for " + onto_name + ": "
                + str(results["total"]/total_rows))
    plot_ax.set_title("Métrica " + translator(name))
    df.head()

def plot_df_ss(plot_ax, df, name, col_names):
    pallete = make_pallete(col_names)
    #bottom_val = min([df[np.isfinite(df[col])][col].min() for col in col_names])
    #top_val = max([df[np.isfinite(df[col])][col].max() for col in col_names])

    print("\tCalculating bins")
    mins = []
    maxes = []
    for col in col_names:
        print("\t\tFor col " + col)
        finite_col = df[col].dropna()
        local_min, local_max = dd.compute(finite_col.min(),
                                        finite_col.max())
        mins.append(local_min)
        maxes.append(local_max)
    bottom_val = min(mins)
    top_val = max(maxes)
    my_bins = np.linspace(bottom_val, top_val, 200)

    print("\tPlotting lines")
    for col in col_names:
        print("\t\tPlotting " + col)
        add_col(plot_ax, df, col, translator(col), my_bins, pallete[col])
    plot_ax.set_title("Histograma da métrica " + translator(name))

def plot_df_corr(plot_ax, df_paths):
    names = [translator(get_filename(corr_file).split(".")[0]) for corr_file in df_paths]
    pallete = make_pallete(names)

    print("Calculating bins")
    mins = []
    maxes = []
    for df_path in df_paths:
        df = dd.read_parquet(df_path)
        print("\tFor df " + df_path)
        finite_col = df['corr'].dropna()
        local_min, local_max = dd.compute(finite_col.min(),
                                        finite_col.max())
        mins.append(local_min)
        maxes.append(local_max)
    bottom_val = min(mins)
    top_val = max(maxes)
    my_bins = np.linspace(bottom_val, top_val, 200)

    print("Calculating histogram for metrics")
    for i in range(len(df_paths)):
        name = translator(names[i])
        print("\tPlotting " + name)
        df_path = df_paths[i]
        df = dd.read_parquet(df_path)
        
        add_col(plot_ax, df, 'corr', name,
                my_bins, pallete[names[i]])

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

def get_frequencies_2d(file_ss, index_x, index_y, converter, my_bins):
    min_cells = max(index_x+1, index_y+1)
    #frequencies = np.array([])
    line_limit = 10000
    failed_lines = 0
    with open(file_ss, 'r') as stream:
        values_x = []
        values_y = []
        header_line = stream.readline().rstrip("\n").split("\t")
        line_chunk, more_lines = get_lines(stream, min_len=min_cells)
        while more_lines:
            for x_str, y_str in [(cells[index_x], cells[index_y]) 
                                for cells in line_chunk]:
                try:
                    value_x, value_y = (converter(x_str), converter(y_str))
                    if np.isfinite(value_x) and np.isfinite(value_y):
                        values_x.append(value_x)
                        values_y.append(value_y)
                except ValueError:
                    failed_lines += 1
            line_chunk, more_lines = get_lines(stream, min_len=min_cells)
        frequencies, bins_x, bins_y = np.histogram2d(values_x, values_y,
                                bins=[my_bins[index_x], my_bins[index_y]])
        coefficients = np.polyfit(values_y, values_x, 5)
        poly = np.poly1d(coefficients)
        y_subset = np.sort(np.random.choice(values_y, 50000, replace=False))
        new_x = poly(y_subset)
        return frequencies, new_x, y_subset
    return None

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
    ss_col_names = ["mf", "bp", "cc"]
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
            for i in tqdm([1,2,3])]
    x_values = [((bin_edges[:-1] + bin_edges[1:]) / 2) 
                for bin_edges in [bins[1], bins[2], bins[3]]]
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
    cols = header[4:]
    print("\t\tGetting frequencies for columns:")
    print("\t\t\t"+str(cols))
    col_indexes = list(range(4,4+len(cols)))
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

    bins = [np.linspace(mins[i], maxes[i], 50) for i in range(len(mins))]
    progress_bar = tqdm(total=len(cols)*3)
    fig, ax = plt.subplots(len(cols), 3, figsize=(12, 4*len(cols)))
    for i in range(len(cols)):
        metric_plots = ax[i]
        col = cols[i]
        for j in range(len(ss_col_names)):
            freqs_2d, predicted_ss, corr_values = get_frequencies_2d(
                                        file_ss, j+1, col_indexes[i],
                                        float, bins)
            #print(str(freqs_2d))
            edges_x = bins[j+1]
            edges_y = bins[col_indexes[i]]
            metric_plots[j].pcolormesh(edges_x, edges_y, freqs_2d,
                            norm=colors.SymLogNorm(linthresh = 0.03,
                                                    vmin=freqs_2d.min(), 
                                                    vmax=freqs_2d.max()))
            metric_plots[j].plot(predicted_ss, corr_values, 
                                color='white', linewidth=4)
            progress_bar.update(1)
            if i == 0:
                metric_plots[j].set_xlabel("Similaridade Semântica ("
                                +translator(ss_col_names[j])+")")
                metric_plots[j].xaxis.set_label_position('top') 
        metric_plots[0].set_ylabel(translator(cols[i]))
    #fig.suptitle("Heatmaps das Métricas de Correlação Para Cada Ontologia")
    fig.tight_layout()
    fig.savefig(output_dir+"/corr_vs_ss.png", bbox_inches='tight')


