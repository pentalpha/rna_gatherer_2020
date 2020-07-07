import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import dask
import dask.dataframe as dd
import pandas as pd
from multiprocessing import Pool, TimeoutError
import multiprocessing
from dask.diagnostics import ProgressBar

processes = max(2, multiprocessing.cpu_count()-1)

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

def legend_without_duplicate_labels(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    fig.legend(*zip(*unique), loc='center left', bbox_to_anchor=(1, 0.5))

if __name__ == "__main__":
    files_ss = sys.argv[1].split(",")
    files_corr = sys.argv[2].split(",")
    dict_file = sys.argv[3]
    annotation = sys.argv[4]
    output_dir = sys.argv[5]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    ss_col_names = ["mf", "bp", "cc"]
    
    names_dict = {}
    print("Reading names dict")
    with open(dict_file,'r') as in_stream:
        for line in in_stream:
            cells = line.rstrip("\n").split("\t")
            names_dict[cells[0]] = cells[1]
    
    n_terms = {}
    print("Reading gene annotation")
    with open(annotation,'r') as in_stream:
        terms_by_gene = {}
        for line in in_stream:
            if not(line.startswith("!") or line.startswith("#")):
                cells = line.rstrip("\n").split("\t")
                if len(cells) >= 5:
                    gene_name = cells[1]
                    go_term = cells[4]
                    if not gene_name in terms_by_gene:
                        terms_by_gene[gene_name] = set()
                    terms_by_gene[gene_name].add(go_term)
                else:
                    print("Invalid line:" + str(cells))
        for gene_name, term_set in terms_by_gene.items():
            if gene_name in names_dict:
                n_terms[names_dict[gene_name]] = len(term_set)
        for gene_name, _id in names_dict.items():
            if not _id in n_terms:
                n_terms[_id] = 0
        print(str(list(n_terms.items())[0:30]))

    '''print("Plot semantic similarities NaN frequency")
    fig, ax = plt.subplots(len(files_ss), 
                        figsize=(5, len(files_ss)*4))
    if len(files_ss) == 1:
        ax = [ax]
    for i in range(len(files_ss)):
        sub_plt = ax[i]
        ss_file = files_ss[i]
        df = dd.read_parquet(ss_file)
        name = get_filename(ss_file).split(".")[0]
        print("\tPlotting " + name)
        plot_df_ss_nan(sub_plt, df, name, n_terms)
        sub_plt.set_xlabel("Diferença entre número de termos nas anotações")
        sub_plt.set_ylabel("Porcentagem (%) de similaridades semânticas nulas.")
    legend_without_duplicate_labels(fig, ax[0])
    fig.tight_layout()
    fig.savefig(output_dir+"/semantic_similarity_nan_hist.png", bbox_inches='tight')'''

    print("Plot semantic similarities")
    fig, ax = plt.subplots(len(files_ss), 
                        figsize=(5, len(files_ss)*4))
    if len(files_ss) == 1:
        ax = [ax]
    for i in range(len(files_ss)):
        sub_plt = ax[i]
        ss_file = files_ss[i]
        df = dd.read_parquet(ss_file, columns=ss_col_names)
        name = get_filename(ss_file).split(".")[0]
        print("\tPlotting " + name)
        plot_df_ss(sub_plt, df, name, ss_col_names)
        sub_plt.set_yscale('symlog')
    legend_without_duplicate_labels(fig, ax[0])
    fig.tight_layout()
    fig.savefig(output_dir+"/semantic_similarity_hist.png", bbox_inches='tight')

    print("Plotting correlations")
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_df_corr(ax, files_corr)
    ax.set_yscale('symlog')
    legend_without_duplicate_labels(fig, ax)
    ax.set_title("Distribution of Correlation Coefficient Values")
    fig.tight_layout()
    fig.savefig(output_dir+"/correlations_hist.png", bbox_inches='tight')

