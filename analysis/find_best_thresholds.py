import sys
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
#import matplotlib.ticker as ticker
import os
import gzip
import io
#import seaborn as sns
#import statistics as stc

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

def get_filename(full_path):
    last = full_path.split("/")[-1]
    #name = ".".join(last.split(".")[:-1])
    return last

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

def load_df_values(file_ss, ss_col_names):
    min_cells = len(ss_col_names)+6
    corr_cols = [[] for a in range(6)]
    ss_cols = [[] for a in range(len(ss_col_names))]
    
    failed_lines = 0
    aprox_lines_to_read = int((25*1024*1024*1024)/195)
    stream = None
    if file_ss.endswith(".tsv"):
        stream = open(file_ss, 'r')
    elif file_ss.endswith(".gzip") or file_ss.endswith(".gz"):
        print("Opening gzip file")
        stream = gzip.open(file_ss, 'rt')
        #stream = io.BufferedReader(gz)
    progress_bar = tqdm(total=aprox_lines_to_read)
    header_line = stream.readline().rstrip("\n").split("\t")[1:]
    ss_indexes = []
    for ss_name in ss_col_names:
        if ss_name in header_line:
            ss_indexes.append(header_line.index(ss_name))
    print("ss column names:"+str(ss_col_names))
    print("ss col indexes:"+str(ss_indexes))
    raw_line = stream.readline()
    progress_bar.update(2)
    while raw_line:
        try:
            np_values = np.fromstring(raw_line.rstrip("\n"),
                                    sep='\t',dtype=float)[1:]
            if min_cells <= len(np_values):
                for i in range(len(ss_indexes)):
                    value = np.float(np_values[ss_indexes[i]])
                    if not np.isfinite(value):
                        value = np.nan
                    ss_cols[i].append(value)
                corr_col_i = 0
                for value in np_values[-6:]:
                    if not np.isfinite(value):
                        corr_cols[corr_col_i].append(np.nan)
                    corr_cols[corr_col_i].append(value)
                    corr_col_i += 1
            else:
                failed_lines += 1
        except ValueError:
            failed_lines += 1
        raw_line = stream.readline()
        progress_bar.update(1)
    print("Failed lines: " + str(failed_lines))
    stream.close()
    progress_bar.close()
    return ss_cols, corr_cols, header_line[-6:]

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (not np.isfinite(x[i]) or np.isnan(x[i])) or x[i] < min_cor:
            f.add(i)
    return f

def lists_to_feature_list(lists, corr_filter):
    features = []
    for i in range(len(lists[0])):
        if not i in corr_filter:
            features.append([col[i] for col in lists])
    return features

def test_thresholds(feature_list, smaller=False, min_val=0.0):
    n_ss = len(feature_list[0])-1
    #ss_totals = [0.0,0.0,0.0]
    ss_totals = [0.0 for x in range(n_ss)]
    #ss_n = [0,0,0]
    ss_n = [0 for x in range(n_ss)]
    #ths = [[], [], []]
    ths = [[] for x in range(n_ss)]
    #ss_avgs = [[], [], []]
    ss_avgs = [[] for x in range(n_ss)]
    
    for i in range(len(feature_list)):
        if i != len(feature_list)-1:
            feature_values = feature_list[i]
            corr_value = feature_values[-1]
            if corr_value >= min_val:
                for ss_index in range(len(ss_totals)):
                    ss_val = feature_values[ss_index]
                    if not np.isnan(ss_val):
                        ss_totals[ss_index] += ss_val
                        ss_n[ss_index] += 1
                        cond = corr_value > feature_list[i+1][-1]
                        if smaller: 
                            cond = corr_value < feature_list[i+1][-1]
                        if cond:
                            new_avg = ss_totals[ss_index] / ss_n[ss_index]
                            ss_avgs[ss_index].append(new_avg)
                            ths[ss_index].append(corr_value)
    #return (ths[0], ss_avgs[0]), (ths[1], ss_avgs[1]), (ths[2], ss_avgs[2])
    return [(ths[i], ss_avgs[i]) for i in range(len(ths))]

'''def test_ths(feature_list, greater_values, min_v):
    if greater_values:
        return test_threasholds(feature_list, min_val=min_v)
    else:
        return test_threasholds_smaller(feature_list, min_val=min_v)'''

def get_confs(ths, conf_levels):
    corr_values, ss_values = ths
    confs = []
    for conf_level in conf_levels:
        threshold = None
        for index in range(len(corr_values))[::-1]:
            if ss_values[index] >= conf_level:
                threshold = corr_values[index]
                break
        confs.append((threshold, conf_level))
    return confs

def get_not_none(pairs, min_val=0.0):
    x = []
    y = []

    for a, b in pairs:
        if a != None and b != None and a >= min_val:
            x.append(a)
            y.append(b)
    
    return x, y

def get_hist(vec, n_bins=10):
    #hist_vec, hist_edges = np.histogram(vec, bins=n_bins)
    plt.hist(vec, bins=n_bins, histtype='stepfilled')
    plt.show()

# %%
if __name__ == "__main__":
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]
    max_ss = float(sys.argv[3])
    ss_col_names = sys.argv[4].split(',')
    #output_dir,file_ss,max_ss = ('test_05','sample05.tsv',1.0)
    '''output_dir,file_ss,max_ss,ss_col_names = ("mgi/test_02_1",
                                 "mgi/sample02.tsv",
                                 1.0,
                                 ['mf','avg'])'''
    if not os.path.exists(output_dir):
        print("creating output dir")
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    #ss_col_names = ["mf", "bp", "cc", "avg", "max"]
    ss_cols, corr_cols, corr_col_names = load_df_values(file_ss, ss_col_names)
    #column_values, header, ids = load_df_values(file_ss, ss_col_names)
    
    def normalize(column, name, is_corr = False):
        if name == "SPR" or name == "PRS":
            print("Normalizing negative values for " + name)
            column = [abs(x) for x in column]
        if is_corr:
            print("Rounding values for " + name)
            column = [round(x, 7) for x in column]
        return column
    
    for i in range(len(ss_cols)):
        col_name = ss_col_names[i]
        ss_cols[i] = normalize(ss_cols[i], col_name)
    
    for i in range(len(corr_cols)):
        col_name = corr_col_names[i]
        corr_cols[i] = normalize(corr_cols[i], col_name, is_corr=True)

# %%
    print("Creating filters...")
    corr_filters = [filter_column(x, 0.0) for x in tqdm(corr_cols)]
    all_corrs_filter = set.union(*corr_filters)
    
# %%
    def use_smaller_values(corr_name):
        if corr_name != "SOB" and corr_name != "FSH":
            return False
        else:
            return True
    i = 0
    ths_by_metric = {}
    for corr_col in tqdm(corr_cols):
        corr_filter = corr_filters[i]
        corr_name = corr_col_names[i]
        print("Testing threasholds for " + corr_name)
        print("\tCreating list with values.")
        feature_list = lists_to_feature_list(ss_cols + [corr_col], corr_filter)
        print("\tSorting list")
        feature_list.sort(key=lambda features: features[-1],
                          reverse=(not use_smaller_values(corr_name)))
        print("\tFinding thresholds")
        ths_by_metric[corr_name] = test_thresholds(feature_list, 
                                                   use_smaller_values(corr_name))
        i += 1

# %%
    print("Finding confidence levels")
    conf_levels = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
    conf_levels += [0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]
    conf_levels += [0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75]
    conf_levels += [0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.99]
    #conf_levels = [round(x, 2) for x in np.arange(0.1, 1.05, 0.05)]
    confs_by_metric = {corr_name: [get_confs(onto_th, conf_levels) for onto_th in onto_ths]
                        for corr_name, onto_ths in tqdm(ths_by_metric.items())}
# %%
    print("Plotting")
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    subplots = []
    name_i = 0
    for i in range(len(axes)):
        for j in range(len(axes[i])):
            sub_plt = axes[i][j]
            subplots.append((corr_col_names[name_i], sub_plt))
            name_i += 1
    zord = 0
    colors = ["darkred", "navy", "forestgreen", "black", "orange"]
    for corr_name, sub_plt in tqdm(subplots):
        ss_ths = ths_by_metric[corr_name]

        sub_plt.yaxis.set_major_locator(MultipleLocator(0.1))
        sub_plt.set_ylim(0,0.65)
        
        sub_plt.set_xlim(-0.05,1.05)
        sub_plt.xaxis.set_major_locator(MultipleLocator(0.2))
        
        if corr_name == "MIC":
            sub_plt.yaxis.set_major_locator(MultipleLocator(0.05))
            sub_plt.set_ylim(0,0.25)
        elif corr_name == "SOB":
            sub_plt.xaxis.set_major_locator(MultipleLocator(2))
            sub_plt.set_xlim(-0.05,15)
        elif corr_name == "DC":
            sub_plt.yaxis.set_major_locator(MultipleLocator(0.1))
            sub_plt.set_ylim(0,1.05)
        elif corr_name== "FSH":
            sub_plt.xaxis.set_major_locator(MultipleLocator(0.25))
            sub_plt.set_xlim(-0.05,1.6)

        sub_plt.grid(True, which="major", axis="y", linestyle="--")

        for i in range(len(ss_ths)):
            x = ss_ths[i][0]
            y = ss_ths[i][1]
            sub_plt.plot(x, y, color=colors[i], 
                         zorder=zord, label=ss_col_names[i].upper())
            zord += 1

        sub_plt.set_title(translator(corr_name))
        sub_plt.set_ylabel("Similaridade Semântica Média")
        sub_plt.set_xlabel("Threshold de Correlação")
        
        if corr_col_names.index(corr_name) == 4:
            #draw legend
            sub_plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                      fancybox=True, shadow=True, ncol=len(ss_ths))
            

    fig.tight_layout()
    fig.savefig(output_dir+"/correlation_thresholds.png", bbox_inches='tight')

    for ss_index in range(len(ss_col_names)):
        ss_name = ss_col_names[ss_index]
        output_table = output_dir+"/confidence_intervals-"+ss_name.upper()+".tsv"
        out_str = "\t".join(["Average Semantic Similarity"] + [str(val) for val in conf_levels])+"\n"
        out_str += "\t".join(["Interval Name"] + [str(i) for i in range(len(conf_levels))])+"\n"
        for corr_name, onto_ths in confs_by_metric.items():
            ths = onto_ths[ss_index]
            out_str += "\t".join([corr_name.upper()] + [str(th) for th, conf_level in ths])+"\n"
        open(output_table,'w').write(out_str)
