import sys
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
#import matplotlib.ticker as ticker
import os
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

def load_df_values(file_ss):
    stream = open(file_ss, 'r')
    header_line = stream.readline().rstrip("\n").split("\t")
    ss_names = header_line[:-6][1:]
    n_ss_measures = len(ss_names)
    single_ss = n_ss_measures == 1
    min_cells = 1+n_ss_measures+6
    columns = []
    for n_cell in range(min_cells):
        columns.append([])
    failed_lines = 0
    aprox_lines_to_read = int((25*1024*1024*1024)/195)
    
    progress_bar = tqdm(total=aprox_lines_to_read)
    
    raw_line = stream.readline()
    progress_bar.update(2)
    while raw_line:
        try:
            np_values = np.fromstring(raw_line.rstrip("\n"),
                                    sep='\t',dtype=float)
            if len(np_values) >= min_cells:
                for i in range(len(np_values)):
                    if not np.isfinite(np_values[i]):
                        np_values[i] = np.nan
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
    return columns[1:], header_line[1:], [int(x) for x in columns[0]], ss_names, single_ss

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (not np.isfinite(x[i]) or np.isnan(x[i])) and x[i] >= min_cor:
            f.add(i)
    return f

def lists_to_feature_list(lists, corr_filter):
    features = []
    for i in range(len(lists[0])):
        if not i in corr_filter:
            features.append([col[i] for col in lists])
    return features

def test_threasholds(feature_list, min_val=-9999, n_ss_cols=1):
    ss_totals = [0.0 for a in range(n_ss_cols)]
    ss_n      =   [0 for a in range(n_ss_cols)]
    ths       =  [[] for a in range(n_ss_cols)]
    ss_avgs   =  [[] for a in range(n_ss_cols)]

    for i in tqdm(range(len(feature_list))):
        if i < len(feature_list)-1:
            feature_values = feature_list[i]
            corr_value = feature_values[-1]
            if corr_value >= min_val:
                for ss_index in range(n_ss_cols):
                    ss_val = feature_values[ss_index]
                    last_total = ss_totals[ss_index]
                    if not np.isnan(ss_val):
                        ss_totals[ss_index] += ss_val
                        ss_n[ss_index] += 1
                        #if corr_value > feature_list[i+1][-1]:
                        new_avg = ss_totals[ss_index] / ss_n[ss_index]
                        ss_avgs[ss_index].append(new_avg)
                        ths[ss_index].append(corr_value)
            else:
                break
    return [(ths[i], ss_avgs[i]) for i in range(len(ss_avgs))]

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

if __name__ == "__main__":
    file_ss = sys.argv[1]
    output_dir = sys.argv[2]
    #max_ss = float(sys.argv[3])
    
    if not os.path.exists(output_dir):
        print("creating output dir")
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    
    '''print("Reading minimum and maximum values for columns")
    unorm_mins, unorm_maxes, header = get_min_max(file_ss)
    
    min_max_path = output_dir+"/"+name+"min_max.tsv"
    with open(min_max_path, 'w') as stream:
        stream.write("\t".join(["Column Name","Min","Max"])+"\n")
        for i in range(len(unorm_mins)):
            stream.write("\t".join([str(x) for x in [header[i], unorm_mins[i], unorm_maxes[i]]]) + "\n")
    
    print(open(min_max_path, 'r').read())

    mins  = [round(x,2) for x in unorm_mins]
    maxes = [round(x,2) for x in unorm_maxes]

    mins = mins[1:]
    maxes = maxes[1:]'''
    column_values, header, ids, ss_col_names, single_ss = load_df_values(file_ss)
    corr_col_names = header[len(ss_col_names):]
    for i in range(len(column_values)):
        column = column_values[i]
        name = header[i]
        '''max_val = max_ss
        if name in corr_col_names:
            max_val = 1.0
            if name == "SOB" or name == "FSH":
                max_val = max(column)
        print("Maximum of " + name + " is " + str(max_val))'''
        if name == "SPR" or name == "PRS":
            print("Normalizing negative values for " + name)
            column_values[i] = [abs(x) for x in column_values[i]]
        if name == "SOB" or name == "FSH":
            print("Value inversion for " + name)
            column_values[i] = [max(0.0, 0.0-x+1.0) for x in column_values[i]]
        '''if max_val != 1.0:
            print("Normalizing values to range 0.0 -> 1.0 for " + name)
            column_values[i] = [min(x/max_val, 1.0) for x in column_values[i]]
        else:
            column_values[i] = [min(x, 1.0) for x in column_values[i]]
        if name in corr_col_names:
            print("Rounding values for " + name)
            column_values[i] = [round(x, 7) for x in column_values[i]]'''

    '''norm_negatives = [mins[i] < 0 for i in range(len(mins))]
    norm_by_max = [maxes[i] > 1.0 for i in range(len(maxes))]
    for i in range(len(column_values)):
        if norm_negatives[i]:
            print("Normalizing negative values of " + header[i])
            column_values[i] = [abs(x) for x in column_values[i]]
        if norm_by_max[i]:
            print("Normalizing values to range 0.0-1.0, for " + header[i])
            column_values[i] = [x/maxes[i] for x in column_values[i]]
        if header[i] in ["SOB","FSH"]:
            column_values[i] = [max(1.0-x, 0) for x in column_values[i]]'''
    print("Creating filters...")
    corr_filters = [filter_column(x, -99999) for x in tqdm(column_values[len(ss_col_names):])]
    all_corrs_filter = set.union(*corr_filters)
    i = 0
    ths_by_metric = {}
    for corr_col in tqdm(column_values[-6:]):
        corr_filter = corr_filters[i]
        corr_name = corr_col_names[i]
        print("Testing threasholds for " + corr_name)
        print("\tCreating list with values.")
        cols_list = [column_values[j] for j in range(len(ss_col_names))] + [corr_col]
        feature_list = lists_to_feature_list(cols_list, corr_filter)
        print("\tSorting list")
        feature_list.sort(key=lambda features: features[-1], reverse=True)
        print("\tFinding thresholds")
        ths_by_metric[corr_name] = test_threasholds(feature_list, min_val=0.0)
        i += 1

    print("Finding confidence levels")
    conf_levels = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
    conf_levels += [0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]
    conf_levels += [0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75]
    conf_levels += [0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.999]
    #conf_levels = [round(x, 2) for x in np.arange(0.1, 1.05, 0.05)]
    confs_by_metric = {corr_name: [get_confs(onto_ths[n], conf_levels) for n in range(len(onto_ths))]
                        for corr_name, onto_ths in tqdm(ths_by_metric.items())}

    print("Plotting")
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    subplots = []
    name_i = 0
    for i in range(len(axes)):
        for j in range(len(axes[i])):
            sub_plt = axes[i][j]
            subplots.append((corr_col_names[name_i], sub_plt))
            name_i += 1
    zord = 0
    colors = ["black", "orange", "darkred", "navy", "forestgreen"]
    for corr_name, sub_plt in tqdm(subplots):
        ss_ths = ths_by_metric[corr_name]
        
        '''mf_confs_raw, bp_confs_raw, cc_confs_raw = confs_by_metric[corr_name]

        mf_confs_x, mf_confs_y = get_not_none(mf_confs_raw, min_val=mf_th[0][-1])
        bp_confs_x, bp_confs_y = get_not_none(bp_confs_raw, min_val=bp_th[0][-1])
        cc_confs_x, cc_confs_y = get_not_none(cc_confs_raw, min_val=cc_th[0][-1])'''
        
        min_y = min([min(y) for x, y in ss_ths])
        min_x = min([min(x) for x, y in ss_ths])
        max_y = max([max(y) for x, y in ss_ths])
        max_x = max([max(x) for x, y in ss_ths])

        y_range = max_y - min_y
        sub_plt.yaxis.set_major_locator(MultipleLocator(y_range/5))
        sub_plt.yaxis.set_minor_locator(MultipleLocator(y_range/10))
        
        x_range = max_x - min_x
        sub_plt.xaxis.set_major_locator(MultipleLocator(x_range/5))
        sub_plt.xaxis.set_minor_locator(MultipleLocator(x_range/10))
        sub_plt.grid(True, which="major", axis="y", linestyle="--")

        for i in range(len(ss_ths)):
            x = ss_ths[i][0]
            y = ss_ths[i][1]
            sub_plt.plot(x, y, color=colors[i], zorder=zord)
            zord += 1

        '''sub_plt.scatter(mf_confs_x, mf_confs_y, marker='>', 
            color="lightcoral", s=40, alpha=0.75, zorder=zord)
        zord += 1
        sub_plt.scatter(bp_confs_x, bp_confs_y, marker='>', 
            color="slateblue", s=40, alpha=0.75, zorder=zord)
        zord += 1
        sub_plt.scatter(cc_confs_x, cc_confs_y, marker='>', 
            color="springgreen", s=40, alpha=0.75, zorder=zord)
        zord += 1'''

        sub_plt.set_title(translator(corr_name))
        sub_plt.set_ylabel("Similaridade Semântica Média")
        sub_plt.set_xlabel("Threshold de Correlação")

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
