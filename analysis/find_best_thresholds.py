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
    min_cells = 3+5+6
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
    return columns[3:], header_line[3:], [int(x) for x in columns[0]]

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

def test_threasholds(feature_list, min_val=0.6):
    ss_totals = [0.0,0.0,0.0]
    ss_n = [0,0,0]
    ths = [[], [], []]
    ss_avgs = [[], [], []]

    for i in range(len(feature_list)):
        if i != len(feature_list)-1:
            feature_values = feature_list[i]
            corr_value = feature_values[-1]
            if corr_value >= min_val:
                for ss_index in [0,1,2]:
                    ss_val = feature_values[ss_index]
                    if not np.isnan(ss_val):
                        ss_totals[ss_index] += ss_val
                        ss_n[ss_index] += 1
                        if corr_value > feature_list[i+1][-1]:
                            new_avg = ss_totals[ss_index] / ss_n[ss_index]
                            ss_avgs[ss_index].append(new_avg)
                            ths[ss_index].append(corr_value)
            else:
                break
    return (ths[0], ss_avgs[0]), (ths[1], ss_avgs[1]), (ths[2], ss_avgs[2])

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
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]
    max_ss = float(sys.argv[3])
    
    if not os.path.exists(output_dir):
        print("creating output dir")
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    ss_col_names = ["mf", "bp", "cc", "avg", "max"]
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
    column_values, header, ids = load_df_values(file_ss)
    corr_col_names = header[5:]
    for i in range(len(column_values)):
        column = column_values[i]
        name = header[i]
        max_val = max_ss
        if name in corr_col_names:
            max_val = 1.0
            if name == "SOB" or name == "FSH":
                max_val = max(column)
        print("Maximum of " + name + " is " + str(max_val))
        if name == "SPR" or name == "PRS":
            print("Normalizing negative values for " + name)
            column_values[i] = [abs(x) for x in column_values[i]]
        if name == "SOB" or name == "FSH":
            print("Value inversion for " + name)
            column_values[i] = [max(max_val-x, 0.0) for x in column_values[i]]
        if max_val != 1.0:
            print("Normalizing values to range 0.0 -> 1.0 for " + name)
            column_values[i] = [min(x/max_val, 1.0) for x in column_values[i]]
        else:
            column_values[i] = [min(x, 1.0) for x in column_values[i]]
        if name in corr_col_names:
            print("Rounding values for " + name)
            column_values[i] = [round(x, 7) for x in column_values[i]]

    corr_col_names = header[5:]

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
    corr_filters = [filter_column(x, 0.0) for x in tqdm(column_values[5:])]
    all_corrs_filter = set.union(*corr_filters)
    i = 0
    ths_by_metric = {}
    for corr_col in tqdm(column_values[5:]):
        corr_filter = corr_filters[i]
        corr_name = corr_col_names[i]
        print("Testing threasholds for " + corr_name)
        print("\tCreating list with values.")
        feature_list = lists_to_feature_list([column_values[0], column_values[1], 
                                            column_values[2], corr_col], corr_filter)
        print("\tSorting list")
        feature_list.sort(key=lambda features: features[-1], reverse=True)
        print("\tFinding thresholds")
        ths_by_metric[corr_name] = test_threasholds(feature_list, min_val=0.0)
        i += 1

    print("Finding confidence levels")
    conf_levels = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
    conf_levels += [0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]
    conf_levels += [0.525, 0.55, 0.575, 0.6, 0.625, 0.65]
    #conf_levels = [round(x, 2) for x in np.arange(0.1, 1.05, 0.05)]
    confs_by_metric = {corr_name: (get_confs(onto_ths[0], conf_levels), 
                                    get_confs(onto_ths[1], conf_levels), 
                                    get_confs(onto_ths[2], conf_levels))
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
    for corr_name, sub_plt in tqdm(subplots):
        mf_th, bp_th, cc_th = ths_by_metric[corr_name]
        '''mf_confs_raw, bp_confs_raw, cc_confs_raw = confs_by_metric[corr_name]

        mf_confs_x, mf_confs_y = get_not_none(mf_confs_raw, min_val=mf_th[0][-1])
        bp_confs_x, bp_confs_y = get_not_none(bp_confs_raw, min_val=bp_th[0][-1])
        cc_confs_x, cc_confs_y = get_not_none(cc_confs_raw, min_val=cc_th[0][-1])'''
        
        sub_plt.yaxis.set_major_locator(MultipleLocator(0.05))
        sub_plt.yaxis.set_minor_locator(MultipleLocator(0.01))
        sub_plt.xaxis.set_major_locator(MultipleLocator(0.1))
        sub_plt.xaxis.set_minor_locator(MultipleLocator(0.05))
        sub_plt.grid(True, which="major", axis="y", linestyle="--")

        sub_plt.plot(mf_th[0], mf_th[1], color="darkred", 
            zorder=zord)
        zord += 1
        sub_plt.plot(bp_th[0], bp_th[1], color="navy", 
            zorder=zord)
        zord += 1
        sub_plt.plot(cc_th[0], cc_th[1], color="forestgreen", 
            zorder=zord)
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

    for ss_index in [0,1,2]:
        ss_name = ss_col_names[ss_index]
        output_table = output_dir+"/confidence_intervals-"+ss_name.upper()+".tsv"
        out_str = "\t".join(["Average Semantic Similarity"] + [str(val) for val in conf_levels])+"\n"
        out_str += "\t".join(["Interval Name"] + [str(i) for i in range(len(conf_levels))])+"\n"
        for corr_name, onto_ths in confs_by_metric.items():
            ths = onto_ths[ss_index]
            out_str += "\t".join([corr_name.upper()] + [str(th) for th, conf_level in ths])+"\n"
        open(output_table,'w').write(out_str)
