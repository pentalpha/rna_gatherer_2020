import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import math
import matplotlib.colors as colors
from sklearn.metrics import mean_squared_error

translate_dict = {"cc": "Componente Celular", "CC": "Componente Celular",
                "mf": "Função Molecular", "MF": "Função Molecular",
                "bp": "Processo Biologico", "BP": "Processo Biologico",
                "max": "Máximo", 
                "avg": "Média",
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

def polyfiter_factory(train_x, train_y, degree):
    coefficients = np.polyfit(train_x, train_y, degree)
    poly = np.poly1d(coefficients)
    return poly

polyfits = {"polyfit(d=3)": lambda x,y: polyfiter_factory(x, y, 3),
            "polyfit(d=5)": lambda x,y: polyfiter_factory(x, y, 5),
            "polyfit(d=9)": lambda x,y: polyfiter_factory(x, y, 9)}

def regression_test(test_x, test_y, func):
    new_y = func(test_x)
    rmse = np.sqrt(mean_squared_error(test_y, new_y))
    return rmse

def get_frequencies_2d(file_ss, index_x, index_y, converter, my_bins, regressors):
    min_cells = max(index_x+1, index_y+1)
    #frequencies = np.array([])
    line_limit = 10000
    failed_lines = 0
    with open(file_ss, 'r') as stream:
        values = []
        header_line = stream.readline().rstrip("\n").split("\t")
        line_chunk, more_lines = get_lines(stream, min_len=min_cells)
        while more_lines:
            for x_str, y_str in [(cells[index_x], cells[index_y]) 
                                for cells in line_chunk]:
                try:
                    new_values = (converter(x_str), converter(y_str))
                    if np.isfinite(new_values[0]) and np.isfinite(new_values[1]):
                        values.append(new_values)
                except ValueError:
                    failed_lines += 1
            line_chunk, more_lines = get_lines(stream, min_len=min_cells)
        frequencies, bins_x, bins_y = np.histogram2d([y for x,y in values], 
                                [x for x,y in values],
                                bins=[my_bins[index_y], my_bins[index_x]])
        subset_indexes = set(np.random.choice(list(range(len(values))), 
                                            int(len(values)*0.6), replace=False))
        x_train = []
        x_test = []
        y_train = []
        y_test = []
        for i in range(len(values)):
            if i in subset_indexes:
                x_train.append(values[i][0])
                y_train.append(values[i][1])
            else:
                x_test.append(values[i][0])
                y_test.append(values[i][1])
        to_plot = np.random.choice(list(range(len(x_test))), 
                                        2000, replace=False)
        points_to_plot = [(x_test[index],y_test[index])
                            for index in to_plot]
        predictions = []
        for reg_name, regressor_maker in regressors.items():
            regressor = regressor_maker(x_train, y_train)
            rmse = regression_test(x_test, y_test, regressor)
            predictions.append((reg_name, rmse, regressor))

        return frequencies, points_to_plot, predictions
    return None

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
    mins = [round(x,1) for x in mins]
    maxes = [round(x,1) for x in maxes]
    
    cols = header[6:]
    col_indexes = list(range(6,6+len(cols)))

    bins = [np.linspace(mins[i], maxes[i], 75) for i in range(len(mins))]
    progress_bar = tqdm(total=len(cols)*5)
    fig, freq_ax = plt.subplots(len(ss_col_names), len(cols),
                                figsize=(3*len(cols), 3*len(ss_col_names)))
    regress_fig, regress_ax = plt.subplots(len(ss_col_names), len(cols),
                                figsize=(4*len(cols), 4*len(ss_col_names)))
    print(str(len(regress_ax)) + " columns of plots")
    print(str(len(regress_ax[0])) + " rows of plots")
    pallete = make_pallete(list(polyfits.keys())+["input"])
    
    for i in range(len(cols)):

        col = cols[i]
        #print("Comparing " + col + " ("+str(col_indexes[i])+") to..")
        for j in range(len(ss_col_names)):
            #print("\t" + ss_col_names[j] + " ("+str(j+1)+").")
            freqs_2d, points_to_plot, regressions = get_frequencies_2d(
                                        file_ss, col_indexes[i], j+1,
                                        float, bins, polyfits)
            #print(str(freqs_2d))
            edges_x = bins[col_indexes[i]]
            edges_y = bins[j+1]
            
            freq_ax[j][i].pcolormesh(edges_x, edges_y, freqs_2d,
                            norm=colors.SymLogNorm(linthresh = 0.03,
                                                    vmin=freqs_2d.min(), 
                                                    vmax=freqs_2d.max()))

            regress_ax[j][i].plot([x for x,y in points_to_plot], 
                                    [y for x,y in points_to_plot], 
                                    '.',
                                    color=pallete["input"],
                                    label="Test Points",
                                    markersize=0.4)
            example_x = np.linspace(edges_x[0], edges_x[-1], 1000)
            for regression_name, rmse, regressor in regressions:
                predicted_y = [regressor(x) for x in example_x]
                regress_ax[j][i].plot(example_x, predicted_y,
                                    color=pallete[regression_name],
                                    label=regression_name+":"+str(rmse)[0:6]+" RMSE",
                                    linewidth=1)

            regress_ax[j][i].legend()
            progress_bar.update(1)

    for i in range(len(cols)):
        freq_ax[0][i].set_xlabel(translator(cols[i]))
        freq_ax[0][i].xaxis.set_label_position('top')

        regress_ax[0][i].set_xlabel(translator(cols[i]))
        regress_ax[0][i].xaxis.set_label_position('top')
        print(str((i,0)) + "xlabel = " + translator(cols[i]))
    for j in range(len(ss_col_names)):
        freq_ax[j][0].set_ylabel(translator(ss_col_names[j]))
        regress_ax[j][0].set_ylabel(translator(ss_col_names[j]))
        print(str((j,0)) + "ylabel = " + translator(ss_col_names[j]))

    fig.tight_layout()
    fig.savefig(output_dir+"/corr_vs_ss.png", bbox_inches='tight')

    regress_fig.tight_layout()
    regress_fig.savefig(output_dir+"/corr_vs_ss_regression.png", bbox_inches='tight')