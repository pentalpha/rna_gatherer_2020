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

def calc_min_max(column_values):
    mins = []
    maxes = []
    for column in column_values:
        mins.append(np.nanmin(column))
        maxes.append(np.nanmax(column))
    return mins, maxes

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

def load_df_values(file_ss):
    min_cells = 1+5+6
    columns = []
    for n_cell in range(min_cells-1):
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
                for i in range(1,len(np_values)):
                    columns[i-1].append(np_values[i])
            else:
                failed_lines += 1
        except ValueError:
            failed_lines += 1
        raw_line = stream.readline()
        progress_bar.update(1)
    print("Failed lines: " + str(failed_lines))
    stream.close()
    return columns, header_line[1:]

def calc_frequencies(column_x, column_y, bin_x, bin_y):
    frequencies, bins_x, bins_y = np.histogram2d(column_y, 
                            column_x, bins=[bin_y, bin_x])
    return frequencies

def make_regression(column_x, column_y, regressors):
    subset_indexes = set(np.random.choice(list(range(len(column_x))), 
                                        int(len(column_x)*0.6), 
                                        replace=False))
    x_train = []
    x_test = []
    y_train = []
    y_test = []
    for i in range(len(column_x)):
        if np.isfinite(column_x[i]) and np.isfinite(column_y[i]):
            if i in subset_indexes:
                x_train.append(column_x[i])
                y_train.append(column_y[i])
            else:
                x_test.append(column_x[i])
                y_test.append(column_y[i])
    to_plot = np.random.choice(list(range(len(x_test))), 
                                    2000, replace=False)
    points_to_plot = [(x_test[index],y_test[index])
                        for index in to_plot]
    predictions = []
    for reg_name, regressor_maker in regressors.items():
        regressor = regressor_maker(x_train, y_train)
        rmse = regression_test(x_test, y_test, regressor)
        predictions.append((reg_name, rmse, regressor))
    return points_to_plot, predictions

if __name__ == "__main__":
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]

    print("Loading data...")
    column_values, header = load_df_values(file_ss)

    print("Calculating min and max values for columns...")
    mins, maxes = calc_min_max(column_values)

    print("Preparing to plot...")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    ss_col_names = ["mf", "bp", "cc", "avg", "max"]
    
    min_max_path = output_dir+"/"+name+"min_max.tsv"
    with open(min_max_path, 'w') as stream:
        stream.write("\t".join(["Column Name","Min","Max"])+"\n")
        for i in range(len(mins)):
            stream.write("\t".join([str(x) for x in [header[i], mins[i],maxes[i]]]) + "\n")
    
    print(open(min_max_path, 'r').read())
    mins = [round(x,1) for x in mins]
    maxes = [round(x,1) for x in maxes]
    
    cols = header[5:]

    bins = [np.linspace(mins[i], maxes[i], 75) for i in range(len(mins))]
    
    fig, freq_ax = plt.subplots(len(ss_col_names), len(cols),
                                figsize=(3*len(cols), 3*len(ss_col_names)))
    regress_fig, regress_ax = plt.subplots(len(ss_col_names), len(cols),
                                figsize=(4*len(cols), 4*len(ss_col_names)))
    print(str(len(regress_ax)) + " columns of plots")
    print(str(len(regress_ax[0])) + " rows of plots")
    pallete = make_pallete(list(polyfits.keys())+["input"])
    
    print("Creating individual plots...")
    progress_bar = tqdm(total=len(cols)*len(ss_col_names))
    for i in range(len(cols)):
        col = cols[i]
        #print("Comparing " + col + " ("+str(col_indexes[i])+") to..")
        for j in range(len(ss_col_names)):
            #print("\t" + ss_col_names[j] + " ("+str(j+1)+").")
            edges_x = bins[len(ss_col_names)+i]
            edges_y = bins[j]
            freqs_2d = calc_frequencies(column_values[len(ss_col_names)+i], 
                                        column_values[j],
                                        edges_x, edges_y)
            
            freq_ax[j][i].pcolormesh(edges_x, edges_y, freqs_2d,
                            norm=colors.SymLogNorm(linthresh = 0.03,
                                                    vmin=freqs_2d.min(), 
                                                    vmax=freqs_2d.max()))

            points_to_plot, predictions = make_regression(
                                            column_values[len(ss_col_names)+i],
                                            column_values[j], 
                                            polyfits)

            regress_ax[j][i].plot([x for x,y in points_to_plot], 
                                    [y for x,y in points_to_plot], 
                                    '.',
                                    color=pallete["input"],
                                    label="Test Points",
                                    markersize=0.4)
            example_x = np.linspace(edges_x[0], edges_x[-1], 1000)
            for regression_name, rmse, regressor in predictions:
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