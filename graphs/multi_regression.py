import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import random
import os
import math
import matplotlib.colors as colors
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression, Ridge
from itertools import combinations
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
import multiprocessing
import dcor
threads = max(2, multiprocessing.cpu_count())

translate_dict = {"cc": "Componente Celular", "CC": "Componente Celular",
                "mf": "Função Molecular", "MF": "Função Molecular",
                "bp": "Processo Biologico", "BP": "Processo Biologico",
                "max": "Símilaridade Semântica Máximo", 
                "avg": "Símilaridade Semântica Média",
                "fsh": "Fisher", "FSH": "Fisher",
                "sob": "Sobolev", "SOB": "Sobolev",
                "mic": "Maximal Information\nCoefficient", 
                "MIC": "Maximal Information\nCoefficient",
                "prs": "Pearson", "PRS": "Pearson",
                "spr": "Spearman", "SPR": "Spearman",
                "dc": "Distance\nCorrelation", "DC": "Distance\nCorrelation"}

full_pallete = ['#ff0000', '#dd0074', '#006600',
                '#000099', '#000000', '#00ff00']

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

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
    '''Find min and max value for all columns'''
    mins = []
    maxes = []
    for column in column_values:
        mins.append(np.nanmin(column))
        maxes.append(np.nanmax(column))
    norm_mins = [round(x,5) for x in mins]
    norm_maxes = [round(x,5) for x in maxes]
    #norm_mins = [math.floor(x) for x in norm_mins]
    #norm_maxes = [math.ceil(x) for x in norm_maxes]
    return mins, maxes, norm_mins, norm_maxes

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
    progress_bar.close()
    return columns, header_line[1:]

def get_combs(ns, values):
    '''Get combinations of several lengths for a list of values'''
    combs = set()
    for n in ns:
        for vals in combinations(values,n):
            new_set = list(set(vals))
            new_set.sort()
            combs.add(",".join([str(x) for x in new_set]))
    result_combs = []
    for comb in combs:
        splited = comb.split(",")
        result_combs.append([int(x) for x in splited])
    return result_combs

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (np.isfinite(x[i]) and not np.isnan(x[i])) and x[i] >= min_cor:
            f.add(i)
    return f

def get_rand_spokesman(row_index, row, filters, min_cor):
    not_null_cols = []
    for i in range(len(row)):
        if row_index in filters[i] and row[i] >= min_cor:
            not_null_cols.append(i)
    if len(not_null_cols) > 0:
        spokesman = random.choice(not_null_cols)
        spokesman_value = row[spokesman]
        return spokesman_value
    else:
        return None

def get_avg_value(row_index, row, filters, min_cor):
    valid_corrs = []
    for x in row:
        if x >= min_cor:
            valid_corrs.append(x)
    if len(valid_corrs) > 0:
        avg = np.nanmean(valid_corrs)
        if avg != np.nan:
            return avg
        else:
            return None
    else:
        return None

def get_first_value(row_index, row, filters, min_cor):
    value = row[0]
    if not math.isnan(value):
        return value
    else:
        return None

def get_norm_filter(columns, filters, min_cor, 
    min_values=300, n_bins=100, spokesman_selector=get_avg_value,
    smaller_frequency=None):
    min_bin_len = int(min_values/n_bins)
    #print("Filtering normal rows")
    bins = []
    min_bin = int(min_cor*n_bins)
    for i in range(n_bins):
        bins.append([])
    #print("Populating frequency table")
    for row_index in range(len(columns[0])):
        row_values = [columns[i][row_index] for i in range(len(columns))]
        spokesman_value = spokesman_selector(row_index, row_values, filters, min_cor)
        if spokesman_value != None:
            bin_i = int(spokesman_value*n_bins)
            if bin_i >= n_bins:
                #print(str(spokesman_value) + " -> " + str(bin_i))
                bin_i = n_bins-1
            bins[bin_i].append(row_index)
    if smaller_frequency == None:
        smaller_frequency_index = len(bins)-1
        for i in range(min_bin, len(bins)):
            l = len(bins[i])
            if l >= min_bin_len:
                if len(bins[smaller_frequency_index]) < min_bin_len:
                    smaller_frequency_index = i
                elif l < len(bins[smaller_frequency_index]):
                    smaller_frequency_index = i
        smaller_frequency = len(bins[smaller_frequency_index])
        #print("Bin " + str(smaller_frequency_index) 
        #    + " has the smaller number of items: " + str(smaller_frequency))
    #print(str([len(b) for b in bins]))
    #print("Randonly choosing from columns")
    for bin_i in range(len(bins)):
        bin_len = len(bins[bin_i])
        if bin_len > smaller_frequency:
            bins[bin_i] = np.random.choice(bins[bin_i], 
                            smaller_frequency, replace=False)
    
    #print("Finishing normalization filter")
    norm_filter = set()
    for b in bins:
        for value in b:
            norm_filter.add(value)
    return norm_filter

def prepare_data_for_regression(columns, valid_rows, set_name):
    '''Filter a set of columns, so that all values are valid'''
    y = []
    xs = []
    for i in valid_rows:
        y.append(columns[0][i])
        xs.append([column[i] for column in columns[1:]])
    '''if len(y) == 0:
        print("Zero values for " + set_name + ". Filter lengths:")
        print(str([len(f) for f in filters]))'''
    return xs, y

def make_svr_rbf(y, x_vars):
    '''Support vector regression with RBF kernel'''
    model = SVR(kernel='rbf', C=100, gamma=0.1, epsilon=.1)
    model.fit(x_vars, y)
    return model, {}

def make_svr_lin(y, x_vars):
    '''Support vector regression with linear kernel'''
    model = SVR(kernel='linear', C=100, gamma='auto')
    model.fit(x_vars, y)
    return model, {}

def make_svr_poly(y, x_vars):
    '''Support vector regression with polynomial kernel'''
    model = SVR(kernel='poly', C=100, gamma='auto', degree=3, epsilon=.1,
               coef0=1)
    model.fit(x_vars, y)
    return model, {}

def make_random_forest(y, x_vars):
    '''Random forest regression'''
    model = RandomForestRegressor(n_jobs = 1)
    model.fit(x_vars, y)
    return model, {}

def make_regression_poly(y, x_vars, d):
    '''Polynomial regression'''
    poly = PolynomialFeatures(degree=d)
    X_ = poly.fit_transform(x_vars)
    model = LinearRegression(n_jobs = threads)
    model.fit(X_, y)
    return model, {"poly": poly}

def make_regression_ridge(y, x_vars, d):
    model = make_pipeline(PolynomialFeatures(), Ridge())
    model.fit(x_vars, y)
    return model, {}

def make_regression_linear(y, x_vars):
    model = LinearRegression(n_jobs = threads)
    model.fit(x_vars, y)
    return model, {}


'''"'''

model_makers = {"LinearRegression": make_regression_linear,
            "PolynomialRegression(d=2)": 
                lambda y, x_vars: make_regression_poly(y, x_vars, 2),
            "PolynomialRegression(d=5)": 
                lambda y, x_vars: make_regression_poly(y, x_vars, 5),
            "PolynomialRidgeRegression(d=2)": 
                lambda y, x_vars: make_regression_ridge(y, x_vars, 2),
            "PolynomialRidgeRegression(d=5)": 
                lambda y, x_vars: make_regression_ridge(y, x_vars, 5),
            "RandomForestRegression": make_random_forest,
            "SupportVectorRegression_RBF": make_svr_rbf,
            "SupportVectorRegression_LIN": make_svr_lin,
            "SupportVectorRegression_POLY": make_svr_poly}

def test_model(result_model, attrs, test_y, test_x_vars):
    '''Tests a trained model and randonly selects test points to plot'''
    result_y = []
    score = 0.0
    rmse = 0.0
    if "poly" in attrs:
        poly_func = attrs["poly"]
        x_ = poly_func.fit_transform(test_x_vars)
        result_y = result_model.predict(x_)
        score = result_model.score(x_, test_y)
        rmse = np.sqrt(mean_squared_error(test_y, result_y))
    else:
        result_y = result_model.predict(test_x_vars)
        score = result_model.score(test_x_vars, test_y)
        rmse = np.sqrt(mean_squared_error(test_y, result_y))
    correlation = dcor.distance_correlation(np.array(test_y), np.array(result_y))
    frequencies, bins_x, bins_y = np.histogram2d(result_y, 
                            test_y, bins=[np.linspace(0.0, 1.0, 100), 
                                        np.linspace(0.0, 1.0, 100)])
    #n_points = min(len(test_y), 10000)
    #point_indexes = np.random.choice(list(range(len(test_y))), n_points, replace=False)
    #points_x = [test_y[i] for i in point_indexes]
    #points_y = [result_y[i] for i in point_indexes]
    return score, rmse, correlation, frequencies

def regress(comb_name, train_x, train_y, test_x, test_y, model_maker, return_dict):
    '''Create a regression model based on certain columns'''
    #progress_bar = tqdm(total=5)
    
    #print("Train samples: " + str([str(len(train_x)), str(len(train_y))]))
    #progress_bar.update(1)
    result_model, attrs = model_maker(train_y, train_x)
    #progress_bar.update(1)
    score, rmse, correlation, frequencies = test_model(result_model, attrs, test_y, test_x)
    #progress_bar.update(1)

    return_dict[comb_name] = (score, rmse, correlation, frequencies, len(train_x))

if __name__ == "__main__":
    #usage: multi_regression.py output_dir data_file min_correlation min_samples sub_perc
    output_dir = sys.argv[1]
    file_ss = sys.argv[2]
    min_cor = float(sys.argv[3])
    min_samples = 1000
    max_rows = int(sys.argv[4])
    n_bins = int(sys.argv[5])

    print("Loading data...")
    column_values, header = load_df_values(file_ss)

    print("Calculating min and max values for columns...")
    unorm_mins, unorm_maxes, mins, maxes = calc_min_max(column_values)
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    name = get_filename(file_ss).split(".")[0]
    ss_col_names = ["mf", "bp", "cc", "avg", "max"]
    
    min_max_path = output_dir+"/"+name+"min_max.tsv"
    with open(min_max_path, 'w') as stream:
        stream.write("\t".join(["Column Name","Not-Normalized Min","Not-Normalized Max"])+"\n")
        for i in range(len(mins)):
            stream.write("\t".join([str(x) 
                        for x in [header[i], unorm_mins[i], unorm_maxes[i]]]) + "\n")
        stream.write("\t".join(["Column Name","Min","Max"])+"\n")
        for i in range(len(mins)):
            stream.write("\t".join([str(x) 
                        for x in [header[i], mins[i],maxes[i]]]) + "\n")
    
    print(open(min_max_path, 'r').read())
    
    '''Normalizing negatives and high values to 0.0 -> 1.0'''
    norm_negatives = [mins[i] < 0 for i in range(len(mins))]
    norm_by_max = [maxes[i] > 1.0 for i in range(len(maxes))]
    for i in range(len(column_values)):
        if norm_negatives[i]:
            print("Normalizing negative values of " + header[i])
            column_values[i] = [abs(x) for x in column_values[i]]
        if norm_by_max[i]:
            print("Normalizing values to range 0.0-1.0, for " + header[i])
            column_values[i] = [x/maxes[i] for x in column_values[i]]
    
    combs = get_combs([1,2,3], [6,7,8,9,10])
    comb_strs = [", ".join([header[i] for i in comb]) for comb in combs]
    print("Combinations of metrics: "
        + str(comb_strs))

    print("Creating filters...")
    column_filters = [filter_column(x, 0.0) for x in tqdm(column_values[:5])] 
    column_filters += [filter_column(x, min_cor) for x in tqdm(column_values[5:])]

    print("Before normalization filter: ")                                           
    print(str([len(column_filter) for column_filter in column_filters[:5]]))
    max_bin_len = math.ceil(max_rows/n_bins)
    for col_i in tqdm(range(5)):
        print("Normalizing " + ss_col_names[col_i])
        norm_filter = get_norm_filter([column_values[col_i]], 
                                    [column_filters[col_i]],
                                    0.0, n_bins=n_bins,
                                    spokesman_selector=get_first_value,
                                    min_values=min_samples,
                                    smaller_frequency=max_bin_len
                                )
        column_filters[col_i] = column_filters[col_i].intersection(norm_filter)
    print("After normalization filter: ")                                           
    print(str([len(column_filter) for column_filter in column_filters[:5]]))
    #print("Separating train and test sets...")
    elements_len = len(column_values[i])
    print("Making regressions...")
    '''Use average and max semantic similarity values'''
    y_indexes = [0,1,2,3,4]
    progress_bar = tqdm(total=len(y_indexes)*len(model_makers.items())*len(combs))
    comb_results = {}
    '''Each y_index is for one semantic similarity value (cc,mf,bp,avg or max)'''
    for y_index in y_indexes:
        '''Each i value is for one combination of metrics'''
        for i in range(len(combs)):
            regression_cols = [column_values[index] for index in [y_index] + combs[i]]
            regression_filters = [column_filters[index] for index in [y_index] + combs[i]]
            regression_master_filter = set.intersection(*regression_filters)
            train_length = int(len(regression_master_filter)*0.7)
            rand_valid_items = np.random.choice(list(regression_master_filter), 
                                train_length, 
                                replace=False)
            train_filter = set(rand_valid_items)
            test_filter = regression_master_filter - train_filter
            
            test_x, test_y = prepare_data_for_regression(regression_cols, 
                                                        test_filter,
                                                        (translator(header[y_index])
                                                        + " " + comb_strs[i]))
            train_x, train_y = prepare_data_for_regression(regression_cols, 
                                                        train_filter,
                                                        (translator(header[y_index]) 
                                                        + " " + comb_strs[i]))
            if len(train_y) >= min_samples and len(test_y) >= min_samples*0.3:
                '''Make and test each model'''
                maker_items = list(model_makers.items())
                maker_item_chunks = chunks(maker_items, threads)
                for maker_chunk in maker_item_chunks:
                    manager = multiprocessing.Manager()
                    return_dict = manager.dict()
                    processes = []
                    last_pid = 0
                    
                    '''Start processes'''
                    for model_name, model_maker in maker_chunk:
                        comb_name = (model_name 
                                    + "\nMétricas: " + comb_strs[i] + "\n" 
                                    + translator(header[y_index]))
                        p = multiprocessing.Process(target=regress, 
                            args=(comb_name, train_x, train_y, test_x, test_y,
                                model_maker, return_dict, ))
                        processes.append(p)
                        p.start()
                    
                    '''Join them'''
                    for p in processes:
                        p.join()
                    progress_bar.update(len(processes))
                    '''Retrieve results'''
                    for key, value in return_dict.items():
                        comb_results[key] = value
                    
                    manager._process.terminate()
                    manager.shutdown()
                    del manager
            else:
                print("Not enough samples: " + str((len(train_y),len(test_y))))
                progress_bar.update(len(model_makers.items()))
    results = list(comb_results.items())
    '''for comb_str, result_items in results:
        score, real_y, predicted_y = result_items
        print("For " + comb_str + ": " + str(score))'''
    print("\n")
    valid_results = []
    for comb_name, result_items in results:
        score, rmse, correlation, frequencies, train_samples = result_items
        if score == -100000:
            print("Too few training samples for " 
                + comb_name 
                + ": " + str(train_samples))
        elif score == -200000:
            print("Too few testing samples for " 
                + comb_name 
                + ": " + str(train_samples))
        else:
            valid_results.append((comb_name, result_items))
    
    def plot_regressions(subplots, sort_by, valid_results, reverse):    
        valid_results.sort(key=lambda x: abs(x[1][sort_by]), reverse=reverse)

        if sort_by == 0:
            print("\nBest scores: " )
        elif sort_by == 1:
            print("\nBest RMSE values: " )
        elif sort_by == 2:
            print("\nBest correlation values: " )

        i = 0
        top_items = valid_results[-8:]
        if reverse == True:
            top_items = valid_results[:8]
            top_items = top_items[::-1]
        for comb_name, result_items in top_items:
            subplot = subplots[i]
            subplot.set_ylim([0.0,1.0])
            score, rmse, dc_corr, frequencies, train_samples = result_items
            subplot.pcolormesh(np.linspace(0.0,1.0,100), 
                                np.linspace(0.0,1.0,100), frequencies,
                                norm=colors.SymLogNorm(linthresh = 0.03,
                                                    vmin=frequencies.min(), 
                                                    vmax=frequencies.max()))
            subplot.plot(np.linspace(0.0,1.0,100), np.linspace(0.0,1.0,100))
            #subplot.plot(points[0], points[1], '.', markersize=0.5)
            subplot.set_title(comb_name + "\nScore=" + str(score)
                            +"\nRMSE="+str(rmse)
                            +"\nCorrelação="+str(dc_corr)
                            +"\nAmostras="+str(train_samples))
            subplot.set_xlabel("Semantic Similarity")
            subplot.set_ylabel("Predicted Semantic Similarity")
            if sort_by == 1:
                print("For " + comb_name.replace("\n", "; ") + ": " + str(rmse))
            elif sort_by == 0:
                print("For " + comb_name.replace("\n", "; ") + ": " + str(score))
            elif sort_by == 2:
                print("For " + comb_name.replace("\n", "; ") + ": " + str(dc_corr))
            i += 1
        
    fig, axes = plt.subplots(3, 8, figsize=(8*4, 20))
    plot_regressions(axes[0], 0, valid_results, False)
    plot_regressions(axes[1], 1, valid_results, True)
    plot_regressions(axes[2], 2, valid_results, False)
    fig.tight_layout()
    out_file = output_dir+"/multi_regression-"+str(min_cor)+".png"
    print(out_file)
    fig.savefig(out_file, bbox_inches='tight')
