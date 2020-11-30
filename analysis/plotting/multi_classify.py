# -*- coding: utf-8 -*-
import sys
import numpy as np
from tqdm import tqdm
import os
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC, SVC, NuSVC
import sklearn
import multiprocessing
from joblib import dump
from sklearn.neighbors import *
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
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

def calc_min_max(column_values):
    '''Find min and max value for all columns'''
    mins = []
    maxes = []
    for column in column_values:
        mins.append(np.nanmin(column))
        maxes.append(np.nanmax(column))
    norm_mins = [round(x,2) for x in mins]
    norm_maxes = [round(x,2) for x in maxes]
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

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (not np.isfinite(x[i]) or np.isnan(x[i])) or x[i] < min_cor:
            f.add(i)
    return f

def get_norm_mask(column, invalid_rows, n_bins=100, subset_len=10000,
                  train_size=0.7):
    bins = []
    for i in range(n_bins):
        bins.append([])
    print("Populating bins")
    for row_index in range(len(column)):
        if not row_index in invalid_rows:
            value = column[row_index]
            bin_i = int(value*n_bins)
            if bin_i >= n_bins:
                bin_i = n_bins-1
            bins[bin_i].append(row_index)
    #print("Original bins\n"+str([len(a) for a in bins]))
    print("Choosing from bins")
    #mask = set()
    items_to_take = subset_len
    items_per_bin = int(round(items_to_take/n_bins,0))
    for i in range(len(bins))[::-1]:
        '''print("items_to_take="+str(items_to_take)
              +", items_per_bin="+str(items_per_bin)
              +", len(bin_i)="+str(len(bins[i])))'''
        if items_per_bin < len(bins[i]):
            bins[i] = np.random.choice(bins[i], 
                            items_per_bin, replace=False)
        #print("took="+str(len(bins[i])))
        n_bins -= 1
        items_to_take -= len(bins[i])
        if n_bins > 0:
            items_per_bin = int(round(items_to_take/n_bins,0))
    print("Subset bins\n"+str([len(a) for a in bins]))
    print("Creating mask from bins")
    train_mask = set()
    test_mask = set()
    for b in bins:
        if len(b) > 0:
            n_to_train = int(round(len(b)*train_size, 0))
            train_selected = set()
            test_selected = set()
            if n_to_train < len(b) and n_to_train > 0:
                train_selected = set(np.random.choice(b,
                                                      n_to_train,
                                                      replace=False))
                for value in b:
                    if not value in train_selected:
                        test_selected.add(value)
            else:
                train_selected = set(b)
            
            for value in train_selected:
                train_mask.add(value)
            for value in test_selected:
                test_mask.add(value)
    #print("Total values selected: " + str(len(mask)))
    return train_mask, test_mask

def ss_to_class(ss_column, min_value):
    return [val >= min_value for val in ss_column]

def get_classify_data(column_values, all_corrs_filter, 
                      min_ss, ss_index, subset_len):
    print("Filtering semantic similarity values")
    ss_filter = filter_column(column_values[ss_index], 0.0)
    print("Filtering correlation values")
    all_filter = all_corrs_filter.union(ss_filter)
    print("Preparing row mask")
    mask = get_norm_mask(column_values[ss_index], all_filter, 
                         n_bins=50, subset_len=subset_len)
    
    print("Preparing data for ML")
    data = prepare_data_for_ml(
        [ss_to_class(column_values[ss_index], min_ss)] + column_values[5:], 
        mask)
    '''print("Evaluating subset")
    above = [0,0,0,0,0,0,0,0,0,0]
    for value in column_values[ss_index]:
        for i in range(len(above)):
            if value >= i/10:
                above[i] += 1
    above_dict = {str(i/10): above[i] for i in range(len(above))}
    print("Number of values above each minimum: ")
    print(str(above_dict))
    trues = 0
    for x in data[1]:
        if x:
            trues += 1
    print("True target values at this subset: " + str(trues))'''
    return data

def try_to_classify(data, iters, train_size):
    print("Splitting train and test sets")
    X_train, X_test, y_train, y_test = train_test_split(data[0], 
                                                        data[1], 
                                                        train_size=train_size)
    clf = MLPClassifier(random_state=0, max_iter=iters)
    print("Fitting")
    clf.fit(X_train, y_train)
    print("Testing")
    res = clf.predict(X_test)
    score = clf.score(X_test, y_test)
    false_positive = 0
    false_negative = 0
    for i in range(len(y_test)):
        original = y_test[i]
        prediction = res[i]
        if prediction and not original:
            false_positive += 1
        elif not prediction and original:
            false_negative += 1
    false_positive_rate = false_positive/len(res)
    false_negative_rate = false_negative/len(res)
    
    return (clf, {"Score": score, "FP": false_positive_rate,
            "FN": false_negative_rate})

def get_conf_model(cols, corr_filter, conf_level, ss_index, 
                   subset_len, train_size, steps):
    data = get_classify_data(cols, corr_filter,
                               conf_level, ss_index, subset_len)
    model, stats = try_to_classify(data, steps, train_size)
    return model, stats
# %%
def prepare_data_for_ml(y_col, x_cols, min_ss, valid_rows):
    '''Filter a set of columns, so that all values are valid'''
    y = []
    xs = []
    for i in valid_rows:
        if i <= len(y_col):
            y.append(y_col[i] >= min_ss)
            xs.append([column[i] for column in x_cols])
        else:
            print(str(i)+" outside of list of len=" +str(len(y_col)))
    return xs, y

def get_wanted_negatives(positives, positive_rate):
        return (positives - (positive_rate * positives)) / positive_rate 
    
def get_usable_fraction(min_val,
                        invalid_rows, 
                        ss_col,
                        positive_rate=0.4,
                        train_part=0.7):
    print("Finding positives")
    positives = []
    negatives = []
    for i in range(len(ss_col)):
        if not i in invalid_rows:
            if ss_col[i] >= min_val:
                positives.append(i)
            else:
                negatives.append(i)
    to_take = int(round(len(positives)*train_part,0))
    train_positives = set(np.random.choice(positives,
                                       to_take,
                                       replace=False))
    test_positives = []
    for x in positives:
        if not x in train_positives:
            test_positives.append(x)
    available_rows = len(ss_col)-len(invalid_rows)
    original_positive_rate = len(positives)/available_rows
    positives = list(train_positives)
    n_negatives = len(negatives)
    wanted_negatives = int(round(
        get_wanted_negatives(len(positives), positive_rate),0))
    if wanted_negatives < n_negatives:
        negatives = np.random.choice(negatives,
                                         wanted_negatives,
                                         replace=False)
    else:
        wanted_positives = int(round(
            get_wanted_negatives(len(negatives), 
                                 1.0-positive_rate),
                                0))
        positives = np.random.choice(positives,
                                     wanted_positives,
                                     replace=False)
    used = (len(positives)+len(negatives))/available_rows
    
    unused_negatives = []
    train_negatives_set = set(negatives)
    for i in range(len(ss_col)):
        if not i in invalid_rows and not i in train_negatives_set:
            if ss_col[i] < min_val:
                unused_negatives.append(i)
    neg_rate = 1.0-original_positive_rate
    unused_to_take = (neg_rate*len(test_positives))/original_positive_rate
    test_negatives = unused_negatives
    if unused_to_take < len(unused_negatives):
        test_negatives = list(np.random.choice(unused_negatives,
                                               int(round(unused_to_take,0)),
                                               replace=False))
    test_positives = list(test_positives)
    
    print("Training positives/negatives: " 
          + str([len(positives),len(negatives)]))
    print("Test positives/negatives: " 
          + str([len(test_positives),len(test_negatives)]))
    return (np.append(positives, negatives),
            np.append(test_positives, test_negatives),
            used)

def get_train_test_simple(onto_filter, onto_col,
                          max_subset_len=100000,
                          train_part=0.7):
    selected = []
    for i in (range(len(onto_col))):
        if not i in onto_filter:
            selected.append(i)
    if len(selected) > max_subset_len:
        selected = list(np.random.choice(selected,
                                    max_subset_len,
                                    replace=False))
    to_train = int(round(train_part*len(selected),0))
    train_set = set(np.random.choice(selected,
                                  to_train,
                                  replace=False))
    test_set = set()
    for i in selected:
        if not i in train_set:
            test_set.add(i)
    return train_set, test_set

def get_train_test(min_val, onto_filter, onto_col,
                   positive_rate=0.5, max_subset_len=100000,
                   train_part=0.7):
    selected_train, selected_test, used = get_usable_fraction(
                                            min_val,
                                            onto_filter, 
                                            onto_col,
                                            positive_rate=positive_rate,
                                            train_part=train_part)
    
    return set(selected_train), set(selected_test)

def safe_div(a,b):
    if b == 0:
        return np.nan
    else:
        return a / b
    

def eval_model(model, test_x, test_y):
    print("Testing")
    res = model.predict(test_x)
    score = model.score(test_x, test_y)
    false_positive = 0
    false_negative = 0
    true_positive = 0
    true_negative = 0
    for i in range(len(test_y)):
        original = test_y[i]
        prediction = res[i]
        if prediction:
            if original:
                true_positive += 1
            else:
                false_positive += 1
        else:
            if original:
                false_negative += 1
            else:
                true_negative += 1
    recall = safe_div(true_positive, true_positive+false_negative)
    precision = safe_div(true_positive, true_positive+false_positive)
    false_positive_rate = safe_div(false_positive, false_positive+true_negative)
    false_negative_rate = safe_div(false_negative, false_negative+true_positive)
    neg_prec = safe_div(true_negative, false_negative+true_negative)
    pos_prec = safe_div(true_positive, false_positive+true_positive)
    return {"Score": score, 
            "Recall": recall*100,
            "Precision": precision*100,
            "FP Rate": false_positive_rate*100,
            "FN Rate": false_negative_rate*100,
            "FP:": false_positive,
            "FN:": false_negative,
            "TP:": true_positive,
            "TN:": true_negative,
            "+ Precision": pos_prec*100,
            "- Precision": neg_prec*100}

def fit_and_eval_ready_sets(clf, train_x, train_y, test_x, test_y):
    clf.fit(train_x, train_y)
    stats = eval_model(clf, test_x, test_y)
    print(str(stats))
    
    return clf, stats

def fit_and_eval(clf, column_values, ss_index, conf_level,
                 onto_filter,
                 positive_rate=0.5, max_subset_len=100000,
                 train_part=0.7):
    
    train, test = get_train_test(conf_level, 
                                 onto_filter, 
                                 column_values[ss_index],
                                 positive_rate=positive_rate,
                                 max_subset_len=max_subset_len,
                                 train_part=train_part)
    print("Train and test sets: " + str([len(train), len(test)]))
    
    
    train_x, train_y = prepare_data_for_ml(column_values[ss_index],
                                       column_values[5:], 
                                       conf_level, train)
    test_x, test_y = prepare_data_for_ml(column_values[ss_index],
                                           column_values[5:], 
                                           conf_level, test)
    
    print("Fitting")
    clf.fit(train_x, train_y)
    stats = eval_model(clf, test_x, test_y)
    print(str(stats))
    
    return clf, stats

# %%
if __name__ == "__main__":
    '''file_ss = sys.argv[1]
    output_dir = sys.argv[2]
    subset_size = int(sys.argv[3])'''
    file_ss = "zfin_resnik.tsv.2m"
    output_dir = "resnik_2m_100k"
    subset_size = 100000
    steps = 100
    if len(sys.argv) > 4:
        steps = int(sys.argv[4])
    
    
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
    print("Creating filters...")
    ss_indexes = [0,1,2,3,4]
    all_corrs_filter = set.union(*[filter_column(x, 0.0) 
                                          for x in tqdm(column_values[5:])])
    onto_filters = [filter_column(column_values[ss_index], 0.0)
                        for ss_index in tqdm(ss_indexes)]
    onto_filters = [f.union(all_corrs_filter) for f in tqdm(onto_filters)]
    
    mins_to_try = [0.2]
    
    classifiers = {}
    max_depths = [None]
    criterions = ["gini"]
    classes_weights = [None, "balanced", 
                       {True: 1, False:3}, {True: 1, False:5}]
    for max_dp in max_depths:
        for crit in criterions:
            for w in classes_weights:
                name1 = ("dtc max_dp="+str(max_dp)
                         +" crit="+str(crit) 
                         + " class_w="+str(w))
                name2 = name1.replace("dtc ", "rfc ")
                classifiers[name2] = RandomForestClassifier(
                    max_depth=max_dp,
                    criterion=crit,
                    class_weight=w,
                    n_jobs=-1)
    # %%
    ss_index = 4
    train_part = 0.7
    positive_rate=0.01
    max_subset_len=300000
    train, test = get_train_test_simple(onto_filters[ss_index], 
                                        column_values[ss_index],
                                        max_subset_len=max_subset_len,
                                        train_part=train_part)
    
    # %%
    steps = 200
    hidden_layer_sizes=[128]
    solver="lbfgs"
    validation_fraction=0.2
    results = []
    min_positives = 500
    prog = tqdm(total=len(mins_to_try)*len(classifiers.items()))
    for min_val in mins_to_try:
        train_x, train_y = prepare_data_for_ml(column_values[ss_index],
                                               column_values[5:],
                                               min_val, train)
        train_positives = 0
        for a in train_y:
            if a:
                train_positives += 1
        if train_positives >= min_positives*train_part:
            test_x, test_y =   prepare_data_for_ml(column_values[ss_index],
                                                   column_values[5:],
                                                   min_val, test)
            for clf_name, clf in classifiers.items():
                print("Training " + clf_name)
                model, stats = fit_and_eval_ready_sets(clf, train_x, train_y, 
                                                test_x, test_y)
                if stats["+ Precision"] > 0 and stats["- Precision"] > 0:
                    results.append((clf_name+"_"+str(min_val), stats))
                prog.update(1)
        else:
            print("Not enough positives for " + str(min_val))
    prog.close()
    #%%
    results.sort(key=lambda x: x[1]["+ Precision"])
    results
    # %%
    train, test = 
    print(str([len(train),len(test)]))
    train_set = ()
    prepare_data_for_ml(,
                                           column_values[5:], 
                                           conf_level, train)
    print("For " + header[ss_index] + " >= "+ str(min_val) + ": "
          + str(len(neg)) + " negatives and " 
          + str(len(pos)) + " positives."
          + " Used " + str(used*100) + "% of available data.")
    
    
    
    test_train_masks = [get_norm_mask(column_values[i], 
                                      onto_filters[i], 
                                      n_bins=50, 
                                      subset_len=400000,
                                      train_size=0.6)
                        for i in tqdm(range(len(onto_filters)))]
    for i in range(len(test_train_masks)):
        print("For ontology " + header[i] + ": " 
              + str(len(test_train_masks[i][0])) + " train items and "
              + str(len(test_train_masks[i][1])) + " test items")
    
    
    
    
    train, test = test_train_masks[0]
    ss_index = 0
    conf_level = 0.4
    
    model, stats = fit_and_eval(column_values, ss_index, conf_level,
                     train, test, 500, 0.2)
    results = []
    
    bar = tqdm(total=len(ss_indexes)*len(mins_to_try))
    for ss_index in ss_indexes:
        for conf_level in mins_to_try:
            model, stats = get_conf_model(column_values, all_corrs_filter, 
                                          conf_level, ss_index, 
                                          subset_size, 0.7, steps)
            print("\t".join([str(x)
                for x in [conf_level, ss_index,                                            
                          stats["Score"], stats["FP"], stats["FN"]]]))
            results.append([conf_level,ss_col_names[ss_index],
                            model,stats])
            bar.update(1)
    
    res_file = output_dir+"/classification_results.tsv"
    
    with open(res_file, 'w') as stream:
        stream.write("Nível de Confiança\tOntologia\tScore"
                     +"\tFalsos Positivos (%)\tFalsos Negativos(%)\n")
        for conf_level, ss_name, model, stats in results:
            line = [str(conf_level), translator(ss_name), str(stats["Score"]), 
                    str(stats["FP"]*100), str(stats["FN"]*100)]
            stream.write("\t".join(line)+"\n")
            model_name = str(conf_level)+"_"+ss_name+".model"
            dump(model, output_dir+"/"+model_name)
    
    
    