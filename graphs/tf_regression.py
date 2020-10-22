# -*- coding: utf-8 -*-
import sys
import numpy as np
from tqdm import tqdm
import os
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import multiprocessing
from joblib import dump
import tensorflow as tf
from datetime import datetime
from packaging import version
from tensorflow import keras

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

def filter_column(x, min_cor):
    f = set()
    for i in range(len(x)):
        if (not np.isfinite(x[i]) or np.isnan(x[i])) or x[i] < min_cor:
            f.add(i)
    return f

def get_norm_mask(column, invalid_rows, 
                  n_bins=100, subset_len=10000,
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
    #print("Subset bins\n"+str([len(a) for a in bins]))
    print("Creating mask from bins")
    mask = []
    for b in bins:
        for value in b:
            mask.append(value)
    train_set = set(np.random.choice(mask, 
                            int(len(mask)*train_size), 
                            replace=False))
    test_set = set()
    for a in mask:
        if not a in train_set:
            test_set.add(a)
    #print("Total values selected: " + str(len(mask)))
    return train_set, test_set

def ss_to_class(ss_column, min_value):
    return [val >= min_value for val in ss_column]

def prepare_data_for_ml_cols(label_vals, feature_vals, feature_names, 
                        valid_rows):
    '''Filter a set of columns, so that all values are valid'''
    y = []
    filtered_x = [[] for a in feature_names]
    for i in valid_rows:
        y.append(label_vals[i])
        for j in range(len(feature_vals)):
            filtered_x[j].append(feature_vals[j][i])
    feature_dict = {feature_names[i]: np.array(filtered_x[i])
                    for i in range(len(feature_names))}
    return feature_dict, y

def prepare_data_for_ml(label_vals, feature_vals, 
                        valid_rows):
    '''Filter a set of columns, so that all values are valid'''
    y = []
    filtered_x = []
    for i in valid_rows:
        y.append(label_vals[i])
        for j in range(len(feature_vals)):
            filtered_x.append([feature_vals[j][i]
                               for j in range(len(feature_vals))])
    return filtered_x, y

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

def from_dataset(ds):
   return ds.make_one_shot_iterator()

if __name__ == "__main__":
    '''file_ss = sys.argv[1]
    output_dir = sys.argv[2]
    subset_size = int(sys.argv[3])'''
    file_ss = "zfin_resnik.tsv.2m"
    output_dir = "regress_2m"
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
    all_corrs_filter = set.union(*[filter_column(x, 0.0) 
                                          for x in tqdm(column_values[5:])])
    ss_col = 0
    mf_filter = filter_column(column_values[ss_col], 0.0).union(all_corrs_filter)
    train_mask, test_mask = get_norm_mask(column_values[ss_col], 
                                    mf_filter, 
                                    n_bins=50, 
                                    subset_len=100000, 
                                    train_size=0.95)
    features_train, labels_train = prepare_data_for_ml(
        column_values[ss_col],
        column_values[5:],
        train_mask)
    features_test, labels_test = prepare_data_for_ml(
        column_values[ss_col],
        column_values[5:],
        test_mask)
    
    from tensorflow.python.keras.models import Sequential
    from tensorflow.python.keras.layers import Flatten
    from tensorflow.python.keras.layers import Dense
    from tensorflow.python.keras.layers import Activation
    from tensorflow.python.keras.layers import Input
    
    model = Sequential([
        Input(shape=(6,)),
        Dense(4, activation='sigmoid'), # dense layer 1
        Dense(2, activation='softmax') # dense layer 2
        #Dense(10, activation='softmax'),  # dense layer 3
    ])
    
    model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
    
    model.fit(features_train, labels_train, epochs=10, 
              batch_size=10, validation_split = 0.2)
    
    
    
    dataset_train = tf.data.Dataset.from_tensor_slices(
        (features_train, labels_train))
    
    '''features_train = tf.constant(X_train)
    features_test = tf.constant(features_test)
    labels_train = tf.constant(y_train)
    labels_test = tf.constant(labels_test)'''
    
    def train_fn():
        return features_train, labels_train
    def test_fn():
        return features_test, labels_test
        
    
    '''dataset_train = tf.data.Dataset.from_tensor_slices(
        (features_train, labels_train))
    dataset_test = tf.data.Dataset.from_tensor_slices(
        (features_test, labels_test))'''
    feature_cols = [tf.feature_column.numeric_column(k)
                    for k in features_train.keys()]
    
    def run_for_arch(arch):
        log_steps = 10
        estimator = tf.estimator.DNNLinearCombinedRegressor(
            linear_feature_columns=[feature_cols],
            config=tf.estimator.RunConfig(
                log_step_count_steps=log_steps))
        
        estimator.train(input_fn=lambda: train_fn(), 
                        steps=200)
    run_for_arch([128,64,32])
    
    small_test = list(test_mask)[:10]
    small_features_test, small_labels_test = prepare_data_for_ml(
        column_values[ss_col],
        column_values[5:],
        header[5:], 
        small_test)
    def small_test_fn():
        return small_features_test, small_labels_test
    
    ev = estimator.evaluate(input_fn=lambda: test_fn())
    def predict_in():
        return features_test
    y = estimator.predict(
        input_fn=lambda: predict_in())
    
    

    
    