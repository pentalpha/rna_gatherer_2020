import sys
#import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import dask
import dask.dataframe as dd

def get_filename(full_path):
    last = full_path.split("/")[-1]
    return last

def add_names_to_set(file_path, names_set):
    with open(file_path, 'r') as stream:
        for line in stream:
            cells = line.split("\t")
            names_set.add(cells[0])
            names_set.add(cells[1])
    return names_set

def replace_with_dict(in_file, out_file, names_dict, header):
    with open(in_file,'r') as in_stream:
        with open(out_file,'w') as out_stream:
            out_stream.write("\t".join(header)+"\n")
            for line in in_stream:
                cells = line.rstrip("\n").replace("\tNone", "\tNaN").split("\t")
                cells[0] = names_dict[cells[0]]
                cells[1] = names_dict[cells[1]]
                out_stream.write("\t".join(cells)+"\n")

if __name__ == "__main__":
    files_ss = sys.argv[1].split(",")
    files_corr = sys.argv[2].split(",")
    all_files = files_ss + files_corr
    output_dir = sys.argv[3]

    '''if len(corr_file_paths) == 1:
        sys.argv[2] == "None"
        corr_file_paths = []'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    dict_file = output_dir + "/gene_ids.tsv"
    names_dict = {}
    if not os.path.exists(dict_file):
        print("Reading files to creage gene ID dict.")
        names_set = set()
        for f in tqdm(all_files):
            names_set = add_names_to_set(f, names_set)
        i = 0
        for name in names_set:
            names_dict[name] = str(i + 1)
            i += 1
        with open(dict_file,'w') as out_stream:
            for name, id in names_dict.items():
                out_stream.write(name+"\t"+str(id)+"\n")
    else:
        print("Reading names dict")
        with open(dict_file,'r') as in_stream:
            for line in in_stream:
                cells = line.rstrip("\n").split("\t")
                names_dict[cells[0]] = cells[1]
    
    print("Creating versions without gene names, but IDs instead")
    def replace_in_files(files, header):
        n_files = []
        for f in tqdm(files):
            new_path = output_dir + "/" + get_filename(f)
            if not os.path.exists(new_path):
                replace_with_dict(f, new_path, names_dict, header)
            else:
                print("Skiping " + f)
            n_files.append(new_path)
        return n_files
    print("\tFor correlations files...")
    new_files_corr = replace_in_files(files_corr, 
                                ["gene_a", "gene_b", "corr"])
    print("\tFor ss files...")
    new_files_ss = replace_in_files(files_ss, 
                                ["gene_a", "gene_b", "mf", "bp", "cc"])
    
    print("Converting new dataframes to parquet format")
    def to_parquet(files, types):
        for f in tqdm(files):
            parquet_file = f.replace(f.split(".")[-1], "parquet")
            if not os.path.exists(parquet_file):
                df = dd.read_csv(f, sep="\t", dtype=types)
                df.to_parquet(parquet_file,
                            engine = "pyarrow")
    print("\tFor correlations files...")
    to_parquet(new_files_corr, {"gene_a": int, "gene_b": int, "corr": float})
    print("\tFor ss files...")
    to_parquet(new_files_ss, {"gene_a": int, "gene_b": int, 
                                "mf": float, "bp": float, "cc": float})