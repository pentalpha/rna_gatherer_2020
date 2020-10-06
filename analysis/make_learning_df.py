import sys
import numpy as np
from tqdm import tqdm
import os
import numpy as np
import math
import pandas as pd

def get_filename(full_path):
    last = full_path.split("/")[-1]
    return last

def get_name(full_path):
    return get_filename(full_path).split(".")[0]

def add_names_to_set(file_path, names_set):
    with open(file_path, 'r') as stream:
        for line in stream:
            cells = line.split("\t")
            names_set.add(cells[0])
            names_set.add(cells[1])
    return names_set

def replace_with_dict(in_file, out_file, names_dict, header): 
    total_lines = 0
    wrote = 0
    pairs_done = set()
    with open(in_file,'r') as in_stream:
        with open(out_file,'w') as out_stream:
            out_stream.write("\t".join(header)+"\n")
            for line in in_stream:
                total_lines += 1
                cells = line.rstrip("\n").split("\t")
                
                not_null = False
                for cell in cells[2:]:
                    if cell != "None" and cells != "NaN":
                        not_null = True
                        break

                if not_null:
                    pair = cells[0]+cells[1]
                    if pair in names_dict:
                        pair_id = names_dict[pair]
                        if not pair_id in pairs_done:
                            pairs_done.add(pair_id)
                            first_part = [str(a) for a in [pair_id]]
                            out_stream.write(
                                "\t".join(first_part + cells[2:]).replace("\tNone", "\tNaN")
                                +"\n")
                            wrote += 1
    print(str(wrote/total_lines) + " of lines mantained for " + in_file)

def get_gene_names(files):
    gene_names = set()
    for f in files:
        with open(f, 'r') as stream:
            stream.readline()
            for line in stream:
                cells = line.rstrip("\n").split("\t")
                if not ("None\tNone\tNone" in line):
                    gene_names.add(cells[0])
                    gene_names.add(cells[1])
    return gene_names

def read_interactions(interactions_file, all_gene_names):
    interactions = set()
    origins = set()
    targets = set()
    with open(interactions_file, 'r') as stream:
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            if cells[0] in all_gene_names and cells[1] in all_gene_names:
                origins.add(cells[0])
                targets.add(cells[1])
                interactions.add((cells[0]+"\t"+cells[1]))
    
    return interactions, origins, targets

if __name__ == "__main__":
    files_corr = sys.argv[1].split(",")
    counts_file = sys.argv[2]
    interactions_file = sys.argv[3]
    output_dir = sys.argv[4]
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    print("Reading counts")
    counts = pd.read_csv(counts_file, sep='\t')
    all_gene_names = set(counts[counts.columns[0]].tolist())
    print(str(len(all_gene_names)))
    print(str(list(all_gene_names)[:4]))
    interactions, origins, targets = read_interactions(interactions_file, all_gene_names)
    print(str(len(interactions)),str(len(origins)),str(len(targets)))
    
    print("Creating interactions DF")
    interactions_df_path = output_dir+"/interaction.tsv"
    if not os.path.exists(interactions_df_path):
        with open(interactions_df_path, 'w') as stream:
            stream.write("gene_a\tgene_b\tinteract\n")
            for a in tqdm(origins):
                for b in targets:
                    new_line = (a+"\t"+b+"\t")
                    if a+"\t"+b in interactions:
                        new_line += "True\n"
                    else:
                        new_line += "False\n"
                    stream.write(new_line)

    print("Creating ID dictionaries")
    pair_to_id = {}
    id_ = 0
    to_log = ""
    origins = list(origins)
    targets = list(targets)
    for i in tqdm(range(len(origins))):
        name_a = origins[i]
        for j in range(len(targets)):
            name_b = targets[j]
            if name_a != name_b:
                pair_name_a = name_a+name_b
                pair_name_b = name_b+name_a
                if not(pair_name_a in pair_to_id) and not(pair_name_b in pair_to_id):
                    pair_to_id[pair_name_a] = id_
                    pair_to_id[pair_name_b] = id_
                    id_ += 1
    
    to_log += (str(id_) + " valid pairs.\n")
    to_log += ("Example pairs: \n")
    to_log += (str(list(pair_to_id.items())[0:10])+"\n")
    print(to_log)
    log = open(output_dir+"/parsing_log.txt", 'w')
    log.write(to_log)
    log.close()

    print("Creating versions without gene names, but IDs instead")
    def replace_in_files(files, header, new_paths=[]):
        n_files = []
        i = 0
        for f in tqdm(files):
            new_path = output_dir + "/" + get_filename(f)
            if len(new_paths) > 0:
                new_path = new_paths[i]
            if not os.path.exists(new_path):
                if len(header) == 1:
                    new_header = header + [get_name(f)]
                    replace_with_dict(f, new_path, pair_to_id, new_header)
                else:
                    replace_with_dict(f, new_path, pair_to_id, header)
            else:
                print("Skiping " + f)
            n_files.append(new_path)
            i += 1
        return n_files
    print("\tFor correlations files...")
    new_files_corr = replace_in_files(files_corr, 
                                ["pair_id"])
    print(str(new_files_corr))
    print("\tFor other files...")
    other_files = replace_in_files([interactions_df_path], 
                                ["pair_id", "interact"],
                                new_paths=[output_dir + "/parsed_interactions.tsv"])
    print("Loading all correlations into memory")
    correlations = {}
    metric_names = []
    expected_size = 0
    for corr_file in tqdm(new_files_corr):
        with open(corr_file, 'r') as stream:
            metric_name = stream.readline().rstrip("\n").split("\t")[-1]
            metric_names.append(metric_name)
            print(str(metric_names))
            line = stream.readline()
            #progress_bar = tqdm(total=id_)
            while line:
                cells = line.rstrip("\n").split("\t")
                id_ = int(cells[0])
                if not id_ in correlations:
                    correlations[id_] = []
                if cells[1] != "NaN":
                    nan_to_add = expected_size - len(correlations[id_])
                    for i in range(nan_to_add):
                        correlations[id_].append(np.nan)
                    correlations[id_].append(float(cells[1]))
                line = stream.readline()
                #progress_bar.update(1)
        expected_size += 1

    print("Making one big tsv")
    for other_file in other_files:
        big_path = other_file.replace(".tsv", "-with_correlations.tsv")
        with open(other_file, 'r') as in_stream:
            with open(big_path, 'w') as out_stream:
                expected_size = 2 + len(metric_names)
                line = in_stream.readline()
                base_header = line.rstrip("\n")
                header = base_header + "\t" + "\t".join(metric_names)
                out_stream.write(header+"\n")
                line = in_stream.readline()
                progress_bar = tqdm(total=id_)
                while line:
                    cells = line.rstrip("\n").split("\t")
                    id_ = int(cells[0])
                    if id_ in correlations:
                        cells += [str(x) for x in correlations[id_]]
                        nan_to_add = "\tnan"*(expected_size - len(cells))
                        out_stream.write("\t".join(cells) + nan_to_add + "\n")
                    line = in_stream.readline()
                    progress_bar.update(1)