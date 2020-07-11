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
    with open(in_file,'r') as in_stream:
        with open(out_file,'w') as out_stream:
            out_stream.write("\t".join(header)+"\n")
            for line in in_stream:
                total_lines += 1
                not_null = False
                cells = line.rstrip("\n").split("\t")
                for cell in cells[2:]:
                    if cell != "None" and cells != "NaN":
                        not_null = True
                        break
                if not_null:
                    pair = cells[0]+cells[1]
                    if pair in names_dict:
                        cells[0] = str(names_dict[pair])
                        out_stream.write(
                            "\t".join(
                                [str(names_dict[pair])] + cells[2:]
                            ).replace("\tNone", "\tNaN")+"\n")
                        wrote += 1
    print(str(wrote/total_lines) + " of lines mantained for " + in_file)

if __name__ == "__main__":
    files_ss = sys.argv[1].split(",")
    files_corr = sys.argv[2].split(",")
    annotation = sys.argv[3]
    output_dir = sys.argv[4]
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print("Reading gene annotation")
    n_terms = {}
    with open(annotation,'r') as in_stream:
        terms_by_gene = {}
        for line in in_stream:
            if not(line.startswith("!") or line.startswith("#")):
                cells = line.rstrip("\n").split("\t")
                if len(cells) >= 5:
                    gene_name = cells[1]
                    go_term = cells[4]
                    if not gene_name in terms_by_gene:
                        terms_by_gene[gene_name] = set()
                    terms_by_gene[gene_name].add(go_term)
                else:
                    print("Invalid line:" + str(cells))
        for gene_name, term_set in terms_by_gene.items():
            n_terms[gene_name] = len(term_set)
    
    useful_genes = []
    for gene_name, n_terms in n_terms.items():
        if n_terms >= 3:
            useful_genes.append(gene_name)

    print(str(len(useful_genes)) + " genes with enough annotation")

    print("Creating ID dictionaries")
    pair_to_id = {}
    id_ = 0
    for i in tqdm(range(len(useful_genes))):
        for j in range(len(useful_genes)):
            pair_name_a = useful_genes[i]+useful_genes[j]
            pair_name_b = useful_genes[j]+useful_genes[i]
            if not(pair_name_a in pair_to_id) and not(pair_name_b in pair_to_id):
                pair_to_id[pair_name_a] = id_
                id_ += 1
    
    print(str(list(pair_to_id.items())[0:10]))

    print("Creating versions without gene names, but IDs instead")
    def replace_in_files(files, header):
        n_files = []
        for f in tqdm(files):
            new_path = output_dir + "/" + get_filename(f)
            if not os.path.exists(new_path):
                if len(header) == 1:
                    new_header = header + [get_name(f)]
                    replace_with_dict(f, new_path, pair_to_id, new_header)
                else:
                    replace_with_dict(f, new_path, pair_to_id, header)
            else:
                print("Skiping " + f)
            n_files.append(new_path)
        return n_files
    print("\tFor correlations files...")
    new_files_corr = replace_in_files(files_corr, 
                                ["pair_id"])
    print("\tFor ss files...")
    new_files_ss = replace_in_files(files_ss, 
                                ["pair_id", "mf", "bp", "cc"])
    '''
    print("Converting new dataframes to parquet format")
    def to_parquet(files, types = None):
        created = []
        for f in tqdm(files):
            dtypes = {"pair_id": int, get_name(f): float}
            if types != None:
                dtypes = types
            parquet_file = f.replace(f.split(".")[-1], "parquet")
            if not os.path.exists(parquet_file):
                df = dd.read_csv(f, sep="\t", dtype=dtypes).set_index("pair_id")
                print(str(df.head()))
                df.to_parquet(parquet_file,
                            engine = "pyarrow")
            created.append(parquet_file)
        return created
    
    print("\tFor correlations files...")
    corr_parquets = to_parquet(new_files_corr)
    print("\tFor ss files...")
    ss_parquets = to_parquet(new_files_ss, types = {"pair_id": int, 
                            "mf": float, "bp": float, "cc": float})

    print("Joining parquets")
    for ss_file in ss_parquets:
        ss_parquet = dd.read_parquet(ss_file, engine="pyarrow")
        temp_ss = ss_file.replace(".parquet", "-with_correlations.parquet")
        for corr_file in tqdm(corr_parquets):
            print("Adding correlations from " + corr_file)
            corr_parquet = dd.read_parquet(corr_file, engine="pyarrow")
            ss_parquet = ss_parquet.merge(corr_parquet, 
                                        left_index=True, right_index=True)
        ss_parquet.to_parquet(temp_ss)
    '''
    print("Loading all correlations into memory")
    correlations = {}
    metric_names = []
    expected_size = 0
    for corr_file in tqdm(new_files_corr):
        with open(corr_file, 'r') as stream:
            metric_name = stream.readline().rstrip("\n").split("\t")[-1]
            metric_names.append(metric_name)
            line = stream.readline()
            progress_bar = tqdm(total=len(useful_genes)*len(useful_genes))
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
                progress_bar.update(1)
        expected_size += 1

    print("Making one big tsv")
    for ss_file in new_files_ss:
        big_path = ss_file.replace(".tsv", "-with_correlations.tsv")
        with open(ss_file, 'r') as in_stream:
            with open(big_path, 'w') as out_stream:
                expected_size = 4 + len(metric_names)
                line = in_stream.readline()
                base_header = line.rstrip("\n")
                header = base_header + "\t" + "\t".join(metric_names)
                out_stream.write(header+"\n")
                line = in_stream.readline()
                progress_bar = tqdm(total=len(useful_genes)*len(useful_genes))
                while line:
                    cells = line.rstrip("\n").split("\t")
                    id_ = int(cells[0])
                    if id_ in correlations:
                        cells += [str(x) for x in correlations[id_]]
                        nan_to_add = "\nan"*(expected_size - len(cells))
                        out_stream.write("\t".join(cells) + nan_to_add + "\n")
                    line = in_stream.readline()
                    progress_bar.update(1)

