import sys
import numpy as np
from tqdm import tqdm
import os
import numpy as np
import math

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

def replace_with_dict(in_file, out_file, names_dict, header, n_terms = None): 
    total_lines = 0
    wrote = 0
    pairs_done = set()
    print_terms = "terms_a" in header
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
                        #pair_id, terms_a, terms_b = str(names_dict[pair])
                        pair_id = names_dict[pair]
                        if not pair_id in pairs_done:
                            pairs_done.add(pair_id)
                            first_part = [str(a) for a in [pair_id]]
                            if print_terms == True:
                                first_part.append(str(n_terms[cells[0]]))
                                first_part.append(str(n_terms[cells[1]]))
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

if __name__ == "__main__":
    files_ss = sys.argv[1].split(",")
    files_corr = sys.argv[2].split(",")
    annotation = sys.argv[3]
    output_dir = sys.argv[4]
    min_terms = None
    max_terms = None
    if len(sys.argv) > 5:
        min_terms = int(sys.argv[5])
        max_terms = int(sys.argv[6])
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    print("Reading gene names")
    gene_names = get_gene_names(files_ss)

    print("Reading gene annotation")
    onto_trans = {"P": "bp", "F": "mf", "C": "cc"}
    n_terms = {}
    terms_by_gene = {}
    with open(annotation,'r') as in_stream:
        for line in in_stream:
            if not(line.startswith("!") or line.startswith("#")):
                cells = line.rstrip("\n").split("\t")
                if len(cells) >= 5:
                    gene_name = cells[1]
                    go_term = cells[4]
                    onto = onto_trans[cells[8]]
                    if not gene_name in terms_by_gene:
                        terms_by_gene[gene_name] = set()
                    terms_by_gene[gene_name].add(go_term)
                else:
                    print("Invalid line:" + str(cells))
        for gene_name, term_set in terms_by_gene.items():
            #for key in onto_trans.values():
            if gene_name in gene_names:
                n_terms[gene_name] = len(term_set)

    n_terms_vec = np.array(list(n_terms.values()))
    mean = n_terms_vec.mean()
    std = n_terms_vec.std()
    if max_terms == None:
        min_terms = math.floor(mean)
        max_terms = math.ceil(mean)
    to_log = ("Mean terms per gene: " + str(mean)+ "\n")
    to_log += ("Standard deviation: " + str(std) + "\n")
    to_log += ("Minimum terms: " + str(min_terms) + "\n")
    to_log += ("Maximum terms: " + str(max_terms) + "\n")

    normally_annotated_genes = []
    all_genes = list(n_terms.keys())
    for gene_name, n_term in n_terms.items():
        if n_term >= min_terms and n_term <= max_terms:
            normally_annotated_genes.append(gene_name)

    to_log += ("Normally annotated genes: " + str(len(normally_annotated_genes)) 
            + " of " + str(len(n_terms_vec)) + "\n")
    print(to_log)
    print(str(len(normally_annotated_genes)*len(all_genes)) + " maximum pairs")
    '''useful_genes = {"bp": set(), "mf": set(), "cc": set()}
    for onto, len_dict in n_terms.items():
        for gene_name, n_terms in len_dict.items():
            if n_terms >= 3:
                useful_genes[onto].add(gene_name)
    all_useful_genes = list(set.union(*[useful_genes["bp"], useful_genes["mf"], useful_genes["cc"]]))'''
    
    print("Creating ID dictionaries")
    pair_to_id = {}
    id_ = 0
    for i in tqdm(range(len(normally_annotated_genes))):
        name_a = normally_annotated_genes[i]
        for j in range(len(all_genes)):
            name_b = all_genes[j]
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
    def replace_in_files(files, header):
        n_files = []
        for f in tqdm(files):
            new_path = output_dir + "/" + get_filename(f)
            if not os.path.exists(new_path):
                if len(header) == 1:
                    new_header = header + [get_name(f)]
                    replace_with_dict(f, new_path, pair_to_id, new_header)
                else:
                    replace_with_dict(f, new_path, pair_to_id, header, n_terms=n_terms)
            else:
                print("Skiping " + f)
            n_files.append(new_path)
        return n_files
    print("\tFor correlations files...")
    new_files_corr = replace_in_files(files_corr, 
                                ["pair_id"])
    print("\tFor ss files...")
    new_files_ss = replace_in_files(files_ss, 
                                ["pair_id", "terms_a", "terms_b", 
                                "mf", "bp", "cc"])
    print("Loading all correlations into memory")
    correlations = {}
    metric_names = []
    expected_size = 0
    for corr_file in tqdm(new_files_corr):
        with open(corr_file, 'r') as stream:
            metric_name = stream.readline().rstrip("\n").split("\t")[-1]
            metric_names.append(metric_name)
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
    for ss_file in new_files_ss:
        big_path = ss_file.replace(".tsv", "-with_correlations.tsv")
        with open(ss_file, 'r') as in_stream:
            with open(big_path, 'w') as out_stream:
                expected_size = 8 + len(metric_names)
                line = in_stream.readline()
                base_header = line.rstrip("\n")
                header = base_header + "\tavg\tmax\t" + "\t".join(metric_names)
                out_stream.write(header+"\n")
                line = in_stream.readline()
                progress_bar = tqdm(total=id_)
                while line:
                    cells = line.rstrip("\n").split("\t")
                    id_ = int(cells[0])
                    if id_ in correlations:
                        '''int_terms = np.fromstring(cells[1]+"\t"+cells[2],
                                                        sep="\t")'''
                        float_ss_values = np.fromstring("\t".join(cells[3:6]),
                                                        sep="\t")
                        max_ss = np.nanmax(float_ss_values)
                        avg_ss = np.nanmean(float_ss_values)
                        cells += [str(avg_ss), str(max_ss)] 
                        cells += [str(x) for x in correlations[id_]]
                        nan_to_add = "\tnan"*(expected_size - len(cells))
                        out_stream.write("\t".join(cells) + nan_to_add + "\n")
                    line = in_stream.readline()
                    progress_bar.update(1)

