import sys
import numpy as np
from tqdm import tqdm
import os
import numpy as np

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

def replace_with_dict(in_file, out_file, names_dict, header, valid_genes=None):
    replace_ss_values = valid_genes != None
    total_lines = 0
    wrote = 0
    with open(in_file,'r') as in_stream:
        with open(out_file,'w') as out_stream:
            out_stream.write("\t".join(header)+"\n")
            for line in in_stream:
                total_lines += 1
                cells = line.rstrip("\n").split("\t")
                if replace_ss_values:
                    mf_valid = ((cells[0] in valid_genes["mf"]) 
                                and (cells[1] in valid_genes["mf"]))
                    bp_valid = ((cells[0] in valid_genes["bp"]) 
                                and (cells[1] in valid_genes["bp"]))
                    cc_valid = ((cells[0] in valid_genes["cc"]) 
                                and (cells[1] in valid_genes["cc"]))
                    cells[2] = cells[2] if mf_valid else "NaN"
                    cells[3] = cells[3] if bp_valid else "NaN"
                    cells[4] = cells[4] if cc_valid else "NaN"
                
                not_null = False
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
    onto_trans = {"P": "bp", "F": "mf", "C": "cc"}
    n_terms = {key: {} for key in onto_trans.values()}
    with open(annotation,'r') as in_stream:
        terms_by_gene = {}
        for line in in_stream:
            if not(line.startswith("!") or line.startswith("#")):
                cells = line.rstrip("\n").split("\t")
                if len(cells) >= 5:
                    gene_name = cells[1]
                    go_term = cells[4]
                    onto = onto_trans[cells[8]]
                    if not gene_name in terms_by_gene:
                        terms_by_gene[gene_name] = {key: set() for key in onto_trans.values()}
                    terms_by_gene[gene_name][onto].add(go_term)
                else:
                    print("Invalid line:" + str(cells))
        for gene_name, term_sets in terms_by_gene.items():
            for key in onto_trans.values():
                n_terms[key][gene_name] = len(term_sets[key])
    
    useful_genes = {"bp": set(), "mf": set(), "cc": set()}
    for onto, len_dict in n_terms.items():
        for gene_name, n_terms in len_dict.items():
            if n_terms >= 3:
                useful_genes[onto].add(gene_name)
    all_useful_genes = list(set.union(*[useful_genes["bp"], useful_genes["mf"], useful_genes["cc"]]))
    
    log = open(output_dir+"/parsing_log.txt", 'w')
    log.write(str(len(all_useful_genes)) + " genes with enough annotation\n")
    log.write("\t" + str(len(useful_genes["bp"])) + " with enough BP terms\n")
    log.write("\t" + str(len(useful_genes["mf"])) + " with enough MF terms\n")
    log.write("\t" + str(len(useful_genes["cc"])) + " with enough CC terms\n")

    print("Creating ID dictionaries")
    pair_to_id = {}
    id_ = 0
    for i in tqdm(range(len(all_useful_genes))):
        for j in range(len(all_useful_genes)):
            name_a = all_useful_genes[i]
            name_b = all_useful_genes[j]
            valid = False
            for onto_name, valid_ids in useful_genes.items():
                if name_a in valid_ids and name_b in valid_ids:
                    valid = True
                    break
            if valid:
                pair_name_a = name_a+name_b
                pair_name_b = name_b+name_a
                if not(pair_name_a in pair_to_id) and not(pair_name_b in pair_to_id):
                    pair_to_id[pair_name_a] = id_
                    id_ += 1
    
    log.write(str(len(pair_to_id.keys())) + " valid pairs.\n")
    log.write("Example pairs: \n")
    log.write(str(list(pair_to_id.items())[0:10])+"\n")
    log.close()

    print("Creating versions without gene names, but IDs instead")
    def replace_in_files(files, header, useful_genes):
        n_files = []
        for f in tqdm(files):
            new_path = output_dir + "/" + get_filename(f)
            if not os.path.exists(new_path):
                if len(header) == 1:
                    new_header = header + [get_name(f)]
                    replace_with_dict(f, new_path, pair_to_id, new_header)
                else:
                    replace_with_dict(f, new_path, pair_to_id, header, valid_genes=useful_genes)
            else:
                print("Skiping " + f)
            n_files.append(new_path)
        return n_files
    print("\tFor correlations files...")
    new_files_corr = replace_in_files(files_corr, 
                                ["pair_id"],
                                None)
    print("\tFor ss files...")
    new_files_ss = replace_in_files(files_ss, 
                                ["pair_id", "mf", "bp", "cc"],
                                useful_genes)
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
                expected_size = 6 + len(metric_names)
                line = in_stream.readline()
                base_header = line.rstrip("\n")
                header = base_header + "\tavg\tmax\t" + "\t".join(metric_names)
                out_stream.write(header+"\n")
                line = in_stream.readline()
                progress_bar = tqdm(total=len(useful_genes)*len(useful_genes))
                while line:
                    cells = line.rstrip("\n").split("\t")
                    id_ = int(cells[0])
                    if id_ in correlations:
                        float_ss_values = np.fromstring("\t".join(cells[1:4]),
                                                        sep="\t")
                        max_ss = np.nanmax(float_ss_values)
                        avg_ss = np.nanmean(float_ss_values)
                        cells += [str(avg_ss), str(max_ss)] 
                        cells += [str(x) for x in correlations[id_]]
                        nan_to_add = "\tnan"*(expected_size - len(cells))
                        out_stream.write("\t".join(cells) + nan_to_add + "\n")
                    line = in_stream.readline()
                    progress_bar.update(1)

