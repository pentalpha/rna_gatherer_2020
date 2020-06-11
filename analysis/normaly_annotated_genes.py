import numpy as np
import sys

'''
Usage:
    python normaly_annotated_genes.py <input_go_associations>.tsv <input_valid_counts>.tsv <output>
'''

input_go = sys.argv[1]
valid_names = sys.argv[2]
output_list = sys.argv[3]

def read_id2go(filepath, valid_ids):
    gos_dict = {}
    go_terms = set()
    all_ids = set()
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                gene_id = cols[0]
                go_str = cols[1]
                aspect = cols[2]
                if gene_id in valid_ids:
                    if not gene_id in gos_dict:
                        gos_dict[gene_id] = set()
                    gos_dict[gene_id].add(go_str)
                    go_terms.add(go_str)
                all_ids.add(gene_id)
    return gos_dict, len(go_terms), len(all_ids)

def read_seq_names(filepath):
    names = set()
    first = True
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            if not first:
                name = raw_line.split("\t")[0].rstrip("\n").lstrip(">")
                names.add(name)
            first = False
    return names

print("Reading valid IDs")
valid_ids = read_seq_names(valid_names)
print("Reading annotation")
id2gos, n_go_terms, n_ids = read_id2go(input_go, valid_ids)

print("Converting data to arrays")
id_go_list = [(_id, gos) for _id, gos in id2gos.items()]
n_gos_list = np.array([len(gos) for _id, gos in id_go_list])
print("Calculating median annotation length")
normal_anno_len = np.median(n_gos_list)
description = ("Median annotation length: " + str(normal_anno_len)
                + "\nValid Annotated genes: " + str(len(id_go_list))
                + "\nGO Terms: " + str(n_go_terms)
                + "\nOriginal genes: " + str(n_ids) + "\n")
print(description)
print("Sorting annotations")
id_go_list.sort(key=lambda x: ((len(x[1]) - normal_anno_len) if len(x[1]) > normal_anno_len else (normal_anno_len - len(x[1]))))
print("Writing results")
id_annolen = "\n".join([_id + "\t" + str(len(gos)) for _id, gos in id_go_list])
output = open(output_list,'w')
output.write(id_annolen)
output.close()

output = open(output_list + ".meta.txt", 'w')
output.write(description)
output.close()