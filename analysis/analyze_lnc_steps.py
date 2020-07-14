import os
import sys
from tqdm import tqdm

aliases = {'non-coding RNA': "ncRNA", 
    'ncRNA': "ncRNA",
    'misc_RNA': "ncRNA",
    'partial mRNA': "mRNA",
}

original_fasta = sys.argv[1]
annotation_dir = sys.argv[2]
samba_only = sys.argv[3]
output = sys.argv[4]

def name_and_type(line):
    type_str = line.split(",")[-1].strip().rstrip("\n")
    if type_str in aliases:
        type_str = aliases[type_str]
    name_str = line.split(" ")[0].lstrip(">")
    return (name_str, type_str)

def compare_to_original(lnc_set, original):
    maybe_lnc = set([name for name, tp in lnc_set])
    
    correct_lnc = original["lncRNA"].intersection(maybe_lnc)
    lnc_notfound = original["lncRNA"] - maybe_lnc

    mRNA_found = original["mRNA"].intersection(maybe_lnc)
    other_ncRNA_found = original["ncRNA"].intersection(maybe_lnc)

    false_positives = (len(mRNA_found) + len(other_ncRNA_found)) / len(maybe_lnc)
    false_negatives = len(lnc_notfound) / len(original["lncRNA"])

    return (len(correct_lnc), len(other_ncRNA_found), len(mRNA_found), 
        round(false_positives*100, 3), round(false_negatives*100, 3))

def read_fasta(path):
    entries = []
    with open(path,'r') as stream:
        last_name = ""
        last_type = ""
        current_line = ""
        for line in stream:
            if line.startswith(">"):
                if last_name != "":
                    if last_type == "ncRNA" and len(current_line) >= 200:
                        last_type = "lncRNA"
                    entries.append((last_name, last_type))
                    current_line = ""
                last_name, last_type = name_and_type(line)
            else:
                if last_name != "":
                    current_line += line.rstrip("\n")
    return entries

def read_samba(path):
    entries = []
    with open(path,'r') as stream:
        stream.readline()
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            if cells[2] == "noncoding":
                entries.append((cells[0],[]))
    return entries

print("Reading original fasta")
original = read_fasta(original_fasta)
original_types = set([type_str for name, type_str in original])
type_sets = {}
for tp in ["ncRNA", "lncRNA", "mRNA"]:
    type_sets[tp] = set()
    for name, type_str in original:
        if ((type_str == tp) or (tp == "ncRNA" 
            and (type_str not in ["lncRNA", "mRNA"]))):
            type_sets[tp].add(name)
    print(tp + ": " + str(len(type_sets[tp])) + " RNAs")

lnc_files = {
            "1. Remoção de RNAs curtos": 
                annotation_dir+"/step_7-filter_small_sequences/long_transcripts.fasta",
            "2. Remover ORFs longos":
                annotation_dir+"/step_8-filter_long_orfs/no_orfs.fasta",
            "3. Testar com RNA Samba":
                annotation_dir+"/step_10-parse_coding_potential/no_coding_potential.fasta",
            "4. Filtrar mRNAs usando NR":
                annotation_dir+"/step_12-read_nr_alignment/not_protein.fasta",
            #"5. Mapear lncRNA no genoma":
            #    annotation_dir+"/step_14-lnc_alignment_parsing/lncRNA_in_genome.fasta"
            }

print("Reading step files")
selected_genes = [(step_name, read_fasta(file_path))
                    for step_name, file_path in tqdm(lnc_files.items())]

selected_genes.append(("Apenas RNA Samba", read_samba(samba_only)))

print("Analysing")
step_stats = [(step_name, compare_to_original(selected, type_sets))
                for step_name, selected in tqdm(selected_genes)]

rows = ["\t".join([step_name]+ [str(stat) for stat in stats])+"\n"
        for step_name, stats in step_stats]

rows = (["Passo\tlncRNA\tOutros ncRNA\tmRNA\t"
        +"Falsos Positivos (%)\tFalsos Negativos (%)\n"]
        + rows)
for row in rows:
    print(row.rstrip("\n"))

with open(output, 'w') as stream:
    stream.writelines(rows)

