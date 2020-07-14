import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import header_to_id
from nexus.util import *
import math
from scipy.stats.stats import pearsonr
from scipy.special import comb
import multiprocessing
import networkx
import obonet
import statsmodels.stats.multitest as multitest

def read_ids2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            go_ids = gos_str.split(";")
            gos_dict[rfam_id] = set(go_ids)
    return gos_dict

def read_rfam2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split()
            rfam_id = cols[0].split(":")[-1]
            if not rfam_id in gos_dict:
                gos_dict[rfam_id] = set()
            go_str = cols[-1]
            gos_dict[rfam_id].add(go_str)
    return gos_dict

def write_id2go(filepath, gos_dict):
    with open(filepath, 'w') as stream:
        for key, gos in gos_dict.items():
            for go in gos:
                stream.write(key + "\t" + go + "\n")

def has_rfam_alignment(df):
    count = 0
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "genbank" in attrs:
            if attrs["genbank"] != "None":
                count += 1
    return count

def has_rfam_id(df):
    count = 0
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count += 1
    return count

def number_of_genbank_ids(df):
    count = set()
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "genbank" in attrs:
            if attrs["genbank"] != "None":
                count.add(attrs["genbank"])
    return len(count)

def get_rfam_ids(df):
    count = set()
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count.add(attrs["rfam"])
    return count

def get_ncrna_type(attrs_str):
    attrs = get_gff_attributes(attrs_str)
    if "type" in attrs:
        return attrs["type"]
    else:
        return "other"

def group_rows(input_rows):
    higher_level = 100
    for row in input_rows:
        level = len(row[0].split(";"))
        if level < higher_level:
            higher_level = level

    head_rows = []
    for row in input_rows:
        level = len(row[0].split(";"))
        if level == higher_level:
            head_rows.append(row)
    if len(head_rows) == len(input_rows):
        input_rows.sort(key=lambda row: row[1], reverse=True)
        return [[row, []] for row in input_rows]
    else:
        row_groups = []
        for head_row in head_rows:
            head_row_id = ";".join(head_row[0].split(";")[0:higher_level])
            group = []
            ids = []
            for row in input_rows:
                if head_row_id in row[0] and head_row_id != row[0]:
                    group.append(row)
                    ids.append(row[0])
            #print("Group of " + head_row_id + ": " + str(ids))
            grouped_group = group_rows(group)
            row_groups.append([head_row, grouped_group])
        row_groups.sort(key=lambda row: row[0][1], reverse=True)
        
        return row_groups

def expand_groups(groups):
    rows = []
    for head_row, sub_rows in groups:
        rows.append(head_row)
        for sub_row, sub_group in sub_rows:
            rows.append(sub_row)
            if len(sub_group) > 0:
                expanded = expand_groups(sub_group)
                rows += expanded
    return rows

def sort_by_genes(input_rows):
    input_rows.sort(key=lambda row: row[0], reverse=False)
    grouped = group_rows(input_rows)
    new_rows = expand_groups(grouped)
    return new_rows


def review_annotations(args, confs, tmpDir, stepDir):
    annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Reviewing annotation sources")
    source_list = annotation["source"].unique().tolist()
    review_rows = []

    print("Retrieving rna_type")
    annotation["rna_type"] = annotation.apply(lambda row: get_ncrna_type(row["attribute"]),
                                                axis = 1)
    print("Counting RNA Types...")
    def count_rna_types(df):
        type_counts = {rna_type: len(type_annotation) 
                    for rna_type, type_annotation in df.groupby(["rna_type"])}
        expanded_type_counts = {}
        for rna_type in type_counts.keys():
            total = 0
            for another_type, counts in type_counts.items():
                if rna_type in another_type:
                    total += counts
            expanded_type_counts[rna_type] = total
        expanded_type_counts["All"] = len(df)
        return expanded_type_counts
        
    counts_by_sources = {"All": count_rna_types(annotation)}
    for source_name, source_annotation in tqdm(annotation.groupby(["source"])):
        counts_by_sources[source_name] = count_rna_types(source_annotation)
    
    print("Counting RFAM families")
    def count_rfam_families(df):
        rfam_ids = {rna_type: get_rfam_ids(type_annotation) 
                    for rna_type, type_annotation in df.groupby(["rna_type"])}
        expanded_rfam_counts = {}
        for rna_type in tqdm(rfam_ids.keys()):
            rfam_total = set()
            for another_type, ids in rfam_ids.items():
                if rna_type in another_type:
                    rfam_total.update(ids)
            expanded_rfam_counts[rna_type] = len(rfam_total)
        expanded_rfam_counts["All"] = len(get_rfam_ids(df))
        return expanded_rfam_counts
    rfam_counts = count_rfam_families(annotation)

    def make_row(type_name, type_annotation):
        new_row = [type_name]
        new_row.append(counts_by_sources["All"][type_name])
        new_row.append(rfam_counts[type_name])
        for source in source_list:
            n = 0
            source_counts = counts_by_sources[source]
            if type_name in source_counts:
                n = source_counts[type_name]
            new_row.append(n)
        return new_row

    print("Grouping by rna_type")
    for rna_type, type_annotation in tqdm(annotation.groupby(["rna_type"])):
        review_rows.append(make_row(rna_type, type_annotation))
    review_rows = sort_by_genes(review_rows)
    review_rows = [make_row("All", annotation)] + review_rows
    print("\n".join([str(row) for row in review_rows]))
    type_review_df = pd.DataFrame(review_rows, 
                            columns=["ncRNA Type", "Total", "Families"] + source_list)
    type_review_df.to_csv(tmpDir+"/type_review.tsv", sep="\t", index = False)
    return True

def write_transcriptome(args, confs, tmpDir, stepDir):
    print("Loading annotation")
    annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
                names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Loading genome: " + args['genome_link'])
    genome_dict = seqListToDict(readSeqsFromFasta(args['genome_link']), header_to_name = header_to_id)
    transcriptome = []
    print("Creating transcriptome file")
    for index, row in annotation.iterrows():
        #print(fasta_header)
        s = genome_dict[str(row["seqname"])] #cant find key PGUA01000001.1 #TODO
        new_header = get_gff_attributes(row["attribute"])["ID"]
        from_seq = int(row["start"])
        to_seq = int(row["end"])
        begin = min(from_seq,to_seq)-1
        up_to = max(from_seq,to_seq)
        new_seq = s[begin:up_to]
        transcriptome.append((new_header, new_seq))
    print("Writing transcriptome")
    writeFastaSeqs(transcriptome, tmpDir + "/transcriptome.fasta")
    return True

def make_id2go(args, confs, tmpDir, stepDir):
    id2go_path = confs["rfam2go"]
    if os.path.exists(id2go_path):
        print("Loading ids2go associations")
        global_ids2go = read_rfam2go(id2go_path)
        print("Loading annotation")
        annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
            names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

        print("Associating IDs to GO terms")
        local_ids2go = {}
        ids = []
        for index, row in annotation.iterrows():
            attrs = get_gff_attributes(row["attribute"])
            ID = attrs["ID"]
            ids.append(ID)
            if "rfam" in attrs:
                rfam_id = attrs["rfam"]
                if rfam_id in global_ids2go:
                    go_list = global_ids2go[rfam_id]
                    local_ids2go[ID] = set(go_list)
        api_annotation = stepDir["get_info"] + "/retrieved_functions.id2go"
        if os.path.exists(api_annotation):
            with open(api_annotation, 'r') as stream:
                for line in stream:
                    cells = line.rstrip("\n").split("\t")
                    if not cells[0] in local_ids2go:
                        local_ids2go[cells[0]] = set()
                    local_ids2go[cells[0]].add(cells[1])
        
        write_id2go(tmpDir + "/id2go.tsv", local_ids2go)
        print("Writing population: " + str(len(ids)) + " ids")
        with open(tmpDir + "/ids.txt", 'w') as stream:
            for ID in ids:
                stream.write(ID + "\n")
        return True
    else:
        print(id2go_path + " does not exist.")
        return False
