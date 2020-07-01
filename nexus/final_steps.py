import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral, header_to_id
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

def number_of_rfam_ids(df):
    count = set()
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count.add(attrs["rfam"])
    return len(count)

def review_annotations(args, confs, tmpDir, stepDir):
    annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Reviewing annotation sources")

    lines = []
    for name, hits in annotation.groupby(["source"]):
        line = {
            "source": name,
            "total": hits.shape[0],
            "with_rfam_id": has_rfam_id(hits),
            "rfam_ids": number_of_rfam_ids(hits)
        }
        lines.append(line)
    lines.append({
        "source": "ALL",
        "total": annotation.shape[0],
        "with_rfam_id": has_rfam_id(annotation),
        "rfam_ids": number_of_rfam_ids(annotation)
    })

    review = pd.DataFrame(lines, columns=["source", "total", "with_rfam_id","rfam_ids"])
    review.to_csv(tmpDir + "/review.tsv", sep="\t", index=False)

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