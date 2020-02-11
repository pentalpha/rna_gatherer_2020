import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.util import runCommand, load_obj, save_obj, write_file, getFilesWith

def read_ids2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            go_ids = gos_str.split(";")
            gos_dict[rfam_id] = go_ids
    return gos_dict

def write_ids2go(filepath, gos_dict):
    with open(filepath, 'w') as stream:
        for key in gos_dict.keys():
            stream.write(key + "\t" + ";".join(gos_dict[key]) + "\n")

def write_transcriptome(args, confs, tmpDir, stepDir):
    print("Loading annotation")
    annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
                names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Loading genome")
    genome_dict = seqListToDict(readSeqsFromFasta(args['genome']))
    transcriptome = []
    print("Creating transcriptome file")
    for index, row in annotation.iterrows():
        #print(fasta_header)
        s = genome_dict[row["seqname"]] #cant find key PGUA01000001.1 #TODO
        new_header = get_gff_attributes(row["attribute"])["ID"]
        from_seq =int(row["start"])-1
        to_seq =int(row["end"])
        new_seq = s[from_seq:to_seq]
        transcriptome.append((new_header, new_seq))
    print("Writing transcriptome")
    writeFastaSeqs(transcriptome, tmpDir + "/transcriptome.fasta")
    return True

def make_ids2go(args, confs, tmpDir, stepDir):
    if "ids2go" in confs:
        ids2go_path = confs["ids2go"]
        if os.path.exists(ids2go_path):
            print("Loading ids2go associations")
            global_ids2go = read_ids2go(ids2go_path)
            print("Loading annotation")
            annotation = pd.read_csv(stepDir["remove_redundancies"] + "/annotation.gff", sep="\t", header=None,
                names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
            local_ids2go = {}
            print("Associating IDs to GO terms")
            ids = []
            for index, row in annotation.iterrows():
                attrs = get_gff_attributes(row["attribute"])
                ID = attrs["ID"]
                ids.append(ID)
                if "rfam" in attrs:
                    rfam_id = attrs["rfam"]
                    if rfam_id in global_ids2go:
                        go_list = global_ids2go[rfam_id]
                        local_ids2go[ID] = go_list
            write_ids2go(tmpDir + "/ids2go.tsv", local_ids2go)
            print("Writing population: " + str(len(ids)) + " ids")
            with open(tmpDir + "/ids.txt", 'w') as stream:
                for ID in ids:
                    stream.write(ID + "\n")
            return True
        else:
            print(ids2go_path + " does not exist.")
            return False
    else:
        print("No RFAM ids2go file specified in config.py")
        return False