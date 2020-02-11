import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.util import runCommand, load_obj, save_obj, write_file, getFilesWith

def get_reference_rfam_ids(args, confs, tmpDir, stepDir):
    if "reference_gff" in args:
        ref = args["reference_gff"]
        annotation = pd.read_csv(ref, sep="\t", header=None, names=["seqname", 
            "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        print("Assigning RFAM id to genes without an RFAM id yet")
        log = tmpDir + "/log.txt"
        out_str = open(log, 'w')
        def update_attribute(row):
            attr_str = row["attribute"]
            attributes = None
            try:
                attributes = get_gff_attributes(attr_str)
            except:
                print("Row could not have attributes parsed:\n\t"+str(row))
            if attributes != None:
                if not("rfam" in attributes):
                    if row["source"] == "RNAcentral":
                        short_id = attributes["ID"].split(".")[0]
                        rfam_id = get_rfam_from_rnacentral(short_id)
                        if rfam_id != None:
                            if len(rfam_id) == 7:
                                attributes["rfam"] = rfam_id
                                out_str.write("Assigned " + rfam_id + " to " + short_id + "\n")
                                print("Assigned " + rfam_id + " to " + short_id + "\n")
                        if not("rfam" in attributes):
                            out_str.write("Could not assign rfam id to " + attributes["ID"]+"\n")
            return get_gff_attributes_str(attributes)
        annotation["attribute"] = annotation.apply(lambda row: update_attribute(row),axis=1)
        annotation.to_csv(tmpDir + "/reference.gff", sep="\t", index=False, header=False)
        out_str.close()
    return True
