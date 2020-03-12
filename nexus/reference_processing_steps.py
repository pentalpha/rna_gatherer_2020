import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges, header_to_id
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.util import runCommand, write_file, getFilesWith
import hashlib
import requests

def get_md5(sequence):
    """
    Calculate md5 for an RNA sequence
    """
    # RNAcentral stores DNA md5 hashes
    sequence = sequence.replace('U','T')
    # get the md5 digest
    m = hashlib.md5()
    m.update(sequence.encode('utf-8'))
    return m.hexdigest()

def get_rnacentral_id_by(argument, value):
    """
    Parse json output and return the RNAcentral id.
    """
    url = 'https://rnacentral.org/api/v1/rna'
    r = requests.get(url, params = {argument: value})
    data = r.json()
    if data['count'] > 0:
        #print("Got response: \n\t" + str(data['results'][0]))
        print(value + " matchs an ID in RNACentral: " + str(data['results'][0]['rnacentral_id']))
        return data['results'][0]['rnacentral_id']
    else:
        print(value + "Does not match an real ID in RNACentral")
        return None

def confirm_rnacentral_id(value):
    """
    Parse json output and return the RNAcentral id.
    """
    url = 'https://rnacentral.org/api/v1/rna/'+value.split("_")[0]
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    data = r.json()
    if 'rnacentral_id' in data:
        print(value + " matchs an real ID in RNACentral")
        return data['rnacentral_id']
    else:
        print(value + " does not match an real ID in RNACentral")
        return None

def retrieve_rnacentral_id(seq_id, seq):
    #print("Trying to retrieve " + seq_id)
    result = None
    if "URS" in seq_id:
        result = confirm_rnacentral_id(seq_id)
    if result == None:
        result = get_rnacentral_id_by("external_id",seq_id)
    if result == None:
        result = get_rnacentral_id_by("md5",get_md5(seq))
    return result

def get_rnacentral_ids(args, confs, tmpDir, stepDir):
    if "reference_fasta" in args:
        seqs = readSeqsFromFasta(args["reference_fasta"])
        id_to_seq = [(header_to_id(header),seq) for header,seq in seqs]
        rnacentral_to_seq = []
        tries = 3
        to_retrieve = set([x for x,y in id_to_seq])
        while tries > 0:
            print("Trying to retrieve " + str(len(to_retrieve)) + " ids.")
            failed = 0
            for seq_id, seq in id_to_seq:
                if seq_id in to_retrieve:
                    result = retrieve_rnacentral_id(seq_id, seq)
                    if result != None:
                        to_retrieve.remove(seq_id)
                        rnacentral_to_seq.append((result,seq))
                    else:
                        failed += 1
            print("\t" + str(failed) + " failed.")
            tries -= 1
            if failed == 0:
                tries = 0
        writeFastaSeqs(rnacentral_to_seq,tmpDir+"/rnacentral_seqs.fasta")
        unindentified = []
        for seq_id, seq in id_to_seq:
            if seq_id in to_retrieve:
                unindentified.append((seq_id,seq))
        writeFastaSeqs(unindentified,tmpDir+"/unknown.fasta")
        return True
    else:
        return True
    
def get_functional_reference(args, confs, tmpDir, stepDir):
    rnacentral_seqs_path = stepDir["get_rnacentral_ids"]+"/rnacentral_seqs.fasta"
    if os.path.exists(rnacentral_seqs_path):
        seqs = readSeqsFromFasta(rnacentral_seqs_path)
        ids = [header_to_id(header) for header,seq in seqs]
    else:
        return True

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
