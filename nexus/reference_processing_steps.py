import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges, header_to_id
from nexus.bioinfo import read_plast_extended
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.bioinfo import minimap_annotation
from nexus.util import runCommand, write_file, getFilesWith, chunks, read_to_list
import hashlib
import requests
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
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
        #print(value + " matchs an ID in RNACentral: " + str(data['results'][0]['rnacentral_id']))
        return data['results'][0]['rnacentral_id']
    else:
        #print(value + "Does not match an real ID in RNACentral")
        return None

def get_rnacentral_json(value):
    url = 'https://rnacentral.org/api/v1/rna/'+value
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    data = r.json()
    if 'rnacentral_id' in data:
        return data
    else:
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
        #print(value + " matchs an real ID in RNACentral")
        return data['rnacentral_id']
    else:
        print(value + " does not match an real ID in RNACentral")
        return None

def retrieve_rnacentral_id(seq_id, seq):
    #print("Trying to retrieve " + seq_id)
    result = None
    try:
        if "URS" in seq_id:
            result = confirm_rnacentral_id(seq_id)
        if result == None:
            result = get_rnacentral_id_by("external_id",seq_id)
        if result == None:
            result = get_rnacentral_id_by("md5",get_md5(seq))
    except Exception:
        print("Could not retrieve information for " + seq_id)
    return result, seq_id, seq

def retrieve_quickgo_annotations(chunk, api_url,taxon_id):
    result = []
    gene_lines = []
    gene_ids = ",".join(chunk)
    for i in range(len(chunk)):
        if "URS" in chunk[i]:
            if not "_"+taxon_id in chunk[i]:
                chunk[i] = chunk[i]+"_"+taxon_id
    gene_ids = ",".join(chunk)
    for aspect in ["biological_process","molecular_function","cellular_component"]:
        requestURL = (api_url+"?selectedFields=geneProductId&selectedFields=goId&geneProductId="
                        +gene_ids+"&taxonId="+taxon_id+"&aspect="+aspect)
        tries = 0
        response_lines = []
        max_tries = 5
        while tries < max_tries:
            if tries > 0:
                #print("Trying again, waiting " + str(tries*tries))
                time.sleep(tries*tries)
            #print("Requesting:\n\t"+requestURL)
            try:
                response = requests.get(requestURL, headers={"Accept":"text/tsv"})
                if response.ok:
                    #print("Got okay response")
                    tries = max_tries
                    text = response.text
                    lines = text.split("\n")[1:]
                    response_lines = [line.split("\t")+[aspect] for line in lines]
                    tries = max_tries
                else:
                    #print("Response not okay.")
                    #print(response.text)
                    try:
                        json_response = json.loads(response.text)
                        msgs = json_response["messages"]
                        #print("Analyzing error msgs")
                        for msg in msgs:
                            if "The 'Gene Product ID' parameter contains in" in msg:
                                invalid_ids = msg.split(": ")[-1].split(", ")
                                #print("Invalid ids: "+ str(invalid_ids))
                                #invalid_ids.append(invalid_id)
                                for id_ in invalid_ids:
                                    #print("Replacing " + id_)
                                    requestURL = requestURL.replace(id_,"")
                                    requestURL = requestURL.replace(",,",",")
                                #print("New request URL:\n\t" + requestURL)
                        tries += 1
                    except exception:
                        tries += 1
            except:
                tries += 1
        gene_lines += response_lines
    added_line = False
    for line in gene_lines:
        valid = True
        valid = valid and len(line) > 1
        if valid:
            valid = valid and line[0] != ""
        if valid:
            result += [line[1:]]
            added_line = True
    '''if added_line:
        #print(str(len(result)) + " annotations retrieved")
    else:
        print("No annotations retrieved for " + gene_ids)'''
    return result

def get_ids_from_annotation(annotation):
    ids = []
    def update_attribute(row):
        attr_str = row["attribute"]
        attributes = None
        try:
            attributes = get_gff_attributes(attr_str)
        except:
            print("Row could not have attributes parsed:\n\t"+str(row))
        if attributes != None:
            if "ID" in attributes:
                short_id = attributes["ID"]
                if "URS" in short_id:
                    short_id = attributes["ID"].split("_")[0]
                ids.append(short_id)
        return 0
    annotation.apply(lambda row: update_attribute(row),axis=1)
    return ids

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
            to_retrieves = [list(chunk) for chunk in chunks(list(to_retrieve), 100)]
            for chunk in tqdm(to_retrieves):
                processes = []
                with ThreadPoolExecutor(max_workers=20) as executor:
                    for seq_id, seq in id_to_seq:
                        if seq_id in chunk:
                            processes.append(executor.submit(retrieve_rnacentral_id, seq_id, seq))
                            
                    for task in as_completed(processes):
                        if task.result() != None:
                            result, seq_id, seq = task.result()
                            #print(result)
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
        if len(unindentified) > 0:
            writeFastaSeqs(unindentified,tmpDir+"/unknown.fasta")
        return True
    else:
        return True
    
def get_functional_reference(args, confs, tmpDir, stepDir):
    rnacentral_seqs_path = stepDir["get_rnacentral_ids"]+"/rnacentral_seqs.fasta"
    unknown_seqs_path = stepDir["get_rnacentral_ids"]+"/unknown.fasta"
    if os.path.exists(rnacentral_seqs_path) or os.path.exists(unknown_seqs_path):
        ids = []
        if os.path.exists(rnacentral_seqs_path):
            seqs = readSeqsFromFasta(rnacentral_seqs_path)
            ids = [header_to_id(header) for header,seq in seqs]
        if os.path.exists(unknown_seqs_path):
            unknown_seqs = readSeqsFromFasta(unknown_seqs_path)
            ids += [header_to_id(header) for header,seq in unknown_seqs]
        annotation_ids = []
        if "reference_gff" in args:
            ref = args["reference_gff"]
            annotation = pd.read_csv(ref, sep="\t", header=None, names=["seqname", 
                "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
            annotation_ids = get_ids_from_annotation(annotation)
            ids += annotation_ids
        id_set = set(ids)
        id_chunks = chunks(list(id_set),50)
        
        unmaped = ""
        maped = set(annotation_ids)
        for id_ in id_set:
            if not id_ in maped:
                unmaped += id_ + "\n"
        if len(unmaped) > 0:
            with open(tmpDir+"/unmaped_ids.txt",'w') as stream:
                stream.write(unmaped)

        results = set()
        annotations_retrieved = False
        if not "taxon_id" in args:
            print("Cannot retrieve functional annotations without taxon id.")
            return True
        for chunk in tqdm(list(id_chunks)):
            annotations = retrieve_quickgo_annotations(chunk,confs["quickgo_api"],args["taxon_id"])
            for id_,go,aspect in annotations:
                results.add((id_,go,aspect))
                annotations_retrieved = True
        if annotations_retrieved:
            with open(tmpDir+"/quickgo_annotations.tsv",'w') as stream:
                for id_,go,aspect in results:
                    stream.write(id_+"\t"+go+"\t"+aspect+"\n")
        return True
    else:
        return True

def map_to_genome(args, confs, tmpDir, stepDir):
    unmaped_path = stepDir["get_functional_reference"]+"/unmaped_ids.txt"
    if os.path.exists(unmaped_path):
        unmaped_ids = set(read_to_list(unmaped_path))
        rnacentral_seqs_path = stepDir["get_rnacentral_ids"]+"/rnacentral_seqs.fasta"
        unknown_seqs_path = stepDir["get_rnacentral_ids"]+"/unknown.fasta"
        possible_seqs = []
        if os.path.exists(rnacentral_seqs_path):
            possible_seqs = readSeqsFromFasta(rnacentral_seqs_path)
        if os.path.exists(unknown_seqs_path):
            possible_seqs += readSeqsFromFasta(unknown_seqs_path)
        unmaped_seqs = []
        for header,seq in possible_seqs:
            id_ = header_to_id(header)
            if "URS" in id_:
                    id_ = id_.split("_")[0] 
            if id_ in unmaped_ids:
                unmaped_seqs.append((id_,seq))
            else:
                print(id_ + " already has a mapping.")
        print("Found " + str(len(unmaped_seqs)) + " RNAs without genome mapping.")
        unmaped_path = tmpDir+"/unmaped.fasta"
        writeFastaSeqs(unmaped_seqs,unmaped_path)
        genome_alignment = tmpDir + "/reference_mapped.paf"
        cmd = " ".join([confs["minimap2"],
            "-x splice:hq -uf -t", str(confs["threads"]),
            args["genome_index"], unmaped_path, 
            ">", genome_alignment])
        code = runCommand(cmd)
        if code != 0:
            return False
        
        output_gff = tmpDir + "/reference_mapped.gff"
        annotated_fasta = tmpDir + "/reference_mapped.fasta"
        mapped_fasta = minimap_annotation(genome_alignment, unmaped_path, output_gff,
                        annotated_fasta, source="reference_mapping", mol_type=None)
        return True
    return True

def get_reference_rfam_ids(args, confs, tmpDir, stepDir):
    new_mappings_path = stepDir["map_to_genome"] + "/reference_mapped.gff"
    if "reference_gff" in args or os.path.exists(new_mappings_path):
        all_mappings = tmpDir + "/no_rfam.gff"
        if "reference_gff" in args:
            print("Assigning RFAM ids to reference mappings.")
            runCommand("cat " + args["reference_gff"] + " >> " + all_mappings)
        if os.path.exists(new_mappings_path):
            print("Assigning RFAM ids to new reference mappings.")
            runCommand("cat " + new_mappings_path + " >> " + all_mappings)
        
        ref = all_mappings
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
                    short_id = attributes["ID"]
                    if "URS" in short_id:
                        short_id = attributes["ID"].split(".")[0].split("_")[0]
                        rfam_id = get_rfam_from_rnacentral(short_id)
                        if rfam_id != None:
                            if len(rfam_id) == 7:
                                attributes["rfam"] = rfam_id
                                out_str.write("Assigned " + rfam_id + " to " + short_id + "\n")
                                #print("Assigned " + rfam_id + " to " + short_id + "\n")
                        if not("rfam" in attributes):
                            out_str.write("Could not assign rfam id to " + attributes["ID"]+"\n")
                            #print("Could not assign rfam id to " + attributes["ID"]+"\n")
            return get_gff_attributes_str(attributes)
        annotation["attribute"] = annotation.apply(lambda row: update_attribute(row),axis=1)
        annotation.to_csv(tmpDir + "/reference.gff", sep="\t", index=False, header=False)
        out_str.close()
    return True

def get_rnacentral_info(args, confs, tmpDir, stepDir):
    ref = stepDir["get_reference_rfam_ids"] + "/reference.gff"
    if os.path.exists(ref):
        annotation = pd.read_csv(ref, sep="\t", header=None, names=["seqname", 
            "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        print("Looking for info over RNACentral API")
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
                    short_id = attributes["ID"]
                    if "URS" in short_id:
                        short_id = attributes["ID"].split(".")[0].split("_")[0]
                        rnacentral_json = get_rnacentral_json(short_id)
                        if rnacentral_json != None:
                            if not "description" in attributes and "description" in rnacentral_json:
                                attributes["description"] = rnacentral_json["description"]
                            if not "type" in attributes and "rna_type" in rnacentral_json:
                                attributes["type"] = rnacentral_json["rna_type"]
                            print("Retrieved info for " + short_id)
            return get_gff_attributes_str(attributes)
        annotation["attribute"] = annotation.apply(lambda row: update_attribute(row),axis=1)
        annotation.to_csv(tmpDir + "/reference.gff", sep="\t", index=False, header=False)
        out_str.close()
    return True
