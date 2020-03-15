#Module: bioinfo.py

from nexus.util import sliceString
from nexus.util import runCommand
import pandas as pd
#import hashlib
import requests
import os

rnacentral2rfam = {}

def readSeqsFromFasta(input_fasta):
    seqs = list()
    with open(input_fasta, "r") as fasta:
        contigName = ""
        seq = ""
        lines = fasta.readlines()
        for line in lines:
            if len(line) > 0:
                if line[0] == '>':
                    if len(seq) > 0:
                        seq = seq.replace("\n","")
                        seqs.append((contigName, seq))
                    contigName = line
                    seq = ""
                else:
                    if len(seq) > 0:
                        seq += "\n"
                    seq += line.rstrip("\n").lstrip("\n")
        seq = seq.replace("\n","")
        contigName = contigName.replace("\n","")
        seqs.append((contigName,seq))
    return seqs

def name_function(header):
    return header.rstrip("\n").lstrip(">")

def shortFastaHeader(seq):
    name = seq[0].split("|")[0]
    name = name.split(" ")[0]
    return (name,seq[1])

def header_to_id(seq):
    seq = seq.lstrip('>').rstrip('\n')
    name = seq.split("|")[0]
    name = name.split(" ")[0]
    return name

def seqListToDict(seqs, header_to_name=name_function):
    x = {}
    for pair in seqs:
        key = header_to_name(pair[0])
        val = pair[1]
        x[key] = val
    return x

def writeFastaSeqs(validSeqs, output_fasta, lineWidth = 60):
    seqs = list()
    for pair in validSeqs:
        header = pair[0]
        if header[-1] != '\n':
            header += "\n"
        if header[0] != '>':
            header = ">" + header
        seq = pair[1]
        slices = sliceString(seq.replace('\n', '').replace('\r', ''), lineWidth)
        contig = ""
        for i in range(0,len(slices)):
            contig += slices[i] + "\n"
        #contig += "\n"
        seqs.append((header,contig))

    with open(output_fasta, "w") as fasta:
        for pair in seqs:
            fasta.write(pair[0])
            fasta.write(pair[1])

def getType(header):
    pos = header.find("type:")
    if(pos <0):
        return "type:undefined"
    else:
        return header[pos:].split(" ")[0]

def writeSeqsWithUniqueHeaders(output_fasta, seqs, base_name):
    newSeqs = list()
    completeSeqs = list()
    partialSeqs = list()
    with open(output_fasta + ".names_table.tsv", "w") as translationTable:
        for i in range(0,len(seqs)):
            newName = ">"+base_name+"_"+str(i+1)
            if base_name == "CDS":
                seqType = getType(seqs[i][0])
                newName += "_" + seqType
                if seqType == "type:complete":
                    completeSeqs.append((newName, seqs[i][1]))
                else:
                    partialSeqs.append((newName, seqs[i][1]))
            newSeqs.append((newName, seqs[i][1]))
            translationTable.write(newName + "\t" + seqs[i][0])
    writeFastaSeqs(newSeqs, output_fasta)
    if len(completeSeqs) > 0:
        print(str((float(len(completeSeqs)) / float(len(newSeqs))) * 100) + "% of CDSs are complete")
        writeFastaSeqs(completeSeqs, output_fasta + ".complete")
    if len(partialSeqs) > 0:
        writeFastaSeqs(partialSeqs, output_fasta + ".partial")

def writeFastaWithUniqueHeaders(input_fasta, base_name="Contig"):
    seqs = readSeqsFromFasta(input_fasta)
    runCommand("cp " + input_fasta + " " + input_fasta + ".old_names")
    writeSeqsWithUniqueHeaders(input_fasta, seqs, base_name=base_name)

def filterSeqs(seqs, ids):
    lower_ids = set()
    for i in ids:
        lower_ids.add(i.lower())
    ids = lower_ids
    validSeqs = list()
    invalidSeqs = list()
    for header, seq in seqs:
        h = header.lower()
        valid = ((h.split()[0].lstrip(">").rstrip("\n") in ids)
            or (h.rstrip("\n") in ids)
            or (h in ids))
        if valid:
            validSeqs.append((header, seq))
        else:
            invalidSeqs.append((header, seq))
    return validSeqs, invalidSeqs

def getFastaHeaders(input_file):
    headers = list()
    #i = 0
    with open(input_file) as f:
        line = f.readline()
        while line != "":
            if line[0] == ">":
                #i += 1
                headers.append(line.lstrip(">").rstrip("\n"))
            line = f.readline()

    #print(
    return headers

def get_best_mapping(df):
    sorted = df.sort_values(["evalue", "bitscore"],
            ascending=[True, False])
    hits = []
    if len(sorted) > 0:
        for i in range(len(sorted)):
            hits.append(sorted.iloc[i])
            break
        return hits
    else:
        return None

def best_hit(df):
    sorted = df.sort_values(["evalue","length","pident"], 
            ascending=[True,False,False])
    if len(sorted) > 0:
        return sorted.iloc[0]
    else:
        return None

def get_strand(xstart,xend,ystart,yend):
    sig_x = "+"
    if xstart > xend:
        sig_x = "-"
    sig_y = "+"
    if ystart > yend:
        sig_y = "-"
    if(sig_x == sig_y):
        return "+"
    else:
        return "-"


def blast(query, db, output_dir, max_evalue = 0.0000000001, threads=8, blast_type="blastn", 
        db_id_sep_char=" ", source="db_name", remaining_fasta="auto_name"):
    import os
    import pandas as pd

    print("Blasting query to DB")
    db_name = os.path.basename(db).split(".")[0]
    if source == "db_name":
        source = db_name
    query_name = os.path.basename(query).split(".")[0]
    search_name = query_name+".to."+db_name
    output = output_dir + "/"+search_name+"_results.tsv"
    #cmd = " ".join([blast_type, "-db", db, "-query", query,
    #    "-out", output, "-outfmt", "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'", "-num_threads", str(threads), "-evalue", str(max_evalue)])
    cmd = " ".join([blast_type, "-d", db, "-i", query, "-e", str(max_evalue), "-a", 
                str(threads), "-outfmt 1", "-o", output])
    code = runCommand(cmd)
    return code == 0


def blast_annotate(query, db, output_dir, max_evalue = 0.0000000001, threads=8, blast_type="blastn", 
        db_id_sep_char=" ", source="db_name", remaining_fasta="auto_name", run_blast=True, alternative_outputdir=None):
    import os
    import pandas as pd

    print("Blasting query to DB")

    create_db = True
    if os.path.exists(db+".nhr"):
        create_db = False
    
    if create_db:
        cmd = " ".join(["makeblastdb -in " + db + " -dbtype nucl"])
        code = runCommand(cmd)
        if code != 0:
            print("Could not create database for given genome.")
            return False, ""

    db_name = os.path.basename(db).split(".")[0]
    if source == "db_name":
        source = db_name
    query_name = os.path.basename(query).split(".")[0]
    search_name = query_name+".to."+db_name
    output = output_dir + "/"+search_name+"_results.tsv"
    cmd = " ".join([blast_type, "-db", db, "-query", query,
        "-out", output, "-outfmt", 
        "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'", 
        "-num_threads", str(threads), "-evalue", str(max_evalue)])
    #cmd = " ".join([blast_type, "-d", db, "-i", query, "-e", str(max_evalue), "-a", 
    #            str(threads), "-outfmt 1", "-o", output])
    if run_blast:
        code = runCommand(cmd)
        if code != 0:
            return False, ""
    else:
        if alternative_outputdir != None:
            output = alternative_outputdir + "/"+search_name+"_results.tsv"

    print("Reading full seq names from DB")
    rnas = getFastaHeaders(db)
    full_names = dict()
    #x = 0
    for entry in rnas:
        parts = entry.split(db_id_sep_char)
        full_names[parts[0]] = " ".join(parts[1:])

    print("Parsing blast output")
    seqs = readSeqsFromFasta(query)
    seq_lens = {}
    for seq in seqs:
        seq_lens[seq[0].rstrip("\n").lstrip(">")] = len(seq[1])
    example_sequence_names = list(seq_lens.keys())[:5]
    print("Some sequence keys: " + str(example_sequence_names))
    blast_df = pd.read_csv(output, sep='\t', header=None, index_col=False,
            names=["qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    blast_df = blast_df.astype({"pident": 'float32', "length": 'int32', "mismatch": 'int32',
                "gapopen": 'int32', "qstart": 'int32', "qend": 'int32',
                "sstart": 'int32', "send": 'int32', "evalue": 'float64', "bitscore": 'float64'})
    print(str(blast_df.head()))
    print("Calculating coverage")
    blast_df["coverage"] = blast_df.apply(lambda row: row["length"] / seq_lens[row['qseqid']], axis=1)
    best_hits = dict()
    print(str(len(blast_df)) + " alignments")
    min_coverage = 0.99
    min_pid = 99
    print("Filtering blast results")
    blast_df = blast_df[blast_df["pident"] > min_pid]
    print(str(len(blast_df)) + " alignments filtered by pident")
    blast_df = blast_df[blast_df["coverage"] > min_coverage]
    print(str(len(blast_df)) + " alignments filtered by coverage")

    print("Choosing best hits")
    unique = 0
    for name, hits in blast_df.groupby(["qseqid"]):
        hit = get_best_mapping(hits)
        if hit != None:
            for i in range(len(hit)):
                best_hits[name+"."+str(i)] = hit[i]
            unique += 1
    print(str(unique) + " transcripts with genome mapping.")
    print(str(len(best_hits.keys())) + " total mappings.")

    print("Writing gff file about seqs identified")
    rows = []
    for name in best_hits:
        hit = best_hits[name]
        '''print(str(hit))
        print(type(str(hit)))
        print(hi)
        print(type(hit["sstart"]))
        print(type(hit["sstart"].item()))'''
        int_sstart = int(hit["sstart"])
        int_send = int(hit["send"])
        start = min(int_sstart, int_send)
        end = max(int_sstart, int_send)
        #full_name = "."
        #if hit["sseqid"] in full_names:
        #    full_name = full_names[hit["sseqid"]]
        row = {"seqname": hit["sseqid"], "source": source,
        "feature": "transcript", "start": str(start),
        "end":str(end), "score": ".",
        "strand": get_strand(int(hit["qstart"]),int(hit["qend"]),int(hit["sstart"]),int(hit["send"])),
        "frame": ".",
        "attribute":"ID="+name+";evalue="+str(hit["evalue"])
        +";coverage="+str(hit["coverage"])+";pident="+str(hit["pident"])}
        rows.append(row)
    gff = pd.DataFrame(rows, columns = ["seqname", "source",
        "feature", "start", "end", "score", "strand",
        "frame", "attribute"])
    gff_name = output_dir+"/"+search_name+"_found.gff"
    gff.to_csv(gff_name, sep="\t", index=False, header = False)
    print(str(len(seqs)) + " transcripts analyzed.")
    print(str(len(gff)) + " mappings detected and annotated on " + gff_name)

    '''print("Writing fasta files.")
    known = set([raw_name.split(".")[0] for raw_name in gff["seqname"].unique().tolist()])
    print("Some known sequences: " + str(list(known)[:5]))
    print("Some sequences: " + str([x[0] for x in seqs[:5]]))
    knownSeqs, unknownSeqs = filterSeqs(seqs, known)

    fasta_name = output_dir + "/" + remaining_fasta
    if remaining_fasta=="auto_name":
        fasta_name = output_dir+"/"+search_name+"_missing.fasta"
    writeFastaSeqs(unknownSeqs, fasta_name)
    writeFastaSeqs(knownSeqs, gff_name.rstrip("gff")+"fasta")
    print(str(len(gff)) + " detected and annotated on " + gff_name)
    print(str(len(unknownSeqs)) + " unknown seqs remaining on " + fasta_name)'''
    return True, gff_name

def read_plast_extended(path):
    df = pd.read_csv(path, sep="\t", header=None)
    col_names = ["query ID", "subject ID", "percent identities", "alignment length", "nb. misses", "nb. gaps", "query begin", "query end", "subject begin", "subject end", "e-value", "bit score", "getNbPositives", "query length", "query frame", "query translated", "query coverage", "nb. gaps in query", "subject length", "subject frame", "subject translated", "subject coverage", "nb. gaps in subject"]
    df.columns = col_names
    return df

def get_subject_aligned(df):
    subj_aligned = dict()
   #subj_aligned_inverse = dict()
    for index, row in df.iterrows():
        seq_name = row["subject ID"]
        start = int(row["subject begin"])-1
        end = int(row["subject end"])-1
        #if start > end:
        if not(seq_name in subj_aligned):
            subj_aligned[seq_name] = []
        if start > end:
            subj_aligned[seq_name].append((end,start))
        else:
            subj_aligned[seq_name].append((start, end))
        #else:
        #    if not(seq_name in subj_aligned_inverse):
        #        subj_aligned_inverse[seq_name] = []
        #    subj_aligned_inverse[seq_name].append((end,start))
    return subj_aligned

def has_intersection(range1, range2):
    values = [range1[0], range1[1], range2[0], range2[1]]
    values.sort()
    #print(str(values))
    unique1 = (values[0] in range1) and (values[1] in range1)
    unique2 = (values[0] in range2) and (values[1] in range2)
    return (not(unique1) 
            or not(unique2))

def extend_ranges(range1, range2):
    values = [range1[0], range1[1], range2[0], range2[1]]
    return (min(values), max(values))

def count_intersections(ranges):
    n = 0
    for i in range(len(ranges)):
        for j in range(len(ranges)):
            if i != j:
                if(has_intersection(ranges[i], ranges[j])):
                    n += 1
    return n

def cluster_ranges(ranges):
    if len(ranges) == 1:
        return ranges
    intersect_found = True
    len1 = len(ranges)
    while(intersect_found):
        intersect_found = False
        for i in range(len(ranges)):
            range1 = ranges[i]
            for j in range(len(ranges)):
                range2 = ranges[j]
                if i != j:
                    intersect_found = has_intersection(range1, range2)
                    break
            if(intersect_found):
                break
        if(intersect_found):
            ranges.remove(range1)
            ranges.remove(range2)
            ranges.append(extend_ranges(range1, range2))
    len2 = len(ranges)
    #if len1 != len2:
    #    print(str(len1 - len2))
    return ranges

def cluster_all_ranges(aligned):
    print("\nClustering " + str(len(aligned)) + " alignments:")
    #lens1 = [len(x) for x in aligned.values()]
    for key in aligned.keys():
        aligned[key] = cluster_ranges(aligned[key])
    #lens2 = [len(x) for x in aligned.values()]
    #diff_lens = [lens1[i] - lens2[i] for i in range(len(lens1))]
    #diff_lens.sort()
    #lens2.sort()
    #print(str(lens2))

def get_gff_attributes(attrs):
    try:
        parts = attrs.split(";")
    except:
        print("Could not split the following attributes:")
        print(str(attrs))
        raise ValueError("Attribute spliting error")
    values = {}
    for part in parts:
        ab = part.split("=")
        if len(ab) == 2:
            name = ab[0]
            val = ab[1].lstrip("'").rstrip("'")
            values[name] = val
    return values

def get_gff_attributes_str(attrs):
    x = ""
    for key in attrs.keys():
        value = attrs[key]
        x += key+"="+value+";"
    x = x.rstrip(";")
    return x

def load_rnacentral2rfam():
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with open(global_data+"/rnacentral2rfam.tsv",'r') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = "URS"+cells[0]
            rfam = "RF"+cells[1]
            rnacentral2rfam[rnacentral] = rfam

def get_rfam_from_rnacentral(id_, print_response=False):
    if len(rnacentral2rfam.keys()) == 0:
        load_rnacentral2rfam()

    if id_ in rnacentral2rfam:
        return rnacentral2rfam[id_]
    elif id_.split("_")[0] in rnacentral2rfam:
        return rnacentral2rfam[id_.split("_"[0])]
    else:
        return None