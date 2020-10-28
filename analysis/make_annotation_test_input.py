import sys
from tqdm import tqdm
import pandas as pd

def compare_smaller_count(counts1, counts2):
    pass

def anno_segment_compare(item1, item2):
    seq_name1, start1, end1, counts1 = item1
    seq_name2, start2, end2, counts2 = item2

    if len(counts1.keys()) > len(counts2.keys()):
        return 1
    elif len(counts1.keys()) == len(counts2.keys()):
        mins1 = list(counts1.values())
        if len(mins1) == 0:
            return 0
        mins1.sort()
        mins2 = list(counts2.values())
        mins2.sort()

        #comparison = 0
        for i in range(min(len(mins1),len(mins2))):
            if mins1[i] > mins2[i]:
                return 1
            elif mins1[i] < mins2[i]:
                return -1
        return 0

        '''smaller_count1 = min(counts1.values()) if len(counts1.values()) > 0 else 0
        smaller_count2 = min(counts2.values()) if len(counts2.values()) > 0 else 0
        if smaller_count1 >= smaller_count2:
            return 1
        elif smaller_count1 == smaller_count2:
            return 0
        else:
            return -1'''
    else:
        return -1

def anno_segment_bigger(segment1, segment2):
    comparison = anno_segment_compare(segment1, segment2)
    if comparison == 1:
        return True
    else:
        return False

def print_segment(segment):
    seq_name, start, end, counts = segment
    print(seq_name+"["+str(start)+":"+str(end)+"] = "
        + ";".join([str(key)+"="+str(val) for key, val in counts.items()]))

def segment_to_seq(all_seqs, segment):
    seq_name, start, end, counts = segment
    if seq_name in all_seqs:
        seq = all_seqs[seq_name][start:end]
        return seq_name, seq
    else:
        return seq_name, None

def get_segment(annotations, seq_name, start, window):
    end = start+window
    segment_annos = annotations[(annotations['start'] >= start) & (annotations['end'] <= end)]
    counts = {source: len(annos) for source, annos in segment_annos.groupby(['source'])}
    return seq_name, start, end, counts

def get_gff_attributes(attrs):
    try:
        parts = attrs.split(";")
    except:
        print("Could not split the following attributes:")
        print(str(attrs))
        raise ValueError("Attribute spliting error")
    last_named_part = 0
    for i in range(1,len(parts)):
        if "=" in parts[i]:
            last_named_part = i
        else:
            parts[last_named_part] += ";"+parts[i]
            parts[i] = ""
    values = {}
    for part in parts:
        if len(part) > 0:
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

def read_gff(gff_path, feature = "transcript", get_types = True):
    print("Loading " + gff_path)
    df = pd.read_csv(gff_path, sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", 
        "score", "strand", "frame", "attribute"])
    #print("\tFiltering features")
    df = df[df["feature"] == feature]
    #print("\tUpdating attributes")
    if get_types:
        #df["attribute"] = df.apply(lambda row: row["attribute"], axis=1)
        df["rna_type"]  = df.apply(
                            lambda row: get_gff_attributes(row["attribute"])["type"],
                            axis = 1)
        df["ID"]  = df.apply(lambda row: get_gff_attributes(row["attribute"])["ID"],
                            axis = 1)
    chr_groups = {name: annos.sort_values(by=['start', 'end']) for name, annos in df.groupby(["seqname"])}
    return chr_groups

def get_ids_from_segments(segments, source_name, annotations):
    ids = set()
    for seg in segments:
        seq_name, start, end, counts = seg
        if source_name in counts:
            if seq_name in annotations:
                annos = annotations[seq_name]
                source_annos = annos[annos['source'] == source_name]
                selected_ann = source_annos[(source_annos['start'] >= start) & (source_annos['end'] <= end)]
                ids_selected = selected_ann['ID'].tolist()
                ids.update(ids_selected)
    return ids

def get_best_segment(seq_name, annotations, step, window):
    current_start = 0
    end = max(annotations['end'].tolist())
    max_iters = end/step
    #progress_bar = tqdm(total=max_iters)
    best_segment = get_segment(annotations, seq_name, current_start, window)
    current_start += step
    #progress_bar.update(1)
    while current_start+window <= end:
        segment = get_segment(annotations, seq_name, current_start, window)
        if anno_segment_bigger(segment, best_segment):
            best_segment = segment
            #print_segment(best_segment)
        current_start += step
        #progress_bar.update(1)
    #progress_bar.close()
    return best_segment


def sliceString(string, max_len):
    slices = []
    for x in range(0,int(len(string)/max_len)+1):
            slices.append(string[x*max_len:(x+1)*max_len])
    return slices

def seqListToDict(seqs):
    x = {}
    for pair in seqs:
        key = pair[0].rstrip("\n").lstrip(">")
        val = pair[1]
        x[key] = val
    return x

def shortFastaHeader(seq):
    name = seq[0].split("|")[0]
    name = name.split(" ")[0]
    return (name,seq[1])

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
        contig += "\n"
        seqs.append((header,contig))

    with open(output_fasta, "w") as fasta:
        for pair in seqs:
            fasta.write(pair[0])
            fasta.write(pair[1])

def readSeqsFromFasta(input_fasta):
    seqs = dict()
    with open(input_fasta, "r") as fasta:
        contigName = ""
        seq = ""
        lines = fasta.readlines()
        for line in lines:
            if len(line) > 0:
                if line[0] == '>':
                    if len(seq) > 0:
                        seq = seq.replace("\n","")
                        contigName = contigName.replace("\n","").replace(">","").split()[0]
                        #print("loaded " + contigName)
                        seqs[contigName] = seq
                    contigName = line
                    seq = ""
                else:
                    if len(seq) > 0:
                        seq += "\n"
                    seq += line.rstrip("\n").lstrip("\n")
        seq = seq.replace("\n","")
        contigName = contigName.replace("\n","").replace(">","").split()[0]
        #print("loaded " + contigName)
        seqs[contigName] = seq
    return seqs

def readSeqsFromFasta_list(input_fasta):
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
                        contigName = contigName.replace("\n","").replace(">","").split()[0]
                        seqs.append((contigName, seq))
                    contigName = line
                    seq = ""
                else:
                    if len(seq) > 0:
                        seq += "\n"
                    seq += line.rstrip("\n").lstrip("\n")
        seq = seq.replace("\n","")
        contigName = contigName.replace("\n","").replace(">","").split()[0]
        seqs.append((contigName,seq))
    return seqs

def filterSeqs(seqs, ids):
    lower_ids = set()
    for i in ids:
        lower_ids.add(i.lower())
    ids = lower_ids
    validSeqs = list()
    invalidSeqs = list()
    for header, seq in tqdm(seqs):
        h = header.lower()
        h = h.lstrip(">").rstrip("\n").split(" ")[0].split("|")[0]
        valid = False
        for a in ids:
            if a in h:
                valid = True
        if not valid:
            valid = ((h in ids)
                or (h.rstrip("\n") in ids)
                or (header.lower() in ids))
        if valid:
            validSeqs.append((header, seq))
        else:
            invalidSeqs.append((header, seq))
    return validSeqs, invalidSeqs

def filter_fasta(valid_ids, fasta_path, output_path):
    seqs = readSeqsFromFasta_list(fasta_path)
    valid, invalid = filterSeqs(seqs, valid_ids)
    writeFastaSeqs(valid, output_path)

def filter_gff(valid_ids, gff_path, output_path):
    with open(output_path, 'w') as out_stream:
        with open(gff_path, 'r') as in_stream:
            for line in in_stream:
                last_cell = line.rstrip("\n").split("\t")[-1]
                found = False
                for a in valid_ids:
                    if a in last_cell:
                        found = True
                if found:
                    out_stream.write(line)

annotation_path = sys.argv[1]
genome_path = sys.argv[2]
transcripts_fasta_path = sys.argv[3]
reference_fasta_path = sys.argv[4]
db_fasta_path = sys.argv[5]
reference_gff_path = sys.argv[6]
step = 100000
window = 7000000
contigs_to_get = 3
output = sys.argv[7]

gff_data = read_gff(annotation_path)

segments = []
import functools
for chr_name, annotations in tqdm(gff_data.items()):
    new_segment = get_best_segment(chr_name, annotations, step, window)
    segments.append(new_segment)
    segments.sort(key=functools.cmp_to_key(anno_segment_compare), reverse=True)

best_segments = segments[:contigs_to_get]
for segment in best_segments:
    print_segment(segment)

rnasamba_ids = get_ids_from_segments(best_segments, 
    "rnasamba", gff_data)
reference_mapping_ids = get_ids_from_segments(best_segments, 
    "reference_mapping", gff_data)
db_alignment_ids = get_ids_from_segments(best_segments, 
    "db_alignment", gff_data)
reference_ids = get_ids_from_segments(best_segments, 
    "reference", gff_data)

print(rnasamba_ids)
filter_fasta([".".join(x.split(".")[:-1]) for x in rnasamba_ids], 
            transcripts_fasta_path, output+"-transcripts.fasta")
print(reference_mapping_ids)
filter_fasta([".".join(x.split(".")[:-1]) for x in reference_mapping_ids], 
            reference_fasta_path, output+"-reference.fasta")
print(db_alignment_ids)
filter_fasta([".".join(x.split(".")[:-1]) for x in db_alignment_ids], 
            db_fasta_path, output+"-db.fasta")
print(reference_ids)
filter_gff([".".join(x.split(".")[:-1]) for x in reference_ids], 
            reference_gff_path, output+"-reference.gff")

print("Loading genome")
all_seqs = readSeqsFromFasta(genome_path)
seqs = {}
for segment in best_segments:
    seq_name, seq = segment_to_seq(all_seqs, segment)
    if seq != None:
        seqs[seq_name] = seq
        print_segment(segment)
    else:
        print(seq_name + " not found in genome.")
with open(output + "-genome.fasta", 'w') as stream:
    for seq_name, seq in seqs.items():
        stream.write(">"+seq_name+"\n"+seq+"\n")
