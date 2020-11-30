
from gatherer.bioinfo import get_gff_attributes, get_gff_attributes_str
from gatherer.bioinfo import header_to_id
from gatherer.util import *
from gatherer.rna_type import *
import sys
import pandas as pd

def matching_coords(a, b, max_outside):
    diff_end = max(0, a[1] - b[1])
    diff_start = max(0, b[0] - a[0])
    return diff_end + diff_start <= max_outside

def sorted_match(my_coords, ref_coords):
    max_out = [(my_coord[1] - my_coord[0])*0.05 for my_coord in my_coords]

    starts = [(my_coords[i][0], True, i) for i in range(len(my_coords))]
    starts += [(ref_coords[i][0], False, i) for i in range(len(ref_coords))]

    ends = [(my_coords[i][1], True, i) for i in range(len(my_coords))]
    ends += [(ref_coords[i][1], False, i) for i in range(len(ref_coords))]

    edges = [(True, (coord, not_ref, index)) for coord, not_ref, index in starts]
    edges += [(False, (coord, not_ref, index)) for coord, not_ref, index in ends]

    edges.sort(key=lambda x: x[1][0])

    possible_matchs = {}
    keeping_track_of = set()
    maybe_ending = set()

    for is_start, values in edges:
        
        coord, not_ref, index = values
        if not_ref:
            if is_start:
                keeping_track_of.add(index)
                possible_matchs[index] = set()
            elif index in keeping_track_of:
                maybe_ending.add(index)
                #keeping_track_of.remove(index)
        else:
            remove_by_distance = []
            for fading_index in maybe_ending:
                if (coord - my_coords[fading_index][1]) > max_out[fading_index]:
                    remove_by_distance.append(fading_index)
            for to_remove in remove_by_distance:
                maybe_ending.remove(to_remove)
                keeping_track_of.remove(to_remove)
            
            for my_coord in keeping_track_of:
                possible_matchs[my_coord].add(index)

    matched = set()
    for my_index, ref_indexes in possible_matchs.items():
        start, end = my_coords[my_index]
        max_outside = max_out[my_index]
        for ref_index in ref_indexes:
            other_start, other_end = ref_coords[ref_index]
            result = matching_coords((start, end), 
                                    (other_start, other_end), 
                                    max_outside)
            if result:
                matched.add(my_index)
                break
    
    return matched

def update_attrs(attr_str):
    attrs = get_gff_attributes(attr_str)
    if "family" in attrs:
        attrs["rfam"] = attrs["family"]
    
    if not "rfam" in attrs:
        new_rfam = get_rfam_from_rnacentral(attrs["ID"])
        if new_rfam != None:
            attrs["rfam"] = new_rfam

    if "rfam" in attrs:
        new_type = get_rna_type(attrs["rfam"])
        attrs["type"] = new_type
    else:
        if not "type" in attrs:
            attrs["type"] = "other"
        #Replace 'misc_rna' with 'other'
        if attrs["type"] in unkown_rna:
            attrs["type"] = "other"
        attrs["type"] = get_full_type(attrs["type"])

    return get_gff_attributes_str(attrs)

def read_gff(gff_path, feature = "transcript", get_types = False):
    print("Loading " + gff_path)
    df = pd.read_csv(gff_path, sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", 
        "score", "strand", "frame", "attribute"])
    #print("\tFiltering features")
    df = df[df["feature"] == feature]
    #print("\tUpdating attributes")
    if get_types:
        df["attribute"] = df.apply(lambda row: row["attribute"], axis=1)
        df["rna_type"]  = df.apply(
                            lambda row: get_gff_attributes(row["attribute"])["type"],
                            axis = 1)
    return df

def get_coords(df):
    coords_by_chr = {}
    chr_groups = df.groupby(["seqname"])
    for chr_name, group in chr_groups:
        coords = []
        for index, row in group.iterrows():
            x = int(row["start"])
            y = int(row["end"])
            coords.append((min(x,y), max(x,y)))
        coords.sort(key=lambda x: x[0])
        coords_by_chr[str(chr_name)] = coords
    return coords_by_chr

annotation_path = sys.argv[1]
reference_path = sys.argv[2]
output_path = sys.argv[3]

annotation_df = read_gff(annotation_path, get_types=True)
reference_transcripts = read_gff(reference_path)
reference_exons = read_gff(reference_path, feature="noncoding_exon")

print("Extracting coordinates")
transcripts_by_chr = get_coords(reference_transcripts)
exons_by_chr = get_coords(reference_exons)

print("Matching coordinates for each type")
def match_groups(my_coords, ref_coords):
    count = 0
    matched = set()
    for chr_name, my_chr in my_coords.items():
        #print("\tChr: " + chr_name)
        if chr_name in ref_coords:
            ref_chr = ref_coords[chr_name]
            local_matchs = sorted_match(my_chr, ref_chr)
            for local_match in local_matchs:
                matched.add(str(chr_name)+"_"+str(local_match))
    return matched

def all_names(my_coords):
    matched = set()
    for chr_name, my_chr in my_coords:
        for i in range(len(my_chr)):
            matched.add(str(chr_name)+"_"+str(i))

raw_matchs_by_type = {}
n_types = len(annotation_df["rna_type"].unique().tolist())
progress_bar = tqdm(total=n_types)
for rna_type, rna_group in annotation_df.groupby(["rna_type"]):
    #print("Matching for " + rna_type)
    my_coords = get_coords(rna_group)
    matched_transcripts = match_groups(my_coords, transcripts_by_chr)
    matched_exons = match_groups(my_coords, exons_by_chr)
    all_matched = set.union(matched_exons, matched_transcripts)
    raw_matchs_by_type[rna_type] = {
        "matchs": len(all_matched),
        "total": len(rna_group)}
    progress_bar.update(1)
progress_bar.close()

print("RNA Type\tExon Match Rate\tTranscript Match Rate")
for rna_type, counts_dict in raw_matchs_by_type.items():
    #exons = counts_dict["exons"]
    transcripts = counts_dict["matchs"]
    total = counts_dict["total"]
    #print(rna_type+"\t"+str(exons/total)+"\t"+str(transcripts/total))
    print(rna_type+"\t"+str(total)+"\t"+str(transcripts/total))

expanded_matchs_by_type = {}
for rna_type in raw_matchs_by_type.keys():
    expanded_matchs = {"total": 0, "matchs": 0}
    for other_type, stats in raw_matchs_by_type.items():
        if rna_type in other_type:
            expanded_matchs["total"] += stats["total"]
            expanded_matchs["matchs"] += stats["matchs"]
    expanded_matchs_by_type[rna_type] = expanded_matchs

raw_rows = []
for rna_type, stats in expanded_matchs_by_type.items():
    total = stats["total"]
    transcripts = stats["matchs"]
    perc = (transcripts/total)* 100.0
    raw_rows.append([rna_type, total, transcripts, perc])

raw_rows.sort(key=lambda x: x[0], reverse=True)

rows = "\n".join(["\t".join([str(x) for x in row]) 
                            for row in raw_rows])+"\n"
with open(output_path, 'w') as stream:
    stream.write("Tipo de ncRNA\tTotal\tCompatíveis com a Referência\tCompatíveis(%)\n")
    stream.write(rows)
print(rows)