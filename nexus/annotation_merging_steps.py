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
from nexus.util import runCommand, write_file, getFilesWith

def filter_non_transcripts(gff_path, gff_output_path):
    ref_annotation = pd.read_csv(gff_path, sep="\t", header=None, names=["seqname", 
            "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    ref_annotation = ref_annotation[ref_annotation["feature"] == "transcript"]
    ref_annotation_file = gff_output_path
    ref_annotation.to_csv(ref_annotation_file, sep="\t", index=False, header=False)

def run_gffcompare(args, confs, tmpDir, stepDir):
    rfam_mappings = stepDir["parse_infernal"] + "/rfam_annotation_genome.gff"
    tRNA_mappings = stepDir["parse_trna"] + "/tRNAs.gff"
    novel_mappings = tmpDir + "/novel_mappings.gff"
    filter_non_transcripts(rfam_mappings, novel_mappings)
    gffs = [novel_mappings,tRNA_mappings]

    lnc_mappings = stepDir["lnc_alignment_parsing"] + "/lncRNA_annotation.gff"
    if os.path.exists(lnc_mappings):
        gffs = gffs + [lnc_mappings]

    if "reference_gff" in args:
        ref = stepDir["get_reference_rfam_ids"] + "/reference.gff"
        ref_output = tmpDir + "/reference.gff"
        filter_non_transcripts(ref, ref_output)
        gffs = [os.path.abspath(ref)] + gffs

    cmd = " ".join(["cd ", tmpDir, " && ", confs["gffcompare"], "-o gffcmp"] + gffs)
    ret = runCommand(cmd)
    ret2 = runCommand(" ".join(["cat"]+gffs+[">", tmpDir+"/all_mappings.gff"]))
    return ret == 0 and ret2 == 0

def best_id_in_source(ids, hits, source):
    best = ids[0]
    for i in ids:
        best_hit = hits[best]
        hit = hits[i]
        best_len = int(best_hit['end']) - int(best_hit['start'])
        current_len = int(hit['end']) - int(hit['start'])
        if current_len > best_len:
            best = i
    '''if source == "RNAcentral":
    elif source == "cmscan":
    elif source == "LGC":'''
    return best

def best_id(ids, hits):
    id_by_source = {}
    for id_str in ids:
        hit = hits[id_str]
        source = hit["source"]
        if not source in id_by_source:
            id_by_source[source] = list()
        id_by_source[source].append(id_str)
    best_by_source = {}
    for source in id_by_source.keys():
        best_by_source[source] = best_id_in_source(id_by_source[source], hits, source)
    if "RNAcentral" in best_by_source:
        return best_by_source["RNAcentral"]
    elif "tRNAscan-SE" in best_by_source:
        return best_by_source["tRNAscan-SE"]
    elif "cmscan" in best_by_source:
        return best_by_source["cmscan"]
    elif "LGC" in best_by_source:
        return best_by_source["LGC"]
    else:
        print("Error: no known source in " + str(id_by_source))
        return None

def remove_redundancies(args, confs, tmpDir, stepDir):
    print("Reading raw annotation")
    raw_annotation = pd.read_csv(stepDir["run_gffcompare"] + "/all_mappings.gff", sep="\t", header=None)
    raw_annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand",
                                "frame", "attribute"]
    hits = {}
    for index, row in raw_annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        hits[attrs["ID"]] = row

    redundant_ids = []
    loci_file = stepDir["run_gffcompare"] + "/gffcmp.loci"
    print("Reading loci file: " + loci_file)
    loci = pd.read_csv(loci_file, sep="\t")
    new_gff_lines = []
    novel_ids = 0
    not_found_desc = 0
    with open(loci_file, "r") as loci_stream:
        for raw_line in loci_stream:
            line = raw_line.rstrip("\n")
            cells = line.split("\t")
            ids_cells = cells[2:len(cells)]
            ids = []
            for cell in ids_cells:
                for subcell in cell.split(","):
                    if subcell != "-":
                        ids.append(subcell)
            redundant_ids.append(ids)

    print("Calculating best ids")
    non_redundant_ids = [best_id(id_list, hits) for id_list in redundant_ids]
    def is_filtered(row, valid_ids):
        attrs = get_gff_attributes(row["attribute"])
        ID = attrs["ID"]
        filtered = not ID in valid_ids
        return filtered
    raw_annotation["filtered"] = raw_annotation.apply(
        lambda row: is_filtered(row,non_redundant_ids),
        axis=1
    )

    print("All annotations: " + str(len(raw_annotation)))
    raw_annotation = raw_annotation[raw_annotation["filtered"] == False]
    print("Without redundancies: " + str(len(raw_annotation)))
    #annotation.is_copy = False
    raw_annotation = raw_annotation.drop('filtered', axis=1)
    print("Writing final annotation")
    raw_annotation.to_csv(tmpDir + "/annotation.gff", sep="\t", index=False, header=False)

    return True

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
            #"with_rfam_alignment": has_rfam_alignment(hits),
            #"genbank_ids": number_of_genbank_ids(hits),
            "with_rfam_id": has_rfam_id(hits),
            "rfam_ids": number_of_rfam_ids(hits)
        }
        lines.append(line)
    lines.append({
        "source": "ALL",
        "total": annotation.shape[0],
        #"with_rfam_alignment": has_rfam_alignment(annotation),
        #"genbank_ids": number_of_genbank_ids(annotation),
        "with_rfam_id": has_rfam_id(annotation),
        "rfam_ids": number_of_rfam_ids(annotation)
    })

    review = pd.DataFrame(lines, columns=["source", "total", "with_rfam_id","rfam_ids"])
    review.to_csv(tmpDir + "/review.tsv", sep="\t", index=False)

    return True

'''def get_best_id(ids, scores):
    best = 0
    for index in range(0, len(ids)):
        if ids[index] in scores:
            if scores[ids[index]] > scores[ids[best]]:
                best = index
        else:
            return (False, ids[index])
    return (True, ids[best])'''

'''def get_new_genes(args, confs, tmpDir, stepDir):
    print("Reading novel annotation")
    novel_annotation = pd.read_csv(stepDir["parse_infernal"] + "/rfam_annotation_genome.gff", sep="\t", header=None)
    novel_annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand",
                                "frame", "attribute"]
    #novel_annotation = novel_annotation[novel_annotation["feature"] != "noncoding_exon"]
    descriptions = {}
    for index, row in novel_annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        descriptions[attrs["ID"]] = row["attribute"]
        #print("descriptions["+attrs["ID"]+"]="+row["attribute"])
    
    #return False
    
    print("Reading novel gene scores")
    scores = {}
    with open(stepDir["parse_infernal"] + "/scores.tsv", 'r') as scores_file:
        for raw_line in scores_file:
            line = raw_line.rstrip("\n").split("\t")
            scores[line[0]] = float(line[1])

    loci_file = stepDir["run_gffcompare"] + "/gffcmp.loci"
    print("Reading loci file: " + loci_file)
    loci = pd.read_csv(loci_file, sep="\t")
    new_gff_lines = []
    novel_ids = 0
    not_found_desc = 0
    with open(loci_file, "r") as loci_stream:
        for raw_line in loci_stream:
            line = raw_line.rstrip("\n")
            cells = line.split("\t")
            loc_str = cells[1]
            loc_str_halfs = loc_str.split("[")
            sequence_id = loc_str_halfs[0]
            strand = loc_str_halfs[1][0]
            range_str = loc_str_halfs[1][2:len(loc_str_halfs[1])]
            start_end = range_str.split("-")
            start = start_end[0]
            end = start_end[1]
            ids_cells = cells[2:len(cells)]
            ids = []
            for cell in ids_cells:
                for subcell in cell.split(","):
                    if subcell != "-":
                        ids.append(subcell)

            novel, id = get_best_id(ids, scores)
            if novel:
                if id in descriptions:
                    novel_ids += 1
                    gff_line = "\t".join([sequence_id, "cmscan", "transcript",
                        start, end, ".", strand, ".", descriptions[id]])
                    new_gff_lines.append(gff_line)
                else:
                    not_found_desc += 1
                    print("Could not find " + id + " in " + stepDir["parse_infernal"] + "/rfam_annotation_genome.gff")
            
    
    with open(tmpDir + "/novel_annotation.gff", "w") as stream:
        for line in new_gff_lines:
            stream.write(line + "\n")
    
    print("After merging with reference annotation/removing redundancy, " + str(novel_ids) + " were left.")
    print("And could not find descriptions for " + str(not_found_desc))

    return True

def make_complete_annotation(args, confs, tmpDir, stepDir):
    print("Loading data")
    annotation_file = tmpDir + "/annotation.gff"
    lnc_mappings = stepDir["lnc_alignment_parsing"] + "/long_transcripts_LncADeep.to.gigas_genome-short_found.gff"
    infernal_mappings = stepDir["get_new_genes"] + "/novel_annotation.gff"
    if "reference_gff" in args:
        ref_annotation = pd.read_csv(args["reference_gff"], sep="\t", header=None, names=["seqname", 
            "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        ref_annotation = ref_annotation[ref_annotation["feature"] == "transcript"]
        ref_annotation_file = tmpDir + "/ref_annotation.gff"
        ref_annotation.to_csv(ref_annotation_file, sep="\t", index=False, header=False)

        code = runCommand("cat " + ref_annotation_file + " "
            + infernal_mappings + " " + lnc_mappings
            + " > " + annotation_file)
        if code != 0:
            return False
    else:
        code = runCommand("cat " + infernal_mappings + " " + lnc_mappings
            + " > " + annotation_file)
        if code != 0:
            return False

    annotation = pd.read_csv(annotation_file, sep="\t", header=None, names=["seqname", 
        "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

    print("Reading rnacentral2rfam")
    mapping_file = args["rnacentral2rfam"]
    rnacentral2rfam = {}
    with open(mapping_file, 'r') as stream:
        for line in stream:
            cells = line.split("\t")
            rnacentral2rfam[cells[0]] = cells[2]

    print("Assigning RFAM id to genes without an RFAM id yet")
    def update_attribute(row):
        attr_str = row["attribute"]
        attributes = get_gff_attributes(attr_str)
        if not("rfam" in attributes):
            if row["source"] == "RNAcentral":
                short_id = attributes["ID"].split(".")[0]
                rfam_id = get_rfam_from_rnacentral(short_id)
                if rfam_id != None:
                    if len(rfam_id) == 7:
                        attributes["rfam"] = rfam_id
                        print("Assigned " + rfam_id + " to " + short_id)
                if not("rfam" in attributes):
                    print("Could not assign rfam id to " + attributes["ID"])
        return get_gff_attributes_str(attributes)
    annotation["attribute"] = annotation.apply(lambda row: update_attribute(row),axis=1)
    print("Assigning RFAM id to genes without an RFAM id yet (again)")
    annotation["attribute"] = annotation.apply(lambda row: update_attribute(row),axis=1)
    #print("Assigned RFAM id to " + str(rfam_assigns)+ " transcripts")

    print("Reading genome")
    genome_dict = seqListToDict(readSeqsFromFasta(args['genome_link']))
    transcriptome = []
    transcriptome_gff = []
    print("Creating transcriptome file")
    for index, row in annotation.iterrows():
        #print(fasta_header)
        s = genome_dict[row["seqname"]] #cant find key PGUA01000001.1 #TODO
        new_header = get_gff_attributes(row["attribute"])["ID"]
        from_seq =int(row["start"])-1
        to_seq =int(row["end"])
        new_seq = s[from_seq:to_seq]
        transcriptome.append((new_header, new_seq))
        gff_line = {"seqname": new_header, "source": row["source"], "feature": row["feature"], 
                    "start": str(1), "end": str(len(new_seq)), "score": row["score"], 
                    "strand": row["strand"], "frame": row["frame"], "attribute": row["attribute"]}
        transcriptome_gff.append(gff_line)
    writeFastaSeqs(transcriptome, tmpDir + "/transcriptome.fasta")

    print("Writing updated annotation")
    gff_local = pd.DataFrame(transcriptome_gff, columns=["seqname", "source", "feature", 
        "start", "end", "score", "strand", "frame", "attribute"])
    gff_local.to_csv(tmpDir + "/transcriptome_annotation.gff", sep="\t", index=False, header=False)
    annotation.to_csv(annotation_file, sep="\t", index=False, header=False)

    return True'''