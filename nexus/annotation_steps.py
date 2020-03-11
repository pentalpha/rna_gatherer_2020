import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges, shortFastaHeader
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.util import runCommand, write_file, getFilesWith

def infernal(fasta, cmscan, rfam, threads):
    output_name = fasta.rstrip(".fasta") + ".tsv"
    output_name2 = fasta.rstrip(".fasta") + ".out"
    new_fasta = fasta + "_running"
    runCommand("mv " + fasta + " " + new_fasta)
    cmd = (cmscan + " -o " + output_name2 + " --tblout " + output_name + " -E 0.01 --acc --cpu "
        + str(threads) + " " + rfam + " " + new_fasta)
    runCommand("rm " + output_name)
    runCommand("rm " + fasta.rstrip(".fasta") + ".out")
    code = runCommand(cmd)
    if(code != 0):
        runCommand("mv " + new_fasta + " " + fasta)
    else:
        runCommand("mv " + new_fasta + " " + fasta + "_done")

def check_progress(output_dir):
    total = 0
    done = 0
    for f in getFilesWith(output_dir, ".fasta"):
        if "_done" in f:
            done += 1
        total += 1
    complete = (done == total)
    return complete, total, done

def print_progress(output_dir):
    complete, total, done = check_progress(output_dir)
    print("Infernal processing: " + str(done) + "/" + str(total))
    return complete

def start(output_dir, cmscan, rfam, threads):
    while(True):
        fasta = ""
        for f in getFilesWith(output_dir, ".fasta"):
            if not(("_running" in f) or ("_done" in f)):
                fasta = f
                break
        if fasta == "":
            break
        infernal(fasta, cmscan, rfam, threads)
        print_progress(output_dir)
        
def erase_comments(path):
    tmp = path+".tmp"
    runCommand("grep -v '^#' " + path + " > " + tmp)
    runCommand("mv " + tmp + " " + path)

def get_infernal_output(output_dir, output_file):
    paths = getFilesWith(output_dir, ".tsv")
    paths2 = getFilesWith(output_dir, ".out")
    if len(paths) > 0:
        if len(paths2) > 0:
            results2 = " ".join(paths2)
            runCommand("cat " + results2 + " > " + output_file.rstrip(".tsv")+".out")
            erase_comments(output_file.rstrip(".tsv")+".out")
        results = " ".join(paths)
        runCommand("cat " + results + " > " + output_file + ".tsv")
        erase_comments(output_file)
        return True
    else:
        print("No results ready yet")
        return False

def read_infernal_output(infernal_path):
    infernal_tsv = infernal_path
    full_names = {}
    lines = []
    print("Reading lines from infernal tabular output")
    with open(infernal_tsv) as f:
        line = f.readline()
        to_print = 0
        while line != "":
            if line[0] != "#":
                elements = [x for x in line.split() if x != ""]
                cols = {"rna_name": elements[0],"rfam":elements[1], "seqname": elements[2], 
                        "qstart": elements[7], "qend": elements[8], "strand": elements[9],
                        "evalue": elements[15]}
                if to_print > 0:
                    print("Line:\n\t"+line)
                    print("\t"+str(cols))
                to_print -= 1
                lines.append(cols)
            line = f.readline()
    #print(str(lines))

    print("Creating gff file about seqs identified")
    rows = []
    name_used = {}
    for hit in lines:
        base_name = hit['seqname']+"_"+hit['rna_name']
        if not base_name in name_used:
            name_used[base_name] = 0
        unique_name = base_name+"_"+str(name_used[base_name])
        name_used[base_name] += 1
        start = min(int(hit["qstart"]),int(hit["qend"]))
        end = max(int(hit["qstart"]),int(hit["qend"]))
        row = {"seqname": hit['seqname'], "source": "cmscan",
        "feature": "transcript", "start": start,
        "end":end, "score": hit["evalue"],
        "strand": hit["strand"],
        "frame": ".", "attribute":"ID="+unique_name+";name="+hit['rna_name']
            +";evalue="+hit['evalue']+";rfam="+hit['rfam']}
        rows.append(row)
    gff = pd.DataFrame(rows, columns = ["seqname", "source",
        "feature", "start", "end", "score", "strand", 
        "frame", "attribute"])
    #gff_name = tmpDir+"/infernal_found.gff"
    #gff.to_csv(gff_name, sep="\t", index=False, header = False)

    #print("Writing fasta files.")
    #known = set(gff["seqname"].unique().tolist())
    return gff

def split_genome(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    if not "genome" in args:
        print("Cannot run annotation without a path to a genome.")
        return False
    fasta_path = os.path.abspath(args["genome"])
    n = 100
    #creating shortcut to genome fasta
    runCommand("ln -s " + fasta_path + " " + args["genome_link"])
    fasta_path = args["genome_link"]
    print("Reading input fasta")
    seqs = readSeqsFromFasta(fasta_path)
    total_length = sum([len(entry[1]) for entry in seqs])
    print("Total length:" + str(total_length))
    max_length = int(total_length / n)
    
    current_length = 0
    part = []
    parts = []
    print("Spliting parts of fasta")
    cumulative = 0
    parts_done = 0
    for seq in seqs:
        if n > parts_done:
            max_length = int((total_length-cumulative) / (n-parts_done))
        if ((current_length >= max_length)):
            parts.append(part)
            parts_done += 1
            part = []
            current_length = 0
        part.append(shortFastaHeader(seq))
        cumulative += len(seq[1])
        current_length += len(seq[1])
    if len(part) > 0:
        parts.append(part)
    
    file_names = [output_dir + "/" + str(i) + ".fasta" for i in range(len(parts))]
    runCommand("mkdir " + output_dir)
    print("Writing fasta files")
    for i in range(len(parts)):
        writeFastaSeqs(parts[i], file_names[i])
    return True

def run_infernal(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    if print_progress(output_dir):
        return True
    else:
        start(output_dir, confs["cmscan"], confs["rfam_cm"], confs["threads"])
        return print_progress(output_dir)

def merge_infernal_outs(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    output_file = tmpDir + "/infernal"
    success = get_infernal_output(output_dir,output_file)
    return success

def parse_infernal(args, confs, tmpDir, stepDir):
    infernal_tsv = stepDir["merge_infernal_outs"] + "/infernal.tsv"
    gff = read_infernal_output(infernal_tsv)
    print("Writing .gff file")
    gff.to_csv(tmpDir + "/rfam_annotation_genome.gff", sep="\t", index=False, header=False)
    return True

def analyze(args, confs, tmpDir, stepDir):
    print("Analyzing cmscan annotation")
    annotation = pd.read_csv(stepDir["parse_infernal"] + "/rfam_annotation_genome.gff", sep="\t", header=None)
    annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    rfams = dict()
    n_families = 0
    ids = set()
    exons = 0
    for index, row in annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        rfam = attrs["rfam"]
        if row["feature"] == "noncoding_exon":
            exons += 1
        if rfam in rfams:
            rfams[rfam] += 1
        else:
            rfams[rfam] = 1
            n_families += 1
        ids.add(attrs["transcript_id"])
    with open(tmpDir + "/stats.txt", 'w') as stream: 
        stream.write(str(len(ids)) + " transcript_id\n")
        stream.write(str(exons) + " exons\n")
        stream.write(str(n_families) + " families\n")

    pairs = []
    with open(tmpDir + "/family_sizes.tsv", 'w') as stream:
        for key in rfams.keys():
            pairs.append((key, rfams[key]))
            stream.write(key + "\t" + str(rfams[key]) + "\n")
    pairs.sort(key=lambda x: x[1])
    sorted_keys = [x[0] for x in pairs]
    sorted_vals = [x[1] for x in pairs]
    import matplotlib as mpl
    mpl.rcParams['ytick.labelsize'] = 6
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, len(sorted_keys)/10))
    #plt.xticks(fontsize=8)
    plt.barh([4*x for x in range(len(sorted_keys))], sorted_vals, height=3.0,align='center', tick_label=sorted_keys)
    for i, v in enumerate(sorted_vals):
        plt.text(v, 4*i, " " + str(v), fontsize=6, color='blue', va='center', fontweight='bold')
    #lt.xticks(range(len(sorted_keys)), sorted_keys)
    plt.savefig(args, confs, tmpDir + "/family_sizes.png", dpi = 300, format='png')
    
    return True

'''def run_infernal(args, confs, tmpDir, stepDir):
    infernal_tsv = os.path.abspath(args["infernal_output"])
    #infernal_out = infernal_tsv.rstrip(".tsv")+".out"
    genome_fasta = os.path.abspath(args["genome"])

    runCommand("cp " + infernal_tsv + " " + tmpDir + "/infernal.tsv")
    #runCommand("cp " + infernal_out + " " + tmpDir + "/new_input.out")
    runCommand("cp " + genome_fasta + " " + tmpDir + "/input.fasta")

    ''''''inputFasta = tmpDir + "/input.fasta"
    cmd = " ".join([confs["structRNAfinder"] + " -i", inputFasta, 
        "-d" , confs["rfam_cm"], "-t", 
        tmpDir + "/infernal.tsv -c", str(confs['threads'])])
    if args['best_hits'] == 'False':
        cmd += " -r"
        print("Reporting all hits")
    else:
        print("Retorping best hits")
    runCommand("cd " + tmpDir + " && " + cmd)''''''
    return True

def parse_infernal(args, confs, tmpDir, stepDir):
    #input_file_name = inputFasta.split("/")[-1]
    #input_no_extension = input_file_name.split(".")
    #input_no_extension = ".".join(input_no_extension[:-1])
    #input_no_extension = ".".join(inputFasta.split("/")[-1].split(".")[:-1])
    infernal_tsv = stepDir["merge_infernal_outs"] + "/infernal.tsv"
    gff = read_infernal_output(infernal_tsv)
    ''''''base_name = stepDir["run_infernal"] + "/input_filtered."
    bed_path = base_name + "bed"
    tab_path = base_name + "tab"

    print("Reading .bed results")
    bed = pd.read_csv(bed_path, sep="\t", header=None)
    bed.columns = ["query", "start", "end", "id", "score", "strand"]
    print("Reading .tab results")
    tab = pd.read_csv(tab_path, sep="\t", header=None, encoding='ISO-8859-1')
    tab = tab.drop([0],axis=0)
    tab.columns = ["Contigs","Family","ID","Score","E-value","from_seq","to_seq","type","domain","description"]
    print("Loaded "+str(len(tab))+ " lines from .tab file.")
    #print(tab.head())

    print("Reading fasta")
    seqs = readSeqsFromFasta(stepDir["run_infernal"] + "/input.fasta")
    print("Creating dict from fasta")
    seqs_dict = seqListToDict(seqs)
    n_genes_in_seq = {}
    for key in seqs_dict.keys():
        n_genes_in_seq[key] = 0
    print("Loaded "+str(len(seqs_dict))+" sequences from input fasta.")
    transcriptome = []
    seq = {}
    n = 0
    families = set()
    print("Creating transcriptome file")
    for index, row in tab.iterrows():
        #name = row["Contigs"]+"-"+str(row["from_seq"])+"-"+str(row["to_seq"])
        split = row[0].split("_")
        fasta_header = "_".join(split[:-1])
        if len(split) == 1:
            fasta_header = row[0]
        #print(fasta_header)
        if not fasta_header in seqs_dict:
            fasta_header = row[0]
            if not fasta_header in seqs_dict:
                print(fasta_header + " from " + str(row) + " not found in seqs_dict")
                pass
        s = seqs_dict[fasta_header]
        #new_header = row[0] + " " +  str(row["Family"]) + " " + str(row["ID"])
        new_header = fasta_header + "-ncRNA_" + str(n_genes_in_seq[fasta_header] + 1)
        n_genes_in_seq[fasta_header] += 1
        families.add(row["ID"])
        from_seq =int(row["from_seq"])
        to_seq =int(row["to_seq"])
        start = min(from_seq,to_seq) - 1
        end = max(from_seq,to_seq)
        new_seq = s[start:end]
        seq[row[0]] = new_seq
        transcriptome.append((new_header, new_seq))
        n += 1
    print("Families: " + str(len(families)))
    print("Indexed " + str(n) + " entries from .tab file")
    print("Indexing .bed results")
    query_strand = []
    subject_id = []
    i = 0
    for index, row in bed.iterrows():
        #name = row["query"]+"-"+str(row["start"])+"-"+str(row["end"])
        query_strand.append(row["strand"])
        subject_id.append(row["id"])

    print("Creating .gff lines")
    lines = []
    lines2 = []
    scores = {}
    for index, row in tab.iterrows():
        #name = row["Contigs"]+"-"+str(row["from_seq"])+"-"+str(row["to_seq"])
        split = row[0].split("_")
        fasta_header = "_".join(split[:-1])
        if len(split) == 1:
            fasta_header = row[0]
        sequence = seq[row[0]]
        row_type = row["type"]
        if type(row_type) != str:
            if np.isnan(row_type):
                row_type = "other"
            else:
                print(str(row))
                print("Invalid entry 'type': not an str or NaN")
                return False

        type_words = row_type.split(";")
        tp = []
        for word in type_words:
            if len(word) > 1:
                tp.append(word)
        if len(tp) >= 2:
            tp = tp[1]
        else:
            tp = tp[0]
        tp = tp.replace(" ", "")
        if tp == "":
            tp = row_type
        from_seq = int(row["from_seq"])
        to_seq = int(row["to_seq"])
        start = min(from_seq, to_seq)
        end = max(from_seq, to_seq)
        scores[str(row[0])] = row["Score"]
        line = {"seqname":  fasta_header,
                "source": "cmscan", "feature": "transcript",
                "start": str(start),
                "end": str(end),
                "score": ".",
                "strand": str(query_strand[i]),
                "frame": ".",
                "attribute": "rfam=" + str(row["ID"])
                + ";family=" + str(row["Family"])
                + ";transcript_id=" + str(row[0])
                + ";type=" + tp
                + ";evalue=" + str(row["E-value"])
                + ";ID="+str(row[0])}
        gene_line1 = line.copy()
        gene_line1["feature"] = "noncoding_exon"
        gene_line1["attribute"] = gene_line1["attribute"] + ":ncRNA_exon1;Parent="+str(row[0])
        line2 = {"seqname": row[0],
                "source": "cmscan", "feature": "transcript",
                "start": str(1),
                "end": str(len(sequence)),
                "score": ".",
                "strand": str(query_strand[i]),
                "frame": ".",
                "attribute": "rfam=" + str(row["ID"])
                + ";family=" + str(row["Family"])
                + ";transcript_id=" + str(row[0])
                + ";type=" + tp
                + ";evalue=" + str(row["E-value"])
                + ";ID="+str(row[0])}
        gene_line2 = line2.copy()
        gene_line2["feature"] = "noncoding_exon"
        lines.append(line)
        lines.append(gene_line1)
        lines2.append(line2)
        lines2.append(gene_line2)
        i += 1;'''
'''
    print("Writing .gff file")
    ''''''gff = pd.DataFrame(lines, columns=["seqname", "source", "feature", 
        "start", "end","score", "strand", "frame", "attribute"])
    gff_local = pd.DataFrame(lines2, columns=["seqname", "source", "feature", 
        "start", "end","score", "strand", "frame", "attribute"])''''''
    gff.to_csv(tmpDir + "/rfam_annotation_genome.gff", sep="\t", index=False, header=False)
    ''''''gff_local.to_csv(tmpDir + "/rfam_annotation_ncRNA.gff", sep="\t", index=False, header=False)
    print("Writing scores.tsv")
    with open(tmpDir + "/scores.tsv", "w") as scores_file:
        for key in scores.keys():
            scores_file.write(key + "\t" + str(scores[key]) + "\n")
    print("Writing .fasta file")
    writeFastaSeqs(transcriptome, tmpDir + "/ncRNA.fasta")'''
    #identified, not_found = filterSeqs(seqs, [s.split('_')[0] for s in query_strand])
    #writeFastaSeqs(not_found, tmpDir + "/no_rfam.fasta")
'''
    return True
'''

'''def gene_ontology(args, confs, tmpDir, stepDir):
    rfam2go_path = args["rfam2go"]
    rfam2go = {}
    with open(rfam2go_path) as stream:
        for line in stream:
            x = line.rstrip("\n")
            rfam_and_gos = x.split("\t")
            rfam2go[rfam_and_gos[0]] = rfam_and_gos[1].split(";")
    annotation = pd.read_csv(stepDir["parse_infernal"] + "/rfam_annotation.gff", sep="\t")'''
    #rfams = set(list(rfam2go.keys()))
    #print(str(len(rfams)))
'''gene_to_go = ""
    genes = ""
    for index, row in annotation.iterrows():
        seqname = row["seqname"]
        attrs = get_gff_attributes(row["attribute"])
        rfam = attrs["id"]
        if rfam in rfam2go:
            gene_to_go += (seqname + "\t"
                    + ";".join(rfam2go[rfam]) + "\n");
            genes += seqname + "\n"
    write_file(gene_to_go, tmpDir + "/gigas_association.tsv")
    write_file(genes, tmpDir + "/gigas_genes.tsv")'''
'''#print(str(rfam2go))
    rfam_list = ""
    for family in families:
        if family in rfam2go:
            rfam_list += family + "\n"
#            gos = rfam2go[family]
#            #print(str(gos))
#            for go in gos:
#                go_list += go + "\n"
    with open(args, confs, tmpDir + "/rfam_file.txt", "w") as stream:
        stream.write(rfam_list)
'''
'''
 example find_enrichment command:
        python ~/software/goatools/scripts/find_enrichment.py rfam_file.txt \
            ~/data/go_rfam_population.tsv ~/data/go_association.tsv \
            --pval=0.05 --method=fdr_bh --pval_field=fdr_bh \
            --outfile=result-no_duplicates.tsv --obo ~/data/go-basic.obo \
            --goslim ~/data/goslim_generic.obo --taxid=113544'''
'''    return True'''


'''def make_family_fastas(args, confs, tmpDir, stepDir):
    print("Reading annotation")
    transcript_fasta = readSeqsFromFasta(stepDir["parse_infernal"] + "/ncRNA.fasta")
    annotation = pd.read_csv(stepDir["make_complete_annotation"] + "/annotation.gff", sep="\t", header=None)
    annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand",
                        "frame", "attribute"]

    print("Grouping transcripts by family")
    families = {}
    for index, row in annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            family = attrs["rfam"]
            if not(family in families):
                families[family] = []
            families[family].append(attrs["ID"])

    print("Filtering seqs and writing group fastas")
    for key in families.keys():
        ids = families[key]
        valid, invalid = filterSeqs(transcript_fasta, ids)
        writeFastaSeqs(valid, tmpDir + "/" + key + ".fasta")

    return True
'''
'''
def search_novel_names(args, confs, tmpDir, stepDir):
    print("Aligning family fastas to RFAM")
    input_fastas = getFilesWith(stepDir["make_family_fastas"], ".fasta")
    print("Running alignmens:")
    for i in tqdm(range(len(input_fastas))):
        input_fasta = input_fastas[i]
        rfam_id = input_fasta.split("/")[-1].split(".")[0]
        rfam_fasta = confs["rfam_fastas"] + "/" + rfam_id + ".fa"
        output = tmpDir + "/" + rfam_id + ".tsv"
        if not os.path.exists(rfam_fasta):
            print(rfam_fasta + " not found, ending pipeline.")
            return False
        code = runCommand(" ".join(["plast -p plastn", "-d", rfam_fasta, "-i", input_fasta, "-e 0.001",
                        "-a", str(threads), "-outfmt 1", "-bargraph", "-o", output, " > /dev/null"]), print_cmd=False)
        if code != 0:
            return False
    return True
'''
'''
def annotate_novel_names(args, confs, tmpDir, stepDir):
    ''''''print("Retrieving RFAM gene names")
    family_fastas = getFilesWith(stepDir["make_family_fastas"], ".fasta")
    used_families = set()
    for fasta in family_fastas:
        used_families.add(os.path.split(fasta)[1].rstrip(".fasta"))
    rfam_fastas = getFilesWith(args["rfam_fastas"], ".fa")
    seq_extended_name = {}
    for i in tqdm(range(len(rfam_fastas))):
        fasta = rfam_fastas[i]
        rfam_id = os.path.split(fasta)[1].rstrip(".fa")
        if rfam_id in used_families:
            headers = getFastaHeaders(fasta)
            for header in headers:
                words = header.split(" ")
                id = words[0]
                name = " ".join(words[1:])
                seq_extended_name[id] = name
    
    test = {}
    limit = 7
    for k in seq_extended_name.keys():
        test[k] = seq_extended_name[k]
        limit -= 1
        if limit < 0:
            break
    print(str(test))'''
'''
    print("Retrieving alignments to RFAM")
    alignment_files = getFilesWith(stepDir["search_novel_names"], ".tsv")
    alignments_file = tmpDir + "/alignments.tsv"
    code = runCommand("cat " + " ".join(alignment_files) + " > " + alignments_file, print_cmd=False)
    if code != 0:
        return False

    print("Reading annotation")
    annotation = pd.read_csv(stepDir["make_complete_annotation"] + "/annotation.gff", sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Reading alignments")
    alignments = pd.read_csv(alignments_file, sep="\t", header=None, names=["qseqid", "sseqid", "pident", 
        "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    genbank_accn = {}
    n_aligned = len(alignments["qseqid"].unique())
    n_annotated = len(annotation.index)

    processed = 0
    print_interval = n_aligned / 7
    last_print = 0
    print("Retrieving best hits from alignments")
    for name, hits in alignments.groupby(["qseqid"]):
        hit = best_hit(hits)
        genbank_accn[name] = hit["sseqid"]
        processed += 1
        if processed - last_print >= print_interval:
            print("Processed " + str(processed)+"/"+str(n_aligned))
            last_print = processed

    print("Assigning genbank accn to novel genes")
    def update_attribute(attr_str):
        attributes = get_gff_attributes(attr_str)
        if attributes["ID"] in genbank_accn:
            attributes["genbank"] = genbank_accn[attributes["ID"]]
        return get_gff_attributes_str(attributes)

    annotation["attribute"] = annotation.apply(lambda row: update_attribute(row["attribute"]),axis=1)

    print("Writing updated annotation")
    annotation.to_csv(args, confs, tmpDir + "/annotation.gff", header=False, sep="\t", index=False)
    print(str(n_aligned) + "/" + str(n_annotated) + " successfully aligned to RFAM sequences.")
    return True
'''