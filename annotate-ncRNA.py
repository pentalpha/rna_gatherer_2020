import os
import sys
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from nexus.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from nexus.bioinfo import cluster_all_ranges
from nexus.bioinfo import read_plast_extended, best_hit
from nexus.bioinfo import get_gff_attributes, get_gff_attributes_str
from nexus.bioinfo import get_rfam_from_rnacentral
from nexus.util import runCommand, write_file, getFilesWith

from config import configs

from nexus.pipeline import Pipeline

from nexus.annotation_steps import *
from nexus.reference_processing_steps import *
from nexus.annotation_merging_steps import *
from nexus.lnc_steps import *
from nexus.counting_steps import *

#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--infernal-output", required=False,
        help="Lazy-infernal output")
    ap.add_argument("-g", "--genome", required=False,
        help="Input genome")
    ap.add_argument("-o", "--output", required=True,
        help="Directory where the results of each pipeline step will be stored.")
    ap.add_argument("-sf", "--start-from", required=False,
        help="In case you already runned the entire pipeline, restart from a specific step.")
    ap.add_argument("-st", "--stop-at", required=False,
        help="Only run the pipeline until a specific step.")
    '''ap.add_argument("-bh", "--best-hits", required=False,
            help="Defalt: False. If True, structRNAfinder wont report all the results")'''
    ap.add_argument("-gff", "--reference-gff", required=False,
        help="Reference annotation for the given genome.")
    ap.add_argument("-tr", "--transcriptome", required=False,
        help="Fasta file with transcripts that could be lncRNA molecules.")
    '''ap.add_argument("-ct", "--coding-transcriptome", required=False,
        help="Fasta file with nucleotide sequences for coding genes.")
    ap.add_argument("-fq", "--fastq-directory", required=False,
        help="Directory containing sequencing fastq files (fq, fq.gz , fastq or fastq.gz).")
    ap.add_argument("-go", "--coding-gene-ontology", required=False,
        help="id2gos file corresponding to the sequences in the coding transcriptome.")
    ap.add_argument("-K", "--k-min-coexpressions", required=False,
        default=1, help=("The minimum number of ncRNAs a Coding Gene must be coexpressed with."
                        +" Increasing the value improves accuracy of functional assignments, but"
                        +" may restrict the results. Default: 1."))'''
    return vars(ap.parse_args())

#parsing arguments
cmdArgs = getArgs()
outputdir = os.path.abspath(cmdArgs["output"])
argsfile = outputdir + "/args.json"
args = {}
if os.path.exists(argsfile):
    with open(argsfile, "r") as input_stream:
        content = "\n".join(input_stream.readlines())
        args = eval(content)

for arg in cmdArgs:
    if cmdArgs[arg] is not None:
        args[arg] = cmdArgs[arg]

if not "genome" in args:
    print("No --genome defined, cannot run.")
    quit()

if not "best_hits" in args:
    args["best_hits"] = "False"

#if not os.path.isfile(inputFasta):
#    sys.exit("Input does not exist!")
if not("max_mem" in confs):
    confs["max_mem"] = 2000
maxMem = confs["max_mem"]
if not("threads" in confs):
    confs["threads"] = 1
threads = confs["threads"]

args["data_dir"] = outputdir + "/data"
if not os.path.exists(args["data_dir"]):
    runCommand("mkdir " + args["data_dir"])
args["genome_link"] = args["data_dir"] + "/genome.fasta"

#plast_cmd = [args['plast'], "-p plastx", "-d", NR, 
#        '-i', query_fasta, "-e", "0.0001", "-a", str(threads), 
#        "-outfmt 1", "-bargraph", "-o", tmpDir + "/plast.tsv"]

if __name__ == '__main__':
    stepFuncs = [("split_genome", split_genome),
                ("run_infernal", run_infernal),
                ("merge_infernal_outs", merge_infernal_outs),
                ("parse_infernal", parse_infernal),
                ("get_reference_rfam_ids", get_reference_rfam_ids),
                ("filter_small_sequences", filter_small_sequences),
                ("filter_long_orfs", filter_long_orfs),
                ("test_coding_potential", test_coding_potential),
                ("nr_alignment", nr_alignment),
                ("read_nr_alignment", read_nr_alignment),
                ("align_to_dbs", align_to_dbs),
                ("lnc_alignment", lnc_alignment),
                ("lnc_alignment_parsing", lnc_alignment_parsing),
                #("analyze", analyze),
                ("run_gffcompare", run_gffcompare),
                ("remove_redundancies", remove_redundancies),
                #("make_complete_annotation", make_complete_annotation),
                #("make_family_fastas", make_family_fastas),
                #("search_novel_names", search_novel_names),
                #("annotate_novel_names", annotate_novel_names),
                ("review_annotations", review_annotations),
                ("write_transcriptome", write_transcriptome),
                ("make_ids2go", make_ids2go)]

    pipe = Pipeline(args, confs, stepFuncs, args["output"])
    if pipe.ready:
        print("Arguments file is " + argsfile)
        with open(argsfile, "w") as output_stream:
            args_str = str(args)
            print("Writing args\n"+args_str)
            output_stream.write(args_str)
        if ("start_from" in args) and ("stop_at" in args):
            pipe.run(start_from=args["start_from"], stop_at=args["stop_at"])
        elif "start_from" in args:
            pipe.run(start_from=args["start_from"])
        elif "stop_at" in args:
            pipe.run(stop_at=args["stop_at"]) 
        else:
            pipe.run()

