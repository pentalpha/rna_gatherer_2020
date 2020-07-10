import os
import sys
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from nexus.util import runCommand

from config import *

from nexus.pipeline import Pipeline

from nexus.annotation_steps import *
from nexus.tRNA_scan_se import *
from nexus.reference_processing_steps import *
from nexus.annotation_merging_steps import *
from nexus.lnc_steps import *
from nexus.final_steps import *
from nexus.alignment_steps import *

optional_files = ["non_redundant", "rna_dbs"]
require_files(optional_files, mandatory=False)
mandatory_files = ["go_obo", "rfam_cm"]
require_files(mandatory_files)

#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-g", "--genome", required=False,
        help="Input genome")
    ap.add_argument("-o", "--output", required=True,
        help="Directory where the results of each pipeline step will be stored.")
    ap.add_argument("-sf", "--start-from", required=False,
        help="In case you already runned the entire pipeline, restart from a specific step.")
    ap.add_argument("-st", "--stop-at", required=False,
        help="Only run the pipeline until a specific step.")
    ap.add_argument("-gff", "--reference-gff", required=False,
        help=("Reference annotation for the given genome. Such reference annotations are"
            +" available for several species at " + confs["coordinates_ftp"]))
    ap.add_argument("-ref", "--reference-fasta", required=False,
        help=("Reference sequences of ncRNA for this genome. RNA Nexus will try to find "
            +" their genomic mapping and also retrieve information about them from RNACentral and QuickGO."))
    ap.add_argument("-tx", "--taxon-id", required=False,
        help="Taxon ID for the species, required to retrieve annotations from QuickGo")
    ap.add_argument("-tr", "--transcriptome", required=False,
        help="Fasta file with transcripts that could be lncRNA molecules.")
    ap.add_argument("-ol", "--orf-length", required=False, default = 150,
        help = "Minimum ORF length for a sequence be considered as coding (TransDecoder.LongOrfs).")
    ap.add_argument("-db", "--use-dbs", required=False, default="True",
        help = ("Map transcripts from ncRNA databases, when available, to genome."
            +" True (default) / False."))  
    return vars(ap.parse_args())

#parsing arguments
cmdArgs = getArgs()
outputdir = os.path.abspath(cmdArgs["output"])
if not os.path.exists(outputdir):
    runCommand("mkdir " + outputdir)
argsfile = outputdir + "/args.json"
args = {}
if os.path.exists(argsfile):
    with open(argsfile, "r") as input_stream:
        content = "\n".join(input_stream.readlines())
        args = eval(content)

for arg in cmdArgs:
    if cmdArgs[arg] is not None:
        args[arg] = cmdArgs[arg]

if not "best_hits" in args:
    args["best_hits"] = "False"

#if not os.path.isfile(inputFasta):
#    sys.exit("Input does not exist!")
maxMem = confs["max_mem"]
if not("threads" in confs):
    confs["threads"] = 1
threads = confs["threads"]

args["data_dir"] = outputdir + "/data"
if not os.path.exists(args["data_dir"]):
    runCommand("mkdir " + args["data_dir"])

args["genome_link"] = args["data_dir"] + "/genome.fasta"
if "genome" in args:
    runCommand("ln -s " + args["genome"] + " " + args["genome_link"])
else:
    if not os.path.exists(args["genome_link"]):
        print("No path to genome in arguments given now or previously,"
              +" please specify a genome to use.")
        quit()
index_path = args["genome_link"] + ".mmi"
if not os.path.exists(index_path):
    print("Indexing genome for future mappings...")
    cmd = " ".join([confs["minimap2"], "-x splice:hq -uf -d", 
         index_path, args["genome_link"]])
    code = runCommand(cmd)
    if code != 0:
        print("Could not find or create minimap2 index for genome.")
        quit()
args["genome_index"] = index_path

global_data = os.path.dirname(os.path.realpath(__file__)) + "/data"
confs["rfam2go"] = global_data + "/rfam2go"
if cmdArgs["use_dbs"] == "False":
    confs["rna_dbs"] = {}

if __name__ == '__main__':
    stepFuncs = [("split_genome", split_genome),
                ("run_infernal", run_infernal),
                ("merge_infernal_outs", merge_infernal_outs),
                ("parse_infernal", parse_infernal),
                ("run_trnascan", run_trnascan),
                ("parse_trna", parse_trna),
                ("filter_small_sequences", filter_small_sequences),
                ("filter_long_orfs", filter_long_orfs),
                ("test_coding_potential", test_coding_potential),
                ("parse_coding_potential", parse_coding_potential),
                ("nr_alignment", nr_alignment),
                ("read_nr_alignment", read_nr_alignment),
                ("lnc_alignment_minimap", lnc_alignment_minimap),
                ("lnc_alignment_parsing", lnc_alignment_parsing),
                ("ncrna_alignment_minimap", ncrna_alignment_minimap),
                ("ncrna_alignment_parsing", ncrna_alignment_parsing),
                ("prepare_ref_annotation", prepare_ref_annotation),
                ("map_to_genome", map_to_genome),
                ("get_info", get_info),
                ("run_gffcompare", run_gffcompare),
                ("remove_redundancies", remove_redundancies),
                ("review_annotations", review_annotations),
                ("write_transcriptome", write_transcriptome),
                ("make_id2go", make_id2go)]

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

