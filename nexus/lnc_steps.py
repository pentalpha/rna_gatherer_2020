from nexus.bioinfo import *
from nexus.util import *
import pandas as pd
import os
from tqdm import tqdm

def filter_small_sequences(args, confs, tmpDir, stepDir):
    #tmpDir = tempDir[0]
    if "transcriptome" in args:
        seqs = readSeqsFromFasta(args["transcriptome"])
        longSeqs = list()

        smallestLength = 9999999
        filtered = list()
        for entry in seqs:
            if len(entry[1]) < smallestLength:
                smallestLength = len(entry[1])
            if len(entry[1]) >= 200:
                longSeqs.append(entry)
            else:
                filtered.append(entry)
        print(str(len(filtered)) + " transcripts were too short to be a lncRNA. Smallest: " + str(smallestLength))
        print(str(len(longSeqs)) + " long enough to be lncRNA.")
        print("Writing them down")
        output_fasta = tmpDir + "/long_transcripts.fasta"
        writeFastaSeqs(longSeqs, output_fasta)
    return True

def filter_long_orfs(args, confs, tmpDir, stepDir):
    #tmpDir = tempDir[1]
    longOrfs = confs["long_orfs"]
    long_transcripts = stepDir["filter_small_sequences"] + "/long_transcripts.fasta"
    if(os.path.exists(long_transcripts)):
        print("Finding ORFs with TransDecoder.LongORFs")
        longOrfsCmd = longOrfs + " -t " + long_transcripts
        code = runCommand("cd " + tmpDir + " && " + longOrfsCmd)
        if code != 0:
            return False
        orfsFile = tmpDir + "/long_transcripts.fasta.transdecoder_dir/longest_orfs.pep"
        orfs = readSeqsFromFasta(orfsFile)
        longOrfs = set()
        for entry in orfs:
            header = entry[0]
            seqLen = len(entry[1])
            name = header.rstrip("\n").lstrip(">").split()[0].split("|")[0].split(".p")[0]
            if seqLen > 100:
                longOrfs.add(name)
        print("Some longOrfs names: " + str(list(longOrfs)[:6]))
        #filtering out transcripts with long orfs
        seqs = readSeqsFromFasta(long_transcripts)
        noLongOrf = list()
        filtered = list()
        filtered, noLongOrf = filterSeqs(seqs, longOrfs)

        print(str(len(filtered)) + " transcripts filtered out")
        writeFastaSeqs(noLongOrf, tmpDir + "/no_orfs.fasta")
        return True
    else:
        return True

def test_coding_potential(args, confs, tmpDir, stepDir):
    long_transcripts_path = stepDir["filter_long_orfs"] + "/no_orfs.fasta"
    if os.path.exists(long_transcripts_path):
        output = tmpDir + "/lgc.csv"
        cmd = " ".join(["python", confs["lgc"], long_transcripts_path, output])
        code = runCommand(cmd)
        if code != 0:
            return False
        noncoding = set()
        with open(output, 'r') as in_stream:
            for raw_line in in_stream.readlines():
                if not raw_line.startswith("#"):
                    cells = raw_line.rstrip("\n").split()
                    if cells[4] == "Non-coding":
                        noncoding.add(cells[0])
        seqs = readSeqsFromFasta(long_transcripts_path)
        noncoding, coding = filterSeqs(seqs, noncoding)

        print(str(len(coding)) + " transcripts with coding potential filtered out")
        writeFastaSeqs(noncoding, tmpDir + "/non_coding.fasta")
    return True

def nr_alignment(args, confs, tmpDir, stepDir):
    query_fasta = stepDir["test_coding_potential"] + "/non_coding.fasta"
    if os.path.exists(query_fasta):
        db = confs["non_redundant"]
        print("Starting diamond to blastx transcripts against NR proteins DB")
        diamond_cmd = [confs["diamond"],"blastx","--db",db,'-q',query_fasta,
                "--out",tmpDir + "/blast.tsv","--outfmt","6",
                "--threads",str(confs['threads']),"--evalue","0.0001"]
        #print(tmpDir)
        #print(str(diamond_cmd))
        code = runCommand("cd " + tmpDir + " && " + " ".join(diamond_cmd))
        return code == 0
    else:
        print("Skipping NR alignment")
        return True

def read_nr_alignment(args, confs, tmpDir, stepDir):
    blastFile = stepDir["nr_alignment"] + "/blast.tsv"
    if os.path.exists(blastFile):
        query_fasta = stepDir["test_coding_potential"] + "/non_coding.fasta"

        print("Filtering according to blastx results")

        seqs = readSeqsFromFasta(query_fasta)
        blast_df = pd.read_csv(blastFile,sep='\t',header=None)
        ids = set(blast_df[0].tolist())
        invalidSeqs, validSeqs = filterSeqs(seqs, ids)

        print(str(len(invalidSeqs)) + " transcripts with protein match filtered")
        print(str(len(validSeqs)) + " transcripts remaining")
        writeFastaSeqs(validSeqs, tmpDir + "/not_protein.fasta")
    else:
        print("Skiping NR alignment parsing")
    return True

def align_to_dbs(args, confs, tmpDir, stepDir):
    fasta_query = stepDir["read_nr_alignment"] + "/not_protein.fasta"
    if os.path.exists(fasta_query):
        dbs_dir = confs["ncrna_dbs_directory"]
        fastas_set = set(getFilesWith(dbs_dir, ".fasta"))
        fastas_set.update(set(getFilesWith(dbs_dir, ".fa")))
        fastas_list = [fasta for fasta in fastas_set]
        sorted_fastas = sorted(fastas_list, key = os.path.getsize)
        print("small ncRNA databases: " + str(sorted_fastas))
        
        for i in tqdm(range(len(sorted_fastas))):
            db = sorted_fastas[i]
            db_name = os.path.basename(db).split(".")[0]
            output = tmpDir + "/"+db_name+".tsv"
            plast = confs["plast"]
            cmd = " ".join([plast,"-p", "plastn", "-d", db, "-i", fasta_query, "-e", "1e-15", "-a", 
                    str(confs["threads"]), "-outfmt 1", "-o", output])
            code = runCommand(cmd)
            if code != 0:
                return False
            seqs = readSeqsFromFasta(fasta_query)
            try:
                blast_df = pd.read_csv(output,sep='\t',header=None)
                ids = set(blast_df[0].tolist())
                invalidSeqs, validSeqs = filterSeqs(seqs, ids)
                print("Searched " + str(len(seqs)) + " on " + db_name 
                    + " and " + str(len(invalidSeqs)) + " were filtered.")
                fasta_query = tmpDir+"/no_"+db_name+".fasta"
                writeFastaSeqs(validSeqs, fasta_query)
            except pd.io.common.EmptyDataError:
                print("Searched " + str(len(seqs)) + " on " + db_name 
                    + " and 0 were filtered.")
                fasta_query = tmpDir+"/no_"+db_name+".fasta"
                writeFastaSeqs(seqs, fasta_query)
        runCommand("cp " + fasta_query + " " + tmpDir+"/lncRNA_only.fasta")
        return True
    return True

def lnc_alignment(args, confs, tmpDir, stepDir):
    lnc_fasta = stepDir["align_to_dbs"] + "/lncRNA_only.fasta"
    if os.path.exists(lnc_fasta):
        code = blast(lnc_fasta, args["genome_link"], tmpDir, 
            threads=confs["threads"], blast_type=(confs["plast"] + " -p plastn"),
            source="LncADeep")
        return code
    else:
        return True

def lnc_alignment_parsing(args, confs, tmpDir, stepDir):
    lnc_fasta = stepDir["align_to_dbs"] + "/lncRNA_only.fasta"
    if os.path.exists(lnc_fasta):
        success = blast_annotate(lnc_fasta, args["genome_link"], tmpDir, 
            threads=confs["threads"], blast_type=(confs["plast"] + " -p plastn"),
            source="LGC", run_blast=False, alternative_outputdir=stepDir["lnc_alignment"])
        if not success:
            return False
        annotation = tmpDir + "/lncRNA_only.to.gigas_genome-short_found.gff"
        runCommand("mv " + annotation + " " + tmpDir + "/lncRNA_annotation.gff")
        annotation = tmpDir + "/lncRNA_annotation.gff"

        print("Writing report")

        fastas = ([args["transcriptome"],
            stepDir["filter_small_sequences"] + "/long_transcripts.fasta",
            stepDir["filter_long_orfs"] + "/no_orfs.fasta",
            stepDir["test_coding_potential"] + "/non_coding.fasta",
            stepDir["read_nr_alignment"] + "/not_protein.fasta"]
        + sorted(getFilesWith(stepDir["align_to_dbs"], ".fasta"), key = os.path.getsize, reverse=True))

        
        with open(tmpDir + "/report.txt", 'w') as stream:
            for i in tqdm(range(len(fastas))):
                fasta = fastas[i]
                seqs = readSeqsFromFasta(fasta)
                stream.write(fasta+": "+str(len(seqs))+" sequences.\n")
            annotation_file = open(annotation, 'r')
            stream.write(annotation + ": " + str(len(annotation_file.readlines()))+" mapped to genome.\n")
            annotation_file.close()
        return True
    else:
        return True