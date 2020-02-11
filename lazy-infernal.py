import nexus.bioinfo as bioinfo
import nexus.util as util
import sys
import os

'''
split file: python lazy-infernal.py prepare <input fasta> <directory> <number of parts>
start processing: python lazy-infernal.py start <directory> <rfam model> <threads>
get current results: python lazy-infernal.py get <directory>
'''

operation = sys.argv[1]

def prepare(fasta_path, output_dir, n):
    print("Reading input fasta")
    seqs = bioinfo.readSeqsFromFasta(fasta_path)
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
        part.append(bioinfo.shortFastaHeader(seq))
        cumulative += len(seq[1])
        current_length += len(seq[1])
    if len(part) > 0:
        parts.append(part)
    
    file_names = [output_dir + "/" + str(i) + ".fasta" for i in range(len(parts))]
    util.runCommand("mkdir " + output_dir)
    print("Writing fasta files")
    for i in range(len(parts)):
        bioinfo.writeFastaSeqs(parts[i], file_names[i])

def getFilesWith(directory, name_part):
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(directory):
        for file in f:
            if name_part in file:
                files.append(os.path.join(r, file))
    return files

def infernal(fasta, rfam, threads):
    output_name = fasta.rstrip(".fasta") + ".tsv"
    new_fasta = fasta + "_running"
    util.runCommand("mv " + fasta + " " + new_fasta)
    cmd = ("cmscan -o " + fasta.rstrip(".fasta") + ".out --tblout " + output_name + " -E 0.01 --acc --cpu "
        + str(threads) + " " + rfam + " " + new_fasta)
    util.runCommand("rm " + output_name)
    util.runCommand("rm " + fasta.rstrip(".fasta") + ".out")
    code = util.runCommand(cmd)
    if(code != 0):
        util.runCommand("mv " + new_fasta + " " + fasta)
    else:
        util.runCommand("mv " + new_fasta + " " + fasta + "_done")


def start(output_dir, rfam, threads):
    while(True):
        fasta = ""
        for f in getFilesWith(output_dir, ".fasta"):
            if not(("_running" in f) or ("_done" in f)):
                fasta = f
                break
        if fasta == "":
            break
        infernal(fasta, rfam, threads)

def erase_comments(path):
    tmp = path+".tmp"
    util.runCommand("grep -v '^#' " + path + " > " + tmp)
    util.runCommand("mv " + tmp + " " + path)

def get(output_dir, output_file):
    paths = getFilesWith(output_dir, ".tsv")
    paths2 = getFilesWith(output_dir, ".out")
    if len(paths) > 0:
        if len(paths2) > 0:
            results2 = " ".join(paths2)
            util.runCommand("cat " + results2 + " > " + output_file.rstrip(".tsv")+".out")
            erase_comments(output_file.rstrip(".tsv")+".out")
        results = " ".join(paths)
        util.runCommand("cat " + results + " > " + output_file)
        erase_comments(output_file)
    else:
        print("No results ready yet")


if(operation == "prepare"):
    prepare(sys.argv[2], sys.argv[3], int(sys.argv[4]))
elif(operation == "start"):
    start(sys.argv[2], sys.argv[3], int(sys.argv[4]))
elif(operation == "get"):
    get(sys.argv[2], sys.argv[3])
else:
    print("Unknown operation")
