import sys
from ete3 import NCBITaxa

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

ncbi = NCBITaxa()

def is_contaminant_taxid(raw_line):
    name = raw_line.rstrip("\n").lstrip(">").split(" ").split("|")[0]
    taxid = int(name.split("_")[-1])
    lineage = ncbi.get_lineage(taxid)
    #bacteria or vector taxids
    if 2 in lineage or 29278 in lineage:
        return True
    else:
        return False

with open(input_fasta, "r") as fasta:
    with open(output_fasta,"w") as contaminants:
        write_seq = False
        for line in fasta:
            if len(line) > 0:
                if line[0] == '>':
                    write_seq = is_contaminant_taxid(line)
                if write_seq:
                    contaminants.write(line)