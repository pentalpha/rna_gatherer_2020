import sys
from util import runCommand
import os

rfam_path = "Rfam.seed"
if not os.path.exists(rfam_path):
    seed_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz"
    code = run_command("wget " + seed_ftp)
    if code != 0:
        print("Could not download Rfam family informations.")
        quit()

    code = run_command("gzip -d Rfam.seed.gz")
    if code != 0:
        print("Could not extract Rfam family informations from gzip")
        quit()

rfam2rnacentral_path = "rfam.tsv"
if not os.path.exists(rfam2rnacentral_path):
    rfam2rnacentral_ftp = ("ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/"
                        +"current_release/id_mapping/database_mappings/rfam.tsv")
    code = run_command("wget " + rfam2rnacentral_ftp)
    if code != 0:
        print("Could not download rnacentral rfam mappings.")
        quit()

output_dir = sys.argv[1]

lines_writen = set()
family_type = dict()

print("Making rnacentral2rfam from " + rfam2rnacentral_path)
with open(rfam2rnacentral_path,'r') as input_stream:
    with open(output_dir+"/rnacentral2rfam.tsv", 'w') as output_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = cells[0]
            rfam = cells[2]
            tp = cells[4]
            new_line = rnacentral[3:]+"\t"+rfam[2:]+"\n"
            if not new_line in lines_writen:
                output_stream.write(rnacentral[3:]+"\t"+rfam[2:]+"\n")
                lines_writen.add(new_line)

print("Making rfam2type from " + rfam_path)
rfam_types = []
with open(rfam_path, 'r', encoding='Windows-1252') as input_stream:
    current_rfam = ""
    ACs = []
    TPs = []
    looking_for_type = False
    for raw_line in input_stream.readlines():
        line = raw_line.strip()
        if looking_for_type:
            if "#=GF TP" in line:
                TPs.append(line.replace(" ", "").split("#=GFTP")[-1])
                looking_for_type = False
        elif "#=GF AC" in line:
            ACs.append(line.replace(" ", "").split("#=GFAC")[-1])
            looking_for_type = True
    for i in range(len(ACs)):
        rfam_id = ACs[i]
        full_type = TPs[i].rstrip(";").lstrip(";").split(";")
        ncRNA_type = full_type[0]
        if ncRNA_type == "Gene":
            if len(full_type) >= 3:
                ncRNA_type = full_type[2]
            elif len(full_type) == 2:
                ncRNA_type = full_type[1]
            else:
                ncRNA_type = "other"
        rfam_types.append((rfam_id, ncRNA_type, "; ".join(full_type)))

with open(output_dir+"/rfam2type.tsv",'w') as output_stream:
    for rfam_id, rna_type, full_type in rfam_types:
        output_stream.write(rfam_id+"\t"+rna_type+"\t"+full_type+"\n")