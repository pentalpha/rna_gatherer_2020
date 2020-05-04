import sys

'''ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/rfam.tsv'''
rfam2rnacentral_path = sys.argv[1]
'''ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz'''
family_path = sys.argv[2]
output_dir = sys.argv[3]

lines_writen = set()
family_type = dict()

with open(rfam2rnacentral_path,'r') as input_stream:
    with open(output_dir+"/rnacentral2rfam.tsv", 'w') as output_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = cells[0]
            rfam = cells[2]
            tp = cells[4]
            family_type[rfam] = tp
            new_line = rnacentral[3:]+"\t"+rfam[2:]+"\n"
            if not new_line in lines_writen:
                output_stream.write(rnacentral[3:]+"\t"+rfam[2:]+"\n")
                lines_writen.add(new_line)

with open(family_path, 'r') as input_stream:
    for line in input_stream.readlines():
        cells = line.rstrip("\n").split()
        rfam = cells[0]
        name = cells[1]
        if not rfam in family_type:
            family_type[rfam] = name
        elif family_type[rfam] == "other" or family_type[rfam] == "misc_rna":
            family_type[rfam] = name

with open(output_dir+"/rfam2type.tsv",'w') as output_stream:
    for rfam_id, rna_type in family_type.items():
        output_stream.write(rfam_id+"\t"+rna_type+"\n")