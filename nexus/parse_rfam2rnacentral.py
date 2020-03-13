import sys

rfam2rnacentral_path = sys.argv[1]
output_dir = sys.argv[2]

lines_writen = set()

with open(rfam2rnacentral_path,'r') as input_stream:
    with open(output_dir+"/rnacentral2rfam.tsv",'w') as output_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = cells[0]
            rfam = cells[2]
            new_line = rnacentral[3:]+"\t"+rfam[2:]+"\n"
            if not new_line in lines_writen:
                output_stream.write(rnacentral[3:]+"\t"+rfam[2:]+"\n")
                lines_writen.add(new_line)