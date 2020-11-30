import sys

'''
Joins the "bests" table for each ontology into a single file that 
    lists the best metric combinations for each confidence level and ontology.
usage: make_metric_table.py BP_FILE MF_FILE CC_FILE <output_file>
'''

onto_files = {"BP": sys.argv[1], "MF": sys.argv[2], "CC": sys.argv[3]}
output_path = sys.argv[4]

onto_tables = {}

for onto, onto_file in onto_files.items():
    with open(onto_file, 'r') as stream:
        table = []
        first = True
        for raw_line in stream.readlines():
            if not first:
                cells = raw_line.rstrip('\n').split("\t")
                table.append([cells[1],cells[2]])
            first = False
        onto_tables[onto] = table

output = []
for i in range(len(onto_tables["BP"])):
    conf_level = onto_tables["BP"][i][0]
    bp_metrics = onto_tables["BP"][i][1]
    mf_metrics = onto_tables["MF"][i][1]
    cc_metrics = onto_tables["CC"][i][1]
    output.append([conf_level,bp_metrics,mf_metrics,cc_metrics])

with open(output_path, 'w') as stream:
    stream.write("\n".join(["\t".join(table_line) for table_line in output]) + "\n")