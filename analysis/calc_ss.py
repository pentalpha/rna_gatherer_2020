import sys
import os
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

id_chunk_len = 1000
pair_chunk_len = 300000
'''ontologies = ["molecular_function",
    "biological_process",
    "cellular_component"]'''
ontologies = ["molecular_function"]
threads = max(2, multiprocessing.cpu_count()-1)

'''
usage:
python calc_ss.py gene_names_file previous_sim_file 
    fs_metric_name annotation_gaf ontology_go new_sim_file
'''
gene_names = sys.argv[1]
calculated_sims_file = sys.argv[2]
fs_metric = sys.argv[3]
gaf_file = sys.argv[4]
ontology_file = sys.argv[5]
output_sims = sys.argv[-1]


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def runCommand(cmd, print_cmd=True):
    if print_cmd:
        print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def append_sims(sims, output_file):
    if len(sims) > 0:
        with open(output_file, "a+") as out_stream:
            for sim in sims:
                new_line = sim + "\n"
                out_stream.write(new_line)

def calc_thread(pid, some_pairs):
    #create input file
    temp_input_path = pid+".tsv.temp"
    with open(temp_input_path, "w+") as temp_file:
        for pair in some_pairs:
            temp_file.write(pair+"\n")
    
    temp_outs = []
    for onto in ontologies:
        temp_out = fs_metric + "." + temp_input_path + "." + onto
        cmd = " ".join(["fastsemsim", "--ontology_file", ontology_file,
        "--ac_type gaf-2.0", "--ac_file", gaf_file, "--tss", fs_metric, 
        "--query_mode pairs --query_ss_type obj --query_input ac",
        "--query_file", temp_input_path, "-o GeneOntology", 
        "--root", onto, "--task SS", "--output_file", temp_out])
        code = runCommand(cmd,print_cmd=False)
        if code != 0:
            print("The command:\t"+cmd+"\n\tresulted in an error:"+str(code))
            return None
        temp_outs.append(temp_out)
    
    local_new_sims = {}
    for index in range(len(temp_outs)):
        out = temp_outs[index]
        added_value = False
        first = True
        with open(out, 'r') as temp_stream:
            for raw_line in temp_stream:
                if not first:
                    cells = raw_line.rstrip("\n").split("\t")
                    local_pair = cells[0]+"\t"+cells[1]
                    if not local_pair in local_new_sims:
                        local_new_sims[local_pair] = [None,None,None]
                    if len(cells) >= 3:
                        if len(cells[2]) > 0:
                            try:
                                local_new_sims[local_pair][index] = float(cells[2])
                                added_value = True
                            except ValueError:
                                print("Could not convert third column to float:")
                                print(str(cells))
                                return None
                else:
                    first = False
        os.remove(out)
        if not added_value:
            print("The command:\t"+cmd+"\n\tresulted in an error:")
            print("No values to be added from " + out)
            return None
    
    lines_produced = [pair_str + "\t" + "\t".join([str(x) for x in sims]) 
                    for pair_str, sims in local_new_sims.items()]
    os.remove(temp_input_path)
    return lines_produced

def calc_sims(pair_lists):
    processes = []
    result_lines = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        i = 1
        for some_pairs in pair_lists:
            processes.append(executor.submit(calc_thread,
                                            str(i), some_pairs))
            i += 1
        print("Created " + str(i-1) + " workers.")
        for task in as_completed(processes):
            new_lines = task.result()
            if new_lines != None:
                result_lines += new_lines
                print("batch finished")
            else:
                print("batch failed")
    return result_lines

print("Loading IDs")
ids = []
with open(gene_names,'r') as stream:
    for raw_line in stream:
        ids.append(raw_line.rstrip("\n").split("\t")[0])
id_chunks = list(chunks(ids, id_chunk_len))

print("Loading previous sims")
calculated_pairs = set()
calculated_sims = []
if os.path.exists(calculated_sims_file):
    with open(calculated_sims_file,'r') as stream:
        for raw_line in stream:
            cells = raw_line.rstrip("\n").split("\t")
            calculated_pairs.add(cells[0]+"\t"+cells[1])
            calculated_pairs.add(cells[1]+"\t"+cells[0])
            calculated_sims.append(cells)
append_sims(calculated_sims, output_sims)

print("Creating pairs")
chunk_number = 0
for id_chunk in tqdm(id_chunks):
    print("Current chunk:" + str(chunk_number))
    pairs_to_calc = []
    for a in id_chunk:
        for b in ids:
            if b != a:
                new_pair = a+"\t"+b
                if not new_pair in calculated_pairs:
                    pairs_to_calc.append(new_pair)

    print("Separating batchs of " + str(pair_chunk_len))
    pair_chunks = list(chunks(pairs_to_calc, pair_chunk_len))
    print("Created " + str(len(pair_chunks)) + " batchs.")
    chunk_number += 1

    new_sims = calc_sims(pair_chunks)

    print("Writing current results")
    append_sims(new_sims, output_sims)