import sys
from random import seed
from random import random
import os
from tqdm import tqdm

print("Starting up")
seed()
gene_names = sys.argv[1]
files = sys.argv[2].split(",")
chance = float(sys.argv[3])
output_dir = sys.argv[4]

new_files = {f: output_dir + "/" + f.split("/")[-1] for f in files}

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

names = set()
with open(gene_names, 'r') as in_stream:
    for line in in_stream:
        names.add(line.split("\t")[-1].rstrip("\n"))

print("Creating pairs")
allowed_pairs = set()
total = 0
for i in tqdm(names):
    for j in names:
        if random() <= chance:
            allowed_pairs.add(i+"\t"+j)
        total += 1

print("Allowed " + str(len(allowed_pairs)) + " out of " + str(total) + " possible pairs.")

print("Writing subsets")
for in_path, out_path in tqdm(new_files.items()):
    print(in_path + " -> " + out_path)
    with open(in_path, 'r') as in_stream:
        first = True
        with open(out_path, 'w') as out_stream:
            for line in in_stream:
                if not first:
                    cells = line.rstrip("\n").split("\t")
                    pair_name = cells[0]+"\t"+cells[1]
                    if pair_name in allowed_pairs:
                        out_stream.write(line)
                else:
                    out_stream.write(line)
                    first = False