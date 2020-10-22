from nexus.util import get_subdirs, file_name
import os

def load_confidence(intervals_file):
    #print("Reading " + intervals_file)
    with open(intervals_file, 'r') as stream:
        lines = [line.rstrip("\n").split("\t") for line in stream.readlines()]
        th_lines = lines[2:]
        confidences = [{} for i in range(len(th_lines[0])-1)]
        for cells in th_lines:
            metric = cells[0]
            for i in range(1,len(cells)):
                if cells[i] == "None":
                    confidences[i-1][metric] = None
                else:
                    confidences[i-1][metric] = float(cells[i])
        #for i in range(len(confidences)):
        #    print("Confidence "+ str(i)+": " + str(confidences[i]))
        return confidences
    return []

def get_available_species():
    confs_dir = (os.path.dirname(os.path.realpath(__file__)) 
                    + "/../data/confidence_levels/")
    species_dirs = get_subdirs(confs_dir)
    return [file_name(f) for f in species_dirs]

def load_confidence_levels(species):
    confs_dir = (os.path.dirname(os.path.realpath(__file__)) 
                    + "/../data/confidence_levels/")
    species_dir = confs_dir+"/"+species
    confs = {"MF": load_confidence(species_dir+"/confidence_intervals-MF.tsv"),
            "BP": load_confidence(species_dir+"/confidence_intervals-BP.tsv"),
            "CC": load_confidence(species_dir+"/confidence_intervals-CC.tsv")}
    return confs

def geometric_filter(value, th):
    if value <= th:
        return value
    else:
        return None

def normal_filter(value, th):
    if value >= th:
        return value
    else:
        return None

def geometric_pass(value, th):
    return value <= th

def normal_pass(value, th):
    return value >= th

def compare_to_th(value, th, metric):
    if metric != "SOB" and metric != "FSH":
        return normal_pass(value, th)
    else:
        return geometric_pass(value, th)
