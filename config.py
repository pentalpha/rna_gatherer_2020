#Module: config.py
#Rename this file as "config.py" and customize it
import multiprocessing
import os
import json
import psutil
from nexus.util import runCommand

project_dir = os.path.dirname(os.path.realpath(__file__))
config_json = project_dir + "/config.json"
if not os.path.exists(config_json):
    print("config.json file not found, creating a new one.")
    runCommand("cp " + project_dir+"/config.dummy.json " + config_json)

config_str = open(config_json, 'r').read()
config_str = config_str.replace("PROJECT_DIR", project_dir)
configs = json.loads(config_str)
if not "threads" in configs:
    configs["threads"] = str(max(2, multiprocessing.cpu_count()-1))
if not "max_mem" in configs:
    configs["max_mem"] = str(int((psutil.virtual_memory().total/1024/1024)*0.8))

print("Using " + configs["threads"] + " threads and " 
    + configs["max_mem"] + " MB of RAM.")

def missing_from_config(possibly_missing):
    missing = []
    for x in possibly_missing:
        if not x in configs:
            missing.append(x)
        else:
            if configs[x] == "":
                missing.append(x)
            elif x == "rna_dbs":
                if len(configs[x]) == 0:
                    missing.append(x)
    return missing

def require_files(possibly_missing, mandatory = True):
    missing = missing_from_config(possibly_missing)
    if len(missing) > 0:
        if mandatory:
            print("The following mandatory files are missing: " + str(missing))
            print("Insert their locations in the configuration file:\n\t" + config_json)
            quit()
        else:
            print("The following optional files are missing: " + str(missing))
            print("If you consider their data necessary, insert their locations in the " 
                + "configuration file:\n\t" + config_json)