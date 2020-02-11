import subprocess
import os
from tqdm import tqdm

def runCommand(cmd, print_cmd=True):
        if print_cmd:
                print("\t> " + cmd)
        process = subprocess.call(cmd, shell=True)
        return process

def sliceString(string, max_len):
        slices = []
        for x in range(0,int(len(string)/max_len)+1):
                slices.append(string[x*max_len:(x+1)*max_len])
        return slices

def write_file(content, path):
        with open(path, 'w') as stream:
                stream.write(content)

def getFilesWith(directory, name_part):
        files = []
        # r=root, d=directories, f = files
        for r, d, f in os.walk(directory):
                for file in f:
                        if name_part in file:
                                files.append(os.path.join(r, file))
        return files

def rm_last_part(path, c):
        parts = path.split(c)
        if len(parts) > 1:
                return c.join(parts[:-1])
        else:
                return parts[0]

def replace_last(original, to_repl, repl):
        k = original.rfind(to_repl)
        new_string = original[:k] + repl + original[k+len(to_repl):]
        return new_string

def file_name(path):
        no_file_extension = rm_last_part(path, ".")
        name = no_file_extension.split("/")[-1]
        return name

def get_subdirs(dir):
    subs = os.listdir(dir)
    subdirs = []
    for sub in subs:
        path = dir + "/" + sub
        if os.path.isdir(path):
            subdirs.append(path)
    return subdirs

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_iterator(my_list, show = False):
    if show:
        return tqdm(range(len(my_list)))
    else:
        return range(len(my_list))