import requests
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import time
import obonet

#def my_obo_parser(obo_path):

def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
                yield lst[i:i + n]

def get_gene_annotation(gene_name, 
            base_url = "https://rnacentral.org/api/v1/rna/"):
    url = base_url+gene_name.split("_")[0]+"/go-annotations/"+gene_name.split("_")[1]

    r = requests.get(url, params = {"format": "json"})
    if r.status_code == 200:
        data = r.json()
        associations = []
        for result in data:
            associations.append((result["rna_id"], result["go_term_id"]))
        return gene_name, associations, r.status_code
    else:
        return gene_name, None, r.status_code

def go_from_rnacentral(id_list, tries = 3):
    associations = set()
    while tries > 0:
        print("Trying to retrieve " + str(len(id_list)) + " ids.")
        failed = 0
        to_retrieves = [list(chunk) for chunk in chunks(list(id_list), 100)]
        too_many = False
        for chunk in tqdm(to_retrieves):
            processes = []
            with ThreadPoolExecutor(max_workers=25) as executor:
                for seq_id in id_list:
                    if seq_id in chunk:
                        processes.append(
                            executor.submit(get_gene_annotation, seq_id))

                for task in as_completed(processes):
                    gene_name, new_associations, response = task.result()
                    if new_associations != None:
                        if len(new_associations) > 0:
                            associations.update(new_associations)
                            #print(str(len(new_associations)) + " annotations for " + gene_id)
                        else:
                            print("No annotations for " + gene_name)
                        id_list.remove(gene_name)
                    else:
                        print(gene_name + " failed with " + str(response))
                        if response == 429:
                            too_many = True
                        failed += 1
        
        print("\t" + str(failed) + " failed.")
        if too_many:
            time.sleep(3)
        else:
            tries -= 1
        if failed == 0:
            tries = 0
    return associations

if __name__ == "__main__":
    urs_list_path = sys.argv[1]
    obo_path = sys.argv[2]
    taxid = ""
    if len(sys.argv) == 4:
        taxid = "_"+sys.argv[3]
        print(taxid)
    
    urs = [x.rstrip("\n")+taxid for x in open(urs_list_path, 'r').readlines()]
    go_associations = list(go_from_rnacentral(urs))
    go_associations.sort()

    print("Processing go.obo")
    obo = obonet.read_obo(obo_path)
    obo_nodes = obo.nodes(data=True)
    id_to_aspect = {id_: data.get('namespace') for id_, data in obo_nodes}
    for ID in list(id_to_aspect.keys()):
        if "alt_id" in obo_nodes[ID]:
            for alt_id in obo_nodes[ID]["alt_id"]:
                id_to_aspect[alt_id] = id_to_aspect[ID]
    
    print("Writing associations")
    with open (urs_list_path+".id2go", 'w') as stream:
        for gene, go in tqdm(go_associations):
            stream.write(gene + "\t" + go + "\t" + id_to_aspect[go] + "\n")