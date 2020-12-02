from matplotlib.cm import get_cmap
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator
import random
import sys

translate_dict = {"cc": "Componente Celular", "CC": "Componente Celular",
                "mf": "Função Molecular", "MF": "Função Molecular",
                "bp": "Processo Biologico", "BP": "Processo Biologico",
                "max": "Similaridade Semântica Máxima", 
                "avg": "Similaridade Semântiac Média",
                "fsh": "Fisher", "FSH": "Fisher",
                "sob": "Sobolev", "SOB": "Sobolev",
                "mic": "Maximal Information\nCoefficient", 
                "MIC": "Maximal Information\nCoefficient",
                "prs": "Pearson", "PRS": "Pearson",
                "spr": "Spearman", "SPR": "Spearman",
                "dc": "Distance\nCorrelation", "DC": "Distance\nCorrelation"}

translate_dict_eng = {"cc": "Cellular Component", "CC": "Cellular Component",
                "mf": "Molecular Function", "MF": "Molecular Function",
                "bp": "Biological Process", "BP": "Biological Process",
                "max": "Maximum Semantic Similarity", 
                "avg": "Average Semantic Similarity",
                "fsh": "Fisher", "FSH": "Fisher",
                "sob": "Sobolev", "SOB": "Sobolev",
                "mic": "Maximal Information\nCoefficient", 
                "MIC": "Maximal Information\nCoefficient",
                "prs": "Pearson", "PRS": "Pearson",
                "spr": "Spearman", "SPR": "Spearman",
                "dc": "Distance\nCorrelation", "DC": "Distance\nCorrelation"}

def translator(name):
    name.replace("_"," ")
    words = name.split()
    new_words = [(translate_dict_eng[word] if word in translate_dict_eng else word)
                    for word in words]
    new_name = " ".join(new_words)
    return new_name

def get_col(lines, col_i, conversor=None):
    col = [line[col_i] for line in lines]
    if conversor != None:
        col = [conversor(x) for x in col]
    return col

def comb_sort(comb_str):
    parts = comb_str.split('-')
    parts.sort()
    return '-'.join(parts)

def read_df(df_path):
    lines = [line.rstrip("\n").split("\t")[1:] for line in open(df_path,'r').readlines()[1:]]
    lines.sort(key=lambda x: int(x[0]))
    confs_list = get_col(lines, 0, int)
    comb_list = [comb_sort(x) for x in get_col(lines, 1)]
    confs_list = get_col(lines, 0, int)
    completion_list = get_col(lines, 4, float)
    path_perc_list = get_col(lines, 5, float)
    
    return confs_list, comb_list, completion_list, path_perc_list


confs_dir = sys.argv[1]
#confs_dir = "/home/pitagoras/main/experiments/predictions/mgi_simgic_tpm-exon"
ontos = ["MF", "BP", "CC"]
colors_a = [get_cmap('tab20').colors[((i+1)*2)-1] for i in range(10)]
colors_b = [get_cmap('tab20').colors[i*2] for i in range(10)]
colors = colors_a + colors_b
random.shuffle(colors)
best_dfs = [confs_dir+"/"+ont+"-best_predictions.tsv" for ont in ontos]
dfs = [read_df(p) for p in best_dfs]
all_combs = set()

for confs_list, comb_list, completion_list, path_perc_list in dfs:
    all_combs.update(comb_list)
comb_freq = {comb_name: 0 for comb_name in all_combs}
for confs_list, comb_list, completion_list, path_perc_list in dfs:
    for comb in comb_list:
        comb_freq[comb] += 1
all_combs = list(all_combs)
all_combs.sort()
comb_colors = {all_combs[i]: colors[i] for i in range(len(all_combs))}
all_combs.sort(key=lambda x: comb_freq[x], reverse=True)
output_plot = confs_dir + "/bests_plot.png"
#%%
fig, axes = plt.subplots(3, 1, figsize=(7, 8))
y_start = 0.0
for onto_name, df, ax in zip(ontos, dfs, axes):
    confs_list, comb_list, completion_list, path_perc_list = df
    vertical_lines = []
    vertical_sections = []
    section_comb = []
    last_border = -5
    for i in range(1,len(comb_list)):
        if comb_list[i] != comb_list[i-1]:
            border = (confs_list[i]+confs_list[i-1])/2
            vertical_lines.append(border)
            vertical_sections.append((last_border, border))
            section_comb.append(comb_list[i-1])
            last_border = border
    vertical_sections.append((last_border, confs_list[-1]+10))
    section_comb.append(comb_list[-1])
    point_colors = [comb_colors[comb] for comb in comb_list]
    for sect, comb in zip(vertical_sections, section_comb):
        x1, x2 = sect
        rect = patches.Rectangle((x1,y_start),x2-x1,103,
                                 linewidth=0,
                                 facecolor=comb_colors[comb],
                                 fill=True)
        ax.add_patch(rect)
    #ax.vlines(vertical_lines, 0.0, 103.0)
    ax.plot(confs_list, path_perc_list, 
            color='red', linewidth=1, marker='o', alpha=0.75)
    '''ax.plot(confs_list, completion_list, 
            color='blue', linewidth=1, marker='o', alpha=0.4)'''
    
    for sect, comb in zip(vertical_sections, section_comb):
        x1, x2 = sect
        ax.annotate(comb, 
                    ((max(0.0,x1)+min(x2,confs_list[-1]))/2, 
                     y_start+3),
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    rotation='vertical',
                    fontsize='medium')
    
    ax.set_ylim(y_start,103.0)
    ax.set_xlim(-0.5,confs_list[-1]+0.5)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.set_title(translator(onto_name))

#axes[-2].set_ylabel("Associações relacionadas\nà uma referência (%)")
#axes[-1].set_xlabel("Nível de Confiança")
axes[-2].set_ylabel("Associations related\nto a reference (Q2)")
axes[-1].set_xlabel("Confidence Level")

'''leg_patches = [plt.plot([],[], marker="s", ms=10, ls="", 
                     mec=None, color=comb_colors[all_combs[i]], 
                     label="{:s}".format(all_combs[i]))[0]  
           for i in range(len(all_combs))]
axes[-2].legend(handles=leg_patches, bbox_to_anchor=(1.5, 0.5), 
           loc='center', ncol=1, facecolor="plum", numpoints=1,
           fontsize=9)'''

fig.tight_layout()
fig.savefig(confs_dir+"/bests.png", bbox_inches='tight')
#%%