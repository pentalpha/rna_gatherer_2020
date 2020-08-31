import pandas as pd
import sys

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import row
from bokeh.plotting import ColumnDataSource
from bokeh.models import HoverTool, Legend, LegendItem
from bokeh.models import Span, Label
from bokeh.models import ContinuousColorMapper, ColorBar
from bokeh.palettes import inferno,Inferno,Cividis,cividis
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import math

def find_aspect(file_name):
    name = file_name.split("/")[-1].upper()
    mf_count = name.count("MF")
    cc_count = name.count("CC")
    bp_count = name.count("BP")
    if mf_count > cc_count and mf_count > bp_count:
        return "Função Molecular"
    elif cc_count > bp_count:
        return "Componente Celular"
    else:
        return "Processo Biológico"

input_df_paths = [sys.argv[1], sys.argv[2], sys.argv[3]]
output = sys.argv[4]

dfs = {find_aspect(input_df): pd.read_csv(input_df,sep=",",index_col=0) for input_df in input_df_paths}

has_path_percs = []
for name, df in dfs.items():
    print(str(df.head()))
    has_path_percs += [int(x) for x in df['confidenceLevel'].tolist()]
min_path_perc = min(has_path_percs)
max_path_perc = max(has_path_percs)
has_path_range = max_path_perc - min_path_perc
print(str(min_path_perc))
print(str(max_path_perc))
print(str(has_path_range))
def convert_perc_to_position(perc, min_perc, perc_range):
    if perc == min_perc:
        return 1
    else:
        relative_perc = perc - min_perc
        percentil = (relative_perc * 255) / perc_range
        return int(percentil)
pallete = 'cividis'
pallete_func = cividis
inferno100 = pallete_func(255)
#print(str(inferno100))
colormap = {str(i): inferno100[convert_perc_to_position(i,min_path_perc,has_path_range)-1] 
            for i in has_path_percs}
print(str(colormap))

f, ax = plt.subplots(4,figsize=(5,14),gridspec_kw={'height_ratios': [10, 10, 10, 0.8]})
current_axis = 0
box = None
y_values = set()
x_values = set()
for name, df in dfs.items():
    axis = ax[current_axis]
    color_series = pd.Series(name='Cor')
    df['Cor'] = df.apply(lambda row: colormap[str(row['confidenceLevel'])],axis=1)

    '''print("Invalid colors:")
    for color_str in df['Cor'].tolist():
        if color_str[0] != '#' or len(color_str) != 7:
            print(color_str)'''

    inferno_color_mapper = ContinuousColorMapper(palette=inferno100, 
                                            low = min_path_perc*100, 
                                            high = max_path_perc*100)

    not_best = df[df['best']=='NOT_BEST']
    best = df[df['best']=='BEST']

    size_scalar = [40 if best_str == "BEST" else 10 for best_str in df['best'].tolist()]
    xs = []
    ys = []
    colors = []

    xs_best = []
    ys_best = []
    colors_best = []

    xs_others = []
    ys_others = []
    colors_others = []

    print(name)
    axis.set_title(name)
    for index, row in df.iterrows():
        x_value = row['Completude']
        y_value = row['hasPath']
        color_value = colormap[str(row['confidenceLevel'])]
        if row['best'] != "BEST":
            if row['confidenceLevel'] >= 0:
                xs.append(x_value)
                ys.append(y_value)
                colors.append(color_value)
            else:
                xs_others.append(x_value)
                ys_others.append(y_value)
                colors_others.append(color_value)
        else:
            print("\t".join([str(x) for x in 
                [x_value, y_value, row['confidenceLevel'], row['metrics'], row['Cor']]]))
            xs_best.append(x_value)
            ys_best.append(y_value)
            colors_best.append(color_value)
        y_values.add(y_value)
        x_values.add(x_value)

    #f, ax = plt.subplots(2,figsize=(6, 6))

    xs.reverse()
    ys.reverse()
    colors.reverse()
    axis.scatter(xs, ys, s=100, c=colors, alpha=0.75, marker='o', label = "Predição")
    
    '''axis.scatter(xs_best, ys_best, s=200, c='black', alpha=1.0, 
        marker='D')'''
    axis.scatter(xs_best, ys_best, s=100, c=colors_best, alpha=0.9, 
        marker='D', label = 'Predição de Melhor Qualidade', edgecolors= ["#000000"]*len(xs_best))
    
    '''axis.scatter(xs_others, ys_others, s=200, c='black', alpha=1.0, 
        marker='X')'''
    axis.scatter(xs_others, ys_others, s=100, c=colors_others, alpha=1.0, 
        marker='X', label = 'Outras Ferramentas', edgecolors= ["#000000"]*len(xs_others))
    
    #ax.colorbar()
    axis.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    axis.xaxis.set_major_locator(ticker.MultipleLocator(10))
    axis.set_axisbelow(True)
    axis.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    axis.yaxis.set_major_locator(ticker.MultipleLocator(5))

    axis.grid(axis='x')
    axis.grid(axis='y')
    axis.set_ylabel('Associações relacionadas à uma referência (%)')
    axis.set_xlabel('Referências presentes na predição (%)')

    #box = axis.get_position()
    
    #axis.set_position([box.x0+0.015, box.y0-0.33, 0.8, 0.75])
    box = axis.get_position()
    
    
    #axis.legend(loc='upper right', labelspacing=1.2, shadow=True, borderpad=1.2)
    #bbox_to_anchor=(1, 0.5)
    #plt.savefig(name.split("/")[-1].split("-")[0] + "-predictions_comparison.png")

    current_axis+=1

y_min = math.floor(min(y_values)/20)*20
x_max = math.ceil(max(x_values)/20)*20

for axis in [ax[0],ax[1],ax[2]]:
    axis.set_ylim(y_min,101)
    axis.set_xlim(0,x_max)

cmap = plt.cm.get_cmap(pallete)
colors = cmap(np.arange(cmap.N))
ax[current_axis].imshow([colors], extent=[min_path_perc, max_path_perc, 0, 0.3])
plt.setp(ax[current_axis].get_yticklabels(), visible=False)
ax[current_axis].tick_params(axis='y', which='both', length=0)
ax[current_axis].xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax[current_axis].xaxis.set_major_locator(ticker.MultipleLocator(1))
box1 = ax[current_axis].get_position()
ax[current_axis].set_position([box.x0, box.y0-0.21, 0.8, 0.2])
ax[current_axis].set_xlabel('Nível de Confiança')

#circle = mlines.Line2D([], [], color='blue', marker='*',
#                          markersize=15, label='Blue stars')
texts = ['Predição', 'Outras Ferramentas', 'Predição de Melhor Qualidade']
markers = ['o','D','X']
edgec = ["#DADAFF","#000000","#000000"]

patches = [plt.scatter([],[], marker=markers[i], color="#DADAFF", s=100,
            label="{:s}".format(texts[i]), edgecolors=edgec[i])  for i in range(len(texts)) ]
legend = ax[current_axis].legend(handles=patches, ncol=1, facecolor="white", 
            numpoints=1, #bbox_to_anchor=(0., 1.02, 1., .102), 
            loc='lower center', bbox_to_anchor=(0.0, -5.5),
            mode="expand", borderaxespad=0.)

#plt.legend(handles=[blue_line])

f.tight_layout()
f.savefig(output)