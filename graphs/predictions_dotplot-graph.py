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

input_df_paths = [sys.argv[1], sys.argv[2], sys.argv[3]]

dfs = {input_df: pd.read_csv(input_df,sep="\t",index_col=0) for input_df in input_df_paths}

has_path_percs = []
for name, df in dfs.items():
    print(str(df.head()))
    has_path_percs += df['hasPath'].tolist()
min_path_perc = min(has_path_percs)
max_path_perc = max(has_path_percs)
has_path_range = max_path_perc - min_path_perc

def convert_perc_to_position(perc, min_perc, perc_range):
    relative_perc = perc - min_perc
    percentil = (relative_perc * 255) / perc_range
    return int(percentil)
pallete = 'cividis'
pallete_func = cividis
inferno100 = pallete_func(255)
print(str(inferno100))
colormap = {str(i): inferno100[convert_perc_to_position(i,min_path_perc,has_path_range)-1] 
            for i in has_path_percs}

last = 3
for name, df in dfs.items():
    color_series = pd.Series(name='Cor')
    df['Cor'] = df.apply(lambda row: colormap[str(row['hasPath'])],axis=1)

    print("Invalid colors:")
    for color_str in df['Cor'].tolist():
        if color_str[0] != '#' or len(color_str) != 7:
            print(color_str)

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

    for index, row in df.iterrows():
        if row['best'] != "BEST":
            if row['confidenceLevel'] >= 0:
                xs.append(row['Completude'])
                ys.append(row['Incompletude'])
                colors.append(colormap[str(row['hasPath'])])
            else:
                xs_others.append(row['Completude'])
                ys_others.append(row['Incompletude'])
                colors_others.append(colormap[str(row['hasPath'])])
        else:
            xs_best.append(row['Completude'])
            ys_best.append(row['Incompletude'])
            colors_best.append(colormap[str(row['hasPath'])])

    f, ax = plt.subplots(2,figsize=(6, 6))

    ax[0].set_title("Comparação de todas as predições"
            +" com a anotação referência")
    xs.reverse()
    ys.reverse()
    colors.reverse()
    ax[0].scatter(xs, ys, s=100, c=colors, alpha=0.8, marker='o', label = "Predição")
    
    '''ax[0].scatter(xs_best, ys_best, s=200, c='black', alpha=1.0, 
        marker='D')'''
    ax[0].scatter(xs_best, ys_best, s=100, c=colors_best, alpha=0.8, 
        marker='D', label = 'Predição de Melhor Qualidade', edgecolors= ["#000000"]*len(xs_best))
    
    '''ax[0].scatter(xs_others, ys_others, s=200, c='black', alpha=1.0, 
        marker='X')'''
    ax[0].scatter(xs_others, ys_others, s=100, c=colors_others, alpha=1.0, 
        marker='X', label = 'Outras Ferramentas', edgecolors= ["#000000"]*len(xs_others))
    
    #ax.colorbar()
    ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(20))
    ax[0].set_axisbelow(True)
    ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20))

    cmap = plt.cm.get_cmap(pallete)
    colors = cmap(np.arange(cmap.N))
    ax[1].imshow([colors], extent=[min_path_perc, max_path_perc, 0, 0.3])
    plt.setp(ax[1].get_yticklabels(), visible=False)
    ax[1].tick_params(axis='y', which='both', length=0)
    ax[1].xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax[1].xaxis.set_major_locator(ticker.MultipleLocator(1))

    ax[0].grid(axis='x')
    ax[0].grid(axis='y')
    ax[0].set_ylabel('Referências ausentes da predição  (%)')
    ax[0].set_xlabel('Referências presentes na predição (%)')

    box = ax[0].get_position()
    
    ax[0].set_position([box.x0+0.015, box.y0-0.33, 0.8, 0.75])
    box = ax[0].get_position()
    
    box1 = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0-0.21, 0.8, 0.2])
    ax[1].set_xlabel('Associações relacionadas à uma referência (%)')
    ax[0].legend(loc='upper right', labelspacing=1.2, shadow=True, borderpad=1.2)
    #bbox_to_anchor=(1, 0.5)
    plt.savefig(name.split("/")[-1].split("-")[0] + "-predictions_comparison.png")

    last -= 1
