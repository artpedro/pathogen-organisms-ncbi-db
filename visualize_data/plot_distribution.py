import os
import pandas as pd
import json as js
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
# species distribution

def readGroupsNames():
    with open(os.path.normpath("data/groups/groups_name_id.json"),"r") as data:
        groups = js.load(data)
        return groups

def write_pdf(fname, figures, resolution=(1366, 768),path='.'):
    file_path = os.path.normpath(path + '/' + fname)
    doc = PdfPages(file_path)
    for fig in figures:
        fig_width, fig_height = fig.get_size_inches()
        dpi = fig.get_dpi()

        target_width, target_height = resolution
        width_ratio = target_width / (fig_width * dpi)
        height_ratio = target_height / (fig_height * dpi)
        ratio = min(width_ratio, height_ratio)

        new_fig_width = fig_width * ratio
        new_fig_height = fig_height * ratio
        fig.set_size_inches(new_fig_width, new_fig_height)

        fig.savefig(doc, format='pdf', dpi=dpi)
    doc.close()

def plotSingleDistribution(name,key='species'):
    metadata = f'data/groups_info/{name}/{name}_metadata'
    print(name)
    if os.path.exists(metadata):
    
        with open(metadata,'r') as file:    
            data = js.load(file)
            species = data[key]

            # Armazenando a informação em ordem crescente
            sorted_species = sorted(species.items(), key=lambda x: x[1])
            names = [item[0].replace(' ','\n') for item in sorted_species]
            count = [item[1] for item in sorted_species]
            for idx,x in enumerate(count):
                if sum(count)/100 > x:
                    del names[idx]
                    del count[idx]
                    

            # Normalizar os valores de contagem entre 0 e 1
            normalized_counts = np.array(count) / max(count)

            # Criação do gráfico de barras verticais com cores proporcionais aos valores
            colors = cm.Blues(normalized_counts)
            
            fig = plt.figure()
            plt.barh(names,count, color=colors, edgecolor='black')
            for i in range(len(count)):
                plt.text(count[i], i, str(count[i]), ha='left', va='center')
            plt.grid(True, axis='x', linestyle='--', alpha=0.9)

            plt.gca().set_axisbelow(True)
            
            plt.xlim(0, max(count))
            
            x_points = [i for i in set(count)]
            
            if max(count) <= 50:
                for i in range(0,max(count),10):
                    x_points.append(i)
            if max(count) > 50 and max(count) < 100:
                for i in range(0,max(count),15):
                    x_points.append(i)
            if max(count) >= 100 and max(count) <= 500:
                for i in range(0,max(count),25):
                    x_points.append(i)
            if max(count) >= 500:
                for i in range(0,max(count),100):
                    x_points.append(i)
            
            x_points = sorted(x_points)

            thrs_points = []
            for idx,x in enumerate(x_points):
                if idx == 0:
                    thrs_points.append(x)
                else:
                    if x_points[idx-1] + 5 >= x:
                        continue
                    else:
                        thrs_points.append(x)
        
            x_points = [i for i in sorted(set(thrs_points))]
            
            plt.xticks(x_points,x_points)
            plt.xlabel('Count')
            plt.ylabel(key)
            plt.title(f'{key} distribution')
        return fig        
    else:
        print('Sem metadata disponível')
        return None
for group in readGroupsNames():
    species = [plotSingleDistribution(group,'species')]
    hosts = [plotSingleDistribution(group,'hosts')]
    if species[0]:
        write_pdf(f"sdis_{group}.pdf",[plotSingleDistribution(group,'species')],path=f'data/groups_info/{group}')
    if hosts[0]:    
        write_pdf(f"hdis_{group}.pdf",[plotSingleDistribution(group,'hosts')],path=f'data/groups_info/{group}')
'''
# host distribution
with open('data/groups_info/Bacillus_cereus_group/Bacillus_cereus_group_metadata','r') as file:
    a = js.load(file)
    species = a['hosts']
    sorted_data = sorted(species.items(), key=lambda x: x[1])
    print(sorted_data)
    names = [item[0] for item in sorted_data]
    count = [item[1] for item in sorted_data]
    # Normalizar os valores de contagem entre 0 e 1
    normalized_counts = np.array(count) / max(count)

    # Criação do gráfico de barras verticais com cores proporcionais aos valores
    colors = cm.Blues(normalized_counts)  # Utiliza o color map 'viridis' para definir as cores
    
    plt.barh(names,count, color=colors, edgecolor='black')
    plt.grid(True, axis='x', linestyle='--', alpha=0.9)
    plt.gca().set_axisbelow(True)
    plt.xlim(0, max(count))
    x_points = [i for i in set(count)]
    step = int(max(count)/10)
    print(step)
    for i in range(0,max(count),step):
        x_points.append(i)
    plt.xticks(x_points,x_points)
    plt.xlabel('Count')
    plt.ylabel('Species')
    plt.title('Species distribution')
    # Displaying the plot
    plt.show()
'''

