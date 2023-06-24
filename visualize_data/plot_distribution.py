import pandas as pd
import json as js
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

'''
# species distribution
with open('data/groups_info/Bacillus_cereus_group/Bacillus_cereus_group_metadata','r') as file:
    a = js.load(file)
    species = a['species']
    
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
    for i in range(0,max(count),5):
        x_points.append(i)
    plt.xticks(x_points,x_points)
    plt.xlabel('Count')
    plt.ylabel('Species')
    plt.title('Species distribution')
    


    # Displaying the plot
    #plt.show()
'''

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