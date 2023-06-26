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

def plotSingleDistribution(name='',key='species',all=False):
    print(name)
    metadata = f'data/groups_info/{name}/{name}_metadata'
    if all:
        metadata = f'data/metadata.json'
    # Verificando se existem dados válidos
    if os.path.exists(metadata):
    
        # Acessando as informações
        with open(metadata,'r') as file:    
            data = js.load(file)
            species = data[key]
            if len(species) == 0:
                return False
            # Filtrando valores muito baixos
            values = list(species.values())
            sum_data = sum(values)
            filtered_species = {k:v for k,v in species.items() if v > sum_data/100}

            # Armazenando a informação em ordem crescente
            sorted_species = sorted(filtered_species.items(), key=lambda x: x[1])

            # Separando os nomes e as contagens
            names = [item[0].replace(' ','\n') for item in sorted_species]
            count = [item[1] for item in sorted_species]
            

            # Normalizar os valores de contagem entre 0 e 1
            normalized_counts = np.array(count) / max(count)

            # Criação da paleta de cores proporcionais aos valores
            colors = cm.Blues(normalized_counts)
            
            # Iniciando a imagem
            fig = plt.figure()

            # Iniciando gráfico
            plt.barh(names,count, color=colors, edgecolor='black')
            
            # Mostrando tamanho do lado das barras
            for i in range(len(count)):
                plt.text(count[i], i, " "+str(count[i]), ha='left', va='center')

            # Gerando grade
            plt.grid(True, axis='x', linestyle='--', alpha=0.9)

            plt.gca().set_axisbelow(True)
            
            # Limites do gráfico
            plt.xlim(0, max(count))
            
            # Pontos para serem plotados em X
            x_points = []
            
            # Selecionando pontos para adicionar baseado na escala
            if max(count) <= 25:
                for i in range(0,max(count),5):
                    x_points.append(i)
            if max(count) <= 50 and max(count)>25:
                for i in range(0,max(count),10):
                    x_points.append(i)
            if max(count) > 50 and max(count) < 100:
                for i in range(0,max(count),15):
                    x_points.append(i)
            if max(count) >= 100 and max(count) <= 500:
                for i in range(0,max(count),25):
                    x_points.append(i)
            if max(count) >= 500 and max(count) <1000:
                for i in range(0,max(count),100):
                    x_points.append(i)
            if max(count) >= 1000:
                for i in range(0,max(count),500):
                    x_points.append(i)
            '''
            # Elimina pontos muito próximos no eixo x
            x_points = sorted(x_points)
            print(x_points)
            last_x = max(x_points)
            thrs_points = []

            for idx,x in enumerate(x_points):
                if idx == 0:
                    thrs_points.append(x)
                else:
                    if last_x > 500:
                        if x_points[idx-1] + 20 >= x:
                            continue
                        else:
                            thrs_points.append(x)
                    
                    elif last_x > 100:
                        if x_points[idx-1] + 5 >= x:
                            continue
                        else:
                            thrs_points.append(x)
                    else:
                        if x_points[idx-1] + 3 >=x:
                            continue
                        else:
                            thrs_points.append(x)

            x_points = [i for i in sorted(set(thrs_points))]
            '''
            plt.xticks(x_points,x_points)
            
            # Labels
            plt.xlabel('Count')
            plt.ylabel(key)
            plt.title(f'{key} distribution')
        return fig        
    else:
        print('Sem metadata disponível')
        return None

def writeAll():
    for group in readGroupsNames():
        species = [plotSingleDistribution(group,'species')]
        hosts = [plotSingleDistribution(group,'hosts')]

        if species:
            write_pdf(f"sdis_{group}.pdf",species,path=f'data/groups_info/{group}')
        if hosts:    
            write_pdf(f"hdis_{group}.pdf",hosts,path=f'data/groups_info/{group}')

def writeSingle(group):
    species = [plotSingleDistribution(group,'species')]
    hosts = [plotSingleDistribution(group,'hosts')]
    print(hosts)
    if species:
        write_pdf(f"sdis_{group}.pdf",species,path=f'data/groups_info/{group}')
    if hosts:    
        write_pdf(f"hdis_{group}.pdf",hosts,path=f'data/groups_info/{group}')

def writeGeneral():
    species = [plotSingleDistribution(all=True)]
    host = [plotSingleDistribution(key='hosts',all=True)]
    if species:
        write_pdf(f"sdis_all.pdf",species,path="data/")
    if host:    
        write_pdf(f"hdis_all.pdf",host,path="data/")
writeGeneral()
writeAll()