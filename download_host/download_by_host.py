import os
import sys
id_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../extract_ids'))
sys.path.append(id_path)
import ids_mining


def get_biosamples_info():
    with open('extract_names/species_names.txt','r') as names:
        all_names = [i.rstrip('\n') for i in names.readlines()]
    for name in all_names:
        with open(f'data/ids/{name}/{name}_biosample_ids.txt','r') as biosamples_ids:
            ids = biosamples_ids.readlines()
            ids = ids[1:]
            for id in ids:
                biosample = id.rstrip('\n')
                info = os.system(f'efetch -db biosample -id {biosample} -format docsum | xtract -pattern DocumentSummary -element Attribute > download_host/{biosample}_info_temp.txt')
                title = os.system(f'efetch -db biosample -id {biosample} -format docsum | xtract -pattern DocumentSummary -element Attribute@attribute_name > download_host/{biosample}_title_temp.txt')
                with open(f'download_host/{biosample}_title_temp.txt','r') as titles:
                    with open(f'download_host/{biosample}_info_temp.txt','r') as infos:
                        all_titles = [i.strip("\n") for i in titles.readline().split(sep="\t")]
                        all_infos = [i.strip("\n") for i in infos.readline().split(sep="\t")]
                        biosample_info = {}
                        for title,info in zip(all_titles,all_infos):
                            biosample_info[title] = info
                        os.remove(f'download_host/{biosample}_title_temp.txt')
                        os.remove(f'download_host/{biosample}_info_temp.txt')
                print(biosample_info)

                    


                
get_biosamples_info()


