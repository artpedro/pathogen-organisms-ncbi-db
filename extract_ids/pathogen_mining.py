import json
import os
from tqdm import tqdm
import urllib
import pandas as pd

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)

class Group():
    def __init__(self,name):
        # nome do grupo
        self.name = name
        
        self.info_path = f"data/groups_info/{self.name}"
        self.tsv_path = self.info_path + '/tsv'

        # buscando o id do grupo
        with open('data/groups/groups_name_id.json') as log:
            data = json.load(log)
            self.pat_id = data[self.name]
        self.tsv_file = f'{self.tsv_path}/{self.name}|{self.pat_id}.tsv'

        # url para extrair as informações em tsv
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'

        # criar repositório para os dados
        if not os.path.exists(self.tsv_path):
            os.mkdir(self.info_path)
            os.mkdir(self.tsv_path)
    
    def checkDataUpdate(self):
        '''
        Verifica se a versão do .tsv salvo precisa ser atualizada para uma nova
        True: Precisa ser atualizada
        False: Não precisa ser atualizada
        '''
        tsv_log = os.listdir(self.tsv_path)
        if tsv_log == []:
            return True
        tsv_log = tsv_log[0]
        tsv_version = tsv_log.split('|')[1].rstrip('.tsv')
        if tsv_version != self.pat_id:
            os.remove(self.tsv_path+'/'+tsv_log)
            return True
        else:
            return False

    def getPatData(self):
        download_url(self.table_url,self.tsv_file)

    def filterTsv(self):
        fields = ['scientific_name','strain','host','asm_level','asm_acc','biosample_acc']
        filtered_df = pd.DataFrame(columns=fields)
        # fix read
        df = pd.read_csv(self.tsv_file,header=1)
        print(df)


a = Group('Edwardsiella_tarda')
a.filterTsv()

'''
with open("data/groups/groups_name_id.json","r") as file:
    names = json.load(file)
    for name in names.keys():
        a = Group(name)
        if a.checkDataUpdate():
            a.getPatData()          '''