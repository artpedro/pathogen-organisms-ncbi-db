import json
import os
import requests
from tqdm import tqdm
from pathlib import Path
import urllib

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
        
        # url para extrair as informações em tsv
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'

        # criar repositório para os dados
        if not os.path.exists(self.info_path):
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
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'
        download_url(self.table_url,f'{self.tsv_path}/{self.name}|{self.pat_id}.tsv')

with open("data/groups/groups_name_id.json","r") as file:
    names = json.load(file)
    for name in names.keys():
        a = Group(name)
        if a.checkDataUpdate():
            a.getPatData()          