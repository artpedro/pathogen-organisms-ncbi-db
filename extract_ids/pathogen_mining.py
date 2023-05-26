import json as js
import os
from tqdm import tqdm
import requests
import urllib
import pandas as pd
import sys
from tqdm import tqdm


def download(url: str, fname: str):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    # Can also replace 'file' with a io.BytesIO object
    with open(fname, 'wb') as file, tqdm(
        desc=fname,
        total=total,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in resp.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)

def download_url(url, file):
    with open(file, "wb") as f:
        print("Downloading %s" % file)
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )    
                sys.stdout.flush()


'''
class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlopen(url, filename=output_path, reporthook=t.update_to)
'''
class Group():
    def __init__(self,name):
        # nome do grupo
        self.name = name
        
        self.info_path = f"data/groups_info/{self.name}"
        self.tsv_path = self.info_path + '/tsv'
        self.filtered_json_path = self.info_path + f'/{self.name}_suitable.json'

        # buscando o id do grupo
        with open('data/groups/groups_name_id.json') as log:
            data = js.load(log)
            self.pat_id = data[self.name]
        self.tsv_file = f'{self.tsv_path}/{self.name}|{self.pat_id}.tsv'
        
        if os.path.exists(self.info_path + f'/{self.name}_suitable.json'):
            with open(self.filtered_json_path,'r') as log:
                self.filtered_json = js.load(log)
                

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
        try:
            download(self.table_url,self.tsv_file)
        except:
            download_url(self.table_url,self.tsv_file)

    def getFilteredTsv(self):
        '''
        Lê os arquivos tabulares de cada grupo patogênico e gera um .json com as informações
        importantes dos registros que possuem genoma completo e informação sobre o hospedeiro
        '''

        # campos importantes
        fields = ['scientific_name','strain','host','asm_level','asm_acc','biosample_acc']
        
        # lendo arquivo .tsv em um dataframe
        df = pd.read_csv(self.tsv_file,sep='\t',usecols=fields)
        
        filtered_df = pd.DataFrame(columns=fields)
        
        # filtrando as linhas dispensáveis do dataframe
        for i,row in df.iterrows():
            if row[fields[3]] == 'Complete Genome':
                filtered_df.loc[len(filtered_df)] = row
        self.filtered_df = filtered_df[filtered_df[fields[2]].notnull()]
        
    
        # armazenando o dataframe filtrado em um .json 
        with open(self.filtered_json_path, 'w') as file:
            json_file = self.filtered_df.to_json(orient="records",indent=1)
            file.write(json_file)

    def readFilteredTsv(self):
        '''
        Extrai informações do .tsv filtrado
        '''
        # informações coletadas
        self.species = []
        self.hosts = []
        self.count = 0
        
        if hasattr(self,'filtered_df'):
            for i,row in self.filtered_df.iterrows():
                if row['host'] not in self.hosts:
                    self.hosts.append(row['host'])
                if row['scientific_name'] not in self.species:
                    self.species.append(row['scientific_name'])
            self.count = len(self.filtered_df.index)
            return True
        elif hasattr(self,'filtered_json'):
            for i in self.filtered_json:
                if i['host'] not in self.hosts:
                    self.hosts.append(i['host'])
                if i['scientific_name'] not in self.species:
                    self.species.append(i['scientific_name'])
            return True
        else:
            print('Sem informações salvas')
            return False

        
        


with open("data/groups/groups_name_id.json","r") as file:
    names = js.load(file)
    total_count = 0
    all_count = {}
    all_hosts = set()

    for name in names.keys():
        a = Group(name)
        if a.checkDataUpdate():       
            a.getPatData()
        if not a.readFilteredTsv():
            a.getFilteredTsv()
            a.readFilteredTsv()    
        total_count = total_count + a.count
        all_count[a.name] = [a.count, a.species, a.hosts]
        all_hosts.update(a.hosts)
    for i in all_count:
        print(f'{i}: {all_count[i][0]} entradas\n\nEspécies presentes: {all_count[i][1]}\n\nHospedeiros presentes: {all_count[i][2]}\n\n-----------------------------------------------------------------\n')
    print(f'total: {total_count} entradas')

    print(all_hosts)
    print(len(all_hosts))

