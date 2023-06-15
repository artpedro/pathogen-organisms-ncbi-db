import json as js
import os
from tqdm import tqdm
import requests
import pandas as pd
import sys
from tqdm import tqdm
import time
from groups_name_id import *


def download(url: str, fname: str):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    descr = fname.split('/')[-1]
    # Can also replace 'file' with a io.BytesIO object
    with open(fname, 'wb') as file, tqdm(
        desc=descr,
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

class Group():
    def __init__(self,name):
        # nome do grupo
        self.name = name
        
        self.info_path = os.path.normpath(f"data/groups_info/{self.name}")
        self.tsv_path = os.path.normpath(self.info_path + '/tsv')
        self.usable_json_path = os.path.normpath(self.info_path + f'/{self.name}_usable.json')
        self.filtered_json_path = os.path.normpath(self.info_path + f'/{self.name}_filtered.json')

        # buscando o id do grupo
        with open(os.path.normpath('data/groups/groups_name_id.json')) as log:
            data = js.load(log)
            self.pat_id = data[self.name]
        self.tsv_file = os.path.normpath(f'{self.tsv_path}/{self.name}_{self.pat_id}.tsv')
        
        if os.path.exists(os.path.normpath(self.info_path + f'/{self.name}_filtered.json')):
            with open(self.filtered_json_path,'r') as log:
                self.filtered_json = js.load(log)
                self.filtered_df = pd.DataFrame(self.filtered_json)
            with open(self.usable_json_path,'r') as log:
                self.usable_json = js.load(log)
                self.usable_df = pd.DataFrame(self.usable_json)
                

        # url para extrair as informações em tsv
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'

        # criar repositório para os dados
        if not os.path.exists(self.tsv_path):
            os.makedirs(os.path.normpath(self.info_path))
            os.makedirs(os.path.normpath(self.tsv_path))
    
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
        tsv_version = tsv_log.split('_')[1].rstrip('.tsv')
        if tsv_version != self.pat_id:
            os.remove(os.path.normpath(self.tsv_path+'/'+tsv_log))
            return True
        else:
            return False

    def getPatData(self):
        try:
            download(self.table_url,self.tsv_file)
            self.checkData()        
        except:
            download_url(self.table_url,self.tsv_file)
            self.checkData()

    def checkData(self):
        with open(self.tsv_file,'r') as tsv:
            if tsv.readline()[0] != '#':
                print('Erro no download')
                if hasattr(self,"retry_download"):
                    print(self.retry_download)
                    self.retry_download = self.retry_download + 1
                    if self.retry_download > 2:
                        self.pat_id = getSinglePathogenId(self.name)
                        refreshSingleGroup(self.name)
                        self.tsv_file = os.path.normpath(f'{self.tsv_path}/{self.name}_{self.pat_id}.tsv')
                        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'
                    self.getPatData()
                else:
                    self.retry_download = 0
                    self.getPatData()
            else:
                print("Download bem sucessido")

    def getUsableTsv(self):
        '''
        Lê os arquivos tabulares de cada grupo patogênico e gera um .json com as informações
        importantes dos registros que possuem genoma completo e informação sobre o hospedeiro

        Retorna falso se não houver nenhum
        '''

        # campos importantes
        fields = ['scientific_name','strain','host','asm_level','asm_acc','biosample_acc']
        
        # lendo arquivo .tsv em um dataframe
        df = pd.read_csv(self.tsv_file,sep='\t',usecols=fields)
        
    
        # filtrando as linhas dispensáveis do dataframe
        self.usable_df = df.loc[df['asm_level'] == 'Complete Genome']
        if self.usable_df.empty:
            return False
        self.usable_df = self.usable_df[self.usable_df['host'].notnull()]
        if self.usable_df.empty:
            return False
        self.usable_df = self.usable_df[self.usable_df['asm_acc'].notnull()]
        if self.usable_df.empty:
            return False

        # normalizando a coluna host
        self.usable_df['host'] = self.usable_df['host'].apply(lambda x:x.upper())
        
    
        # armazenando o dataframe com os registros disponíveis antes da filtragem

        with open(self.usable_json_path,'w') as file:
            self.usable_json = self.usable_df.to_json(orient="records",indent=1)
            file.write(self.usable_json)
        return True
    def filterTsv(self):    
        human = 'HOMO|HUMAN|HOMO|SAPIENS'
        chicken = 'GALLUS|CHICKEN'
        bovine = 'TAURUS|BOVINAE|COW|BOVINE'
        fish = 'AQUA|FISH|CARP'

        self.filtered_df = self.usable_df

        self.filtered_df.loc[self.filtered_df['host'].str.contains(human),'host'] = 'HOMO SAPIENS'
        self.filtered_df.loc[self.filtered_df['host'].str.contains(chicken),'host'] = 'CHICKEN'
        self.filtered_df.loc[self.filtered_df['host'].str.contains(bovine),'host'] = 'BOVINE'
        self.filtered_df.loc[self.filtered_df['host'].str.contains(fish),'host'] = 'FISH'
    
        # armazenando o dataframe filtrado em um .json 
        with open(self.filtered_json_path, 'w') as file:
            self.filtered_json = self.filtered_df.to_json(orient="records",indent=1)
            file.write(self.filtered_json)

    def readFilteredTsv(self,filter=True):
        '''
        Extrai informações do .json gerado pelo getUsableTsv() ou filterTsv() (filter = True).
        As informações são redundantes.
        '''
        # informações coletadas
        if filter:
            if hasattr(self,'filtered_df'):
                df = self.filtered_df
        elif hasattr(self,'usable_df'):
                df = self.usable_df
        else:
            print('Sem registros salvos')
            return False
                    
        self.species = df['scientific_name'].tolist()
        self.hosts = df['host'].tolist()
        self.strains = df['strain'].tolist()
        self.count = df.shape[0]

        # salvando contagem dos dados
        self.hosts_dic = {k:self.hosts.count(k) for k in set(self.hosts)}
        self.species_dic = {k:self.species.count(k) for k in set(self.species)}
        self.strains_dic = {k:self.strains.count(k) for k in set(self.strains) if self.strains.count(k) > 1}
        
        # filtragem inicial - terminar
        filtered_species = {k:[" ".join(j) for j in [i.split()[:2] for i in self.species]].count(k) for k in set([" ".join(j) for j in [i.split()[:2] for i in self.species]])}     
        self.species_fdic = filtered_species
        self.subsp_dic = {k:[i.split()[3] for i in self.species if len(i.split())>3].count(k) for k in set([i.split()[3] for i in self.species if len(i.split())>3])}
    def makeGroupMetadata(self):
        '''
        Reune informações do readFilteredTsv() e armazena em um arquivo para posterior representação gráfica
        '''
        metadata = {'group':self.name,'count':self.count,'species':self.species_fdic,'subsp.':self.subsp_dic,'scientific_name':self.species_dic,'strain':self.strains_dic, 'hosts':self.hosts_dic}
        
        with open(os.path.normpath(self.info_path + f'/{self.name}_metadata'),'w') as log:
            content = js.dumps(metadata,indent=1)
            log.write(content)
        


def readGroupsNames():
    with open(os.path.normpath("data/groups/groups_name_id.json"),"r") as data:
        groups = js.load(data)
        return groups
    
def mountExample(name="Edwardsiella_tarda"):
    '''
    Organiza e extrai as informações de somente um grupo patogênico.
    default = "Edwardsiella_tarda" 
    '''
    obj = Group(name)
    obj.getPatData()
    obj.getUsableTsv()
    print(f'Grupo {name} montado')
    return obj

def updateExample(name="Edwardsiella_tarda"):
    '''
    Verifica se é necessário atualizar somente o grupo patogênico em "name"
    e encaminha para atualização.
    default = "Edwardsiella_tarda"
    '''
    obj = Group(name)
    if obj.checkDataUpdate():
        obj.getPatData()
        obj.getUsableTsv()
        print("Informações atualizadas")
    else:
        print("As informações do registro não precisam ser atualizadas")


def mountData():
    '''
    Organiza e extrai as informações de todos os grupos patogênicos.
    '''
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        if obj.checkDataUpdate():    
            obj.getPatData()
        obj.getUsableTsv()


def updateData():
    '''
    Verifica se o banco precisa ser atualizado se baseando nos registros
    do groups_name_id.json e realiza as atualizações necessárias
    '''
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        if obj.checkDataUpdate():
            obj.getPatData()
            obj.getUsableTsv()
        else:
            obj.getUsableTsv()

def readAllData():
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        if obj.getUsableTsv():
            obj.filterTsv()
            obj.readFilteredTsv()
            obj.makeGroupMetadata()
            print(group,' lido')
        
def readSingleData(group="Edwardsiella_tarda"):
    obj = Group(group)
    if obj.getUsableTsv():
        obj.filterTsv()
        obj.readFilteredTsv()
        obj.makeGroupMetadata()
    else:
        print('Sem registros úteis')
    
# Terminar função
def generalMetadata():
    groups = readGroupsNames()
    metadata = {}
    for group in groups:
        with open(os.path.normpath(f'data/groups_info/{group}/{group}_metadata'),'r') as log:
            pass        


'''
start = time.time()
mountData()
end = time.time()
minutos = int((end - start) // 60)
segundos = (end - start) % 60
print(f'mount runtime: {minutos}:{segundos}')
'''
start = time.time()
readAllData()
end = time.time()
minutos = int((end - start) // 60)
segundos = (end - start) % 60
print(f'read runtime: {minutos}:{segundos}')
# readAllData() = 0:50
# new readAllData() = 0:05




'''
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
            a.getUsableTsv()
            a.readFilteredTsv()    
        total_count = total_count + a.count
        all_count[a.name] = [a.count, a.species, a.hosts]
        all_hosts.update(a.hosts)
    for i in all_count:
        print(f'{i}: {all_count[i][0]} entradas\n\nEspécies presentes: {all_count[i][1]}\n\nHospedeiros presentes: {all_count[i][2]}\n\n-----------------------------------------------------------------\n')
    print(f'total: {total_count} entradas')

    print(all_hosts)
    print(len(all_hosts))
'''
