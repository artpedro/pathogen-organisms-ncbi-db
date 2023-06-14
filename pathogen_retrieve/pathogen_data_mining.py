import json as js
import os
from tqdm import tqdm
import requests
import urllib
import pandas as pd
import sys
from tqdm import tqdm
import time


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

class Group():
    def __init__(self,name):
        # nome do grupo
        self.name = name
        
        self.info_path = os.path.normpath(f"data/groups_info/{self.name}")
        self.tsv_path = os.path.normpath(self.info_path + '/tsv')
        self.filtered_json_path = os.path.normpath(self.info_path + f'/{self.name}_suitable.json')

        # buscando o id do grupo
        with open(os.path.normpath('data/groups/groups_name_id.json')) as log:
            data = js.load(log)
            self.pat_id = data[self.name]
        self.tsv_file = os.path.normpath(f'{self.tsv_path}/{self.name}_{self.pat_id}.tsv')
        
        if os.path.exists(os.path.normpath(self.info_path + f'/{self.name}_suitable.json')):
            with open(self.filtered_json_path,'r') as log:
                self.filtered_json = js.load(log)
                

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
        Extrai informações do .json gerado pelo getFilteredTsv() filtrado.
        As informações são redundantes.
        '''
        # informações coletadas
        self.species = []
        self.hosts = []
        self.strains = []
        self.count = 0
        
        # testar isso
        if hasattr(self,'filtered_df'):
            self.hosts = [row['host'] for i,row in self.filtered_df.iterrows()]
            self.species = [row['scientific_name'] for i,row in self.filtered_df.iterrows()]
            self.strains = [row['strain'] for i,row in self.filtered_df.iterrows()]
            self.count = len(self.species)
        elif hasattr(self,'filtered_json'):
            self.hosts = [row['host'] for row in self.filtered_json]
            self.species = [row['scientific_name'] for row in self.filtered_json]
            self.strains = [row['strain'] for row in self.filtered_json]
            self.count = len(self.species)
        elif self.filtered_json:
            self.hosts = [row['host'] for row in self.filtered_json]
            self.species = [row['scientific_name'] for row in self.filtered_json]
            self.strains = [row['strain'] for row in self.filtered_json]
            self.count = len(self.species)         
        else:
            print('Sem informações salvas')
            return False
        # salvando contagem dos dados
        self.hosts_dic = {k:self.hosts.count(k) for k in set(self.hosts)}
        self.species_dic = {k:self.species.count(k) for k in set(self.species)}
        self.strains_dic = {k:self.strains.count(k) for k in set(self.strains) if self.strains.count(k) > 1}
        
        # filtragem inicial - terminar
        filtered_species = {k:[" ".join(j) for j in [i.split()[:2] for i in self.species]].count(k) for k in set([" ".join(j) for j in [i.split()[:2] for i in self.species]])}
                            
        self.species_dic = filtered_species
            
        deleted_hosts = []
        if 'Chicken' not in self.hosts_dic.keys():
            self.hosts_dic['Chicken'] = 0
        if 'Homo sapiens' not in self.hosts_dic.keys():
            self.hosts_dic['Homo sapiens'] = 0
        for host in self.hosts_dic:
            if host == 'Homo sapiens':
                continue
            if host == 'Chicken':
                continue
            elif ('HOMO' in host.upper().split()) or ('HUMAN' in host.upper().split()):
                self.hosts_dic['Homo sapiens'] = self.hosts_dic['Homo sapiens'] + self.hosts_dic[host] 
                deleted_hosts.append(host)
            elif ('GALLUS' in host.upper().split()) or ('HEN' in host.upper().split()):
                if 'Chicken' in self.hosts_dic.keys():
                    self.hosts_dic['Chicken'] = self.hosts_dic['Chicken'] + self.hosts_dic[host]
                deleted_hosts.append(host)
            
        for i in deleted_hosts:
            del self.hosts_dic[i]
        
    def makeGroupMetadata(self):
        '''
        Reune informações do readFilteredTsv() e armazena em um arquivo para posterior representação gráfica
        '''
        metadata = {'group':self.name,'count':self.count,'species':self.species_dic,'strain':self.strains_dic, 'hosts':self.hosts_dic}

        with open(os.path.normpath(self.info_path + f'/{self.name}_metadata'),'w') as log:
            content = js.dumps(metadata,indent=1)
            log.write(content)
        


def readGroupsNames():
    with open(os.path.normpath("data/groups/groups_name_id.json"),"r") as data:
        groups = js.load(data)
        return groups
    
def mountExample(name="Edwardsiellaa_tarda"):
    '''
    Organiza e extrai as informações de somente um grupo patogênico.
    default = "Edwardsiella_tarda" 
    '''
    obj = Group(name)
    obj.getPatData()
    obj.getFilteredTsv()
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
        obj.getFilteredTsv()
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
        obj.getFilteredTsv()

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
            obj.getFilteredTsv()
        else:
            obj.getFilteredTsv()

def readAllData():
    # terminar essa função
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        obj.readFilteredTsv()
        obj.makeGroupMetadata()
        print(group,' lido')
        

        #print(
#            f'''Grupo {group}:
#
#Quantidade de registros: {obj.count} registros
#
#Espécies: {obj.species}
#Tipos: {obj.species_dic}
#
#Hospedeiros: {obj.hosts}
#Tipos: {obj.hosts_dic}
#
#Strains: {obj.strains}
#Tipos: {obj.strains_dic}
#'''
#             )
        
def readSingleData(group="Edwardsiella_tarda"):
    obj = Group(group)
    obj.readFilteredTsv()
    obj.makeGroupMetadata()
    print(
        f'''Grupo {group}:

Quantidade de registros: {obj.count} registros

Espécies: {obj.species}
Tipos: {obj.species_dic}

Hospedeiros: {obj.hosts}
Tipos: {obj.hosts_dic}

Strains: {obj.strains}
Tipos: {obj.strains_dic}
'''
        )
# Terminar função
def generalMetadata():
    groups = readGroupsNames()
    metadata = {}
    for group in groups:
        with open(os.path.normpath(f'data/groups_info/{group}/{group}_metadata'),'r') as log:
            pass        
        
start = time.time()
mountData()
readAllData()
end = time.time()
minutos = int((end - start) // 60)
segundos = (end - start) % 60
print(f'runtime: {minutos}:{segundos}')


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
'''
