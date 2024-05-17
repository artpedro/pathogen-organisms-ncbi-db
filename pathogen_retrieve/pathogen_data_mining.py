import json as js
import os
from tqdm import tqdm
import requests
import pandas as pd
import sys
from tqdm import tqdm
import time
from groups_name_id import *
import warnings

warnings.filterwarnings('ignore')

# funções utilizadas para baixar os registros do ncbi
def download(url: str, fname: str):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    descr = fname.split('/')[-1]
    
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

# classe principal para extração e manipulação dos grupos patogênicos presentes no banco de dados pathogen
class Group():
    def __init__(self,name,from_log=True):
        # nome do grupo
        self.name = name
        
        # caminho geral para armazenamento das informações
        self.info_path = os.path.normpath(f"data/groups_info/{self.name}")

        # caminho para o arquivo tabular baixado no banco pathogen
        self.tsv_path = os.path.normpath(self.info_path + '/tsv')

        # caminho para o json com os registros de assemblies com
        #   1- informação de host disponível
        #   2- assembly de genoma completo
        #   3- disponibilidade no banco de dados assembly
        
        self.usable_json_path = os.path.normpath(self.info_path + f'/{self.name}_usable.json')
        
        # caminho para o json com registros de assemblies após filtragens
        #   1- agrupamento de hospedeiros semelhantes
        #   2- eliminação de hospedeiros irrelevantes*
        self.filtered_json_path = os.path.normpath(self.info_path + f'/{self.name}_filtered.json')

        # buscando o pat-id do grupo
        with open(os.path.normpath('data/groups/groups_name_id.json')) as log:
            data = js.load(log)
            
            # armazenando o pat-id do grupo
            self.pat_id = data[self.name]

        # caminho para o arquivo tabular extraído do pathogen + nome correto do arquivo
        self.tsv_file = os.path.normpath(f'{self.tsv_path}/{self.name}_{self.pat_id}.tsv')
        
        # procurando se já existem informações filtradas ou semi-filtradas
        # e armazenado-as
        if from_log:    
            if os.path.exists(self.filtered_json_path):
                with open(self.filtered_json_path,'r') as log:
                    self.filtered_json = js.load(log)

                    # informações filtradas em dataframe
                    self.filtered_df = pd.DataFrame(self.filtered_json)
            
            if os.path.exists(self.usable_json_path):
                with open(self.usable_json_path,'r') as log:
                    self.usable_json = js.load(log)

                    # informações brutas em dataframe
                    self.usable_df = pd.DataFrame(self.usable_json)
                    

        # url para extrair as informações tabulares do banco pathogen
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'

        # criar repositório para os dados
        if not os.path.exists(self.info_path):
            os.makedirs(os.path.normpath(self.info_path))
            os.makedirs(os.path.normpath(self.tsv_path))
    
    def checkDataUpdate(self):
        '''
        Encaminha o arquivo para ser atualizado
        True: Precisa ser atualizada
        False: Não precisa ser atualizada
        '''
        # verifica o que está presente na pasta com informações
        usable_log = os.listdir(self.info_path)

        # se só tiver a pasta tsv ou não houver nada, encaminha para a atualização
        if usable_log == ['tsv'] or usable_log == []:
            print(f'{self.name}: Nenhum registro encontrado')
            return True
        
        # se encontrar um usable.json, verifica a tag do arquivo
        # se a tag mostrar uma versão anterior do pathogen, encaminha para a atualização
        
        if f"{self.name}_usable.json" in usable_log:
            with open(self.usable_json_path,'r') as file:
                data = js.load(file)
                if data[0][self.name] == self.pat_id:
                    print(f'{self.name}: Registro já atualizado')
                    return False
                else:
                    return True
        # caso chegue até aqui, também sera encaminhado a atualização
        else:
            return True

    def getPatData(self):
        '''
        Baixa o arquivo tabular contendo todas as informações
        acerca do grupo presentes do ncbi pathogen
        '''
        try:
            if not os.path.exists(self.tsv_path):
                os.makedirs(os.path.normpath(self.tsv_path))
            download(self.table_url,self.tsv_file)
            self.checkData()        
        except:
            if not os.path.exists(self.tsv_path):
                os.makedirs(os.path.normpath(self.tsv_path))
            download_url(self.table_url,self.tsv_file)
            self.checkData()

    def checkData(self):
        '''
        Verifica se o download do arquivo tabular foi realizado com sucesso,
        caso ocorra algum erro, a atualização é forçada
        '''
        # lendo o arquivo tabular
        with open(self.tsv_file,'r') as tsv:
            # geralmente o arquivo começa com '#'
            if tsv.readline()[0] != '#':
                print(f'{self.name}: Erro no download')
                
                # o atributo "retry_download" controla a quantidade de tentativa de downloads
                if hasattr(self,"retry_download"):
                    self.retry_download = self.retry_download + 1
                    if self.retry_download > 2:
                        # atualiza o pat-id e os atributos dependentes dele forçando a atualização
                        self.pat_id = getSinglePathogenId(self.name)
                        new_id = refreshSingleGroup(self.name)
                        self.pat_id = new_id
                        self.tsv_file = os.path.normpath(f'{self.tsv_path}/{self.name}_{self.pat_id}.tsv')
                        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'
                    # tenta baixar de novo
                    self.getPatData()
                else:
                    # inicializa o retry_download e tenta baixar de novo
                    self.retry_download = 0
                    self.getPatData()
            else:
                print(f"{self.name}: Download bem sucessido")

    def getUsableTsv(self):
        '''
        Lê os arquivos tabulares de cada grupo patogênico e gera um .json com as informações
        importantes dos registros que possuem genoma completo e informação sobre o hospedeiro

        Retorna falso se não houver nenhuma
        '''

        # certificando a versão do registro

        # verificação se existe um usable_json desatualizado
        # se estiver desatualizado, é encaminhado para o download e a
        # recursivadade é iniciada
        if os.path.exists(self.usable_json_path):
            if hasattr(self,'iter'):
                print(f'{self.name}: Gerando novo registro')
            else:  
                with open(self.usable_json_path,'r') as file:
                    data = json.load(file)
                    if data[0][self.name] == self.pat_id:
                        print(f'{self.name}: Informações já atualizadas')
                    else:
                        print(f"{self.name}: Informações desatualizadas")
                        self.getPatData()
                        self.iter = 1
                        self.getUsableTsv()
                        print(f'{self.name}: Informações atualizadas')
                        return True
        # caso não exista registros, eles são extraídos
        else:
            self.getPatData()

        # campos importantes no tsv
        fields = ['scientific_name','strain','host','host_disease','asm_level','asm_acc','biosample_acc','collection_date','geo_loc_name']
        # lendo arquivo .tsv em um dataframe
        df = pd.read_csv(self.tsv_file,sep='\t',usecols=fields)

        # filtrando as linhas dispensáveis do dataframe

        # somente genoma completo
        self.usable_df = df.loc[df['asm_level'] == 'Complete Genome']
        if self.usable_df.empty:
            print(f'{self.name}: Sem informações válidas')
            os.remove(self.tsv_file)
            return False
        # somente com informações de host disponíveis
        self.usable_df = self.usable_df[self.usable_df['host'].notnull()]
        if self.usable_df.empty:
            print(f'{self.name}: Sem informações válidas')
            os.remove(self.tsv_file)
            return False
        # somente com assembly id disponível
        self.usable_df = self.usable_df[self.usable_df['asm_acc'].notnull()]
        if self.usable_df.empty:
            print(f'{self.name}: Sem informações válidas')
            os.remove(self.tsv_file)
            return False

        # padronizando a coluna host
        self.usable_df['host'] = self.usable_df['host'].apply(lambda x:x.upper())
        self.usable_df['host_disease'] = self.usable_df['host_disease'].apply(lambda x:x.upper() if isinstance(x,str) else x)
        self.usable_df['geo_loc_name'] = self.usable_df['geo_loc_name'].apply(lambda x:x.upper().split(sep=':')[0] if isinstance(x,str) else x)
            
        # armazenando o dataframe com os registros disponíveis e uma tag indicando o grupo e o pat-id

        with open(self.usable_json_path,'w') as file:
            self.usable_json = self.usable_df.to_json(orient="records",indent=1)
            self.usable_json = js.loads(self.usable_json)
            self.tag = {self.name:self.pat_id,
                        'date':time.ctime()}
            self.usable_json.insert(0,self.tag)
            self.usable_json = js.dumps(self.usable_json,indent=1)
            file.write(self.usable_json)

        # remove o tsv_file do cache
        os.remove(self.tsv_file)
        return True
        
    
    def filterTsv(self):    
        '''
        Lê o UsableTsv de cada grupo patogênico e gera um .json filtrado, garantindo somente registros válidos,
        além de padronizar os Host em grupos específicos

        Retorna falso se não houver nenhuma informação para filtrar ou se a filtragem resultar em um json vazio
        '''

        # strings de filtragem PROBLEMA: peixes/insetos/plantas e variados
        filters = {
        'HUMAN':       'HOMO|HUMAN|SAPIENS',
        'CHICKEN':     'GALLUS|CHICKEN|POULTRY|BROILER|AVIAN|DUCK|HEN', #AVIAN DUCK ANATIDAE TURKEY
        'BOVINE':      'TAURUS|BOVINAE|COW|BOVINE|CATTLE|BEEF|CALF|MILK', #MILK
        'FISH':        'AQUA|FISH|CARP|SALMON',
        'HORSE':       'EQUUS|CABALLUS|EQUINE',
        'SHEEP':       'OVINE|SHEEP|OVIS',
        'PIG':         'PIG|SWINE|PORCINE|PORK|PIGLET|BOAR',
        'ENVIRONMENT': 'ENVIRONMENT|SOIL|HOG|SCROFA|ENVIRONMENTAL|WATER', # WATER
        'DOG':         'CANIS|LUPUS|CANINE|DOG'
        }

        if hasattr(self,'usable_df'):
            # Aproveitando cache
            self.filtered_df = self.usable_df

        else:
            if os.path.exists(self.usable_json_path):
                self.usable_df = pd.read_json(self.usable_json_path)
                self.filtered_df = self.usable_df

                # remove as colunas geradas pela tag
                self.filtered_df.drop(columns=[f"{self.name}"],inplace=True)
                self.filtered_df.drop(0,inplace=True)
            else:
                print(f'{self.name}: Sem informações para filtrar')
                return False
        
        # Segunda verificação de integridade
        def check_empty(df):
            if df.empty:
                print(f'{self.name}: Sem informações válidas')
                #os.remove(self.tsv_file)
                return True
            
        # somente genoma completo
        self.filtered_df = self.filtered_df.loc[self.filtered_df['asm_level'] == 'Complete Genome']
        if check_empty(self.filtered_df):
            return False
        
        # somente com informações de host disponíveis
        self.filtered_df = self.filtered_df[self.usable_df['host'].notnull()]
        if check_empty(self.filtered_df):
            return False
        
        # somente com assembly id disponível
        self.filtered_df = self.filtered_df[self.filtered_df['asm_acc'].notnull()]
        if check_empty(self.filtered_df):
            return False
        
        # padronizando a coluna host
        self.filtered_df['host'] = self.filtered_df['host'].apply(lambda x:x.upper())
        self.filtered_df['host_disease'] = self.filtered_df['host_disease'].apply(lambda x:x.upper() if isinstance(x,str) else x)
        self.filtered_df['geo_loc_name'] = self.filtered_df['geo_loc_name'].apply(lambda x:x.upper().split(sep=':')[0] if isinstance(x,str) else x)

        # FILTRAGEM HOST
        
        def apply_filters(df,filters=filters):
            # Para cada filtro
            for key, filter in filters.items():
                # Para cada palavra no filtro
                for var in filter.split(sep="|"):
                    # Aplica o filtro
                    df.loc[df['host'].str.contains(var),'host'] = key
            return df
        
        self.filtered_df = apply_filters(self.filtered_df)            

        # armazenando o dataframe filtrado em um .json 
        with open(self.filtered_json_path, 'w') as file:
            self.filtered_json = self.filtered_df.to_json(orient="records",indent=1)
            self.filtered_json = js.loads(self.filtered_json)
            for entry in self.filtered_json:
                entry['fasta'] = os.path.normpath(f"data/groups_info/{self.name}/{entry['asm_acc']}")

            # tag para com a data de atualização e o pat_id
            self.tag = {self.name:self.pat_id,
                        'date':time.ctime()}
            
            self.filtered_json.insert(0,self.tag)
            
            self.filtered_json = js.dumps(self.filtered_json,indent=1)
            
            file.write(self.filtered_json)
        return True

    def readFilteredTsv(self,filter=True):
        '''
        Extrai informações do .json gerado pelo getUsableTsv() ou filterTsv() (filter = True).
        As informações são redundantes.
        '''
        # coletando informações
        if filter:
            if hasattr(self,'filtered_df'):
                df = self.filtered_df
        elif hasattr(self,'usable_df'):
                df = self.usable_df
        else:
            print(f'{self.name}: Sem registros salvos')
            return False
                    
        self.species = df['scientific_name'].tolist()
        self.hosts = df['host'].tolist()
        self.strains = df['strain'].tolist()
        if not df['host_disease'].empty:
            self.host_diseases = df['host_disease'].tolist()
        else:
            self.host_diseases = []
        self.geo = df['geo_loc_name'].tolist()
        self.count = df.shape[0]

        # salvando contagem dos dados
        self.hosts_dic = {k:self.hosts.count(k) for k in set(self.hosts)}
        self.species_dic = {k:self.species.count(k) for k in set(self.species)}
        self.strains_dic = {k:self.strains.count(k) for k in set(self.strains) if self.strains.count(k) > 1}
        self.host_diseases_dic = {k:self.host_diseases.count(k) for k in set(self.host_diseases)}
        self.geo_dic = {k:self.geo.count(k) for k in set(self.geo)}

        # filtragem inicial dos nomes de espécie
        self.species_fdic = {k:[" ".join(j) for j in [i.split()[:2] for i in self.species]].count(k) for k in set([" ".join(j) for j in [i.split()[:2] for i in self.species]])}     
        self.subsp_dic = {" ".join(k.split()[:4]):[j for j in [i.split()[3] for i in self.species if ((len(i.split())>3) and (i.split()[2] == "subsp."))]].count(k.split()[3]) for k in [i for i in self.species_dic if ((len(i.split())>3) and (i.split()[2] == "subsp."))]}
    
    def makeGroupMetadata(self):
        '''
        Reune informações do readFilteredTsv() e armazena em um arquivo para posterior representação gráfica
        '''
        metadata = {'group':self.name,'count':self.count,'species':self.species_fdic,'subsp.':self.subsp_dic,'scientific_name':self.species_dic,'strain':self.strains_dic,'location':self.geo_dic, 'hosts':self.hosts_dic, 'hosts_diseases':self.host_diseases_dic}
        
        with open(os.path.normpath(self.info_path + f'/{self.name}_metadata'),'w') as log:
            content = js.dumps(metadata,indent=1)
            log.write(content)
    
        
def readGroupsNames():
    '''
    Lê o nome dos patógenos presentes no Pathogen no NCBI
    '''
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
    obj.filterTsv()
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
        obj.filterTsv()
        print("Informações atualizadas")
    else:
        print("As informações do registro não precisam ser atualizadas")

def resetExample(name="Edwardsiella_tarda"):
    '''
    Reseta somente um grupo
    '''
    obj = Group(name,from_log=False)
    obj.getPatData()
    obj.getUsableTsv()
    obj.filterTsv()
    print("Informações atualizadas")
def resetData():
    '''
    Reconstroi o banco
    '''
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group,from_log=False)
        obj.getPatData()
        obj.getUsableTsv()
        obj.filterTsv()

def updateData():
    '''
    Verifica se o banco precisa ser atualizado se baseando nos registros
    do groups_name_id.json e realiza as atualizações necessárias
    '''
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        if obj.checkDataUpdate():
            obj.getUsableTsv()
            obj.filterTsv()

def readAllData():
    '''
    Lê as informações disponíveis em todos os registros e armazena em metadata
    '''
    groups = readGroupsNames()
    for group in groups:
        obj = Group(group)
        if obj.filterTsv():
            obj.readFilteredTsv()
            obj.makeGroupMetadata()
            print(group,' lido')
        else:
            print('Sem registros válidos')
        
def readSingleData(group="Edwardsiella_tarda"):
    '''
    Lê as informações disponíveis do grupo patogênico em "group".
    default = "Edwardsiella_tarda"
    '''
    obj = Group(group)
    if obj.getUsableTsv():
        obj.filterTsv()
        obj.readFilteredTsv()
        obj.makeGroupMetadata()
    else:
        print('Sem registros úteis')
    
def generalMetadata():
    '''
    Gera um arquivo metadata com informações de todas os registros com metadata
    '''
    groups = readGroupsNames()
    metadata = {'count':0,'species':{},'hosts':{},'subsp.':{},'location':{},'hosts_diseases':{}}
    for group in groups:
        print(group)
        if os.path.exists(os.path.normpath(f'data/groups_info/{group}/{group}_metadata')):
            with open(os.path.normpath(f'data/groups_info/{group}/{group}_metadata'),'r') as log:
                data = js.load(log)
                metadata['count'] += data['count']
                
                for sp in data['species']:
                    metadata['species'][sp] = data['species'][sp]
                for sb in data['subsp.']:
                    metadata['subsp.'][sb] = data['subsp.'][sb]
                for host in data['hosts']:
                    if host in metadata['hosts']:
                        metadata['hosts'][host] += data['hosts'][host]
                    else:
                        metadata['hosts'][host] = data['hosts'][host]
                for disease in data['hosts_diseases']:
                    if disease in metadata['hosts_diseases']:
                        metadata['hosts_diseases'][disease] += data['hosts_diseases'][disease]
                    else:
                        metadata['hosts_diseases'][disease] = data['hosts_diseases'][disease]
                
                for loc in data['location']:
                    if loc in metadata['location']:
                        metadata['location'][loc] += data['location'][loc]
                    else:
                        metadata['location'][loc] = data['location'][loc]
    with open('data/metadata.json','w') as file:
        info = js.dumps(metadata,indent=1)
        file.write(info)
            
# PRECISA-SE ARRUMAR O FORMATO DO JSON FILTERED

if __name__ == "__main__":
    print("Refreshing data...")
    start = time.time()
    refresh()
    resetData()
    end = time.time()
    minutos = int((end - start) // 60)
    segundos = (end - start) % 60
    print(f'mount runtime: {minutos}:{segundos}')
    print("Reading data...")
    start = time.time()
    readAllData()
    end = time.time()
    minutos = int((end - start) // 60)
    segundos = int((end - start) % 60)
    #print(f'read runtime: {minutos}:{segundos}')
    # readAllData() = 0:50
    # new readAllData() = 0:05

    print("Generating metadata...")
    start = time.time()
    generalMetadata()
    end = time.time()
    minutos = int((end - start) // 60)
    segundos = int((end - start) % 60)
    #print(f'generalmetadata runtime: {minutos}:{segundos}')
