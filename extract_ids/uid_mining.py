from Bio import Entrez
import os
import sys
import subprocess
from time import sleep
from datetime import datetime
import json

# email para contato caso dê problema no entrez
Entrez.email = "arturpedromartins@gmail.com"

class Species():
    def __init__(self,name):
        # salvando o nome da especie
        self.name = name
        self.db = "assembly"

        print(f'\n\nEspécie {self.name} iniciada\n\n')
        
        self.path = f'data/ids/{self.name}'
        if not os.path.exists(self.path):
            os.mkdir(self.path)        
        
        ''' 
        # extraindo informações se já existir um log
        if os.path.exists(f'extract_ids/logs/{self.name}_log_ids.txt'):
            with open(f'extract_ids/logs/{self.name}_log_ids.txt', 'r') as log:
                self.ids = [i.split() for i in log.readlines() if len(i.split()) == 2]
                self.assembly_ids = [i[0] for i in self.ids]
                self.biosample_ids = [i[1] for i in self.ids]
        else:
            print('\n\nEsta espécie não possui registros baixados')    
       '''

    def extract_ids(self):
        ''' 
            o método nao recebe argumentos

            .extract_ids() é responsável por receber a lista de UIDs presente no atributo "record"
            e adquirir o AssemblyAcession e o BioSampleID, para posterior download e classificação
            do patogênico. extract_ids() não retorna nada, mas armazena os dados em atributos do objeto.
        '''
        # handle do entrez esearch
        handle = Entrez.esearch(db=self.db, term=f"(\"{self.name}\"[Organism] OR {self.name}[All Fields]) AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])",idtype="acc",retmax = '100000')
        sleep(0.37)
        # dicionário contendo informações da query acima
    
        self.record = Entrez.read(handle)
        handle.close()
        # armazenando os uids em uma variável
        self.new_uids = self.record['IdList']

        # iniciando lista com os acessos no Assembly e o identificador único no BioSample
        self.asm_idlist=[]
        self.bsm_uidlist=[]
        print("\nExtraindo os ids dos registros de genoma completo encontrados da espécie\n")
        
        handle = Entrez.esummary(db=self.db, id=self.new_uids)
        record = Entrez.read(handle)
        summary = record['DocumentSummarySet']['DocumentSummary']
        for uid in summary:
            self.asm_idlist.append(uid['AssemblyAccession'])
            self.bsm_uidlist.append(uid['BioSampleId'])
        return True
    def write_ids_files(self):
        # escreve os ids salvos na memória em arquivos dentro de data/ids
        with open(self.path+f'/{self.name}_uids.json','w') as file:
            data = {'label':
                     {'name':f'{self.name}',
                      'count':f'{len(self.new_uids)}'
                      ,'date':f'{datetime.now()}'
                      },
                      'uids':[i for i in self.new_uids]
                      }
            content = json.dumps(data,indent=1)
            file.write(content)
        with open(self.path+f'/{self.name}_ids.json','w') as file:
            # a estrutura do arquivo será em json com uma chave sendo a "etiqueta" do registro
            # e a outra sendo as entradas, no formato: "AssemblyID":"BioSampleUID"

            data = {'label':
                     {'name':f'{self.name}',
                      'count':f'{len(self.new_uids)}'
                      ,'date':f'{datetime.now()}'
                      },
                      'entries':{}
                      }
            for i,j in zip(self.asm_idlist,self.bsm_uidlist):
                data['entries'][i]=j        
            content = json.dumps(data,indent=1)
            file.write(content)       

        print(f'\n\nIDs de {self.name} salvos em {self.path}\n\n')

    def check_new_entries(self):
        # verifica se a quantidade de entradas da espécie está desatualizada com o ncbi
        # true = novas entradas que precisam ser extraídas
        # false = sem novas entradas
        # OBS: os ids precisam estar salvos na memória pelo extract_ids()
        # salvando o output do esearch

        new_count = len(self.new_uids)
        print(f'Quantidade de registros encontrados: {new_count}')
        try:
            with open(self.path+f'/{self.name}_uids.json', 'r') as log:
                
                # extraindo a quantidade de acessos em log
                data = json.load(log)
                log_count = int(data['label']['count'])
                print(f'\nQuantidade de registros armazenados: {log_count}')
            if int(new_count) > log_count:
                print(f'\n\n\nNOVAS {int(new_count)-int(log_count)} ENTRADAS DETECTADAS PARA ' +f'{self.name}\n\n\n'.upper())
                return True
            print("Nenhuma nova entrada detectada")
            return False
        except:
            print('Nenhum registro presente')
            return True

args = sys.argv
if args[1] == "mount":
    with open('extract_names/species_names.txt','r') as names:
        all_names = [i.rstrip('\n') for i in names.readlines()]

        for name in all_names:
            sleep(0.37)
            obj = Species(name)
            if obj.extract_ids():
                if obj.check_new_entries():
                    obj.write_ids_files()        
