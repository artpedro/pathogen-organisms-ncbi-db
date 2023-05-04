import os
import sys
import subprocess
from datetime import datetime
class Species():
    def __init__(self,name):
        # salvando o nome da especie
        self.name = name

        # salvando o comando do esearch
        # obs: somente assemblies de genoma completo
        self.esearch = ["esearch", "-db", "assembly" ,"-query" , f"'\"{self.name}\"[Organism] AND latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter]'"]

        # salvando o output do esearch
        #self.esearch_output = subprocess.check_output(self.esearch, text=True)
        
        # salvando a quantidade de resultados da pesquisa
        #self.count =  int(self.esearch_output.split('</Count>')[0].split('<Count>')[1])
        
        print(f'\n\nEspécie {self.name} iniciada\n\n')
        
        # extraindo informações se já existir um log
        if os.path.exists(f'extract_ids/logs/{self.name}_log_ids.txt'):
            with open(f'extract_ids/logs/{self.name}_log_ids.txt', 'r') as log:
                self.ids = [i.split() for i in log.readlines() if len(i.split()) == 2]
                self.assembly_ids = [i[0] for i in self.ids]
                self.biosample_ids = [i[1] for i in self.ids]
        else:
            print('\n\nEsta espécie não possui registros baixados')    
       

    def extract_ids(self):
        # extraindo os ids e colocando em um arquivo temporário
        os.system(f"esearch -db assembly -query '\"{self.name}\"[Organism] AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])' | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > data/{self.name}_temp_ids.txt")
        self.ids = []
        #lendo o arquivo temporário e extraindo os IDs
        with open(f"data/{self.name}_temp_ids.txt","r") as temp:
            lines = temp.readlines()
            for line in lines:
                line = line.rstrip("\n")
                id = line.split('\t')
                self.ids.append(id)

        # ids de assembly salvos no objeto
        self.assembly_ids= [i[0] for i in self.ids]
        # ids de biosample salvos no objeto
        self.biosample_ids= [i[1] for i in self.ids]
        self.count = len(self.ids)

        # removendo o arquivo temporário
        os.remove(f'data/{self.name}_temp_ids.txt')
        print(f'\n\nIDs de {self.name} salvos na memória\n\n')
        self.log_ids()

    def log_ids(self):
        #aqui eu vou pegar todos os IDs e atualizar o log com possíveis novos 

        # arquivo em logs com o id assembly e o id biosample
        with open(f'extract_ids/logs/{self.name}_log_ids.txt','w') as log:
            log.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.ids:
                log.write(f'{i[0]} {i[1]}\n')
        print(f'\n\nLOG de {self.name} gerado\n\n')
    
    def write_ids_files(self):
        # escreve os ids salvos na memória em arquivos dentro de data/ids

        # caminho único de cada espécie
        self.path = f'data/ids/{self.name}'
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        
        # salvando os dois ids em um arquivo
        with open(self.path+f'/{self.name}_ids.txt','w') as ids:
            ids.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.ids:    
                ids.write(f'{i[0]} {i[1]}\n')

        # salvando os ids de assembly
        with open(self.path+f'/{self.name}_assembly_ids.txt', 'w') as assembly:
            assembly.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.assembly_ids:
                assembly.write(f'{i}\n')

        # salvando os ids de biosample        
        with open(self.path+f'/{self.name}_biosample_ids.txt','w') as biosample:
            biosample.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.biosample_ids:
                biosample.write(f'{i}\n')
        print(f'\n\nIDs de {self.name} salvos em {self.path}\n\n')
        self.log_ids()      
    
    def check_new_entries(self):
        # verifica se a quantidade de entradas da espécie está desatualizada com o ncbi
        # true = novas entradas que precisam ser extraídas
        # false = sem novas entradas

        # salvando o output do esearch
        esearch_out = subprocess.check_output(self.esearch,text=True)
        print(esearch_out)
        new_count =  esearch_out.split('</Count>')[0].split('<Count>')[1]
        print(new_count)
        with open(f'extract_ids/logs/{self.name}_log_ids.txt', 'r') as log:
            # extraindo a quantidade de acessos em log
            label = log.readline()
            log_count = int(label.split(sep='|')[1].lstrip(' número de entradas: '))
            print(log_count)
        if new_count != log_count:
            print(f'\n\n\nNOVAS {int(new_count)-int(log_count)} ENTRADAS DETECTADAS PARA' +f'{self.name}\n\n\n'.upper())
            return True
        print("Nenhuma nova entrada detectada")

        return False
    
    def update_ids(self):
        with open(f'extract_ids/logs/{self.name}_log_ids.txt','r') as log:
            # lendo o arquivo log para saber quando foi feita a última verificação
            label = log.readline()
            last_date = label.split(sep='|')[2].lstrip(' date: ')
            last_date = last_date.split()[0].replace('-','/')
            print("Última vez atualizado: "+last_date)

            # pesquisando pelo esearch por entradas posteriores a data do log
            new_esearch = ["esearch", "-db", "assembly" ,"-query" , f"((\"{last_date}\"[AsmUpdateDate] : \"3000\"[AsmUpdateDate]) AND \"{self.name}\"[Organism]) AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])"]
            new_esearch_str = " ".join(new_esearch)
            out_new_esearch =  subprocess.check_output(new_esearch,text=True)
            print(f'out: {out_new_esearch}')
            new_count = out_new_esearch.split('</Count>')[0].split('<Count>')[1]

            if int(new_count) != 0:
                print(f'Novas {new_count} entradas detectadas')
                os.system(f"{new_esearch_str} | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > data/{self.name}_temp_ids.txt")
                with open(f'data/{self.name}_temp_ids.txt','r') as temp:
                    new_ids = temp.readlines()
                    new_ids = [i.split() for i in new_ids]
                    new_asm = [i[0] for i in new_ids]
                    new_bio = [i[1] for i in new_ids]
                    self.ids = self.ids + new_ids
                    self.assembly_ids = self.assembly_ids + new_asm
                    self.biosample_ids = self.biosample_ids + new_bio
                    with open(f'{self.path}/{self.name}_ids.txt','a') as ids:
                        for i in new_ids:
                            ids.write(f'{i[0]} {i[1]}\n')
                    with open(f'{self.path}/{self.name}_assembly_ids.txt', 'a') as assembly:
                        for i in new_ids:
                            assembly.write(f'{i[0]}\n')
                    with open(self.path+f'/{self.name}_biosample_ids.txt','a') as biosample:
                        for i in self.biosample_ids:
                            biosample.write(f'{i[1]}\n')
                self.log_ids()
            else:
                print(f'Nenhum novo registro encontrado para {self.name}')

args = sys.argv
if args[1] == "mount":
    with open('extract_names/species_names.txt','r') as names:
        all_names = [i.rstrip('\n') for i in names.readlines()]
        all_species = {}
        for name in all_names:
            all_species[name] = Species(name)
            if all_species[name].check_new_entries():
                all_species[name].update_ids()
        
    
#ed.check_new_entries()         
#ed.log_ids()
#ed.update_ids()


        