import os
import subprocess
from datetime import datetime
class Species():
    def __init__(self,name):
        # salvando o nome da especie
        self.name = name

        # salvando o comando do esearch
        self.esearch = ["esearch", "-db", "assembly" ,"-query" , f"'\"{self.name}\"[Organism] AND latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter] AND (latest[filter] AND all[filter] NOT anomalous[filter])'"]

        # salvando o output do esearch
        self.esearch_output = subprocess.check_output(self.esearch, text=True)
        
        self.count =  self.esearch_output.split('</Count>')[0].split('<Count>')[1]

        print(f'\n\nEspécie {self.name} iniciada\n\n')

    def log_ids(self):
        #aqui eu vou pegar todos os IDs e atualizar o log com possíveis novos 
        if not os.path.exists(f'logs/{self.name}_log_ids.txt'):
            with open(f'logs/{self.name}_log_ids.txt','w') as log:
                log.write(f'{self.name} | número de assemblies: {self.count} | date: {datetime.now()}')
                for i,j in self.assembly_ids,self.biosample_ids:
                    log.write(f'{i},{j}')
            print(f'\n\nLOG de {self.name} gerado\n\n')
        else:
            self.refresh_ids()

    
    def extract_ids(self):
        # extraindo os ids e colocando em um arquivo temporário
        os.system(f"esearch -db assembly -query '\"{self.name}\"[Organism] AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])' | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > {self.name}_temp_ids.txt")
        self.ids = []
        #lendo o arquivo temporário e extraindo os IDs
        with open("temp_ids","r") as temp:
            lines = temp.readlines()
            for line in lines:
                line = line.rstrip("\n")
                id = line.split('\t')
                self.ids.append(id)
        self.assembly_ids= [i[0] for i in self.ids]
        self.biosample_ids= [i[1] for i in self.ids]
        os.system(f'rm {self.name}_temp_ids.txt')
        print(f'\n\nIDs de {self.name} salvos na memória\n\n')

    def write_ids_files(self):
        path = f'ids/{self.name}'
        with open(path+f'/{self.name}_ids.txt','w') as ids:
            ids.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i,j in self.assembly_ids,self.biosample_ids:    
                ids.write(f'{i},{j}\n')
        with open(path+f'/{self.name}_assembly_ids.txt', 'w') as assembly:
            assembly.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.assembly_ids:
                assembly.write(f'{i}\n')
        with open(path+f'/{self.name}_biosample_ids.txt','w') as biosample:
            biosample.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.biosample_ids:
                biosample.write(f'{i}\n')
        print(f'\n\nIDs de {self.name} salvos em {path}\n\n')        
    #nao esquecer de limpar os arquivos temporários
    def check_new_entries(self):
        # salvando o output do esearch
        esearch_out = subprocess.check_output(self.esearch,text=True)
        new_count =  esearch_out.split('</Count>')[0].split('<Count>')[1]
        if new_count != self.count:
            print(f'\n\n\nNOVAS {int(new_count)-int(self.count)} ENTRADAS DETECTADAS PARA' +f'{self.name}\n\n\n'.upper())
        return True
    def update_ids(self):
        with open(f'logs/{self.name}_log_ids.txt','r') as log:
            label = log.readline()
            last_date = label.split('| date: ')[1]
            last_date = last_date.split()[0].replace('-','/')
            new_esearch = ["esearch", "-db", "assembly" ,"-query" , f"((\"{last_date}\"[AsmUpdateDate] : \"3000\"[AsmUpdateDate]) AND \"{self.name}\"[Organism]) AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])"]
            out_new_esearch =  subprocess.check_output(new_esearch,text=True).split('</Count>')[0].split('<Count>')[1]
            new_count = out_new_esearch.split('</Count>')[0].split('<Count>')[1]

            print(f'Novas {new_count} entradas detectadas')
            os.system(f"echo {out_new_esearch}  | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > {self.name}_temp_ids.txt")
            #ler arquivo e adicionar aos atributos corretos

            

            



        