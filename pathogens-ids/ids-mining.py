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
        with open(f'pathogens-ids/logs/{self.name}_log_ids.txt','w') as log:
            log.write(f'{self.name} | número de assemblies: {self.count} | date: {datetime.now()}\n')
            for i in self.ids:
                log.write(f'{i[0]} {i[1]}\n')
        print(f'\n\nLOG de {self.name} gerado\n\n')
    
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
        self.assembly_ids= [i[0] for i in self.ids]
        self.biosample_ids= [i[1] for i in self.ids]
        os.remove(f'data/{self.name}_temp_ids.txt')
        print(f'\n\nIDs de {self.name} salvos na memória\n\n')

    def write_ids_files(self):
        self.path = f'data/ids/{self.name}'
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        with open(self.path+f'/{self.name}_ids.txt','w') as ids:
            ids.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.ids:    
                ids.write(f'{i[0]} {i[1]}\n')
        with open(self.path+f'/{self.name}_assembly_ids.txt', 'w') as assembly:
            assembly.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.assembly_ids:
                assembly.write(f'{i}\n')
        with open(self.path+f'/{self.name}_biosample_ids.txt','w') as biosample:
            biosample.write(f'{self.name} | número de entradas: {self.count} | date: {datetime.now()}\n')
            for i in self.biosample_ids:
                biosample.write(f'{i}\n')
        print(f'\n\nIDs de {self.name} salvos em {self.path}\n\n')        
    
    def check_new_entries(self):
        # salvando o output do esearch
        esearch_out = subprocess.check_output(self.esearch,text=True)
        new_count =  esearch_out.split('</Count>')[0].split('<Count>')[1]
        if new_count != self.count:
            print(f'\n\n\nNOVAS {int(new_count)-int(self.count)} ENTRADAS DETECTADAS PARA' +f'{self.name}\n\n\n'.upper())
            return True
        print("Nenhuma nova entrada detectada")
        return False
    def update_ids(self):
        with open(f'pathogens-ids/logs/{self.name}_log_ids.txt','r') as log:
            label = log.readline()
            last_date = label.split('| date: ')[1]
            last_date = last_date.split()[0].replace('-','/')
            print("Última vez atualizado: "+last_date)
            new_esearch = ["esearch", "-db", "assembly" ,"-query" , f"((\"{last_date}\"[AsmUpdateDate] : \"3000\"[AsmUpdateDate]) AND \"{self.name}\"[Organism]) AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])"]
            new_esearch_str = " ".join(new_esearch)
            out_new_esearch =  subprocess.check_output(new_esearch,text=True)
            print(f'out: {out_new_esearch}')
            new_count = out_new_esearch.split('</Count>')[0].split('<Count>')[1]
            #Arrumar aqui
            print(f'Novas {new_count} entradas detectadas')
            os.system(f"{new_esearch_str} | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > data/{self.name}_temp_ids.txt")
            with open(f'data/{self.name}_temp_ids.txt','r') as temp:
                new_ids = temp.readlines()
                new_ids = [i.split() for i in new_ids]
                with open(f'{self.path}/{self.name}_ids.txt','a') as ids:
                    for i in new_ids:
                        ids.write(f'{i[0]} {i[1]}\n')
                with open(f'{self.path}/{self.name}_assembly_ids.txt', 'a') as assembly:
                    for i in new_ids:
                        assembly.write(f'{i[0]}\n')
                with open(self.path+f'/{self.name}_biosample_ids.txt','a') as biosample:
                    for i in self.biosample_ids:
                        biosample.write(f'{i[1]}\n')


ed = Species('Edwardsiella_tarda')
ed.extract_ids()
print(ed.ids)
print(ed.assembly_ids)
print(ed.biosample_ids)
ed.write_ids_files()
ed.check_new_entries()         
ed.log_ids()
ed.update_ids()



        