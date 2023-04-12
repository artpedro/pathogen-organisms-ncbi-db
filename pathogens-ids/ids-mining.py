import os
import subprocess
class Species():
    def __init__(self,name):
        # salvando o nome da especie
        self.name = name

        # salvando o comando do esearch
        self.esearch = ["esearch", "-db", "assembly" ,"-query" , "'\"Salmonella Enterica\"[Organism] AND latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter] AND (latest[filter] AND all[filter] NOT anomalous[filter])'"]

        # salvando o output do esearch
        self.esearch_output = subprocess.check_output(self.esearch, text=True)
        self.count =  subprocess.check_output(f"echo {self.esearch_output}| xtract -pattern Count -element Count|cat",text=True)
    def extract_ids(self):
        # extraindo os ids e colocando em um arquivo temporário
        os.system("esearch -db assembly -query '\"Salmonella enterica\"[Organism] AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])' | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > temp_ids")
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
    def log_ids(self):
        #aqui eu vou pegar todos os IDs e atualizar o log com possíveis novos  
        pass
    def save_ids(self,path):
        #aqui eu vou salvar os IDs em um txt com uma label
        pass
    #nao esquecer de limpar os arquivos temporários
        



        