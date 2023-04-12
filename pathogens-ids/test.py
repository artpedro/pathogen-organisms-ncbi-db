import os
import subprocess
#command = ["esearch" ,"-db" ,"assembly" ,"-query" ,"'\"Salmonella enterica\"[Organism] AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])'", "|" ,"efetch", "-format", "docsum" ,"|","xtract" "-pattern", "DocumentSummary","-element" ,"AssemblyAccession" ,"-element" ,"BioSampleAccn","|","cat"]
#sub_IDS = subprocess.check_output(command,text=True)
#print(str(type(sub_IDS))+" "+ sub_IDS)
os.system("esearch -db assembly -query '\"Salmonella enterica\"[Organism] AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])' | efetch -format docsum |xtract -pattern DocumentSummary -element AssemblyAccession -element BioSampleAccn > temp_ids")
ids = []
with open("temp_ids","r") as temp:
    lines = temp.readlines()
    for line in lines:
        line = line.rstrip("\n")
        id = line.split('\t')
        ids.append(id)
biosample = [i[1] for i in ids]
print(biosample)

