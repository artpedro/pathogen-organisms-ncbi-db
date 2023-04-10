from bs4 import BeautifulSoup
import os

# Link para a tabela dos organismos patogênicos
url_ncbi_pat = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
os.system(f"wget {url_ncbi_pat} -P pathogens-names")

# Codigo-fonte da pagina

with open("pathogens-names/index.html","r") as html:
    soup = BeautifulSoup(html, 'html.parser')
    tags = soup.find_all('a')
    with open("pathogens-names/filtered.html", "w") as t:
        for a in tags:
            t.write(str(a))
os.system("rm pathogens-names/index.html")

names =[]

with open("pathogens-names/filtered.html") as html:
    soup = BeautifulSoup(html,'html.parser')
    for tag in soup.find_all("a"):
        content = tag.get_text()
        if content[-1] == "/":
            names.append(content[0:-1])
    with open("pathogens-names/species_names","w") as file:
        for name in names:
            file.write(name+"\n")

os.system('rm pathogens-names/filtered.html')