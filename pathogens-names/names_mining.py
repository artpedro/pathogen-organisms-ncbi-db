from bs4 import BeautifulSoup
import os
import sys

# Possíveis inputs: ncbi_pat; add; remove; list

operation = sys.argv[1]
if operation == 'add' or operation == 'remove':
    species = sys.argv[2]

# Atualizar pelo NCBI_Pathogen 
if (operation == "ncbi_pat") or not operation:
    # Link para a tabela dos organismos patogênicos
    url_ncbi_pat = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"

    # Extraindo o html da página
    os.system(f"wget {url_ncbi_pat} -P pathogens-names")

    # Abrindo o html e gerando um .html somente com as tags <a>
    with open("pathogens-names/index.html","r") as html:
        soup = BeautifulSoup(html, 'html.parser')
        tags = soup.find_all('a')
        # Escrevendo um arquivo somente com as tags filtradas
        with open("pathogens-names/filtered.html", "w") as t:
            for a in tags:
                t.write(str(a))
    os.system("rm pathogens-names/index.html")

    # Lista de nomes iniciada
    names =[]
    
    # Filtrando somente os nomes das espécies do .html e alimentando a lista
    with open("pathogens-names/filtered.html") as html:
        soup = BeautifulSoup(html,'html.parser')
        
        # Encontrando somente os <a> que representam repositórios
        for tag in soup.find_all("a"):
            content = tag.get_text()
            if content[-1] == "/":
                names.append(content[0:-1])
        
        # Atualizar lista
        if os.path.exists("pathogens-names/species_names.txt"):
            with open("pathogens-names/species_names.txt","r") as rfile:
                # Nomes já no arquivo
                names_infile = [line.rstrip("\n") for line in rfile.readlines()]
                with open("pathogens-names/species_names.txt","w") as file:
                    new_names = names_infile
                    for name in names:
                        if name not in names_infile:
                            new_names.append(name)
                            print(f"Adicionando {name} ao species_names.txt")      
                    for name in new_names:
                        file.write(name + "\n")
        else:
            with open("pathogens-names/species_names.txt","w") as file:
                for name in names:
                        file.write(name+"\n")
    os.system('rm pathogens-names/filtered.html')

# Adicionar nome manualmente a lista
if operation == "add" and species:
    with open("pathogens-names/species_names.txt", "r") as names:
        all = names.readlines()
        # Checando se o nome já consta na lista
        if f"{species}\n" in all:
            print(f"{species} já está na lista")
        else:    
            all.append(species+"\n")
            with open("pathogens-names/species_names.txt","w") as new_file:
                for name in sorted(all):
                    new_file.write(name)
                print(species+" adicionada ao species_name.txt")

# Remover nome manualmente da lista
if operation == "remove" and species:
    with open("pathogens-names/species_names.txt", "r") as names:
        all = [line.rstrip("\n") for line in names.readlines()]
        if species in all:
            all.remove(species)
            with open("pathogens-names/species_names.txt","w") as new_names:
                for name in all:
                    new_names.write(name+"\n")
                print(species+" removida do species_name.txt")
        else:
            print(species+" não está na lista")