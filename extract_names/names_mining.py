from bs4 import BeautifulSoup
import os
import sys
url_ncbi_pat = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
raw_html = "extract_names/index.html"
filtered_html = "extract_names/filtered.html"
species_names = "extract_names/species_names.txt"

def check_names():
    os.system(f'cat {species_names}')

# Acessa e armazena os nomes dos organismos patogênicos presentes no ncbi em uma variável list
def get_ncbi_pat():
    # Lista de nomes iniciada
    names =[]
    # Extraindo o html da página
    os.system(f"wget {url_ncbi_pat} -P extract_names")

    # Abrindo o html e gerando um .html somente com as tags <a>
    with open(raw_html,"r") as raw:
        soup = BeautifulSoup(raw, 'html.parser')
        tags = soup.find_all('a')
        # Escrevendo um arquivo somente com as tags filtradas
        with open(filtered_html, "w") as filtered:
            for a in tags:
                filtered.write(str(a))
    os.remove(raw_html)

    # Filtrando somente os nomes das espécies do .html e alimentando a lista
    with open(filtered_html,"r") as html:
        soup = BeautifulSoup(html,'html.parser')
        
        # Encontrando somente os <a> que representam repositórios
        for tag in soup.find_all("a"):
            content = tag.get_text()
            if content[-1] == "/":
                names.append(content[0:-1])
    os.remove(filtered_html)
    return names
# Recebe uma lista e adiciona os nomes dessa lista ao arquivo species_names.txt
def write_ncbi_pat(species):
    species_names = "extract_names/species_names.txt"
    with open(species_names,'w') as file:
        file.write('\n'.join(species))

# Verifica o arquivo species_names.txt e atualiza com possíveis novos nomes presentes no ncbi            
def refresh_ncbi_pat():
    species_names = "extract_names/species_names.txt"
    with open(species_names,"r") as rfile:
        # Nomes já no arquivo
        names_infile = [line.rstrip("\n") for line in rfile.readlines()]
        new_names = names_infile
        print('Verificando novos nomes')
        for name in get_ncbi_pat():
            if name not in names_infile:
                new_names.append(name)
                print(f"Adicionando {name} ao species_names.txt")      
        write_ncbi_pat(new_names)

# pela mudança de extração dos dados, essas funções perdem o sentido
'''
# Adicionar nome manualmente a lista
def add_name(name):
    species_names = "extract_names/species_names.txt"
    with open(species_names, "r") as names:
        all = [line.rstrip("\n") for line in names.readlines()]
        # Checando se o nome já consta na lista
        if name in all:
            print(f"{name} já está na lista")
        else:    
            all.append(name)
            with open(species_names,"w") as new_file:
                new_file.write("\n".join(sorted(all)))
                print(name+" adicionada ao species_name.txt")

# Remover nome manualmente da lista
def remove_name(name):
    with open(species_names, "r") as names:
        all = [line.rstrip("\n") for line in names.readlines()]
        if name in all:
            all.remove(name)
            with open(species_names,"w") as new_names:
                new_names.write("\n".join(all))
                print(name+" removida do species_name.txt")
        else:
            print(name+" não está na lista")
'''   

ncbi_names = get_ncbi_pat()

if os.path.exists("extract_names/species_names.txt"):
    refresh_ncbi_pat()
else:
    write_ncbi_pat(ncbi_names)


'''
        # handle ncbi_pat case
    elif sys.argv[1] == "add":
        if sys.argv[2]:
            name = sys.argv[2]
            add_name(name)
        else:
            print("Uso: python3 script.py <check|ncbi_pat|add|remove> <name>")
    elif sys.argv[1] == "remove":
        if sys.argv[2]:
            name = sys.argv[2]
            remove_name(name)
        else:
            print("Uso: python3 script.py <check|ncbi_pat|add|remove> <name>")
    '''      
