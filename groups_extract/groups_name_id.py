from bs4 import BeautifulSoup
import urllib3
import json
import os

url_ncbi_pat = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
group_names = "data/groups/groups_name_id.json"
http = urllib3.PoolManager()

def getHtml(url):
    '''
    Recebe um url e devolve o arquivo HTML recebido pelo metódo GET no protocolo HTTP
    '''

    resp = http.request('GET',url)
    return resp.data.decode('utf-8')

def getPathogenGroups():
    '''
    Devolve uma lista com os nomes dos grupos de espécies presentes no banco de dados
    Pathogen do NCBI
    '''
    names =[]

    html = getHtml(url_ncbi_pat)
    soup = BeautifulSoup(html, 'html.parser')
    for tag in soup.find_all('a'):
        if tag.get_text()[-1] == '/':
            names.append(tag.get_text().rstrip('/'))
    return names

def getPathogenId(names):
    '''
    Recebe a lista de nomes gerada pelo getPathogenGroups() ou refreshGroups() e devolve um dicionário com
    o nome do grupo seguido do ID da sua versão mais atual
    '''

    if names == []:
        return
    name_id = {}

    for name in names:
        url_group = url_ncbi_pat + f'/{name}/latest_kmer/'
        html = getHtml(url_group)
        soup = BeautifulSoup(html,'html.parser')

        contents = [i.get_text() for i in soup.find_all('a') if i.get_text()[-1] != "/"]
        name_id[name] = contents[1].rstrip('.final.descriptor.xml')
    
    return name_id
        
def writeGroups(name_id):
    '''
    Recebe o dicionário gerado pelo getPathogenId() e armazena-o em um arquivo .json
    no repositório data/
    '''
    with open(group_names,'w') as file:
        groups_json = json.dumps(name_id,indent=1)
        file.write(groups_json)
   
def RefreshGroups():
    '''
    
    '''
    if os.path.exists(group_names):
        with open(group_names,'r') as log:
            log_info = json.load(log)
            return [log_info.keys()]
    else:
        print('Sem informações prévias para atualizar\nConsidere usar a função extract')
        return []

    '''species_names = "extract_names/species_names.txt"
    with open(species_names,"r") as rfile:
        # Nomes já no arquivo
        names_infile = [line.rstrip("\n") for line in rfile.readlines()]
        new_names = names_infile
        print('Verificando novos nomes')
        for name in getPathogenGroups():
            if name not in names_infile:
                new_names.append(name)
                print(f"Adicionando {name} ao species_names.txt")      
        write_ncbi_pat(new_names)'''
    pass

#ncbi_names = getPathogenGroups()
writeGroups(getPathogenId(getPathogenGroups()))