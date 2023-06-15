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
    names = [tag.get_text().rstrip('/') for tag in soup.find_all('a') if tag.get_text()[-1] == '/']
    if "BioProject_Hierarchy" in names:
        names.remove("BioProject_Hierarchy")
    return names

def getSinglePathogenId(name):
    url_group = url_ncbi_pat + f'/{name}/latest_kmer/Metadata/'
    html = getHtml(url_group)
    soup = BeautifulSoup(html,'html.parser')
    contents = [i.get_text() for i in soup.find_all('a') if i.get_text()[-1] != "/"]
    return contents[1].rstrip('.metadata.tsv')

def getPathogenId(names):
    '''
    Recebe a lista de nomes gerada pelo getPathogenGroups() ou refreshGroups() e devolve um dicionário com
    o nome do grupo seguido do ID da sua versão mais atual
    '''

    name_id = {}

    for name in names:
        url_group = url_ncbi_pat + f'/{name}/latest_kmer/Metadata/'
        html = getHtml(url_group)
        soup = BeautifulSoup(html,'html.parser')

        contents = [i.get_text() for i in soup.find_all('a') if i.get_text()[-1] != "/"]
        name_id[name] = contents[1].rstrip('.metadata.tsv')
    
    return name_id
        
def writeGroups(name_id):
    '''
    Recebe o dicionário gerado pelo getPathogenId() e armazena-o em um arquivo .json
    no repositório data/
    '''
    with open(group_names,'w') as file:
        groups_json = json.dumps(name_id,indent=1)
        file.write(groups_json)
        print('\nRegistros atualizados')
        
def refreshGroups():
    '''
    Verifica se existe dados sobre os grupos patogênicos e os encaminha para atualizar
    '''
    if os.path.exists(group_names):
        with open(group_names,'r') as log:
            log_info = json.load(log)
            return list(log_info.keys())
    else:
        print('Sem informações prévias para atualizar\nConsidere usar a função extract')
        return False

def refreshSingleGroup(name):
    '''
    Verifica se existe dados sobre os grupos patogênicos e atualiza somente um
    '''
    if os.path.exists(group_names):
        with open(group_names,'r') as log:
            log_info = json.load(log)
            old = log_info[name]
            log_info[name] = getSinglePathogenId(name)
            new = log_info[name]
            if old != new:
                with open(group_names,'w') as file:
                    file.write(json.dumps(log_info,indent=1))
                    print(name,": ",old," ---> ",new)
            else:
                print('Versão mais atual')
# principais

def refresh():
    names = refreshGroups()
    new_ids = getPathogenId(names)
    with open(group_names,'r') as log:
        log_info = json.load(log)
        for i in names:
            if new_ids[i] == log_info[i]:
                print(i,': ',log_info[i])
            else:
                print(i,': ',log_info[i],'--->',new_ids[i])
    writeGroups(new_ids)
def mount():
    print('Extraindo os nomes dos grupos patogênicos do banco Pathogen do NCBI')
    groups = getPathogenGroups()
    print('\nRecuperando o PathogenID mais recente dos grupos de patógenos')
    data = getPathogenId(groups)
    print('\nArmazenando informações extraídas na pasta /data/groups')
    writeGroups(data)
    
