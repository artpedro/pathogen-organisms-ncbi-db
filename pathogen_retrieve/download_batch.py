import json as js
import os
from tqdm import tqdm
import requests
import pandas as pd
import sys
import time
from pathogen_data_mining import Group
from bs4 import BeautifulSoup

def time_func(func):
    start = time.time()
    func()
    end = time.time()
    minutos = int((end - start) // 60)
    segundos = int((end - start) % 60)
    print(f'Runtime: {minutos}:{segundos}')

def download(url: str, fname: str):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    descr = fname.split('/')[-1]
    # Can also replace 'file' with a io.BytesIO object
    with open(fname, 'wb') as file, tqdm(
        desc=descr,
        total=total,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in resp.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)

def download_url(url, file):
    with open(file, "wb") as f:
        print("Downloading %s" % file)
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )    
                sys.stdout.flush()

def readGroupsNames():
    '''
    Lê o nome dos patógenos presentes no Pathogen no NCBI
    '''
    with open(os.path.normpath("data/groups/groups_name_id.json"),"r") as data:
        groups = js.load(data)
        return groups

def download_batch(batch,path):
    '''
    Recebe uma lista de Assembly IDs e faz o download da ultima versão do arquivo fasta compactado
    no caminho informado por path
    '''
    ftp_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'

    if not os.path.exists(path):
        os.mkdir(path)
    for entry in batch:

        # preciso adicionar um if para não baixar assemblies já baixados, talvez uma nova função
        page = f'{entry[0:3]}/{entry[4:7]}/{entry[7:10]}/{entry[10:13]}'
        genome_url = ftp_url + page
        download(genome_url,'temp')
        with open('temp') as html:
            soup = BeautifulSoup(html, 'html.parser')
            names = [tag.get_text().rstrip('/') for tag in soup.find_all('a') if tag.get_text()[-1] == '/']
            full = f'{genome_url}/{names[-1]}/{names[-1]}_genomic.fna.gz'
        download(full,f'{path}/{entry}.fna.gz')

def list_ids(name):
    '''
    Devolve uma lista com todos os Assembly IDs em data correspondente ao grupo informado em "name"
    '''
    ids = []
    with open(f'data/groups_info/{name}/{name}_filtered.json','r') as data:
        all_data = js.load(data)
        for line in all_data:
            if "asm_acc" in line:
                ids.append(line["asm_acc"])
    return ids

def download_group(name='Edwardsiella_tarda'):
    '''
    Baixa todos os arquivos fasta de um grupo patogênico.
    '''
    download_batch(list_ids(name),f'data/fasta/{name}')

def download_all():
    names = readGroupsNames()
    for name in names:
        download_group(name)

