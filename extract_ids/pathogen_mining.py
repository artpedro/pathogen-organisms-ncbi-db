import json
import os
import requests
import shutil
from tqdm.auto import tqdm

import wget
from wget import download

# Problemas na hora de extrair arquivos chunked
# solução:
'''from tqdm import tqdm
import os
import sys
from pathlib import Path
import requests

# This example script downloads python program for mac.

# Home directory of Mac, pathlib.Path module make this easy.
home_path = Path.home()
# This is the sub directory under home directory.
sub_path = "tmp"
# The download link of python.
url = "https://www.python.org/ftp/python/3.8.5/python-3.8.5-macosx10.9.pkg"

# The header of the dl link has a Content-Length which is in bytes.
# The bytes is in string hence has to convert to integer.
filesize = int(requests.head(url).headers["Content-Length"])

# os.path.basename returns python-3.8.5-macosx10.9.pkg,
# without this module I will have to manually split the url by "/"
# then get the last index with -1.
# Example:
# url.split("/")[-1]
filename = os.path.basename(url)

# make the sub directory, exists_ok=True will not have exception if the sub dir does not exists.
# the dir will be created if not exists.
os.makedirs(os.path.join(home_path, sub_path), exist_ok=True)

# The absolute path to download the python program to.
dl_path = os.path.join(home_path, sub_path, filename)
chunk_size = 1024

# Use the requests.get with stream enable, with iter_content by chunk size,
# the contents will be written to the dl_path.
# tqdm tracks the progress by progress.update(datasize)
with requests.get(url, stream=True) as r, open(dl_path, "wb") as f, tqdm(
        unit="B",  # unit string to be displayed.
        unit_scale=True,  # let tqdm to determine the scale in kilo, mega..etc.
        unit_divisor=1024,  # is used when unit_scale is true
        total=filesize,  # the total iteration.
        file=sys.stdout,  # default goes to stderr, this is the display on console.
        desc=filename  # prefix to be displayed on progress bar.
) as progress:
    for chunk in r.iter_content(chunk_size=chunk_size):
        # download the file chunk by chunk
        datasize = f.write(chunk)
        # on each chunk update the progress bar.
        progress.update(datasize)'''
# implementar

# Downloader class.
def download(url,file):
    with requests.get(url, stream=True) as r:
    
    # check header to get content length, in bytes
        total_length = int(r.headers.get("Content-Length"))

    # implement progress bar via tqdm
        with tqdm.wrapattr(r.raw, "read", total=total_length, desc="")as raw:
    
        # save the output to a file
            with open(file, 'wb')as output:
                shutil.copyfileobj(raw, output)

class Group():
    def __init__(self,name):
        # nome do grupo
        self.name = name
        
        self.info_path = f"data/groups_info/{self.name}"
        self.tsv_path = self.info_path + '/tsv'
        
        # buscando o id do grupo
        with open('data/groups/groups_name_id.json') as log:
            data = json.load(log)
            self.pat_id = data[self.name]
        
        # url para extrair as informações em tsv
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'

        # criar repositório para os dados
        if not os.path.exists(self.info_path):
            os.mkdir(self.info_path)
            os.mkdir(self.tsv_path)
    
    def checkDataUpdate(self):
        '''
        Verifica se a versão do .tsv salvo precisa ser atualizada para uma nova
        True: Precisa ser atualizada
        False: Não precisa ser atualizada
        '''
        tsv_log = os.listdir(self.tsv_path)
        if tsv_log == []:
            return True
        tsv_log = tsv_log[0]
        tsv_version = tsv_log.split('|')[1].rstrip('.tsv')
        if tsv_version != self.pat_id:
            os.remove(self.tsv_path+'/'+tsv_log)
            return True
        else:
            return False

    def getPatData(self):
        self.table_url = f'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{self.name}/latest_kmer/Metadata/{self.pat_id}'+'.metadata.tsv'
        download(self.table_url,f'{self.tsv_path}/{self.name}|{self.pat_id}.tsv')

with open("data/groups/groups_name_id.json","r") as file:
    names = json.load(file)
    for name in names.keys():
        a = Group(name)
        if a.checkDataUpdate():
            a.getPatData()        