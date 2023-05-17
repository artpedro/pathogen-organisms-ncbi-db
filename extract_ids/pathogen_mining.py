import os
import wget
from wget import download


# Downloader class.
class downloader():
    def progressBar(self,current,total):
        '''
        Barra de progresso para indicar o quanto falta para baixar
        '''
        print("Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total))
        
    # Create a downloadfile method
    # Accepting the url and the file storage location
    # Set the location to an empty string by default.

    def downloadFile(self, url, location=""):
        '''
        Método para baixar um arquivo dado um url e um caminho para
        armazenar o arquivo baixado. Por padrão, o caminho é o diretório atual
        '''
         # Baixar o arquivo com uma barra de progresso
        download(url, out = location, bar = self.progressBar)


class Group():
    def __init__(self,name):
        self.name = name
        # decidir como organizar essas informações
        self.id_path = f"data/groups_info/{self.name}"
        
        # criar repositório para os dados
        if not os.path.exists(self.id_path):
            os.mkdir(self.id_path)

    def get_pat_data():
        # create url
        # check new versions
        # fetch info
        pass

a = Group("teste2")