import os
import wget

class Group():
    def __init__(self,name):
        self.name = name

        # decidir como organizar essas informações
        self.id_path = f"/data/groups_info/{self.name}/"
        self.assembly_path = f"/data/assemblies/"
    
    def get_pat_data():
        # create url
        # check new versions
        # fetch info
        pass