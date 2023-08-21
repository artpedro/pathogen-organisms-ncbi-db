import json as js
import os
import pandas as pd
import sys
import time
from pymongo import MongoClient

client = MongoClient("localhost",27017)

def readGroupsNames():
    '''
    Lê o nome dos patógenos presentes no Pathogen no NCBI
    '''
    with open(os.path.normpath("data/groups/groups_name_id.json"),"r") as data:
        groups = js.load(data)
        return groups

def groups_version_to_mongo(client):

    db = client["pathogen_db"]
    col = db["groups_versions"]

    with open('data/groups/groups_name_id.json','r') as json:
        versions = js.load(json)

    data = []
    for group,id in versions.items():
        data.append({'group':group,'pathogen_id':id})
    
    col.insert_many(data)
    

def groups_info_to_mongo(client):
    # em desenvolvimento

    db = client["pathogen_db"]
    col = db["groups_data"]

    for group in readGroupsNames():
        
        if os.path.exists(f'data/groups_info/{group}/{group}_filtered.json'):
            
            with open(f'data/groups_info/{group}/{group}_filtered.json','r') as json:
                info = js.load(json)
            
            data = []
            idx = 0
            for entry in info:
                if idx == 0:
                    idx += 1
                    #ignorando primeiro elemento de tag
                    continue

                # por enquanto !!!
                entry.pop('date')

                for key in entry:
                    if entry[key] == None:
                        entry[key] = key
                print(entry)
'''            
with open(f'data/groups_info/Edwardsiella_tarda/Edwardsiella_tarda_filtered.json','r') as json:
    info = js.load(json)
    data = []
    idx = 0
    for entry in info:
        entry.pop('date')
        if idx == 0:
            idx += 1
            #ignorando primeiro elemento de tag
            continue    
        for key in entry:
            if entry[key] == None:
        print(entry)
print(info)
'''