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

def DELETE_ALL():
    db = client["pathogen_db"]
    col = db["groups_versions"]
    col.delete_many({})

    col = db["groups_info"]
    col.delete_many({})

    col = db["all_groups_data"]
    col.delete_many({})

def groups_version_to_mongo(client):

    # acessando a coleção de versoes no banco de dados
    db = client["pathogen_db"]
    col = db["groups_versions"]
    cur = col.find({})

    # armazenando quais grupos já estão no banco de dados
    log = {}
    for i in cur:
        log[i['group']] = i['pathogen_id']
        print(i)
    
    # recuperando as versões de data
    versions = readGroupsNames()

    # atualizando ou adicionando versões
    for group,id in versions.items():
        if (group not in log) or (log[group] != id):
            entry = {'group':group,'pathogen_id':id}
            col.update_one({'group':group},{"$set":entry},upsert=True)
    
def groups_info_to_mongo(client):
    # em desenvolvimento

    db = client["pathogen_db"]
    col = db["all_groups_data"]
    cur = col.find({})

    # armazenando quais grupos já estão no banco de dados
    log = [i['asm_acc'] for i in cur]

    for group,id in readGroupsNames().items():        
        # checando se existe informação sobre
        if os.path.exists(f'data/groups_info/{group}/{group}_filtered.json'):
            # carregando informação
            with open(f'data/groups_info/{group}/{group}_filtered.json','r') as json:
                info = js.load(json)

            # para cada entrada de patogêno
            for entry in info:
                if len(entry) == 2:
                    continue
                entry['group'] = group

                # se o acesso ao assembly dessa entrada não estiver no banco, adiciona-lo
                if entry['asm_acc'] not in log:
                    col.update_one({'asm_acc':entry['asm_acc']},{'$set':entry},upsert=True)
                    #print(f"{entry['asm_acc']} added to mongo")

def force_info_to_mongo(client):
    # igual a anterior, mas não checa se já está no banco, só substitui

    db = client["pathogen_db"]
    col = db["all_groups_data"]
    
    for group in readGroupsNames():
        
        if os.path.exists(f'data/groups_info/{group}/{group}_filtered.json'):

            with open(f'data/groups_info/{group}/{group}_filtered.json','r') as json:
                info = js.load(json)

            for entry in info:
                if len(entry) == 2:
                    continue
                entry['group'] = group
                col.update_one({'asm_acc':entry['asm_acc']},{'$set':entry},upsert=True)
                print(f"{entry['asm_acc']} added to mongo")
 
# force_info_to_mongo(client)