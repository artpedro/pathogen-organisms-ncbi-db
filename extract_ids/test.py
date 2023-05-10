import json
from Bio import Entrez
search_handle = Entrez.esearch(db="assembly", term="(\"Salmonella\"[Organism] OR Salmonella[All Fields]) AND (latest[filter] AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])", usehistory="y", idtype="acc",retmax='100000')
search_results = Entrez.read(search_handle)
search_handle.close()

idlist = search_results["IdList"]
handle = Entrez.esummary(db='assembly', id=idlist)
print('handle gerado')
record = Entrez.read(handle)
print('record gerado')
with open('record.txt','w') as txt:
    record_json = json.dumps(record)
    txt.write(record_json)
    print(record_json)