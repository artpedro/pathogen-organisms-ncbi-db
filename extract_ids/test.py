from Bio import Entrez
Entrez.email = "arturpedromartins@gmail.com"
handle = Entrez.esearch(db="assembly", term="Salmonella_Enterica[Orgn] and \"complete genome\"[filter]",idtype="acc")
record = Entrez.read(handle)
print(record)