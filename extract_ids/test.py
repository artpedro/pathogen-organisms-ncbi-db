from Bio import Entrez
Entrez.email = "arturpedromartins@gmail.com"
search_handle = Entrez.esearch(db="assembly", term="Salmonella_Enterica[Orgn] and \"complete genome\"[filter]", usehistory="y", idtype="acc")
search_results = Entrez.read(search_handle)
search_handle.close()
print(search_results['IdList'])

webenv = record["WebEnv"]
querykey = record["QueryKey"]
print(webenv,querykey)
fetch_handle = Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        webenv=webenv,
        query_key=querykey,
        idtype="acc",
    )
data = fetch_handle.read()
fetch_handle.close()
print(data)