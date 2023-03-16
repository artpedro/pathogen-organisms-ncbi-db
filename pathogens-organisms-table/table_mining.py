
import time
from selenium import webdriver
import pandas as pd
from bs4 import BeautifulSoup

# Inicia o driver
driver = webdriver.Chrome('usr/local/bin/chromedriver')

# Link para a tabela dos organismos patogênicos
url = "https://www.ncbi.nlm.nih.gov/pathogens/organisms/"

# Abrindo o site
driver.get(url)
time.sleep(5)


# Codigo-fonte da pagina
html = driver.page_source
soup = BeautifulSoup(html, 'html.parser')

# Localizando uma tabela no html extraido
table = soup.find('table', {'class': 'home_table'})

print(str(table))

with open("pathogens-organisms-table/table.html", "w") as html:
    html.write(str(table))


headers = []
rows = []

# Extrai o cabeçalho da tabela

index=0 
for th in table.find_all('th'):
    headers.append(th.text.strip())
    index += 1
    if index == 9:
        break
    print(headers)

# Extrai as linhas da tabela
for tr in table.find_all('tr')[1:]:
    row = []
    species = tr.find_all('th')
    row.append(species[0].text.strip())

    for td in tr.find_all('td'):
        row.append(td.text.strip())
    rows.append(row)

# Cria um dataframe do pandas a partir dos cabeçalhos e linhas extraídas da tabela
df = pd.DataFrame(rows, columns=headers)

df.to_csv('pathogens-organisms-table/organisms_table.csv',sep=',',index=False)