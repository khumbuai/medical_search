uid = 'PMID'


import requests

db = 'pubmed'
query = 'asthma[mesh]+AND+leukotrienes[mesh]+AND+2009[pdat]'
base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

url = base + "esearch.fcgi?db={}&term={}&usehistory=y".format(db,query)

r = requests.get(url)


from Bio import Entrez

Entrez.email = "A.N.Other@example.com"
handle = Entrez.esearch(db="pubmed", term="biopython")
record = Entrez.read(handle)
print(record["IdList"])

id = '30013827'

handle = Entrez.esummary(db="pubmed", id=id)
record = Entrez.read(handle)

#handle = Entrez.efetch(db="pubmed", id=id, rettype="gb", retmode="text")
handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
r = Entrez.read(handle)
print(handle.read())