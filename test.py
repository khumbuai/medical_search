from Bio import Entrez
import time

from urllib.error import HTTPError  # for Python 3
from Bio import Medline

def save_search_items():
    Entrez.email = "bjornjobb@gmail.com"
    search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                term="Opuntia[ORGN]",
                                                reldate=365, datetype="pdat",
                                                sort='Best match',
                                                usehistory="y"))

    count = int(search_results["Count"])
    print("Found %i results" % count)

    print(search_results)

    batch_size = 10
    out_handle = open("recent_orchid_papers.txt", "w")
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        try:
            fetch_handle = Entrez.efetch(db="pubmed", rettype="medline",
                                         retmode="text",retstart=start,
                                         retmax=batch_size,
                                         webenv=search_results["WebEnv"],
                                         query_key=search_results["QueryKey"])
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
            else:
                raise
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()

def transorm_text():

    with open("recent_orchid_papers.txt", "r") as f:
        text = f.readlines()

    entries = Medline.parse(text)

    return entries

save_search_items()

x = transorm_text()

for item in x:
    print(item)


