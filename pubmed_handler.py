from Bio import Entrez
from Bio import Medline
from urllib.error import HTTPError  # for Python 3

from utils import try_except

EMAIL = 'example@example.com'
Entrez.email = EMAIL


def _preprocess_query(query, db):
    handle = Entrez.espell(term=query, db=db, usehistory="y")
    record = Entrez.read(handle)
    return_query = record["CorrectedQuery"] if len(record["CorrectedQuery"] ) != 0 else query
    return return_query

def _postprocess_results(results):
    '''
    Transforms string entries in the search results to ints, if applicable
    :param List[Dict] results:
    :return:
    :rtype List[Dict]:
    '''
    for result in results:
        for key, value in result.items():
            try:
                result[key] = int(value)
            except (ValueError, TypeError):
                pass
    return results


@try_except(return_value=[])
def search_ncbi(query, autocorrect=False, db='pubmed', reldate=365, max_results=500):
    '''
    Searches the ncbi databank for the specific query.
    Note that autocorrect=True uses another request to the Entrez webserver.

    The parameters which can be used for the Entrez methods can be found on
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
    The meaning of the keys in the results can be found on
    https://www.nlm.nih.gov/bsd/mms/medlineelements.html

    :param Str db: Database to be searched
    https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    :param Int reldate: The search returns only those items that are no older than reldate days.
    :return: List of search results, containing (among other informations) the title and abstract of the respective
    papers.
    :rtype List[Dict]:
    '''
    if autocorrect:
        query = _preprocess_query(query, db)

    handler = Entrez.esearch(db=db,
                             term=query,
                             reldate=reldate,
                             datetype="pdat",
                             sort='Best match',
                             usehistory="y")

    search_results = Entrez.read(handler)
    results_as_text = []

    print(search_results)

    # Download in batches, since url may break if max_results is large.
    batch_size = 25
    count = min(max_results, int(search_results["Count"]))
    steps = list(range(0, count, batch_size))
    batch_sizes = [batch_size for _ in steps]
    # Last batch_size is chosen such that the total number of retrieved documents is equal to max_results
    if len(steps) > 0:  # catch empty search results
        batch_sizes[-1] = count - steps[-1]

    for start, batch_size in zip(steps, batch_sizes):
        try:
            fetch_handle = Entrez.efetch(db=db,
                                         rettype="medline",
                                         retmode="text",
                                         retstart=start,
                                         retmax=batch_size,
                                         webenv=search_results["WebEnv"],
                                         query_key=search_results["QueryKey"])
            results_as_text += fetch_handle.read().split('\n')
        except HTTPError as err:
            pass

    results = [dict(result) for result in Medline.parse(results_as_text)]

    results = _postprocess_results(results)
    return results


if __name__ == '__main__':
    results = search_ncbi('lung cancer', db='pubmed', reldate=365, max_results=115)
    for result in results:
        print(result)
