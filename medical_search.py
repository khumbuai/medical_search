from pubmed_handler import search_ncbi
from ranker import rank_by_semantic_relevance

def medical_search(query, topn = 5):

    results = search_ncbi(query, db='pubmed', reldate=365, max_results=115)
    ranked_results = rank_by_semantic_relevance(query,results)
    return ranked_results[:topn]

if __name__ == '__main__':

    QUERY = 'lung cancer china'

    results = medical_search(QUERY)
    for result in results:
        print((result['relevance'],result['TI'],result['AB'][:200]))