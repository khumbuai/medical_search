"""
Script for ranking a list of pubmed results by date, ft, ...
"""

from gensim.models.fasttext import FastText
from utils import preprocess
import numpy as np
import nltk
from nltk.corpus import stopwords

ft_model = FastText.load('models/ft/med_model_dim300_win5_min100.bin')
stops = set(stopwords.words('english'))

def _remove_stopwords(sentence):
    """

    :param sentence: list of words
    :return: list of words
    """
    if isinstance(sentence, list):
        return [word for word in sentence if not word in stops]

    else:
        sentence = sentence.split()
        sentence = [word for word in sentence if not word in stops]
        sentence = ' '.join(sentence)
        return sentence



def similarity(w1, w2):
    """
    helper function to catch oov error
    :param w1: string
    :param w2: string
    :return:
    """
    try:
        sim = ft_model.similarity(w1, w2)
    except KeyError:
        sim = 0
    return sim


def get_abstract_similarity(query, abstract_tokens, topn=10):
    """

    :param query: string
    :param abstract_tokens: list of abstract tokens
    :param topn: max of top similar words to consider
    :return:
    """

    vecs = [similarity(query, word) for word in abstract_tokens]
    topn_mean = np.mean(list(sorted(vecs))[-topn:])
    return topn_mean


def rank_by_semantic_relevance(query, list_of_dictionaries):
    """

    :param list_of_dictionaries: [{'title':bla,'abstract':bla bla,...},{}]
    :return:
    """

    w = [0.5, 0.5]  # weighting between title and abstract content

    query = preprocess(query)  # to be consistent with ft training

    titles = [preprocess(item['TI']) for item in list_of_dictionaries]
    titles = [_remove_stopwords(title) for title in titles]
    title_similarities = np.array([similarity(query, title) for title in titles])

    abstracts = [nltk.word_tokenize(preprocess(item['AB'])) for item in list_of_dictionaries]
    abstracts = [_remove_stopwords(abstract) for abstract in abstracts]
    abstract_similarities = np.array([get_abstract_similarity(query, abstract_tokens) for abstract_tokens in abstracts])

    total_similarities = w[0] * title_similarities + w[1] * abstract_similarities
    ranking = np.argsort(total_similarities)[::-1]  # highest rating first

    # add similarity to list of dictionaries
    for i, item in enumerate(list_of_dictionaries):
        item['relevance'] = total_similarities[i]

    # re-order by relevance rank
    ranked_list_of_dictionaries = [list_of_dictionaries[ind] for ind in ranking]

    return ranked_list_of_dictionaries


if __name__ == '__main__':

    # test case
    a1 = 'Prostate cancer is the development of cancer in the prostate, a gland in the male reproductive system.[6] Most prostate cancers are slow growing; however, some grow relatively quickly.[1][3] The cancer cells may spread from the prostate to other area of the body, particularly the bones and lymph nodes.[7] It may initially cause no symptoms.[1] In later stages, it can lead to difficulty urinating, blood in the urine or pain in the pelvis, back, or when urinating.[2] A disease known as benign prostatic hyperplasia may produce similar symptoms.[1] Other late symptoms may include feeling tired due to low levels of red blood cells.[1]'
    a2 = 'Cancer is a group of diseases involving abnormal cell growth with the potential to invade or spread to other parts of the body.[2][8] These contrast with benign tumors, which do not spread to other parts of the body.[8] Possible signs and symptoms include a lump, abnormal bleeding, prolonged cough, unexplained weight loss and a change in bowel movements.[1] While these symptoms may indicate cancer, they may have other causes.[1] Over 100 types of cancers affect humans.[8]'
    a3 = 'Canids are found on all continents except Antarctica, having arrived independently or accompanied human beings over extended periods of time. Canids vary in size from the 2-m-long (6 ft 7 in) gray wolf to the 24-cm-long (9.4 in) fennec fox. The body forms of canids are similar, typically having long muzzles, upright ears, teeth adapted for cracking bones and slicing flesh, long legs, and bushy tails. They are mostly social animals, living together in family units or small groups and behaving cooperatively. Typically, only the dominant pair in a group breeds, and a litter of young is reared annually in an underground den. Canids communicate by scent signals and by vocalizations. They are very intelligent. One canid, the domestic dog, long ago entered into a partnership with humans and today remains one of the most widely kept domestic animals.'

    list_of_dictionaries = [{'TI': 'i have cancer', 'AB': a1}, {'TI': 'i love dogs', 'AB': a3},
                            {'TI': 'i hate cancers', 'AB': a2}]
    query = 'something about cancer'

    ranked_list_of_dictionaries = rank_by_semantic_relevance(query, list_of_dictionaries)
    print([(item['title'], item['relevance']) for item in ranked_list_of_dictionaries])
