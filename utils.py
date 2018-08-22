import nltk
import re

PUNCTUATION_TOKENS = ['.', '..', '...', ',', ';', ':', '(', ')', '"', '\'', '[', ']', '{', '}', '?', '!', '-', u'–',
                      '+', '*', '--', '\'\'', '``', "'"]
PUNCTUATION = '?.!/;:()&+'

def preprocess(text):

    words = nltk.word_tokenize(text)
    words = [x.lower() for x in words]
    words = [x for x in words if x not in PUNCTUATION_TOKENS]
    words = [re.sub('[' + PUNCTUATION + ']', '', x) for x in words]
    text = ' '.join(words)
    """
    function that uses same preprocessing as while ft training
    :param text:
    :return:
    """

    return text