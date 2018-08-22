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


def try_except(return_value=None):
    '''
    Decorator which catches exceptions and returns return_value in case an exception occurred.
    :param return_value:
    :return:
    '''
    def func_decorator(func):

        def inner(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except:
                return return_value
        return inner

    return func_decorator
