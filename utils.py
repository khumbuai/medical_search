

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
