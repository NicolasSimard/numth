import os, json, atexit

class cache_calls:
    """This decorator caches values of one parameter functions (e.g.
    arithmetic functions). For example, the function

    @cache_calls
    def f(n):
        ...
        return ...

    will have all its values cached, even between sessions. For example,
    if f(5) is computed in a python sesion, it will be stored in the file
    f.cache in the flder func_cache and if f(5) is needed in a future
    session, it will be read from the cached file (instead of being
    re-computed).
    """

    def __init__(self, func):
        """Initialise the decorator with the cached values."""

        if not os.path.exists('func_cache'):
            os.makedirs('func_cache')
        self.cache_file = 'func_cache/' + func.__name__ + '.cache'
        if os.path.isfile(self.cache_file):
            with open(self.cache_file,'r') as f:
                self.cached_calls = json.load(f)
        else:
            self.cached_calls = {}
        self.func = func
        self.__doc__ = func.__doc__
        atexit.register(self.close) # Make sure the close function is called

    def __call__(self,n):
        # When json serializes a dictionnary, it stores the keys as strings,
        # since an integer cannot be the key of a dictionnary in Javascript.
        # We therefore convert the key to a string.
        if str(n) in self.cached_calls:
            return self.cached_calls[str(n)]
        else:
            self.cached_calls[str(n)] = self.func(n)
            return self.cached_calls[str(n)]

    def close(self):
        """Dump the calls in the cache file when session is exited.

        When the python session is exited, the cached calls (stored in the
        dictionnary self.cached_calls) are stored in the cache file for later
        use.
        """

        json.dump(self.cached_calls, open(self.cache_file,'w'))


class cache_list:
    """This decorator caches the list of a function which outputs a list.

    Suppose for example that you create a function fibo_list(n) that returns a
    list of the first n fibonacci numbers. If you call f(100) and then f(101),
    you have to re-compute the first 100 fibonacci numbers twice. What this
    decorator does is that it stores the output of f(100), which is a list, so
    that it can be used when f(101) is called. This can save a lot of time!
    """

    def __init__(self, func):
        """Initialise the decorator with the cached list."""

        if not os.path.exists('func_cache'):
            os.makedirs('func_cache')
        self.cache_file = 'func_cache/' + func.__name__ + '.cache'
        if os.path.isfile(self.cache_file):
            print("Initialising cache...")
            with open(self.cache_file,'r') as f:
                self.cached_list = json.load(f)
        else:
            self.cached_list = []
        self.initial_size = len(self.cached_list)
        self.func = func
        self.__doc__ = func.__doc__
        atexit.register(self.close) # Make sure the close function is called

    def __call__(self, n, known = []):
        if len(known) > len(self.cached_list):
            self.cached_list = known
        if n > len(self.cached_list):
            self.cached_list = self.func(n, self.cached_list)
        return self.cached_list[:n]

    def close(self):
        """Dump the list in the cache file when session is exited.

        When the python session is exited, the cached calls (stored in the
        dictionnary self.cached_calls) are stored in the cache file for later
        use.
        """
        if len(self.cached_list) > self.initial_size:
            print("Saving new cache...")
            json.dump(self.cached_list, open(self.cache_file,'w'))
