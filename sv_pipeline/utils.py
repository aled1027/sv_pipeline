"""A module of utilities"""

import itertools
import time

class Timer(object):
    # from : http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()
        if self.name:
            print('==> [Entering: %s]' % self.name)

    def __exit__(self, type, value, traceback):
        if self.name:
            print('<== [Exiting : %s] ' % self.name, end='')
        print('Elapsed: %s' % (time.time() - self.tstart))

def pairwise(iterable):
        """s -> (s0,s1), (s1,s2), (s2, s3), ...

        https://docs.python.org/dev/library/itertools.html
        """
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

def flatten(listOfLists):
    "Flatten one level of nesting"
    return itertools.chain.from_iterable(listOfLists)




