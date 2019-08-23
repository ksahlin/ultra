from itertools import product

# modified code from https://stackoverflow.com/questions/41085553/segment-tree-implementation-in-python
# to support updades

class RangeQuery(object):
    "Data structure providing efficient range queries."

    def __init__(self, items, fn):
        """Build a RangeQuery object for a sequence of items.

        fn -- function taking two items and returning their query
        result, for example "min" to query the range-minimum. It must be
        associative, commutative, and idempotent.

        """
        # Mapping from (start, step) to reduce(fn, items[start:start + 2**step])
        self._rq = rq = {(i, 0): item for i, item in enumerate(items)}
        print(rq)
        self._fn = fn
        n = len(items)
        print(list(product(range(1, n.bit_length()), range(n))))
        for step, i in product(range(1, n.bit_length()), range(n)):
            j = i + 2 ** (step-1)
            if j < n:
                rq[i, step] = fn(rq[i, step-1], rq[j, step-1])
            else:
                rq[i, step] = rq[i, step-1]

    def query(self, start, stop):
        "Return reduce(fn, items[start:stop])."
        j = (stop - start).bit_length() - 1
        x = self._rq[start, j]
        y = self._rq[stop - 2 ** j, j]
        return self._fn(x, y)

import random
import unittest

class TestRangeQuery(unittest.TestCase):
    def test_range_minmax(self):
        n = 1000
        items = [random.randrange(n) for _ in range(n)]
        rqmin = RangeQuery(items, min)
        rqmax = RangeQuery(items, max)
        for _ in range(n):
            start, stop = sorted(random.sample(range(n + 1), 2))
            self.assertEqual(rqmin.query(start, stop), min(items[start:stop]))
            self.assertEqual(rqmax.query(start, stop), max(items[start:stop]))



n = 10
items = [random.randrange(n) for _ in range(n)]
print(items)


rqmin = RangeQuery(items, min)
start, stop = sorted(random.sample(range(n + 1), 2))
print(start, stop)
print(rqmin.query(start, stop))

rqmax = RangeQuery(items, max)
print(rqmax.query(start, stop))