"""helpers shared by cat_vs_refseq.py and merge_audit.py."""
import bisect
import collections


def merge(iv):
    o = []
    for a, b in sorted(iv):
        if o and a <= o[-1][1]:
            o[-1][1] = max(o[-1][1], b)
        else:
            o.append([a, b])
    return o


class Idx:
    def __init__(self, items):  # items: dict id -> (chrom, strand, exons)
        self.by = collections.defaultdict(list)
        for k, (c, s, ex) in items.items():
            self.by[(c, s)].append((ex[0][0], ex[-1][1], k, ex))
        self.lo = {}
        self.ml = {}
        for k, v in self.by.items():
            v.sort()
            self.lo[k] = [x[0] for x in v]
            self.ml[k] = max(x[1] - x[0] for x in v)

    def hits(self, c, s, ex):
        v = self.by.get((c, s))
        if not v:
            return []
        lo, hi = ex[0][0], ex[-1][1]
        i = bisect.bisect_left(self.lo[(c, s)], lo - self.ml[(c, s)])
        o = []
        while i < len(v) and v[i][0] < hi:
            a, b, k, e2 = v[i]
            i += 1
            if b > lo and share(ex, e2):
                o.append(k)
        return o


def share(x, y):
    i = j = 0
    while i < len(x) and j < len(y):
        if x[i][0] < y[j][1] and y[j][0] < x[i][1]:
            return True
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return False
