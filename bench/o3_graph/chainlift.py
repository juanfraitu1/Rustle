"""chainlift.py — minimal UCSC chain liftover (KEY=o2_presence; copied into the repo 2026-10-18 for hg002truth.py). Chain t = source, q = destination.
Chain(path).lift(chrom, s, e) -> list of (qname, qs, qe, strand) ungapped pieces (0-based half-open, q on + strand)."""
import gzip, bisect, collections
class Chain:
    def __init__(self, path):
        self.blocks = collections.defaultdict(list)  # tname -> [(ts, te, qname, qs_fwd_start, qstrand, qsize)]
        op = gzip.open if path.endswith('.gz') else open
        with op(path, 'rt') as f:
            hdr = None
            for ln in f:
                x = ln.split()
                if not x:
                    continue
                if x[0] == 'chain':
                    hdr = x; tname = x[2]; tpos = int(x[5]); qname = x[7]; qsize = int(x[8]); qstr = x[9]; qpos = int(x[10])
                    continue
                size = int(x[0])
                if qstr == '+':
                    self.blocks[tname].append((tpos, tpos + size, qname, qpos, '+'))
                else:  # q coords on the reverse strand: forward = qsize - (qpos + size)
                    self.blocks[tname].append((tpos, tpos + size, qname, qsize - qpos - size, '-'))
                if len(x) == 3:
                    tpos += size + int(x[1]); qpos += size + int(x[2])
        for t in self.blocks:
            self.blocks[t].sort()
        self.starts = {t: [b[0] for b in v] for t, v in self.blocks.items()}
    def lift(self, chrom, s, e):
        out = []
        if chrom not in self.blocks:
            return out
        B = self.blocks[chrom]; S = self.starts[chrom]
        i = max(0, bisect.bisect_right(S, s) - 1)
        while i < len(B) and B[i][0] < e:
            ts, te, qn, qf, st = B[i]
            a, b = max(s, ts), min(e, te)
            if a < b:
                if st == '+':
                    out.append((qn, qf + (a - ts), qf + (b - ts), '+'))
                else:
                    out.append((qn, qf + (te - b), qf + (te - a), '-'))
            i += 1
        return out
