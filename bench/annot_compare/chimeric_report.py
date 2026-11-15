#!/usr/bin/env python3
"""Row 1708 (docs/PREREG_join_rules_2026-11-14.md): report of chimeric molecules, no precision or recall. Input = one line per
supplementary record (name, chrom, pos, strand, primary-part chrom, pos, strand, mapq), from the job's samtools pass.
Each part is put on the RefSeq gene whose span holds its start (any strand, the shortest such gene); a molecule joins a gene
pair when its two parts fall on two different genes. Candidate = >= 5 molecules here and >= 2 in the witness on that pair.
usage: chimeric_report.py SAMPLE.supp.tsv WITNESS.supp.tsv REFSEQ.gtf[.gz] OUT.tsv"""
import bisect, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fusion_audit as FA


def locator(genes):
    """genes: name -> (chrom, start, end). Returns f(chrom, pos) -> name of the shortest gene holding pos, or None."""
    by = collections.defaultdict(list)
    for n, (c, s, e) in genes.items(): by[c].append((s, e, n))
    mx = {}
    for c, v in by.items(): v.sort(); mx[c] = max(e - s for s, e, _ in v)
    st = {c: [x[0] for x in v] for c, v in by.items()}

    def f(c, p):
        v = by.get(c)
        if not v: return None
        i = bisect.bisect_right(st[c], p); best = None
        for s, e, n in reversed(v[bisect.bisect_left(st[c], p - mx[c]):i]):
            if s <= p < e and (best is None or e - s < best[0]): best = (e - s, n)
        return best and best[1]
    return f


def pairs(path, loc):
    P = collections.defaultdict(set); c = collections.Counter(); meta = {}
    for l in open(path):
        q, c1, p1, s1, c2, p2, s2, _ = l.rstrip('\n').split('\t')
        if not c2: continue
        c['supplementary records'] += 1; a, b = loc(c1, int(p1) - 1), loc(c2, int(p2) - 1)
        d = 'other chromosome' if c1 != c2 else 'same chromosome < 1 Mb' if abs(int(p1) - int(p2)) < 1000000 else 'same chromosome >= 1 Mb'
        if 'chrM' in (c1, c2) and c1 != c2: c['one part on chrM'] += 1
        if a is None or b is None: c['a part outside any gene'] += 1; continue
        if a == b: c['both parts in one gene'] += 1; continue
        k = tuple(sorted((a, b))); P[k].add(q); m = meta.setdefault(k, collections.Counter()); m[d] += 1; m['same strand' if s1 == s2 else 'opposite strand'] += 1
    c['molecules joining two genes'] = len(set().union(*P.values())) if P else 0; c['gene pairs'] = len(P)
    return P, meta, c


def main(sample, witness, rs, out):
    R, rmeta = FA.load(rs, 'gene_id'); genes = {rmeta[g][2].removeprefix('gene-'): (c, ex[0][0], ex[-1][1]) for g, (c, s, ex) in R.items()}
    loc = locator(genes); P, meta, c = pairs(sample, loc); Wp, _, _ = pairs(witness, loc)
    for k in c: print(k, c[k], sep='\t')
    rows = sorted(((len(v), len(Wp.get(k, ())), k) for k, v in P.items()), reverse=True)
    big = [r for r in rows if r[0] >= 5]; cand = [r for r in big if r[1] >= 2]
    print(f'gene pairs with >= 5 molecules {len(big)}; candidates (also >= 2 molecules in the witness) {len(cand)}')
    dom = lambda k: max(('other chromosome', 'same chromosome < 1 Mb', 'same chromosome >= 1 Mb'), key=lambda d: meta[k][d])
    print('candidates by distance:', dict(collections.Counter(dom(k) for _, _, k in cand)), '; involving chrM genes', sum(any(genes[g][0] == 'chrM' for g in k) for _, _, k in cand),
          '; name-stem relatives', sum(FA.stem(k[0]) == FA.stem(k[1]) for _, _, k in cand))
    with open(out, 'w') as fh:
        print('gene_a\tgene_b\tmolecules\twitness_molecules\tcandidate\tdistance\tsame_strand\topposite_strand\tchrom_a\tchrom_b', file=fh)
        for n, w, k in big: print(*k, n, w, int(w >= 2), dom(k), meta[k]['same strand'], meta[k]['opposite strand'], genes[k[0]][0], genes[k[1]][0], sep='\t', file=fh)


def _selfcheck():
    f = locator({'A': ('c', 0, 100), 'B': ('c', 40, 60), 'C': ('d', 5, 9)})
    assert f('c', 50) == 'B' and f('c', 10) == 'A' and f('c', 100) is None and f('d', 5) == 'C' and f('e', 1) is None


if __name__ == '__main__':
    _selfcheck(); main(*sys.argv[1:5])
