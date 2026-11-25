#!/usr/bin/env python3
"""What the member places that are on no listed relative are (register 1956).

    place_classes.py MEMBERS.tsv ASSEMBLY.gtf RELATIONS.tsv HGNC.txt RMSK.bed [LO HI]...   (default bands: 0.9 1.01)

Needs DENOVO_CAT_GFF / DENOVO_CAT_TSV. Classes of a (place, source) row, first match wins: listed relative; relative at one
remove (share a listed relative); same named family (HGNC gene group, or symbol without the pseudogene suffix); repeat
(>= half of the footprint in an interspersed or simple repeat); another gene by biotype; nothing annotated."""
import sys, os, re, bisect, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
REP = {'SINE', 'LINE', 'LTR', 'DNA', 'Retroposon', 'Simple_repeat', 'Low_complexity', 'Satellite'}

def merged(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out

def covered(blocks, starts, ivs):
    """bases of `blocks` inside the merged, sorted intervals `ivs`"""
    n = 0
    for a, b in blocks:
        k = max(0, bisect.bisect_right(starts, a) - 1)
        while k < len(ivs) and ivs[k][0] < b: n += max(0, min(b, ivs[k][1]) - max(a, ivs[k][0])); k += 1
    return n

def stem(name): return re.sub(r'P\d+$', '', name or '')

def klass(src, under, rel, grp, name, bio, rep_share):
    """src = source gene or None; under = genes under the place"""
    if src is None: return 'source on no annotated gene: ' + ('repeat' if rep_share >= 0.5 else 'not a repeat')
    under = under - {src}
    if any(g in rel[src] for g in under): return '1 listed relative'
    if any(rel[src] & rel[g] for g in under): return '2 relative at one remove'
    if any((grp[src] & grp[g]) or (stem(name[src]) and stem(name[src]) == stem(name[g])) for g in under): return '3 same named family'
    if rep_share >= 0.5: return '4 repeat'
    if not under: return '6 nothing annotated'
    b = {bio[g] for g in under}
    return '5 another gene: ' + ('protein-coding' if 'protein_coding' in b else 'pseudogene' if any('pseudo' in x for x in b) else 'long non-coding' if b & {'lncRNA', 'lincRNA'} else 'other')

def main(mem, gtf, relp, hgnc, rmsk, *bands):
    import famtruth as F
    bands = [(float(bands[i]), float(bands[i + 1])) for i in range(0, len(bands), 2)] or [(0.9, 1.01)]
    CAT = F.load_cat(); idx = F.exon_index(CAT); ex = collections.defaultdict(list)
    for g, d in CAT.items():
        for a, b in d['exons']: ex[d['chrom']].append((a, b, g))
    for v in ex.values(): v.sort()
    ST = {c: [x[0] for x in v] for c, v in ex.items()}; MX = {c: max(b - a for a, b, _ in v) for c, v in ex.items()}
    def genes_at(c, blocks):
        o = set(); v = ex.get(c, [])
        for s, e in blocks:
            k = bisect.bisect_left(ST.get(c, []), s - MX.get(c, 0))
            while k < len(v) and v[k][0] < e:
                if v[k][1] > s: o.add(v[k][2])
                k += 1
        return o
    rel = collections.defaultdict(set)
    for l in list(open(relp))[1:]:
        f = l.rstrip('\n').split('\t'); rel[f[0]].add(f[2]); rel[f[2]].add(f[0])
    H = [l.rstrip('\n').split('\t') for l in open(hgnc)]; h = {k: i for i, k in enumerate(H[0])}; by_ens = {}; by_sym = {}
    for r in H[1:]:
        gs = set(filter(None, r[h['gene_group_id']].strip('"').split('|'))); by_ens[r[h['ensembl_gene_id']]] = gs; by_sym[r[h['symbol']]] = gs
    grp = collections.defaultdict(set); name = collections.defaultdict(str); bio = collections.defaultdict(str)
    for g, d in CAT.items():
        name[g] = d['name']; bio[g] = d['biotype']; grp[g] = set(by_sym.get(d['name'], set()))
        for e in d['ens']: grp[g] |= by_ens.get(e, set())
    rp = collections.defaultdict(list)
    for l in open(rmsk):
        f = l.split('\t', 8)
        if f[6] in REP: rp[f[0]].append((int(f[1]), int(f[2])))
    rp = {c: merged(v) for c, v in rp.items()}; rs = {c: [a for a, _ in v] for c, v in rp.items()}
    tx = collections.defaultdict(list); meta = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); tx[g].append((int(f[3]) - 1, int(f[4]))); meta[g] = (f[0], f[6] if f[6] in '+-' else '+')
    pairs = open(os.environ['PLACE_CLASSES_PAIRS'], 'w') if os.environ.get('PLACE_CLASSES_PAIRS') else None   # one line per row: band, class, place, source gene, genes under, repeat share
    src_gene = {}; rows = [l.rstrip('\n').split('\t') for l in open(mem)][1:]
    for lo, hi in bands:
        C = collections.Counter(); P = collections.defaultdict(set); R = collections.Counter(); ex_ = collections.defaultdict(list)
        for p in rows:
            if not lo <= float(p[7]) < hi: continue
            s = p[5]
            if s not in src_gene: src_gene[s] = F.map_copy(idx, meta[s][0], meta[s][1], sorted(set(tx[s]))) if s in meta else None
            bl = [tuple(map(int, x.split('-'))) for x in p[12].split(',')]; share = covered(bl, rs.get(p[1], []), rp.get(p[1], [])) / max(1, sum(b - a for a, b in bl))
            under = genes_at(p[1], bl); k = klass(src_gene[s], under, rel, grp, name, bio, share)
            C[k] += 1; P[k].add(p[0]); R[k] += share >= 0.5
            if pairs: pairs.write(f"{lo}\t{k}\t{p[0]}\t{src_gene[s] or '.'}\t{','.join(sorted(under - {src_gene[s]})) or '.'}\t{share:.3f}\n")
            if len(ex_[k]) < 6: ex_[k].append(f"{name[src_gene[s]] if src_gene[s] else s}>{','.join(sorted(name[g] for g in under - {src_gene[s]})[:2]) or '.'}")
        n = sum(v for k, v in C.items() if not k.startswith('source')); n13 = sum(v for k, v in C.items() if k[0] in '123')
        print(f'## score {lo}-{min(hi, 1.0)}: {sum(C.values())} rows; with an annotated source {n}; classes 1-3: {n13} ({n13 / max(1, n):.3f})')
        print('class\trows\tshare of annotated-source rows\tdistinct places\trows that are mostly repeat\texamples (source>place)')
        for k in sorted(C): print(f"{k}\t{C[k]}\t{(C[k] / max(1, n)) if not k.startswith('source') else float('nan'):.3f}\t{len(P[k])}\t{R[k]}\t{' '.join(ex_[k])}")

def selftest():
    iv = merged([(0, 10), (5, 20), (30, 40)]); assert iv == [[0, 20], [30, 40]] and covered([(15, 35)], [a for a, _ in iv], iv) == 10
    rel = collections.defaultdict(set, {'A': {'B', 'P'}, 'B': {'A'}, 'Q': {'P'}, 'P': {'A', 'Q'}}); e = collections.defaultdict(set); nm = collections.defaultdict(str, {'A': 'RPL21', 'Z': 'RPL21P4', 'X': 'FOO'}); bio = collections.defaultdict(str, {'X': 'lncRNA'})
    assert klass('A', {'B'}, rel, e, nm, bio, 0)[0] == '1' and klass('A', {'Q'}, rel, e, nm, bio, 0)[0] == '2' and klass('A', {'Z'}, rel, e, nm, bio, 0)[0] == '3'
    assert klass('A', {'X'}, rel, e, nm, bio, 0.9)[0] == '4' and klass('A', {'X'}, rel, e, nm, bio, 0) == '5 another gene: long non-coding' and klass('A', {'A'}, rel, e, nm, bio, 0)[0] == '6' and klass(None, set(), rel, e, nm, bio, 0).startswith('source')
    print('place_classes selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else: main(*sys.argv[1:])
