#!/usr/bin/env python3
"""Step 0 (not registered, row 1589): would a PROTEIN edge reach the true sibling pairs the families stage misses, without
admitting old paralogs? (user 2026-11-08: "could we extend to using protein sequences")

  prep   GTF GENOME COPIES OUTDIR  -> OUTDIR/prot.faa (longest ORF >= 80 aa of each de novo transcript that maps to an
                                      expressed protein-coding gene; <= 3 distinct ORFs per gene) + OUTDIR/state.pkl
  score  OUTDIR                     -> reads OUTDIR/allvall.tsv (diamond blastp --outfmt 6 qseqid sseqid pident length qlen
                                      slen evalue bitscore) and prints the step-0 table
A gene pair is protein-linked when any ORF pair aligns with e <= 1e-10, both query and subject covered >= 50% (alignment
length / protein length) and identity >= 50% (stricter strata reported). Stop rules (fixed before any protein result):
GO to a registered arm iff >= 100 missed truth pairs are linked AND the precision of the NEW linked pairs (truth / all
coding pairs not predicted now, same-symbol pairs ignored as in pairscore) >= 0.70.
env: as pairscore (DENOVO_CAT_GFF, DENOVO_CAT_TSV, DENOVO_PARALOGS, PAIRSCORE_SAME_SYMBOL=ignore). --selftest checks the ORF finder.
"""
import sys, os, re, collections, itertools, pickle, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product('TCAG', repeat=3),
         'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
COMP = str.maketrans('ACGTN', 'TGCAN')


def longest_orf(seq, min_aa=80):
    """longest ATG..stop ORF on the given strand (an ORF running off the 3' end counts); protein without the stop or None."""
    best = ''
    for fr in range(3):
        prot = ''.join(CODON.get(seq[i:i + 3], 'X') for i in range(fr, len(seq) - 2, 3))
        for chunk in prot.split('*'):
            m = chunk.find('M')
            if m >= 0 and len(chunk) - m > len(best):
                best = chunk[m:]
    return best if len(best) >= min_aa else None


def prep(gtf, genome, copies, out):
    import famtruth as F
    G = F.genes_from_argv([], 'human'); idx = F.exon_index(G)
    expressed = F.expressed_from_gtf(gtf, idx)
    truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
    fams, _, _ = F.load_catalog(copies, G, idx)
    pred = {(x, y) for v in fams.values() for x, y in itertools.combinations(sorted(v), 2)}
    ok = lambda g: G[g]['biotype'] == 'protein_coding' and g in universe and g in expressed
    tx = collections.defaultdict(list); meta = {}
    for line in open(gtf):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
    want = collections.defaultdict(list)
    for t, ex in tx.items():
        ex.sort(); g = F.map_copy(idx, meta[t][0], meta[t][1], ex)
        if g and ok(g):
            want[meta[t][0]].append((t, g, ex, meta[t][1]))
    import pysam
    fa = pysam.FastaFile(genome); orfs = collections.defaultdict(set)
    for chrom, items in want.items():
        for t, g, ex, strand in items:
            s = ''.join(fa.fetch(chrom, a, b) for a, b in ex).upper()
            if strand == '-':
                s = s.translate(COMP)[::-1]
            p = longest_orf(s)
            if p:
                orfs[g].add(p)
    n = 0
    with open(f'{out}/prot.faa', 'w') as w:
        for g, ps in orfs.items():
            for i, p in enumerate(sorted(ps, key=lambda x: (-len(x), x))[:3]):
                w.write(f'>{g}|{i}\n{p}\n'); n += 1
    pickle.dump(dict(names={g: G[g]['name'] for g in G}, truth=truth, pred=pred,
                     coding_genes={g for g in G if G[g]['biotype'] == 'protein_coding' and g in universe}), open(f'{out}/state.pkl', 'wb'))
    print(f'expressed coding genes {sum(1 for g in G if ok(g))}, with an ORF {len(orfs)}, proteins {n}, truth pairs {len(truth)}')


def score(out):
    st = pickle.load(open(f'{out}/state.pkl', 'rb')); names, truth, pred, cg = st['names'], st['truth'], st['pred'], st['coding_genes']
    same = lambda p: os.environ.get('PAIRSCORE_SAME_SYMBOL') == 'ignore' and names[p[0]] == names[p[1]]
    coding_pred = {p for p in pred if p[0] in cg and p[1] in cg and not same(p)}
    missed = {p for p in truth if p not in coding_pred and not same(p)}
    link = {}
    for line in open(f'{out}/allvall.tsv'):
        q, s, pid, ln, ql, sl, ev, bs = line.split('\t')[:8]
        a, b = q.split('|')[0], s.split('|')[0]
        if a == b or float(ev) > 1e-10 or int(ln) < 0.5 * int(ql) or int(ln) < 0.5 * int(sl):
            continue
        p = (min(a, b), max(a, b)); link[p] = max(link.get(p, 0.0), float(pid))
    h = len(coding_pred & truth)
    print(f'truth {len(truth)}  predicted coding {len(coding_pred)} (hit {h})  base precision {h / len(coding_pred):.3f} '
          f'recall {h / len(truth):.3f} (must match pairscore)  missed truth {len(missed)}')
    print('min_identity\tlinked_pairs\tmissed_truth_reached\tnew_pairs\tnew_truth\tnew_precision\tbase_hits_linked\tbase_false_linked')
    for t in (50, 70, 80, 90, 95):
        L = {p for p, v in link.items() if v >= t and not same(p)}
        new = L - coding_pred
        nt = len(new & truth)
        print(f'{t}\t{len(L)}\t{len(L & missed)}\t{len(new)}\t{nt}\t{nt / max(1, len(new)):.3f}\t'
              f'{len(L & coding_pred & truth)}/{len(coding_pred & truth)}\t{len((L & coding_pred) - truth)}/{len(coding_pred - truth)}')
    L = {p for p, v in link.items() if v >= 50 and not same(p)}; new = L - coding_pred; nt = len(new & truth)
    go = len(L & missed) >= 100 and nt / max(1, len(new)) >= 0.70
    print(f'STEP 0 (identity >= 50): missed truth reached {len(L & missed)} (bar 100), new-pair precision '
          f'{nt / max(1, len(new)):.3f} (bar 0.70) -> {"GO" if go else "STOP"}')
    fam = collections.Counter(re.sub(r'[0-9]+[A-Z]?$', '', names[p[0]]) for p in new - truth)
    print('new non-truth pairs, top symbol stems:', ', '.join(f'{k} {v}' for k, v in fam.most_common(10)))
    fam = collections.Counter(re.sub(r'[0-9]+[A-Z]?$', '', names[p[0]]) for p in L & missed)
    print('missed truth reached, top symbol stems:', ', '.join(f'{k} {v}' for k, v in fam.most_common(10)))
    if os.environ.get('PROTEDGE_SD'):  # POST HOC (after the STOP): keep only protein links whose two genes one SD pair holds
        sd = sd_linked({g for p in L for g in p}, os.environ['PROTEDGE_SD'])
        for t in (50, 80, 90):
            Ls = {p for p, v in link.items() if v >= t and not same(p) and p in sd}; new = Ls - coding_pred; nt = len(new & truth)
            print(f'POST HOC protein AND SD-linked, identity >= {t}: missed truth reached {len(Ls & missed)}, new pairs '
                  f'{len(new)}, new truth {nt}, new precision {nt / max(1, len(new)):.3f}')


def sd_linked(genes, bed):
    """gene pairs whose CAT exons one SEDEF pair holds, one gene per copy (>= half exon bases AND span; family.rs rule)."""
    import bisect, famtruth as F
    G = F.load_cat(); by = collections.defaultdict(list)
    for i, l in enumerate(open(bed)):
        if l.startswith('#'): continue
        f = l.split('\t'); by[f[0]].append((int(f[1]), int(f[2]), i, 0)); by[f[3]].append((int(f[4]), int(f[5]), i, 1))
    for v in by.values(): v.sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}; mx = {c: max(x[1] - x[0] for x in v) for c, v in by.items()}
    H = collections.defaultdict(set)
    for g in genes:
        bl = sorted(G[g]['exons']); c = G[g]['chrom']; lo, hi = bl[0][0], bl[-1][1]
        tot = sum(e - s for s, e in bl); v = by.get(c, [])
        for i in range(bisect.bisect_left(st.get(c, []), lo - mx.get(c, 0)), len(v)):
            s, e, p, side = v[i]
            if s >= hi: break
            ins = sum(max(0, min(y, e) - max(x, s)) for x, y in bl)
            if e > lo and tot and 2 * ins >= tot and 2 * max(0, min(hi, e) - max(lo, s)) >= hi - lo: H[g].add((p, side))
    by_pair = collections.defaultdict(set)
    for g, hs in H.items():
        for p, side in hs: by_pair[(p, side)].add(g)
    return {(min(a, b), max(a, b)) for (p, side), gs in by_pair.items() if side == 0 for a in gs for b in by_pair.get((p, 1), ()) if a != b}


def selftest():
    assert longest_orf('CCATG' + 'GCT' * 85 + 'TAACC') == 'M' + 'A' * 85
    assert longest_orf('ATGTAA' * 50) is None
    assert longest_orf('ATG' + 'AAA' * 100) == 'M' + 'K' * 100          # runs off the 3' end
    print('selftest ok')


if __name__ == '__main__':
    cmd = sys.argv[1]
    selftest() if cmd == '--selftest' else prep(*sys.argv[2:6]) if cmd == 'prep' else score(sys.argv[2])
