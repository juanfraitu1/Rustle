#!/usr/bin/env python3
"""Row 1596 (docs/PREREG_paralog_kmer_expression_2026-11-09.md): paralog-specific k-mer expression on RNA reads.

  rnakmer.py run COPIES.tsv GENOME.fa BAM OUT_PREFIX [--threads N] [--s1e S1E.tsv --cat CAT.tsv]
  rnakmer.py arms COPIES.tsv GENOME.fa BAM OUT_PREFIX [--threads N] [--exon-placement]   row 1598: k {21,25,30} x m {1,2,3},
             chromosome split; --exon-placement = row 1599's unique placement (>= 100 bp of primary aligned bases in copy exons)
  rnakmer.py --selftest
SPEC(copy) = canonical k=30 k-mers of the copy's spliced sequence found in no other copy (any family) and at most once in
the genome (every contig). Read -> copy c of family F iff >= 3 distinct SPEC hits on c and 0 on F's other copies.
Writes OUT.per_copy.tsv, OUT.per_family.tsv, OUT.summary.tsv; verdict on stdout.
"""
import sys, os, csv, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parcn_assembly as P

MINHIT, RESOLVED, FLANK = 3, 100, 1000
G = {}                                    # fork-shared state for workers


def exons(s):
    return sorted(tuple(map(int, x.split('-'))) for x in s.split(',') if x)


def load_copies(path):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    fams = sorted({r['family_id'] for r in rows})
    fi = {f: i for i, f in enumerate(fams)}
    for i, r in enumerate(rows):
        r['ci'], r['fi'], r['ex'] = i, fi[r['family_id']], exons(r['exons'])
        r['s'], r['e'] = int(r['locus_start']), int(r['locus_end'])
    return rows, fams


def spec_candidates(seqs):
    """seqs: list of spliced sequences -> (codes sorted, copy index) of k-mers present in exactly one copy."""
    cs, ix = [], []
    for i, s in enumerate(seqs):
        u = np.unique(P.canonical_kmers(s)); u = u[u != 0]
        cs.append(u); ix.append(np.full(len(u), i, dtype=np.int32))
    c, x = np.concatenate(cs), np.concatenate(ix)
    o = np.argsort(c, kind='stable'); c, x = c[o], x[o]
    u, start, n = np.unique(c, return_index=True, return_counts=True)
    keep = n == 1
    return u[keep], x[start[keep]]


def _count_contig(name):
    import pysam
    seq = pysam.FastaFile(G['genome']).fetch(name).encode()
    t = G['targets']; counts = np.zeros(len(t), dtype=np.int64); ch = 1 << 24
    for s in range(0, max(1, len(seq) - P.K + 1), ch):
        codes = P.canonical_kmers(seq[s:s + ch + P.K - 1]); codes = codes[codes != 0]
        pos = np.searchsorted(t, codes); pos[pos >= len(t)] = 0
        hit = t[pos] == codes
        if hit.any():
            counts += np.bincount(pos[hit], minlength=len(t))
    return counts


def hits_of(seq):
    """{copy index: distinct SPEC hits} of a read sequence."""
    k = np.unique(P.canonical_kmers(seq))
    sc = G['spec']
    if len(k) == 0 or len(sc) == 0:
        return {}
    pos = np.searchsorted(sc, k); pos[pos >= len(sc)] = 0
    h = G['spec_copy'][pos[sc[pos] == k]]
    if len(h) == 0:
        return {}
    u, n = np.unique(h, return_counts=True)
    return dict(zip(u.tolist(), n.tolist()))


def _region(ci):
    """reads with a non-supplementary record overlapping copy ci's span +- FLANK."""
    import pysam
    if 'bam' not in G:
        G['bam'] = pysam.AlignmentFile(G['bampath'])
    r = G['copies'][ci]; out = {}
    for a in G['bam'].fetch(r['chrom'], max(0, r['s'] - FLANK), r['e'] + FLANK):
        if a.is_supplementary or a.is_unmapped:
            continue
        q = a.query_name
        d = out.get(q)
        prim = not a.is_secondary
        if (d is None or (prim and not d[2])) and a.query_sequence is not None:   # primary SEQ preferred (secondaries may be clipped)
            d = out[q] = [G.get('hitfn', hits_of)(a.query_sequence), d[1] if d else None, prim]
        if d is None:
            continue
        if not a.is_secondary and a.reference_start < r['e'] and a.reference_end > r['s']:
            ebp = sum(max(0, min(y, e) - max(x, b)) for x, y in a.get_blocks() for b, e in r['ex'])
            d[1] = (r['ci'], a.mapping_quality, ebp)
    return ci, out


def assign(hits, fam_copies):
    """copy of the family the read is assigned to, or None."""
    h = [(c, hits.get(c, 0)) for c in fam_copies]
    on = [c for c, n in h if n > 0]
    return on[0] if len(on) == 1 and hits[on[0]] >= MINHIT else None


def spearman(x, y):
    from scipy.stats import spearmanr
    return spearmanr(x, y).correlation


def run(copies_tsv, genome, bam, out, threads, s1e=None, cat=None):
    import pysam, multiprocessing as mp, time
    t0 = time.time()
    rows, fams = load_copies(copies_tsv)
    fa = pysam.FastaFile(genome)
    seqs = [''.join(fa.fetch(r['chrom'], a, b) for a, b in r['ex']) for r in rows]
    cand, cand_copy = spec_candidates(seqs)
    print(f'copies {len(rows)}, families {len(fams)}, single-copy k-mers {len(cand)} ({time.time() - t0:.0f} s)', flush=True)
    G.update(genome=genome, targets=cand)
    gc = np.zeros(len(cand), dtype=np.int64)
    contigs = sorted(fa.references, key=lambda c: -fa.get_reference_length(c))
    with mp.Pool(threads) as pool:
        for c in pool.imap_unordered(_count_contig, contigs):
            gc += c
    keep = gc <= 1
    G.update(spec=cand[keep], spec_copy=cand_copy[keep], copies=rows, bampath=bam)
    nspec = np.bincount(cand_copy[keep], minlength=len(rows))
    resolved = nspec >= RESOLVED
    print(f'SPEC k-mers {keep.sum()} (genome count <= 1 of {len(cand)}); resolved copies {resolved.sum()} '
          f'({time.time() - t0:.0f} s)', flush=True)
    del G['targets']
    reads = {}                             # name -> [hits, primary copies [(ci, mapq)]], fetched families
    with mp.Pool(threads) as pool:
        for ci, part in pool.imap_unordered(_region, range(len(rows)), chunksize=4):
            f = rows[ci]['fi']
            for q, (h, pr, prim) in part.items():
                d = reads.get(q)
                if d is None:
                    d = reads[q] = [h, set(), set(), prim]
                elif prim and not d[3]:
                    d[0], d[3] = h, True
                if pr:
                    d[1].add(pr)
                d[2].add(f)
    print(f'reads in scope {len(reads)} ({time.time() - t0:.0f} s)', flush=True)
    fam_copies = collections.defaultdict(list)
    for r in rows:
        fam_copies[r['fi']].append(r['ci'])
    kmer_n, uniq_n, tied_asg = collections.Counter(), collections.Counter(), collections.Counter()
    bucket = collections.Counter(); bucket_asg = collections.Counter()
    agree = reach_n = reach_d = 0; tied_n = tied_a = tied_same = 0
    for h, prim, fetched, _ in reads.values():
        fs = fetched | {rows[c]['fi'] for c in h}
        for f in fs:
            a = assign(h, fam_copies[f])
            if a is not None:
                kmer_n[a] += 1
            pc = [(c, m) for c, m, _ in prim if rows[c]['fi'] == f]
            if not pc:
                b = 'off'
            elif len(pc) == 1 and pc[0][1] == 60:
                b = 'unique'; t = pc[0][0]; uniq_n[t] += 1
                if resolved[t]:
                    reach_d += 1
                    if a is not None:
                        reach_n += 1; agree += a == t
            else:
                b = 'tied'; tied_n += 1
                if a is not None:
                    tied_a += 1; tied_same += a in {c for c, _ in pc}; tied_asg[a] += 1
            bucket[b] += 1; bucket_asg[b] += a is not None
    with open(out + '.per_copy.tsv', 'w') as w:
        w.write('family\ttid\tchrom\tstart\tend\tn_spec\tresolved\tkmer_reads\tunique_reads\ttied_assigned\n')
        for r in rows:
            c = r['ci']
            w.write(f"{r['family_id']}\t{r['tid']}\t{r['chrom']}\t{r['s']}\t{r['e']}\t{nspec[c]}\t{int(resolved[c])}\t"
                    f"{kmer_n[c]}\t{uniq_n[c]}\t{tied_asg[c]}\n")
    s1e_of = collections.defaultdict(list)
    if s1e and cat:
        import rnacn
        gidx = rnacn.by_chrom(rnacn.load_genes(s1e, cat))
        for r in rows:
            for g in rnacn.copy_hits_genes(r['chrom'], r['ex'], gidx):
                s1e_of[r['fi']].append(f"{g['gene']}:{g['parcn']}")
    with open(out + '.per_family.tsv', 'w') as w:
        w.write('family\tcopies\tresolved\tkmer_reads\tunresolved_share\tN_eff\ts1e_genes_parCN\n')
        for f, cs in fam_copies.items():
            tot = sum(kmer_n[c] for c in cs)
            neff = 1 / sum((kmer_n[c] / tot) ** 2 for c in cs) if tot >= 20 else float('nan')
            un = sum(kmer_n[c] for c in cs if not resolved[c]) / tot if tot else float('nan')
            w.write(f"{fams[f]}\t{len(cs)}\t{int(sum(resolved[c] for c in cs))}\t{tot}\t{un:.3f}\t{neff:.3f}\t"
                    f"{','.join(sorted(set(s1e_of[f]))) or '.'}\n")
    rc = [r['ci'] for r in rows if resolved[r['ci']] and uniq_n[r['ci']] >= 10]
    rho = spearman([kmer_n[c] for c in rc], [uniq_n[c] for c in rc]) if len(rc) > 2 else float('nan')
    ag, re_ = agree / reach_n if reach_n else float('nan'), reach_n / reach_d if reach_d else float('nan')
    L = [f'copies\t{len(rows)}', f'families\t{len(fams)}', f'spec_kmers\t{int(keep.sum())}',
         f'resolved_copies\t{int(resolved.sum())}', f'reads_in_scope\t{len(reads)}']
    L += [f'bucket_{b}\t{bucket[b]}\tassigned {bucket_asg[b]}' for b in ('unique', 'tied', 'off')]
    L += [f'agreement\t{agree}/{reach_n}\t{ag:.4f}', f'reach\t{reach_n}/{reach_d}\t{re_:.4f}',
          f'spearman\t{rho:.4f}\tcopies {len(rc)}',
          f'tied_assigned\t{tied_a}/{tied_n}', f'tied_assigned_to_primary_copy\t{tied_same}/{tied_a}']
    if resolved.sum() < 200 or reach_n < 1000:
        v = f'STOP (resolved copies {int(resolved.sum())} < 200 or assigned unique reads {reach_n} < 1000)'
    else:
        v = ('PASS' if ag >= 0.98 and re_ >= 0.70 and rho >= 0.90 else 'FAIL') + \
            f' (agreement {ag:.4f} >= 0.98, reach {re_:.4f} >= 0.70, spearman {rho:.4f} >= 0.90)'
    L.append(f'verdict\t{v}')
    open(out + '.summary.tsv', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L) + f'\n({time.time() - t0:.0f} s)')


KS, MS = (21, 25, 30), (1, 2, 3)


def multi_hits(seq):
    out = {}
    for k in KS:
        P.K = k; G['spec'], G['spec_copy'] = G['specs'][k]
        out[k] = hits_of(seq)
    return out


def half_of(chrom):
    """row 1598 split: 'dev' = odd autosomes + chrX, 'held' = even autosomes, chrY, anything else."""
    c = chrom[3:] if chrom.startswith('chr') else chrom
    return 'dev' if c == 'X' or (c.isdigit() and int(c) % 2 == 1) else 'held'


def assign_m(hits, fam_copies, m):
    on = [c for c in fam_copies if hits.get(c, 0) > 0]
    return on[0] if len(on) == 1 and hits[on[0]] >= m else None


def arms(copies_tsv, genome, bam, out, threads, exon_placement=False):
    import pysam, multiprocessing as mp, time
    t0 = time.time()
    rows, fams = load_copies(copies_tsv)
    fa = pysam.FastaFile(genome)
    seqs = [''.join(fa.fetch(r['chrom'], a, b) for a, b in r['ex']) for r in rows]
    contigs = sorted(fa.references, key=lambda c: -fa.get_reference_length(c))
    specs, resolved = {}, {}
    for k in KS:
        P.K = k
        cand, cand_copy = spec_candidates(seqs)
        G.update(genome=genome, targets=cand)
        gc = np.zeros(len(cand), dtype=np.int64)
        with mp.Pool(threads) as pool:
            for c in pool.imap_unordered(_count_contig, contigs):
                gc += c
        keep = gc <= 1
        specs[k] = (cand[keep], cand_copy[keep])
        resolved[k] = np.bincount(cand_copy[keep], minlength=len(rows)) >= RESOLVED
        print(f'k {k}: SPEC {keep.sum()} of {len(cand)}, resolved copies {resolved[k].sum()} ({time.time() - t0:.0f} s)', flush=True)
    del G['targets']
    G.update(specs=specs, copies=rows, bampath=bam, hitfn=multi_hits, exon_placement=exon_placement)
    reads = {}
    with mp.Pool(threads) as pool:
        for ci, part in pool.imap_unordered(_region, range(len(rows)), chunksize=4):
            f = rows[ci]['fi']
            for q, (h, pr, prim) in part.items():
                d = reads.get(q)
                if d is None:
                    d = reads[q] = [h, set(), set(), prim]
                elif prim and not d[3]:
                    d[0], d[3] = h, True
                if pr:
                    d[1].add(pr)
                d[2].add(f)
    print(f'reads in scope {len(reads)} ({time.time() - t0:.0f} s)', flush=True)
    fam_copies = collections.defaultdict(list)
    for r in rows:
        fam_copies[r['fi']].append(r['ci'])
    half = {f: half_of(rows[cs[0]]['chrom']) for f, cs in fam_copies.items()}   # copies sorted by copy_idx in the table
    S = collections.defaultdict(collections.Counter)                           # (half, k, m) -> counters
    kmer_n = collections.defaultdict(collections.Counter); uniq_n = collections.Counter(); bk = collections.Counter()
    for H, prim, fetched, _ in reads.values():
        fs = set(fetched)
        for k in KS:
            fs |= {rows[c]['fi'] for c in H[k]}
        for f in fs:
            if G.get('exon_placement'):                                  # row 1599: exon overlap >= 100 bp
                on = [x for x in prim if rows[x[0]]['fi'] == f]; ex = [x for x in on if x[2] >= 100]
                if len(ex) == 1 and ex[0][1] == 60 and all(x[2] == 0 for x in on if x is not ex[0]):
                    t, b = ex[0][0], 'unique'
                else:
                    t = None
                    b = ('tied' if ex and (ex[0][1] < 60 or len(ex) >= 2) else 'span_only' if on and not ex
                         else 'off' if not on else 'other')
                bk[(half[f], b)] += 1
            else:
                pc = [(c, m) for c, m, _ in prim if rows[c]['fi'] == f]
                t = pc[0][0] if len(pc) == 1 and pc[0][1] == 60 else None
            if t is not None:
                uniq_n[t] += 1
            for k in KS:
                h = H[k]
                for m in MS:
                    a = assign_m(h, fam_copies[f], m); key = (half[f], k, m); c = S[key]
                    if a is not None:
                        kmer_n[(k, m)][a] += 1
                    if t is not None and resolved[k][t]:
                        c['reach_d'] += 1
                        if a is not None:
                            c['reach_n'] += 1; c['agree'] += a == t
                        else:
                            on = [x for x in fam_copies[f] if h.get(x, 0) > 0]
                            c['none' if not on else 'weak' if on == [t] else 'conflict'] += 1
    from scipy.stats import spearmanr
    res = {}
    with open(out + '.arms.tsv', 'w') as w:
        w.write('half\tk\tm\tresolved\tagreement\treach\tspearman\tn_assigned_unique\tunassigned_none\tweak\tconflict\n')
        for hf in ('dev', 'held'):
            for k in KS:
                for m in MS:
                    c = S[(hf, k, m)]
                    rc = [r['ci'] for r in rows if half[r['fi']] == hf and resolved[k][r['ci']] and uniq_n[r['ci']] >= 10]
                    rho = spearmanr([kmer_n[(k, m)][x] for x in rc], [uniq_n[x] for x in rc]).correlation if len(rc) > 2 else float('nan')
                    nres = sum(1 for r in rows if half[r['fi']] == hf and resolved[k][r['ci']])
                    ag = c['agree'] / c['reach_n'] if c['reach_n'] else float('nan')
                    re_ = c['reach_n'] / c['reach_d'] if c['reach_d'] else float('nan')
                    res[(hf, k, m)] = (ag, re_, rho, nres, c['reach_n'])
                    w.write(f"{hf}\t{k}\t{m}\t{nres}\t{ag:.4f}\t{re_:.4f}\t{rho:.4f}\t{c['reach_n']}\t{c['none']}\t{c['weak']}\t{c['conflict']}\n")
    print(open(out + '.arms.tsv').read())
    if bk:
        print('buckets (read, family) pairs: ' + ', '.join(f'{h}/{b} {n}' for (h, b), n in sorted(bk.items())))
    ok = [(res[('dev', k, m)][1], k, m) for k in KS for m in MS if res[('dev', k, m)][0] >= 0.98]
    if not ok:
        print('VERDICT: FAIL (no dev arm with agreement >= 0.98)'); return
    _, k, m = max(ok)
    ag, re_, rho, nres, n = res[('held', k, m)]
    if nres < 200 or n < 1000:
        v = f'STOP (held-out resolved {nres}, assigned unique {n})'
    else:
        v = ('PASS' if ag >= 0.98 and re_ >= 0.70 and rho >= 0.90 else 'FAIL') + f' (agreement {ag:.4f}, reach {re_:.4f}, spearman {rho:.4f})'
    print(f'dev choice k={k} m={m}; VERDICT held-out: {v} ({time.time() - t0:.0f} s)')


def selftest():
    rng = np.random.default_rng(3)
    base = ''.join(rng.choice(list('ACGT'), 400))
    a = base; b = base[:200] + ('A' if base[200] != 'A' else 'C') + base[201:]      # one SNP at 200
    c = ''.join(rng.choice(list('ACGT'), 400))
    codes, cp = spec_candidates([a, b, c])
    # a and b differ by one SNP: 30 k-mers each specific; c fully specific
    assert np.bincount(cp, minlength=3).tolist() == [30, 30, 371], np.bincount(cp, minlength=3)
    G.update(spec=codes, spec_copy=cp)
    h = hits_of(b[150:260])
    assert h == {1: 30} and assign(h, [0, 1]) == 1 and assign({0: 2}, [0, 1]) is None and assign({0: 5, 1: 1}, [0, 1]) is None
    rc = b[150:260][::-1].translate(str.maketrans('ACGT', 'TGCA'))
    assert hits_of(rc) == {1: 30}                               # canonical: strand does not matter
    assert half_of('chr1') == 'dev' and half_of('chrX') == 'dev' and half_of('chr2') == 'held' and half_of('chrY') == 'held'
    assert half_of('chrM') == 'held' and assign_m({0: 1}, [0, 1], 1) == 0 and assign_m({0: 1}, [0, 1], 2) is None
    print('selftest ok')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == '--selftest':
        selftest()
    elif a[0] == 'arms':
        arms(a[1], a[2], a[3], a[4], int(a[a.index('--threads') + 1]) if '--threads' in a else 8, '--exon-placement' in a)
    else:
        o = lambda k, d=None: a[a.index(k) + 1] if k in a else d
        run(a[1], a[2], a[3], a[4], int(o('--threads', 8)), o('--s1e'), o('--cat'))
