#!/usr/bin/env python3
"""PREREG r1622 (docs/PREREG_guided_assembly_2026-11-10.md): guided-assembly prototype, arm G1 (known-chain rescue).

  guide_proto.py split REF_GTF PREFIX            -> PREFIX.guide.gtf, PREFIX.guide_decoy.gtf, PREFIX.decoys.tsv
  guide_proto.py support BAM GUIDE_GTF OUT_TSV   -> full-chain primary / secondary records per multi-exon guide chain
  guide_proto.py score REF_GTF PREFIX SUPPORT_TSV NAME=GTF ...   (first NAME = our de novo GTF)
  guide_proto.py tile BAM GUIDE_GTF OUT_TSV      (Amendment 6 step 0: per guide chain, the least-covered adjacent-intron pair)
  guide_proto.py floor REF_GTF PREFIX SUPPORT_TSV OURS_GTF STRINGTIE_G_GTF   (Amendments 3-4: arms G1a, G1b, G1c(K))
Gene = gene_id (transcript_id when empty). Hidden = md5(gene) ends in a hex digit < 5; GUIDE_HIDE=tx hashes the transcript.
"""
import sys, os, re, gzip, hashlib, random, collections
TX = os.environ.get('GUIDE_HIDE') == 'tx'
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
def load(path):  # transcript -> (chrom, strand, gene, sorted exons 1-based inclusive)
    T = {}
    for l in op(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]*)"', f[8]).group(1); g = re.search(r'gene_id "([^"]*)"', f[8])
        T.setdefault(t, (f[0], f[6], (g.group(1) if g and g.group(1) else t), []))[3].append((int(f[3]), int(f[4])))
    for v in T.values(): v[3].sort()
    return T
chain = lambda e: tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))  # introns 0-based half-open
hidden = lambda g: int(hashlib.md5(g.encode()).hexdigest()[-1], 16) < 5
def write(w, t, c, s, g, ex):
    w.write(f'{c}\tguide\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{s}\t.\ttranscript_id "{t}"; gene_id "{g}";\n')
    for a, b in ex: w.write(f'{c}\tguide\texon\t{a}\t{b}\t.\t{s}\t.\ttranscript_id "{t}"; gene_id "{g}";\n')
mode = sys.argv[1]
if mode == 'split':
    T = load(sys.argv[2]); pre = sys.argv[3]; rnd = random.Random(1622)
    introns = {(c, i) for c, s, g, e in T.values() for i in chain(e)}
    n = k = 0
    with open(pre + '.guide.gtf', 'w') as w, open(pre + '.guide_decoy.gtf', 'w') as wd, open(pre + '.decoys.tsv', 'w') as dl:
        for t, (c, s, g, e) in sorted(T.items()):
            if hidden(t if TX else g): continue
            write(w, t, c, s, g, e); write(wd, t, c, s, g, e)
            if len(e) >= 3:
                n += 1
                if n % 20: continue
                for _ in range(20):  # move the genomic right end of one intron 4-10 bp; the new intron must be in no RefSeq transcript
                    j = rnd.randrange(1, len(e)); d = rnd.choice([-1, 1]) * rnd.randint(4, 10); a, b = e[j]
                    if not (e[j - 1][1] + 20 < a + d < b - 10): continue
                    e2 = e[:j] + [(a + d, b)] + e[j + 1:]
                    if (c, chain(e2)[j - 1]) in introns: continue
                    write(wd, 'DECOY_' + t, c, s, g, e2); dl.write(f'DECOY_{t}\t{c}\t{s}\t{j}\t{d}\n'); k += 1; break
    print(f'transcripts {len(T)}; hidden genes {len({v[2] for v in T.values() if hidden(v[2])})} of {len({v[2] for v in T.values()})}; decoys {k}', file=sys.stderr)
elif mode == 'support':
    import pysam
    want = collections.defaultdict(set)
    for c, s, g, e in load(sys.argv[3]).values():
        if len(e) > 1: want[c].add(chain(e))
    prim = collections.Counter(); sec = collections.Counter()
    for a in pysam.AlignmentFile(sys.argv[2]):
        if a.is_unmapped or a.is_supplementary or 'N' not in a.cigarstring: continue
        ch = []; p = a.reference_start
        for o, n in a.cigartuples:
            if o == 3: ch.append((p, p + n))
            if o in (0, 2, 3, 7, 8): p += n
        ch = tuple(ch); c = a.reference_name
        if ch in want[c]: (sec if a.is_secondary else prim)[(c, ch)] += 1
    with open(sys.argv[4], 'w') as w:
        for c, v in want.items():
            for ch in v: w.write(f"{c}\t{','.join(f'{a}-{b}' for a, b in ch)}\t{prim[(c, ch)]}\t{sec[(c, ch)]}\n")
elif mode == 'score':
    ref = load(sys.argv[2]); pre = sys.argv[3]
    sup = {}
    for l in open(sys.argv[4]):
        c, ch, p, s = l.rstrip('\n').split('\t'); sup[(c, tuple(tuple(map(int, x.split('-'))) for x in ch.split(',')))] = int(p)
    cs = lambda T: {(c, s, chain(e)) for c, s, g, e in T.values() if len(e) > 1}
    hid = {(c, s, chain(e)) for t, (c, s, g, e) in ref.items() if len(e) > 1 and hidden(t if TX else g)}
    G = load(pre + '.guide_decoy.gtf'); dec = cs({t: v for t, v in G.items() if t.startswith('DECOY_')}); gui = cs({t: v for t, v in G.items() if not t.startswith('DECOY_')})
    hid -= gui  # a hidden gene's chain that a guide gene also carries is not novel
    supported = {x for x in gui if sup.get((x[0], x[2]), 0) >= 1}
    print(f'guide chains {len(gui)} (with >= 1 full-chain primary alignment {len(supported)}); hidden chains {len(hid)}; decoys {len(dec)} (with a full-chain primary alignment {sum(sup.get((c, ch), 0) >= 1 for c, s, ch in dec)})')
    arms = [a.split('=', 1) for a in sys.argv[5:]]; ours = cs(load(arms[0][1]))
    g1 = ours | supported
    rows = [(n, cs(load(p))) for n, p in arms]; rows.insert(1, ('ours G1 (prototype)', g1))
    print('arm\tchains\tguide chains found\tof supported\thidden found\tnon-guide chains\tnon-guide precision\tdecoys emitted')
    for n, S in rows:
        ng = S - gui - dec
        print(n, len(S), len(S & gui), f'{len(S & supported) / max(1, len(supported)):.3f}', len(S & hid), len(ng), f'{len(ng & hid) / max(1, len(ng)):.3f}', len(S & dec), sep='\t')
    print(f'G1: de novo chains missing from G1 {len(ours - g1)}; rescued {len(g1 - ours)} by primary alignments: ' + str(sorted(collections.Counter(min(sup.get((c, ch), 0), 5) for c, s, ch in g1 - ours).items())))
elif mode == 'floor':
    ref = load(sys.argv[2]); cs = lambda T: {(c, s, chain(e)) for c, s, g, e in T.values() if len(e) > 1}
    gui = cs(load(sys.argv[3] + '.guide.gtf'))
    hid = {(c, s, chain(e)) for t, (c, s, g, e) in ref.items() if len(e) > 1 and hidden(t if TX else g)} - gui
    sup = set()
    for l in open(sys.argv[4]):
        c, ch, p, s = l.rstrip('\n').split('\t')
        if int(p): sup.add((c, tuple(tuple(map(int, x.split('-'))) for x in ch.split(','))))
    reads = {}
    for l in op(sys.argv[5]):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'transcript': reads[re.search(r'transcript_id "([^"]*)"', f[8]).group(1)] = int(re.search(r'reads "(\d+)"', f[8]).group(1))
    R = collections.Counter()
    for t, (c, s, g, e) in load(sys.argv[5]).items():
        if len(e) > 1: R[(c, s, chain(e))] = max(R[(c, s, chain(e))], reads.get(t, 0))
    ours = set(R); by = collections.defaultdict(list)
    for c, s, ch in gui:
        for i in ch: by[(c, s, i)].append(ch)
    def explained(k):  # 'frag' = contiguous run of a guide chain, 'ret' = ordered subset of one, else None
        c, s, ch = k; out = None
        if any((c, s, i) not in by for i in ch): return None
        for g in by[(c, s, ch[0])]:
            j = g.index(ch[0])
            if g[j:j + len(ch)] == ch: return 'frag'
            it = iter(g[j:])
            if all(i in it for i in ch): out = 'ret'
        return out
    ng = ours - gui; ex = {k: explained(k) for k in ng}
    cls2 = {k for k in ng if any((k[0], k[1], i) in by for i in k[2]) and any((k[0], k[1], i) not in by for i in k[2])}
    g1 = ours | {k for k in gui if (k[0], k[2]) in sup}
    arms = [('ours de novo', ours), ('G1', g1), ('G1a', g1 - {k for k in ng if ex[k] == 'frag'}), ('G1b', g1 - {k for k in ng if ex[k]})]
    arms += [(f'G1c({K})', arms[3][1] - {k for k in cls2 if R[k] < K}) for K in (3, 5, 10, 20)]
    arms.append(('StringTie -G', cs(load(sys.argv[6]))))
    print(f'split {"transcript" if TX else "gene"}: guide chains {len(gui)}, hidden chains {len(hid)}; ours non-guide {len(ng)}: fragments {sum(v == "frag" for v in ex.values())} (hidden {sum(ex[k] == "frag" for k in ng & hid)}), retained {sum(v == "ret" for v in ex.values())} (hidden {sum(ex[k] == "ret" for k in ng & hid)}), new-intron class {len(cls2)} (hidden {len(cls2 & hid)})')
    print('arm\tchains\tguide found\thidden found\tnon-guide\tnon-guide precision')
    for n, S in arms: print(n, len(S), len(S & gui), len(S & hid), len(S - gui), f'{len((S - gui) & hid) / max(1, len(S - gui)):.3f}', sep='\t')
    # Amendment 5: G1d = G1c(10) with every non-guide chain at reads >= 3; curves for bars R1 / R2
    cov = collections.Counter(); T = load(sys.argv[6]); cv = {}
    for l in op(sys.argv[6]):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'transcript': cv[re.search(r'transcript_id "([^"]*)"', f[8]).group(1)] = float(re.search(r'cov "([^"]*)"', f[8]).group(1))
    for t, (c, s, g, e) in T.items():
        if len(e) > 1: cov[(c, s, chain(e))] = max(cov[(c, s, chain(e))], cv.get(t, 0))
    print('curve\tnon-guide\thidden found\tnon-guide precision')
    if os.environ.get('GUIDE_DUMP'):  # parity file for `utilities guide`: the G1d chain set
        with open(os.environ['GUIDE_DUMP'], 'w') as w:
            for c, st, ch in sorted({k for k in arms[6][1] if k in gui or R[k] >= 3}): w.write(f"{c}\t{st}\t{','.join(f'{a}-{b}' for a, b in ch)}\n")
    for t in (3, 5, 10, 20, 50):
        S = {k for k in arms[6][1] - gui if R[k] >= t}; print(f'G1d reads >= {t}', len(S), len(S & hid), f'{len(S & hid) / max(1, len(S)):.3f}', sep='\t')
    for t in (0, 1, 2, 3, 5, 10):
        S = {k for k in arms[-1][1] - gui if cov[k] >= t}; print(f'StringTie -G cov >= {t}', len(S), len(S & hid), f'{len(S & hid) / max(1, len(S)):.3f}', sep='\t')
elif mode == 'tile':
    import pysam
    def fits(g, j, ch, st, en):  # read chain ch = g[j:j+len(ch)], and the read's ends do not run into the path's neighbouring introns
        return g[j:j + len(ch)] == ch and (j == 0 or st >= g[j - 1][1]) and (j + len(ch) == len(g) or en <= g[j + len(ch)][0])
    G = ((10, 20), (30, 40), (50, 60))
    assert fits(G, 1, ((30, 40), (50, 60)), 22, 70) and not fits(G, 1, ((30, 40), (50, 60)), 15, 70) and not fits(G, 0, ((10, 20), (50, 60)), 0, 70)
    chains = sorted({(c, ch) for c, ch in ((c, chain(e)) for c, s, g, e in load(sys.argv[3]).values()) if ch})
    by = collections.defaultdict(list)
    for n, (c, g) in enumerate(chains):
        for j, i in enumerate(g): by[(c, i)].append((n, j))
    cnt = [[0] * max(1, len(g) - 1) for c, g in chains]
    for a in pysam.AlignmentFile(sys.argv[2]):
        if a.is_unmapped or a.is_supplementary or a.is_secondary or 'N' not in a.cigarstring: continue
        ch = []; p = a.reference_start
        for o, n in a.cigartuples:
            if o == 3: ch.append((p, p + n))
            if o in (0, 2, 3, 7, 8): p += n
        ch = tuple(ch)
        for n, j in by.get((a.reference_name, ch[0]), ()):
            g = chains[n][1]
            if not fits(g, j, ch, a.reference_start, p): continue
            if len(g) == 1: cnt[n][0] += 1
            for k in range(j, j + len(ch) - 1): cnt[n][k] += 1
    with open(sys.argv[4], 'w') as w:
        for (c, g), v in zip(chains, cnt): w.write(f"{c}\t{','.join(f'{a}-{b}' for a, b in g)}\t{min(v)}\t{','.join(map(str, v))}\n")
