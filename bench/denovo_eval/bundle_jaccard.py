#!/usr/bin/env python3
"""Row 1762 (docs/PREREG_bundle_jaccard_2026-11-16.md): read-set Jaccard between de novo bundles (loci) as the family rule.
  samtools view -F 0x804 BAM CHROM | bundle_jaccard.py hits LOCI.gff3 CHROM > CHROM.hits   (molecule, record, loci with >= 20 bp)
  sort -k1,1 *.hits | bundle_jaccard.py pairs OUT          -> OUT.n.tsv (locus, molecules), OUT.pairs.tsv (A, B, shared)
  bundle_jaccard.py arms LOCI.gff3 OUT [SEDEF.bed]         -> OUT.J<t>.s<floor>.copies.tsv (pairscore input) + OUT.arms.tsv
  bundle_jaccard.py mcl LOCI.gff3 OUT                      -> OUT.M<t>.I<inflation>.copies.tsv + OUT.mclarms.tsv (Amendment 3)
  bundle_jaccard.py mcl LOCI.gff3 OUT SEDEF.bed            -> OUT.{SJ,SM,DM}<t>.copies.tsv + OUT.sdarms.tsv (Amendment 4: SD gate / damp)
"""
import sys, re, bisect, collections, itertools, functools

def loci(gff):
    """name -> [chrom, strand, start, end, exons (0-based half-open)], in file order."""
    L = {}
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'gene': L[re.search(r'Name=([^;]+)', f[8]).group(1)] = [f[0], f[6], int(f[3]) - 1, int(f[4]), []]
        elif f[2] == 'exon': L[re.search(r'gene=([^;]+)', f[8]).group(1)][4].append((int(f[3]) - 1, int(f[4])))
    return L

def blocks(pos, cigar):
    out, s, e = [], pos, pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(n)
        if op in 'M=XD': e += n
        elif op == 'N':
            out.append((s, e)); s = e = e + n
    out.append((s, e))
    return out
assert blocks(10, '5S10=100N5X2D3=4I1=') == [(10, 20), (120, 131)]

a = sys.argv[1:]
if a[0] == 'hits':
    L = loci(a[1]); names = list(L); NODES = 'nodes' in a  # `nodes`: ids are bundlenodes (exons in file order), not loci
    STRANDED = 'stranded' in a  # `stranded`: a record counts only on bundles of its own strand (reads mapped -uf: alignment strand = transcript strand)
    ex, k = [], 0
    for i, nm in enumerate(names):
        for s, e in L[nm][4]:
            if L[nm][0] == a[2]: ex.append((s, e, k if NODES else i, L[nm][1]))
            k += 1
    ex.sort()
    st = [x[0] for x in ex]; m = max((x[1] - x[0] for x in ex), default=0)
    for line in sys.stdin:
        f = line.split('\t', 6)
        pos = int(f[3]) - 1; share = collections.Counter(); rs = '-' if int(f[1]) & 16 else '+'
        for s, e in blocks(pos, f[5]):
            i = bisect.bisect_left(st, s - m)
            while i < len(ex) and ex[i][0] < e:
                if ex[i][1] > s and (not STRANDED or ex[i][3] == rs): share[ex[i][2]] += min(ex[i][1], e) - max(ex[i][0], s)
                i += 1
        hit = [str(i) for i, v in share.items() if v >= 20]
        if hit: print(f[0], f"{f[2]}:{f[3]}:{f[1]}", ','.join(hit), sep='\t')
elif a[0] == 'both':  # LOCI.gff3 OUT < sorted NODE hits of one bucket of molecules -> OUT.n.tsv, OUT.pairs.tsv, OUT.nodepairs.tsv (all counts, unfiltered)
    L = loci(a[1]); loc = [i for i, nm in enumerate(L) for _ in L[nm][4]]; out = a[2]
    n = collections.Counter(); sh = collections.Counter(); nsh = collections.Counter()
    def flush(recs):
        for l in {loc[x] for r in recs for x in r}: n[l] += 1
        recs = list(set(recs)); seen = set(); lseen = set()  # identical records add no pair
        for r, q in itertools.combinations(recs, 2):
            rl, ql = {loc[x] for x in r}, {loc[x] for x in q}
            for A in rl - ql:
                for B in ql - rl: lseen.add((A, B) if A < B else (B, A))
            for x in r - q:
                for y in q - r: seen.add((x, y) if x < y else (y, x))
        for k in seen: nsh[k] += 1
        for k in lseen: sh[k] += 1
    cur, recs = None, []
    for line in sys.stdin:
        q, _, ls = line.rstrip('\n').split('\t')
        if q != cur:
            if recs: flush(recs)
            cur, recs = q, []
        recs.append(frozenset(map(int, ls.split(','))))
    if recs: flush(recs)
    for suf, d in (('n', n), ('pairs', sh), ('nodepairs', nsh)):
        with open(f"{out}.{suf}.tsv", 'w') as w:
            for k, c in d.items(): w.write(('\t'.join(map(str, k)) if isinstance(k, tuple) else str(k)) + f"\t{c}\n")
elif a[0] == 'nodepairs':  # stdin sorted node hits -> OUT.nodepairs.tsv (x, y, molecules): x and y on two different records
    sh = collections.Counter()
    def flush(recs):
        seen = set()
        for r, q in itertools.combinations(recs, 2):
            for x in r - q:
                for y in q - r: seen.add((x, y) if x < y else (y, x))
        for k in seen: sh[k] += 1
    cur, recs = None, []
    for line in sys.stdin:
        q, _, ls = line.rstrip('\n').split('\t')
        if q != cur:
            if len(recs) > 1: flush(recs)
            cur, recs = q, []
        recs.append(frozenset(map(int, ls.split(','))))
    if len(recs) > 1: flush(recs)
    with open(a[1] + '.nodepairs.tsv', 'w') as w:
        for (x, y), c in sh.items():
            if c >= 2: w.write(f"{x}\t{y}\t{c}\n")
elif a[0] == 'nodearms':  # LOCI.gff3 OUT -> OUT.N<T>.copies.tsv
    L = loci(a[1]); names = list(L); out = a[2]
    node = [(i, e - s) for i, nm in enumerate(names) for s, e in L[nm][4]]
    link = collections.defaultdict(lambda: (set(), set()))
    for line in open(out + '.nodepairs.tsv'):
        x, y, _ = map(int, line.split()); A, B = node[x][0], node[y][0]
        if A == B: continue
        la, lb = L[names[A]], L[names[B]]
        if la[0] == lb[0] and la[2] < lb[3] and lb[2] < la[3]: continue
        if A > B: A, B, x, y = B, A, y, x
        link[A, B][0].add(x); link[A, B][1].add(y)
    tot = collections.Counter()
    for i, bp in node: tot[i] += bp
    E = []
    for (A, B), (xs, ys) in link.items():
        small, ns = (A, xs) if tot[A] <= tot[B] else (B, ys)
        E.append((A, B, sum(node[x][1] for x in ns) / tot[small]))
    with open(out + '.nodearms.tsv', 'w') as rep:
        rep.write("arm\tedges\tfamilies\tloci_in_families\tlargest\n")
        for T in (0.3, 0.5, 0.7):
            p = {}
            def find(x):
                while p.setdefault(x, x) != x: p[x] = p[p[x]]; x = p[x]
                return x
            ed = [(A, B) for A, B, f in E if f >= T]
            for A, B in ed: p[find(A)] = find(B)
            fam = collections.defaultdict(list)
            for x in p: fam[find(x)].append(x)
            with open(f"{out}.N{T}.copies.tsv", 'w') as w:
                w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
                for k, (r, mem) in enumerate(sorted(fam.items())):
                    for ci, x in enumerate(sorted(mem)):
                        c, sdn, s, e, exs = L[names[x]]
                        w.write(f"BN{k}\t{ci}\t{names[x]}\t{c}\t{s}\t{e}\t{len(exs)}\t{sdn}\t0\t{','.join(f'{u}-{v}' for u, v in sorted(exs))}\n")
            rep.write(f"N{T}\t{len(ed)}\t{len(fam)}\t{len(p)}\t{max(map(len, fam.values()), default=0)}\n")
elif a[0] == 'pairs':
    n = collections.Counter(); sh = collections.Counter()
    def flush(recs):
        allL = set().union(*recs)
        for l in allL: n[l] += 1
        if len(recs) > 1 and len(allL) > 1:
            for A, B in itertools.combinations(sorted(allL), 2):  # one record on A not B, another on B not A
                if any(A in r and B not in r for r in recs) and any(B in r and A not in r for r in recs): sh[A, B] += 1
    cur, recs = None, []
    for line in sys.stdin:
        q, _, ls = line.rstrip('\n').split('\t')
        if q != cur:
            if recs: flush(recs)
            cur, recs = q, []
        recs.append(frozenset(map(int, ls.split(','))))
    if recs: flush(recs)
    with open(a[1] + '.n.tsv', 'w') as w:
        for l, c in sorted(n.items()): w.write(f"{l}\t{c}\n")
    with open(a[1] + '.pairs.tsv', 'w') as w:
        for (A, B), c in sorted(sh.items()): w.write(f"{A}\t{B}\t{c}\n")
elif a[0] == 'arms':
    L = loci(a[1]); names = list(L); out = a[2]
    n = {int(l): int(c) for l, c in (x.split() for x in open(out + '.n.tsv'))}
    sd = collections.defaultdict(list)  # SEDEF pair sides per chrom, for the co-duplicated count
    if len(a) > 3:
        for k, line in enumerate(open(a[3])):
            f = line.split('\t')
            if len(f) > 11 and f[1].isdigit() and f[10].isdigit():  # UCSC SEDEF track: partner in columns 10-12
                sd[f[0]].append((int(f[1]), int(f[2]), k)); sd[f[9]].append((int(f[10]), int(f[11]), k))
    @functools.lru_cache(None)
    def sdset(i):
        c, _, s, e, _ = L[names[i]]
        return {k for x, y, k in sd[c] if x < e and y > s}
    E = []
    for line in open(out + '.pairs.tsv'):
        A, B, s = map(int, line.split())
        la, lb = L[names[A]], L[names[B]]
        if la[0] == lb[0] and la[2] < lb[3] and lb[2] < la[3]: continue  # overlapping bundles are never linked
        E.append((A, B, s, s / (n[A] + n[B] - s)))
    with open(out + '.arms.tsv', 'w') as rep:
        rep.write("arm\tedges\tedges_sedef_held\tfamilies\tloci_in_families\tlargest\n")
        for floor in (1, 2):
            for t in (0.01, 0.05, 0.10, 0.20, 0.30):
                p = {}
                def find(x):
                    while p.setdefault(x, x) != x: p[x] = p[p[x]]; x = p[x]
                    return x
                ed = [(A, B) for A, B, s, j in E if s >= floor and j >= t]
                for A, B in ed: p[find(A)] = find(B)
                held = sum(1 for A, B in ed if sdset(A) & sdset(B)) if sd else 'NA'
                fam = collections.defaultdict(list)
                for x in p: fam[find(x)].append(x)
                arm = f"J{t:.2f}.s{floor}"
                with open(f"{out}.{arm}.copies.tsv", 'w') as w:
                    w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
                    for k, (r, mem) in enumerate(sorted(fam.items())):
                        for ci, x in enumerate(sorted(mem)):
                            c, sdn, s, e, exs = L[names[x]]
                            w.write(f"BJ{k}\t{ci}\t{names[x]}\t{c}\t{s}\t{e}\t{len(exs)}\t{sdn}\t{n.get(x, 0)}\t{','.join(f'{u}-{v}' for u, v in sorted(exs))}\n")
                rep.write(f"{arm}\t{len(ed)}\t{held}\t{len(fam)}\t{len(p)}\t{max(map(len, fam.values()), default=0)}\n")
elif a[0] == 'mcl':  # Amendment 3: MCL (bench/mcl_port.py) per connected component of the J >= t graph, weight = J
    import os, subprocess, multiprocessing
    sys.path.insert(0, os.path.expanduser('~/Rustle/bench'))
    import mcl_port
    mcl_port.subprocess = type('S', (), {'run': staticmethod(lambda *x, **k: subprocess.run(*x, timeout=3000, **k))})  # 50 min rule
    L = loci(a[1]); names = list(L); out = a[2]
    n = {int(l): int(c) for l, c in (x.split() for x in open(out + '.n.tsv'))}
    E = []
    for line in open(out + '.pairs.tsv'):
        A, B, s = map(int, line.split())
        la, lb = L[names[A]], L[names[B]]
        if la[0] == lb[0] and la[2] < lb[3] and lb[2] < la[3]: continue  # overlapping bundles are never linked
        E.append((A, B, s / (n[A] + n[B] - s)))
    def one(job):
        infl, ed = job
        nodes = {x for e in ed for x in e}
        try: cl = mcl_port.mcl({(str(A), str(B)): w for (A, B), w in ed.items()}, inflation=infl)
        except subprocess.TimeoutExpired: return [sorted(nodes)], len(nodes)
        cl = [[int(x) for x in c] for c in cl]
        assert sorted(x for c in cl for x in c) == sorted(nodes), "MCL clusters must partition the component"
        return cl, 0
    SD = len(a) > 3  # Amendment 4: SEDEF.bed (chr s e chr s e ...) -> arms SJ / SM / DM, default span rule (sdbio.py `holds`)
    if SD:
        by = collections.defaultdict(list)
        for k, line in enumerate(x for x in open(a[3]) if not x.startswith('#')):
            f = line.split('\t'); by[f[0]].append((int(f[1]), int(f[2]), k, 0)); by[f[3]].append((int(f[4]), int(f[5]), k, 1))
        def holds(bl, s, e):
            tot = sum(y - x for x, y in bl); ins = sum(max(0, min(y, e) - max(x, s)) for x, y in bl); lo, hi = bl[0][0], bl[-1][1]
            return tot > 0 and 2 * ins >= tot and 2 * max(0, min(hi, e) - max(lo, s)) >= hi - lo
        assert holds([(10, 20), (30, 40)], 0, 100) and not holds([(10, 20), (30, 40)], 0, 24) and not holds([(10, 20), (90, 100)], 0, 30)
        @functools.lru_cache(None)
        def H(i):
            c, _, _, _, exs = L[names[i]]; bl = sorted(exs)
            return frozenset((k, side) for s, e, k, side in by[c] if bl and s < bl[-1][1] and e > bl[0][0] and holds(bl, s, e))
        linked = lambda A, B: any((k, 1 - side) in H(B) for k, side in H(A))
        E = [(A, B, w, linked(A, B)) for A, B, w in E]
    with open(out + ('.sdarms.tsv' if SD else '.mclarms.tsv'), 'w') as rep, multiprocessing.Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 4))) as pool:
        rep.write("arm\tedges\tcomponents\tfamilies\tloci_in_families\tlargest\tunclustered_component_sizes\n")
        def arm(name, ed, infl):  # infl None = connected components
            p = {}
            def find(x):
                while p.setdefault(x, x) != x: p[x] = p[p[x]]; x = p[x]
                return x
            for A, B, _ in ed: p[find(A)] = find(B)
            comp = collections.defaultdict(dict)
            for A, B, w in ed: comp[find(A)][(A, B)] = w
            if infl is None: res = [([sorted({x for e in c for x in e})], 0) for c in comp.values()]
            else: res = pool.map(one, sorted(((infl, c) for c in comp.values()), key=lambda j: -len(j[1])), chunksize=1)
            fams = [c for cl, _ in res for c in cl if len(c) > 1]
            with open(f"{out}.{name}.copies.tsv", 'w') as w:
                w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
                for k, mem in enumerate(fams):
                    for ci, x in enumerate(sorted(mem)):
                        c, sdn, s, e, exs = L[names[x]]
                        w.write(f"BM{k}\t{ci}\t{names[x]}\t{c}\t{s}\t{e}\t{len(exs)}\t{sdn}\t{n.get(x, 0)}\t{','.join(f'{u}-{v}' for u, v in sorted(exs))}\n")
            rep.write(f"{name}\t{len(ed)}\t{len(comp)}\t{len(fams)}\t{sum(map(len, fams))}\t{max(map(len, fams), default=0)}\t{','.join(str(u) for _, u in res if u) or '-'}\n"); rep.flush()
        T = (0.30, 0.20, 0.10, 0.05, 0.01)
        if not SD:
            for t in T:
                for infl in (2.0, 4.0): arm(f"M{t:.2f}.I{infl}", [e for e in E if e[2] >= t], infl)
        else:
            for t in T:  # gate arms first (small graphs), damp arms after
                g = [(A, B, w) for A, B, w, h in E if w >= t and h]
                arm(f"SJ{t:.2f}", g, None); arm(f"SM{t:.2f}", g, 2.0)
            for t in T: arm(f"DM{t:.2f}", [(A, B, w if h else 0.25 * w) for A, B, w, h in E if w >= t], 2.0)
