#!/usr/bin/env python3
"""PREREG r1603 (docs/PREREG_tie_groups_2026-11-10.md): tie groups from the reads' own tie sets, and their score.

  tie_groups.py names GENES_TSV CAT_TSV BAM > names.txt
  samtools view -N names.txt BAM | tie_groups.py build PREFIX [exact|as98|asms98] [LINK_FRACTION] [LOCUS_FLOOR] [RMSK_BED[:MIN_FRACTION]] [unstranded]
  tie_groups.py score PREFIX GENES_TSV CAT_TSV SIGNAL_TSV [full]   (full = r1605 scorer; needs the DENOVO_* truth variables)
"""
import sys, re, collections, statistics

def blocks_of(pos, cigar):  # 0-based half-open aligned blocks, split at N
    out = []; s = e = pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(n)
        if op in 'MD=X': e += n
        elif op == 'N': out.append((s, e)); s = e = e + n
    out.append((s, e)); return out
assert blocks_of(10, '5M3N2M1D2M4S') == [(10, 15), (18, 23)]
ov = lambda a, b, x, y: a < y and x < b
def load_cat(path, only=None):
    G = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'gene_id' and (only is None or f[0] in only):
            G[f[0]] = (f[4], f[7], [tuple(map(int, b.split('-'))) for b in f[9].split(',')])
    return G
def table(path):
    with open(path) as fh:
        h = fh.readline().rstrip('\n').split('\t'); return [dict(zip(h, l.rstrip('\n').split('\t'))) for l in fh]

mode = sys.argv[1]
if mode == 'names':
    import pysam
    genes = {r['gene'] for r in table(sys.argv[2])}; G = load_cat(sys.argv[3], genes); bam = pysam.AlignmentFile(sys.argv[4]); names = set()
    for c, _, ex in G.values():
        for a in bam.fetch(c, ex[0][0], ex[-1][1]):
            if not a.is_supplementary and not a.is_unmapped: names.add(a.query_name)
    sys.stdout.write('\n'.join(sorted(names)) + '\n')

elif mode == 'build':
    pre = sys.argv[2]; rule = sys.argv[3] if len(sys.argv) > 3 else 'exact'  # exact (r1603) | as98 | asms98 (r1604)
    assert rule in ('exact', 'as98', 'asms98'), rule
    frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0  # r1605: a link also needs >= frac x the smaller locus's molecules
    unstranded = len(sys.argv) > 7 and sys.argv[7] == 'unstranded'  # r1608: an unspliced record has no strand; its loci are '.'
    recs = collections.defaultdict(list)  # name -> [(AS, chrom, strand, blocks, ms)]
    for l in sys.stdin:
        f = l.split('\t'); flag = int(f[1])
        if flag & 0x804: continue
        m = re.search(r'\tAS:i:(-?\d+)', l)
        ms = re.search(r'\tms:i:(-?\d+)', l)
        recs[f[0]].append((int(m.group(1)) if m else 0, f[2], '-' if flag & 16 else '+', blocks_of(int(f[3]) - 1, f[5]), int(ms.group(1)) if ms else 0))
    tied = []  # (chrom, strand, start, end, blocks, name)
    for n, rs in recs.items():
        b = max(r[0] for r in rs); bm = max(r[4] for r in rs)
        if rule == 'exact': t = [r for r in rs if r[0] == b]
        else: t = [r for r in rs if r[0] >= 0.98 * b > 0 or (rule == 'asms98' and r[4] >= 0.98 * bm > 0)]
        if len(t) >= 2: tied += [(c, '.' if unstranded and len(bl) == 1 else s, bl[0][0], bl[-1][1], bl, n) for _, c, s, bl, _ in t]
    tied.sort(key=lambda r: r[:4]); loci = []; key = None
    for r in tied:
        if loci and key == r[:2]:
            a, b = seed; o = min(b, r[3]) - max(a, r[2])
            if o >= 0.5 * min(b - a, r[3] - r[2]): loci[-1].append(r); continue
        loci.append([r]); key = r[:2]; seed = (r[2], r[3])
    floor = int(sys.argv[5]) if len(sys.argv) > 5 else 1  # r1606: a locus needs >= floor tied molecules to exist
    loci = [L for L in loci if len({r[5] for r in L}) >= floor]
    def model(L):  # the most common intron chain of the locus; ends = medians over the records with it
        ch = lambda r: tuple((r[4][x][1], r[4][x + 1][0]) for x in range(len(r[4]) - 1))
        chain = collections.Counter(ch(r) for r in L).most_common(1)[0][0]; same = [r for r in L if ch(r) == chain]
        s = int(statistics.median(r[2] for r in same)); e = int(statistics.median(r[3] for r in same))
        ex = [(s, chain[0][0])] + [(chain[x][1], chain[x + 1][0]) for x in range(len(chain) - 1)] + [(chain[-1][1], e)] if chain else [(s, e)]
        if any(a >= b for a, b in ex): ex = [(s, e)] if not chain else [(min(s, chain[0][0] - 1), chain[0][0])] + ex[1:-1] + [(chain[-1][1], max(e, chain[-1][1] + 1))]
        return ex
    models = [model(L) for L in loci]
    if len(sys.argv) > 6:  # r1607: RMSK_BED[:MIN_FRACTION], a locus mostly made of interspersed repeat is written to PREFIX.repeat.tsv and takes no part
        import bisect
        path, _, thr = sys.argv[6].partition(':'); thr = float(thr or 0.5); rm = collections.defaultdict(list)
        for l in open(path):
            f = l.split('\t', 8)
            if f[6] in ('LINE', 'SINE', 'LTR', 'DNA', 'Retroposon', 'RC'): rm[f[0]].append((int(f[1]), int(f[2])))
        starts = {}
        for c, v in rm.items():
            v.sort(); m = [list(v[0])]
            for x, y in v[1:]:
                if x <= m[-1][1]: m[-1][1] = max(m[-1][1], y)
                else: m.append([x, y])
            rm[c] = m; starts[c] = [x for x, _ in m]
        def rep_frac(c, ex):
            cov = 0; m = rm.get(c, [])
            for a, b in ex:
                i = max(bisect.bisect_right(starts.get(c, []), a) - 1, 0)
                while i < len(m) and m[i][0] < b: cov += max(0, min(b, m[i][1]) - max(a, m[i][0])); i += 1
            return cov / max(sum(b - a for a, b in ex), 1)
        fr = [rep_frac(L[0][0], ex) for L, ex in zip(loci, models)]
        with open(pre + '.repeat.tsv', 'w') as w:
            for L, ex, x in zip(loci, models, fr):
                if x >= thr: w.write(f'{L[0][0]}\t{ex[0][0]}\t{ex[-1][1]}\t{L[0][1]}\t{len({r[5] for r in L})}\t{x:.3f}\n')
        keep = [i for i, x in enumerate(fr) if x < thr]; loci = [loci[i] for i in keep]; models = [models[i] for i in keep]
    of = collections.defaultdict(set)
    for i, L in enumerate(loci):
        for r in L: of[r[5]].add(i)
    edge = collections.Counter()
    for s in of.values():
        s = sorted(s)
        for i in range(len(s)):
            for j in range(i + 1, len(s)): edge[(s[i], s[j])] += 1
    par = list(range(len(loci)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    nmol = [len({r[5] for r in L}) for L in loci]
    for (i, j), n in edge.items():
        if n >= 2 and n >= frac * min(nmol[i], nmol[j]): par[find(i)] = find(j)
    groups = collections.defaultdict(list)
    for i in range(len(loci)): groups[find(i)].append(i)
    with open(pre + '.gtf', 'w') as g, open(pre + '.groups.tsv', 'w') as t:
        t.write('group\tn_loci\tn_molecules\tloci\n')
        for k, (_, ids) in enumerate(sorted(groups.items())):
            mols = {r[5] for i in ids for r in loci[i]}; names = []
            for i in ids:
                L = loci[i]; ex = models[i]
                c, st = L[0][0], L[0][1]; a = f'gene_id "TG{k}"; transcript_id "TG{k}.L{i}"; group "TG{k}"; locus_molecules "{len({r[5] for r in L})}"; group_molecules "{len(mols)}"; group_loci "{len(ids)}";'
                g.write(f'{c}\ttie_groups\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{st}\t.\t{a}\n')
                for x, y in ex: g.write(f'{c}\ttie_groups\texon\t{x + 1}\t{y}\t.\t{st}\t.\t{a}\n')
                names.append(f'{c}:{ex[0][0]}-{ex[-1][1]}{st}')
            t.write(f'TG{k}\t{len(ids)}\t{len(mols)}\t{",".join(names)}\n')
    print(f'molecules {len(recs)} tied records {len(tied)} loci {len(loci)} groups {len(groups)}', file=sys.stderr)

elif mode == 'score':
    pre, genes_tsv, cat_tsv, sig_tsv = sys.argv[2:6]
    fam = {r['gene']: r['family'] for r in table(genes_tsv)}; G = load_cat(cat_tsv)
    idx = collections.defaultdict(list)  # (chrom, strand, 100 kb bin) -> genes
    for g, (c, s, ex) in G.items():
        for b in range(ex[0][0] // 100000, ex[-1][1] // 100000 + 1): idx[(c, s, b)].append(g)
    M = collections.defaultdict(lambda: [None, None, None, []])  # tid -> chrom, strand, group, exons
    for l in open(pre + '.gtf'):
        f = l.rstrip('\n').split('\t')
        if f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); m = M[t]
        m[0], m[1], m[2] = f[0], f[6], re.search(r'group "([^"]+)"', f[8]).group(1); m[3].append((int(f[3]) - 1, int(f[4])))
    hits = {}; grp = collections.defaultdict(list)
    for t, (c, s, g, ex) in M.items():
        cand = {x for b in range(ex[0][0] // 100000, ex[-1][1] // 100000 + 1) for st in ('+-' if s == '.' else s) for x in idx.get((c, st, b), [])}
        hits[t] = {x for x in cand if any(ov(a, b, p, q) for a, b in ex for p, q in G[x][2])}; grp[g].append(t)
    reach = collections.defaultdict(list)  # gene -> spans of models reaching it
    for t, hs in hits.items():
        for x in hs: reach[x].append(M[t][3][-1][1] - M[t][3][0][0])
    S = table(sig_tsv); ind = [r for r in S if r['signal'] == 'indistinguishable']; rest = [r for r in S if r['signal'] != 'indistinguishable']
    taf = [r for r in ind if r['name'].startswith('TAF11L')]
    print(f"B1 reach {sum(r['gene'] in reach for r in ind)} of {len(ind)} (bar 58)")
    judged = pure = 0; nocat = 0; impure = []
    for g, ts in grp.items():
        fams = collections.Counter(fam[x] for t in ts for x in hits[t] if x in fam)
        if not fams: continue
        top = fams.most_common(1)[0][0]; judged += 1
        foreign = [t for t in ts if hits[t] and not any(fam.get(x) == top for x in hits[t])]
        pure += not foreign; nocat += sum(not hits[t] for t in ts)
        if foreign: impure.append((top, [hits[t] for t in foreign]))
    print(f"B2 purity {pure} of {judged} = {pure / max(judged, 1):.3f} (bar 0.90); loci with no CAT gene in judged groups {nocat}")
    import os
    if os.environ.get('DENOVO_PARALOGS') and os.environ.get('DENOVO_CAT_GFF'):  # r1604 addendum: older Compara paralogues
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
        import famtruth as F
        ens = {g: d['ens'] for g, d in F.load_cat().items()}; par = collections.defaultdict(set)
        for l in open(os.environ['DENOVO_PARALOGS']):
            f = l.rstrip('\n').split('\t')
            if len(f) > 1 and f[1].startswith('ENSG'): par[f[0]].add(f[1]); par[f[1]].add(f[0])
        fam_ens = collections.defaultdict(set)
        for g, fm in fam.items(): fam_ens[fm] |= ens.get(g, set())
        rel = lambda top, hs: any(par[e] & fam_ens[top] for x in hs for e in ens.get(x, ()))
        older = sum(all(rel(top, hs) for hs in fl) for top, fl in impure)
        print(f"split: pure {pure}, older paralogue {older}, unrelated {len(impure) - older} of {judged}; pure or older = {(pure + older) / max(judged, 1):.3f}")
    if len(sys.argv) > 6 and sys.argv[6] == 'full':  # r1605: full truth families, same-family pseudogenes, grouped reach, group size
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
        import famtruth as F
        GG = F.load_cat(); pairs, _ = F.truth_pairs(GG, F.PARALOGS, set(F.LEVELS['primates']))
        full = {g: k for k, c in enumerate(F.components(pairs)) for g in c}
        pre_ = lambda n: re.match(r'[A-Za-z]*', n).group(0)
        stems = collections.defaultdict(set)
        for g, k in full.items(): stems[k].add(pre_(GG[g]['name']))
        def in_family(x, top):
            if full.get(x) == top: return True
            q = pre_(GG[x]['name']) if x in GG and 'pseudogene' in GG[x]['biotype'] else ''
            return bool(q) and (q in stems[top] or (q.endswith('P') and q[:-1] in stems[top]))
        sd = collections.defaultdict(list)  # r1607: (chrom, 100 kb bin) -> (start, end, partner chrom, start, end), both directions
        if os.environ.get('RUSTLE_SEDEF_BED'):
            for l in open(os.environ['RUSTLE_SEDEF_BED']):
                f = l.split('\t', 7)
                if l[0] == '#' or len(f) < 6 or not f[1].isdigit(): continue
                for c1, a1, b1, c2, a2, b2 in ((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])), (f[3], int(f[4]), int(f[5]), f[0], int(f[1]), int(f[2]))):
                    for bn in range(a1 // 100000, b1 // 100000 + 1): sd[(c1, bn)].append((a1, b1, c2, a2, b2))
        span = lambda t: (M[t][0], M[t][3][0][0], M[t][3][-1][1])
        def codup(t, famloci):  # t and a family locus of the group lie in one SEDEF pair (the same side or the two sides)
            c, a, b = span(t)
            for a1, b1, c2, a2, b2 in {x for bn in range(a // 100000, b // 100000 + 1) for x in sd.get((c, bn), [])}:
                if a < b1 and a1 < b and any((fc == c and fa < b1 and a1 < fb) or (fc == c2 and fa < b2 and a2 < fb) for fc, fa, fb in famloci): return True
            return False
        j2 = p2 = pd = 0
        for g, ts in grp.items():
            if not any(x in fam for t in ts for x in hits[t]): continue
            c = collections.Counter(full[x] for t in ts for x in hits[t] if x in full)
            if not c: continue
            top = c.most_common(1)[0][0]; j2 += 1
            foreign = [t for t in ts if hits[t] and not any(in_family(x, top) for x in hits[t])]
            p2 += not foreign
            famloci = [span(t) for t in ts if any(in_family(x, top) for x in hits[t])]
            pd += bool(foreign) and all(codup(t, famloci) for t in foreign)
        print(f"B2s purity counting co-duplicated loci as family {p2 + pd} of {j2} = {(p2 + pd) / max(j2, 1):.3f} (bar 0.90); co-duplicated-only groups {pd}")
        big = {g for g, ts in grp.items() if len(ts) >= 2}; gr = lambda r: any(M[t][2] in big for t in M if r['gene'] in hits[t])
        print(f"B2f purity (full truth + same-family pseudogenes) {p2} of {j2} = {p2 / max(j2, 1):.3f} (bar 0.90)")
        print(f"B1g grouped reach {sum(gr(r) for r in ind)} of {len(ind)} (bar 58); TAF11L grouped {sum(gr(r) for r in taf)} of {len(taf)}")
        print(f"B5 largest group {max(len(v) for v in grp.values())} loci (bar <= 60)")
    print(f"B3 TAF11L reached by a model <= 5 kb: {sum(any(x <= 5000 for x in reach.get(r['gene'], [])) for r in taf)} of {len(taf)} (bar 20)")
    print(f"descriptive: other hard members reached {sum(r['gene'] in reach for r in rest)} of {len(rest)}; groups {len(grp)}, loci {len(M)}; group sizes", sorted(collections.Counter(min(len(v), 20) for v in grp.values()).items()))
    print('unreached indistinguishable:', sorted(r['name'] for r in ind if r['gene'] not in reach))
