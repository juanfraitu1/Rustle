#!/usr/bin/env python3
"""Multi-copy mode, the rescue of copies no two reads agree on (rows 1863-1867; driver stage `multicopy`).

    multicopy_rescue.py candidates BASE.gtf FLOOR1.gtf CAND.gtf
    multicopy_rescue.py select     BASE.gtf CAND.gtf UNION.fam.copies.tsv OUT.gtf OUT.tsv
    multicopy_rescue.py annotate   MULTICOPY.gtf EQ.loci.tsv OUT.gtf GROUPS.tsv
    multicopy_rescue.py --selftest

candidates: the spliced transcripts of FLOOR1 (the same assembly at a 1-read floor) that lie where BASE has no transcript on that
  strand, written with gene and transcript ids prefixed `MC_` (FLOOR1 and BASE number their loci alike).
select: after the families stage has run on BASE + CAND, keep the candidates whose locus is a copy of a family of >= 2 loci.
  OUT.gtf = BASE + the kept candidates (transcript lines tagged multicopy_rescue "1"); OUT.tsv = one line per kept locus:
  gene_id, family_id, family_copies, base_copies (copies of the family that come from BASE).
annotate (row 1871): read counts that do not multiply over indistinguishable copies. EQ.loci.tsv = tools/eq_classes.py on
  MULTICOPY.gtf (exact band: records whose score equals the molecule's best). Every transcript line gets its LOCUS's
  reads_placed (molecules whose best record is here alone), reads_em (tied molecules shared out, sums to each molecule once),
  reads_max (placed + every tied molecule), tie_group / tie_group_loci (loci joined by >= 2 tied molecules; "." = none) and
  reads_tieadj = the transcript's own `reads` x reads_em / reads_max (an approximation: the split is the locus's, not the
  transcript's). GROUPS.tsv = one line per tie group: its loci, molecules placed in it, tied molecules, and the copies the
  scores allow, min_expressed = max(1, loci with >= 2 placed molecules) to max_expressed = its loci.
Why confined to families: genome-wide the candidates' intron chains are 3 to 7 times less often annotated than the base
assembly's (row 1867), while in multi-copy families they recover copies whose reads each carry a different chain (row 1865)."""
import bisect, collections, re, sys

TID = re.compile(r'transcript_id "([^"]+)"'); GID = re.compile(r'gene_id "([^"]+)"')

def spans(path):
    iv = collections.defaultdict(list)
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'transcript': iv[(f[0], f[6])].append((int(f[3]), int(f[4])))
    out = {}
    for k, v in iv.items():
        v.sort(); m = []
        for a, b in v:
            if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
            else: m.append([a, b])
        out[k] = (m, [a for a, _ in m])
    return out

def occupied(sp, chrom, strand, a, b):
    if (chrom, strand) not in sp: return False
    m, st = sp[(chrom, strand)]; i = bisect.bisect_right(st, b) - 1
    return i >= 0 and m[i][1] >= a

def by_transcript(path):
    tx = collections.OrderedDict()
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8:
            m = TID.search(f[8])
            if m: tx.setdefault(m.group(1), []).append(l)
    return tx

def candidates(base, floor1, out, mono=False):
    """mono: single-exon floor-1 transcripts are candidates too (row 1919); they are tagged at `select`."""
    sp = spans(base); n = 0
    with open(out, 'w') as w:
        for t, L in by_transcript(floor1).items():
            ex = [l.split('\t') for l in L if l.split('\t')[2] == 'exon']
            if len(ex) < 2 and not mono: continue
            a, b = min(int(e[3]) for e in ex), max(int(e[4]) for e in ex)
            if occupied(sp, ex[0][0], ex[0][6], a, b): continue
            if len(ex) < 2 and ex[0][6] not in '+-' and (occupied(sp, ex[0][0], '+', a, b) or occupied(sp, ex[0][0], '-', a, b)): continue   # a strandless single exon under any model
            n += 1
            for l in L: w.write(l.replace('gene_id "', 'gene_id "MC_', 1).replace('transcript_id "', 'transcript_id "MC_', 1))
    return n

def select(base, cand, copies, out, table):
    L = [l.rstrip('\n').split('\t') for l in open(copies)]; h = {c: i for i, c in enumerate(L[0])}
    fam = collections.defaultdict(set); of = {}
    for f in L[1:]:
        g = f[h['gene_id']]; fam[f[h['family_id']]].add(g); of[g] = f[h['family_id']]
    have = {m.group(1) for l in open(cand) for m in [GID.search(l)] if m}
    keep = {g for g, k in of.items() if g in have and len(fam[k]) >= 2}
    nex = collections.Counter(TID.search(l.split('\t')[8]).group(1) for l in open(cand) if l.split('\t')[2] == 'exon'); mono_tx = {t: k == 1 for t, k in nex.items()}
    n = 0
    with open(out, 'w') as w:
        w.write(open(base).read())
        for l in open(cand):
            f = l.split('\t')
            if len(f) > 8 and GID.search(f[8]).group(1) in keep:
                if f[2] == 'transcript': n += 1; l = l.rstrip('\n').rstrip().rstrip(';') + f'; multicopy_rescue "{"mono" if mono_tx.get(TID.search(f[8]).group(1)) else "1"}";\n'
                w.write(l)
    with open(table, 'w') as w:
        w.write('gene_id\tfamily_id\tfamily_copies\tbase_copies\n')
        for g in sorted(keep): w.write(f"{g}\t{of[g]}\t{len(fam[of[g]])}\t{sum(not x.startswith('MC_') for x in fam[of[g]])}\n")
    return len(keep), n

def annotate(gtf, eq, out, groups):
    L = [l.rstrip('\n').split('\t') for l in open(eq)]; h = {c: i for i, c in enumerate(L[0])}
    loc = {f[h['locus']]: f for f in L[1:]}
    num = lambda f, c: float(f[h[c]])
    n = 0
    with open(out, 'w') as w:
        for l in open(gtf):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'transcript':
                g = GID.search(f[8]).group(1); e = loc.get(g)
                if e:
                    mx = num(e, 'max'); r = re.search(r'reads "([0-9.]+)"', f[8]); adj = f"{float(r.group(1)) * num(e, 'em') / mx:.1f}" if r and mx > 0 else '.'
                    l = (l.rstrip('\n').rstrip().rstrip(';') + f'; reads_placed "{e[h["placed"]]}"; reads_em "{e[h["em"]]}"; reads_max "{e[h["max"]]}"; '
                         f'tie_group "{e[h["group"]]}"; tie_group_loci "{e[h["group_loci"]]}"; reads_tieadj "{adj}";\n'); n += 1
            w.write(l)
    gr = collections.defaultdict(list)
    for g, e in loc.items():
        if e[h['group']] != '.': gr[e[h['group']]].append(e)
    with open(groups, 'w') as w:
        w.write('tie_group\tloci\tchroms\tplaced\ttied\tmin_expressed\tmax_expressed\tmembers(locus:placed:em)\n')
        for k, es in sorted(gr.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            placed = sum(num(e, 'placed') for e in es); em = sum(num(e, 'em') for e in es)
            lo = max(1, sum(num(e, 'placed') >= 2 for e in es))
            w.write(f"{k}\t{len(es)}\t{len({e[h['chrom']] for e in es})}\t{placed:.0f}\t{em - placed:.0f}\t{lo}\t{len(es)}\t"
                    + ','.join(f"{e[h['locus']]}:{e[h['placed']]}:{e[h['em']]}" for e in sorted(es, key=lambda e: e[h['locus']])) + '\n')
    return n, len(gr)

def selftest():
    import os, tempfile
    def tx(c, s, g, t, exons, extra=''):
        a = f'gene_id "{g}"; transcript_id "{t}";{extra}'
        return ''.join([f"{c}\tx\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{s}\t.\t{a}\n"] + [f"{c}\tx\texon\t{x}\t{y}\t.\t{s}\t.\t{a}\n" for x, y in exons])
    with tempfile.TemporaryDirectory() as d:
        p = lambda n: os.path.join(d, n)
        open(p('b.gtf'), 'w').write(tx('c1', '+', 'G1', 'G1.1', [(100, 200), (300, 400)]))
        open(p('f.gtf'), 'w').write(tx('c1', '+', 'G1', 'G1.1', [(100, 200), (300, 400)])      # same place as base: not a candidate
                                    + tx('c1', '+', 'G1', 'G1.2', [(150, 200), (300, 450)])    # overlaps base: not a candidate
                                    + tx('c1', '-', 'G2', 'G2.1', [(150, 200), (300, 450)])    # other strand: candidate
                                    + tx('c1', '+', 'G3', 'G3.1', [(1000, 1100), (1200, 1300)])  # empty place: candidate
                                    + tx('c1', '+', 'G4', 'G4.1', [(5000, 5100)])               # unspliced: never
                                    + tx('c2', '+', 'G5', 'G5.1', [(10, 20), (30, 40)]))        # empty contig: candidate
        assert candidates(p('b.gtf'), p('f.gtf'), p('c.gtf')) == 3
        assert candidates(p('b.gtf'), p('f.gtf'), p('cm.gtf'), mono=True) == 4 and 'MC_G4.1' in open(p('cm.gtf')).read()
        open(p('km.tsv'), 'w').write('family_id\tcopy_idx\tgene_id\nF0\t0\tG1\nF0\t1\tMC_G4\n'); assert select(p('b.gtf'), p('cm.gtf'), p('km.tsv'), p('om.gtf'), p('om.tsv')) == (1, 1) and 'multicopy_rescue "mono"' in open(p('om.gtf')).read()
        ids = set(TID.findall(open(p('c.gtf')).read())); assert ids == {'MC_G2.1', 'MC_G3.1', 'MC_G5.1'}, ids
        open(p('k.tsv'), 'w').write('family_id\tcopy_idx\tgene_id\nF0\t0\tG1\nF0\t1\tMC_G3\nF1\t0\tMC_G5\nF2\t0\tMC_G2\nF2\t1\tMC_G9\n')
        assert select(p('b.gtf'), p('c.gtf'), p('k.tsv'), p('o.gtf'), p('o.tsv')) == (2, 2)   # G3 (with a base copy), G2 (family of two candidates); G5 alone in its family
        o = open(p('o.gtf')).read(); assert o.count('multicopy_rescue "1"') == 2 and 'MC_G5' not in o and 'G1.1' in o
        rows = [l.split('\t') for l in open(p('o.tsv')).read().splitlines()[1:]]; assert rows == [['MC_G2', 'F2', '2', '0'], ['MC_G3', 'F0', '2', '1']], rows
        open(p('e.tsv'), 'w').write('locus\tchrom\tplaced\tmax\tem\tgroup\tgroup_loci\nG1\tc1\t0\t30\t10.0\tG7\t3\nMC_G3\tc1\t4\t34\t14.0\tG7\t3\nMC_G2\tc1\t5\t5\t5.0\t.\t0\n')
        open(p('b2.gtf'), 'w').write(tx('c1', '+', 'G1', 'G1.1', [(100, 200), (300, 400)], ' reads "30";') + tx('c1', '+', 'MC_G3', 'MC_G3.1', [(1000, 1100), (1200, 1300)], ' reads "17";') + tx('c1', '-', 'MC_G2', 'MC_G2.1', [(150, 200), (300, 450)], ' reads "5";') + tx('c1', '+', 'G8', 'G8.1', [(9000, 9100), (9200, 9300)], ' reads "2";'))
        assert annotate(p('b2.gtf'), p('e.tsv'), p('a.gtf'), p('g.tsv')) == (3, 1)
        a = open(p('a.gtf')).read(); assert 'reads_tieadj "10.0"' in a and 'reads_tieadj "7.0"' in a and 'reads_tieadj "5.0"' in a and a.count('tie_group "G7"') == 2 and a.count('tie_group') == 6, a
        assert len(a.splitlines()) == len(open(p('b2.gtf')).read().splitlines())
        gl = open(p('g.tsv')).read().splitlines(); assert gl[1].split('\t')[:7] == ['G7', '2', '1', '4', '20', '1', '2'], gl
    print('multicopy_rescue selftest ok')

if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['--selftest']: selftest()
    elif len(a) in (4, 5) and a[0] == 'candidates' and a[4:] in ([], ['--mono']): print(f"{candidates(*a[1:4], mono=len(a) == 5)} {'spliced or single-exon' if len(a) == 5 else 'spliced'} floor-1 transcripts where the assembly has none on that strand")
    elif len(a) == 6 and a[0] == 'select':
        g, n = select(*a[1:]); print(f"{g} rescued loci ({n} transcripts) are copies of a family of >= 2 loci")
    elif len(a) == 5 and a[0] == 'annotate':
        n, g = annotate(*a[1:]); print(f"{n} transcripts annotated with tie-adjusted counts; {g} tie groups")
    else: sys.exit(__doc__)
