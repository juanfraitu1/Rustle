#!/usr/bin/env python3
"""Plant reads at the copies of known gene families, with the copy of origin of every read recorded (O2 family benchmark).

    plant.py ANNOT.gff3.gz GENOME.fa OUTDIR [--scenarios E,D,S] [--reads 60] [--isoforms 2] [--seed 1903]
             [--sub 0.001 --ins 0.002 --dele 0.002] [--full 0.7] [--ragged 20] [--private RATE]
    plant.py --from-run PREFIX.gtf PREFIX.fam.copies.tsv GENOME.fa OUTDIR [same options]

Copies: the protein-coding members of famtruth's MAINSTAYS + YAGS on the annotation (CHM13 CAT/Liftoff; env as famtruth).
With --from-run (driver stage `o2`, any genome, no annotation): the copies of every family of >= 2 loci in the run's own
copy table, and as their isoforms the run's own transcripts at those loci; needs GENOME.fa.fai (samtools faidx), no pysam.
Per copy up to --isoforms annotated transcripts with distinct intron chains (most exons first, then longest); a copy's
reads go 2 : 1 to its first and second isoform.
--private RATE (register 1924; default 0 = the sample is the reference): every copy gets its own variants at RATE per
exonic base (80% substitutions, 10% one-base deletions, 10% one-base insertions), by genome position so its isoforms
share them; the reads come from that sequence. Writes sample.fa (the sample's transcripts) and private.tsv.
Scenarios (per family, drawn with --seed):
  E  every copy --reads reads
  D  one copy 10 x --reads, the others --reads / 3
  S  half of the copies (rounded down) silent, the others --reads   (families of >= 2 copies)
Reads: the transcript's sequence on its strand; a share --full complete, the others cut at the 5' end at a uniform point
leaving >= 500 bases; then 0..--ragged bases (default 20) are taken off each end, so reads of one transcript do not all
share their two coordinates (with --ragged 0, the setting of rows 1903-1911, the full-length reads of a transcript are
coordinate duplicates, which Rustle counts once toward its 2-read floor: row 1916); substitutions, insertions,
deletions per base. No read is a readthrough or a chimera.
Writes OUTDIR/loci.gtf (planted transcripts; gene_id = copy; attributes family, name), plan.tsv (family, copy, name,
transcript, exons, length, reads per scenario), SC.fq, SC.truth.tsv (read, copy, transcript).
--selftest runs the pure parts."""
import sys, os, re, gzip, random, argparse, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def pick_isoforms(txs, k):
    """txs {tid: sorted [(s, e)] 1-based closed} -> up to k tids with distinct intron chains, most exons then longest."""
    out, seen = [], set()
    for t in sorted(txs, key=lambda t: (-len(txs[t]), -sum(e - s + 1 for s, e in txs[t]), t)):
        chain = tuple((a[1], b[0]) for a, b in zip(txs[t], txs[t][1:])) or ('mono', txs[t][0])
        if chain in seen: continue
        seen.add(chain); out.append(t)
        if len(out) == k: break
    return out

def plan_counts(fam, reads, rng):
    """fam {family: [copies]} -> {scenario: {copy: reads}}"""
    p = {'E': {}, 'D': {}, 'S': {}}
    for f in sorted(fam):
        gs = fam[f]
        for g in gs: p['E'][g] = reads
        d = rng.choice(gs)
        for g in gs: p['D'][g] = 10 * reads if g == d else max(1, reads // 3)
        if len(gs) >= 2:
            sil = set(rng.sample(gs, len(gs) // 2))
            for g in gs: p['S'][g] = 0 if g in sil else reads
    return p

def split(n, k):
    """n reads over k isoforms, 2 : 1 for two."""
    if k == 1: return [n]
    a = round(n * 2 / 3); return [a, n - a] + [0] * (k - 2)

def mutate(s, rng, sub, ins, dele, full, ragged=0):
    if rng.random() >= full and len(s) > 500: s = s[rng.randint(0, len(s) - 500):]
    if ragged and len(s) > 4 * ragged + 100:   # ragged ends: no two reads of a transcript share both coordinates (row 1916)
        a, b = rng.randint(0, ragged), rng.randint(0, ragged); s = s[a:len(s) - b]
    o = []
    for b in s:
        x = rng.random()
        if x < sub: o.append(rng.choice([y for y in 'ACGT' if y != b]))
        elif x < sub + ins: o.append(b); o.append(rng.choice('ACGT'))
        elif x < sub + ins + dele: pass
        else: o.append(b)
    return ''.join(o)

class Fai:
    """sequence fetch through GENOME.fa.fai (0-based half-open), without pysam"""
    def __init__(self, fa):
        self.fh = open(fa, 'rb'); self.ix = {}
        if not os.path.exists(fa + '.fai'): sys.exit(f'plant: {fa}.fai is missing (samtools faidx {fa})')
        for l in open(fa + '.fai'):
            f = l.rstrip('\n').split('\t'); self.ix[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    def fetch(self, c, s, e):
        n, off, lb, lw = self.ix[c]; s = max(0, s); e = min(n, e)
        if e <= s: return ''
        a = off + s // lb * lw + s % lb; b = off + (e - 1) // lb * lw + (e - 1) % lb + 1
        self.fh.seek(a); return self.fh.read(b - a).replace(b'\n', b'').replace(b'\r', b'').decode()

def from_run(gtf, copies):
    """-> members {gene_id: family_id} (families of >= 2 loci), tx {gene: {tid: [(s, e)]}}, loc {tid: (chrom, strand)}, names"""
    fam = collections.defaultdict(list)
    hdr = None
    for l in open(copies):
        f = l.rstrip('\n').split('\t')
        if hdr is None: hdr = {k: i for i, k in enumerate(f)}; continue
        g = f[hdr['gene_id']]
        if g not in fam[f[hdr['family_id']]]: fam[f[hdr['family_id']]].append(g)
    members = {g: f for f, gs in fam.items() if len(gs) >= 2 for g in gs}
    tx = collections.defaultdict(lambda: collections.defaultdict(list)); loc = {}
    for l in open(gtf):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        g = re.search(r'gene_id "([^"]+)"', f[8]); t = re.search(r'transcript_id "([^"]+)"', f[8])
        if g and t and g.group(1) in members:
            tx[g.group(1)][t.group(1)].append((int(f[3]), int(f[4]))); loc[t.group(1)] = (f[0], f[6] if f[6] in '+-' else '+')
    return members, tx, loc

def selftest():
    import tempfile
    d = tempfile.mkdtemp(); fa = os.path.join(d, 'g.fa')
    open(fa, 'w').write('>c1\nACGTACGTAC\nGGGGCCCCAA\nTT\n>c2\nAAAA\n'); open(fa + '.fai', 'w').write('c1\t22\t4\t10\t11\nc2\t4\t33\t4\t5\n')
    F = Fai(fa); assert F.fetch('c1', 8, 14) == 'ACGGGG' and F.fetch('c1', 18, 22) == 'AATT' and F.fetch('c2', 0, 4) == 'AAAA' and F.fetch('c1', 20, 99) == 'TT'
    open(os.path.join(d, 'x.gtf'), 'w').write('c1\tx\texon\t1\t5\t.\t+\t.\tgene_id "A"; transcript_id "a1";\nc1\tx\texon\t9\t12\t.\t+\t.\tgene_id "A"; transcript_id "a1";\n'
        'c2\tx\texon\t1\t4\t.\t.\t.\tgene_id "B"; transcript_id "b1";\nc2\tx\texon\t1\t4\t.\t-\t.\tgene_id "C"; transcript_id "c1";\n')
    open(os.path.join(d, 'x.tsv'), 'w').write('family_id\tgene_id\nF0\tA\nF0\tB\nF0\tA\nF1\tC\n')
    m, tx, loc = from_run(os.path.join(d, 'x.gtf'), os.path.join(d, 'x.tsv'))
    assert m == {'A': 'F0', 'B': 'F0'} and tx['A']['a1'] == [(1, 5), (9, 12)] and loc['b1'] == ('c2', '+') and 'C' not in tx
    t = {'a': [(1, 10), (21, 30), (41, 50)], 'b': [(1, 10), (21, 30), (41, 60)], 'c': [(1, 10), (41, 50)], 'd': [(5, 50)]}
    assert pick_isoforms(t, 2) == ['b', 'c'], pick_isoforms(t, 2)     # a and b share a chain; b is longer
    assert pick_isoforms(t, 9) == ['b', 'c', 'd']
    r = random.Random(1); p = plan_counts({'F': ['x', 'y', 'z'], 'G': ['q']}, 60, r)
    assert set(p['E'].values()) == {60} and sorted(p['D'][g] for g in 'xyz') == [20, 20, 600] and p['D']['q'] == 600
    assert sorted(p['S'][g] for g in 'xyz') == [0, 60, 60] and 'q' not in p['S']
    assert split(60, 2) == [40, 20] and split(1, 2) == [1, 0] and split(7, 1) == [7]
    s = 'ACGT' * 300
    assert mutate(s, random.Random(2), 0, 0, 0, 1.0) == s and 500 <= len(mutate(s, random.Random(3), 0, 0, 0, 0.0)) <= len(s)
    r = [mutate(s, random.Random(i), 0, 0, 0, 1.0, 20) for i in range(30)]; assert all(len(s) - 40 <= len(x) <= len(s) and x in s for x in r) and len(set(r)) > 10
    print('plant selftest ok')

def main():
    if '--selftest' in sys.argv: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--from-run', nargs=2, metavar=('GTF', 'COPIES')); ap.add_argument('pos', nargs='+')
    ap.add_argument('--scenarios', default='E,D,S'); ap.add_argument('--reads', type=int, default=60); ap.add_argument('--isoforms', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1903); ap.add_argument('--sub', type=float, default=0.001); ap.add_argument('--ins', type=float, default=0.002)
    ap.add_argument('--dele', type=float, default=0.002); ap.add_argument('--full', type=float, default=0.7); ap.add_argument('--ragged', type=int, default=20)
    ap.add_argument('--private', type=float, default=0.0, metavar='RATE')
    a = ap.parse_args()
    if a.from_run:
        if len(a.pos) != 2: ap.error('--from-run GTF COPIES GENOME.fa OUTDIR')
        a.genome, a.out = a.pos
        members, tx, loc = from_run(*a.from_run); G = {g: {'name': g} for g in members}; FA = Fai(a.genome)
    else:
        if len(a.pos) != 3: ap.error('ANNOT.gff3.gz GENOME.fa OUTDIR')
        a.annot, a.genome, a.out = a.pos
        import pysam, famtruth as F
        G = F.load_cat(); pats = {**F.MAINSTAYS, **F.YAGS}; members = {}
        for name, pat in pats.items():
            rx = re.compile(pat)
            for g, d in G.items():
                if rx.match(d['name']) and d['biotype'] == 'protein_coding': members[g] = name
        tx = collections.defaultdict(lambda: collections.defaultdict(list)); t2g = {}; loc = {}
        for l in gzip.open(a.annot, 'rt'):
            f = l.rstrip('\n').split('\t')
            if len(f) < 9: continue
            at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
            if f[2] == 'transcript' and at['Parent'] in members: t2g[at['ID']] = at['Parent']; loc[at['ID']] = (f[0], f[6])
            elif f[2] == 'exon' and at['Parent'] in t2g: tx[t2g[at['Parent']]][at['Parent']].append((int(f[3]), int(f[4])))
        FA = pysam.FastaFile(a.genome)
    rng = random.Random(a.seed); os.makedirs(a.out, exist_ok=True)
    iso, seqs = {}, {}; priv = {}; prng = random.Random(a.seed * 7919 + 1)      # --private: a stream of its own, the read draws are unchanged
    with open(f'{a.out}/loci.gtf', 'w') as w:
        for g in sorted(tx):
            for t in tx[g]: tx[g][t].sort()
            iso[g] = pick_isoforms(tx[g], a.isoforms)
            for t in iso[g]:
                ex = tx[g][t]; c, st = loc[t]
                s = ''.join(FA.fetch(c, x - 1, y) for x, y in ex).upper()
                if a.private > 0:          # the sample's own variants at this copy, by genome position (shared by its isoforms)
                    pv = priv.setdefault(g, {}); o = []; i = 0
                    for x, y in ex:
                        for pos in range(x, y + 1):
                            if pos not in pv:
                                u = prng.random(); b0 = s[i]
                                pv[pos] = None if u >= a.private else ('S' + prng.choice([z for z in 'ACGT' if z != b0])) if u < a.private * 0.8 else 'D' if u < a.private * 0.9 else 'I' + prng.choice('ACGT')
                            v = pv[pos]
                            o.append(s[i] if v is None else v[1] if v[0] == 'S' else '' if v == 'D' else s[i] + v[1]); i += 1
                    s = ''.join(o)
                seqs[t] = s[::-1].translate(COMP) if st == '-' else s
                at = f'gene_id "{g}"; transcript_id "{t}"; family "{members[g]}"; ref_gene_name "{G[g]["name"]}";'
                w.write(f'{c}\tplant\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{st}\t.\t{at}\n')
                for x, y in ex: w.write(f'{c}\tplant\texon\t{x}\t{y}\t.\t{st}\t.\t{at}\n')
    if a.private > 0:
        with open(f'{a.out}/sample.fa', 'w') as w:
            for t in sorted(seqs): w.write(f'>{t}\n{seqs[t]}\n')
        with open(f'{a.out}/private.tsv', 'w') as w:
            w.write('copy\tposition\tvariant\n')
            for g in sorted(priv):
                for pos in sorted(priv[g]):
                    if priv[g][pos] is not None: w.write(f'{g}\t{pos}\t{priv[g][pos]}\n')
    fam = collections.defaultdict(list)
    for g in sorted(iso): fam[members[g]].append(g)
    plan = plan_counts(fam, a.reads, rng); scs = a.scenarios.split(',')
    per = {sc: {g: dict(zip(iso[g], split(n, len(iso[g])))) for g, n in plan[sc].items()} for sc in scs}
    with open(f'{a.out}/plan.tsv', 'w') as w:
        w.write('family\tcopy\tname\ttranscript\texons\tlength\t' + '\t'.join(scs) + '\n')
        for g in sorted(iso, key=lambda g: (members[g], G[g]['name'], g)):
            for t in iso[g]:
                w.write(f"{members[g]}\t{g}\t{G[g]['name']}\t{t}\t{len(tx[g][t])}\t{len(seqs[t])}\t" + '\t'.join(str(per[sc].get(g, {}).get(t, 'NA')) for sc in scs) + '\n')
    for sc in scs:
        n = 0
        with open(f'{a.out}/{sc}.fq', 'w') as q, open(f'{a.out}/{sc}.truth.tsv', 'w') as tr:
            for g in sorted(per[sc]):
                for t, k in per[sc][g].items():
                    for i in range(k):
                        s = mutate(seqs[t], rng, a.sub, a.ins, a.dele, a.full, a.ragged); nm = f'{sc}|{g}|{t}|{i}'; n += 1
                        q.write(f'@{nm}\n{s}\n+\n{"I" * len(s)}\n'); tr.write(f'{nm}\t{g}\t{t}\n')
        print(sc, 'copies', len(per[sc]), 'silent', sum(1 for g in plan[sc] if plan[sc][g] == 0), 'reads', n)

if __name__ == '__main__':
    main()
