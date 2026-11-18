#!/usr/bin/env python3
"""PREREG r1775 (docs/PREREG_eichler_baseline_2026-11-17.md): read-origin truth that needs no catalogue and no arm, and the
score of placement rules against it.

  neutral_truth.py select BAM OUT.fq          reads with a primary on a truth-family member (env as famtruth, PANEL_SPECIES)
  neutral_truth.py score REF.fa MAT.mmi PAT.mmi OUTDIR MM.bam [NAME=OTHER.bam ...]
        MM.bam / OTHER.bam = the selected reads aligned to REF by minimap2 / another aligner (all records, any order)

Placements = a read's records merged into loci. Origin = the placement whose reference sequence lands, on the
individual's haplotypes, where the read itself has its best (slack 0) hit. Arms: primary, eichler-keep, eichler-move,
and each OTHER aligner's primary.

KNOWN DEFECT, repaired under NT_REPAIR=1 (r1785; opt-in, without it the old test runs) (docs/NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md): the origin test is a span overlap without strand or exon check,
and a placement with no near-perfect haplotype hit is dropped as a candidate; both mislabel reads (chr1, row 1783). Repair
it under its own pre-registered row before scoring another arm."""
import sys, os, re, hashlib, collections, subprocess
import pysam

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CAPS = {'S0': 5000, 'S1': 30000, 'S2': 30000}; MAXLOCI = 12; TH = int(os.environ.get('NT_THREADS', os.cpu_count()))
if os.environ.get('NT_CAP'): CAPS = {k: int(os.environ['NT_CAP']) for k in CAPS}  # smoke test

def selected(name, pct, offset=0):
    """hash sample of the reads: residues offset .. offset + pct - 1 of 100 (another offset = reads not yet looked at)"""
    return (int(hashlib.sha256(name.encode()).hexdigest()[:8], 16) - offset) % 100 < pct

def select(bam_path, out):
    PCT = int(os.environ.get('NT_SELECT_PCT', 100)); OFFSET = int(os.environ.get('NT_SELECT_OFFSET', 0))  # hash sample of the selected reads, for deep samples
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
    import famtruth as F
    G = F.load_genes(os.environ.get('PANEL_SPECIES', 'human'))
    pairs, _ = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['simiiformes']))
    bam = pysam.AlignmentFile(bam_path); seen = set(); n = 0
    with open(out, 'w') as w:
        for f in F.components(pairs):
            for g in f:
                v = G[g]
                if not v['exons'] or v['chrom'] not in bam.references: continue
                for r in bam.fetch(v['chrom'], v['exons'][0][0], v['exons'][-1][1]):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_name in seen: continue
                    if r.is_reverse == (v['strand'] == '+'): continue
                    if not any(a < y and x < b for a, b in r.get_blocks() for x, y in v['exons']): continue
                    seen.add(r.query_name)
                    if not selected(r.query_name, PCT, OFFSET): continue
                    s = r.query_sequence
                    if r.is_reverse: s = s.translate(COMP)[::-1]
                    w.write(f">{r.query_name}\n{s}\n"); n += 1
    print(f"selected {n} reads with a primary on a truth-family member")

def records(path):
    """read -> [(chrom, strand, start, end, blocks, AS, is_primary)], supplementary records skipped."""
    d = collections.defaultdict(list)
    for r in pysam.AlignmentFile(path, check_sq=False):
        if r.is_unmapped or r.is_supplementary: continue
        d[r.query_name].append((r.reference_name, '-' if r.is_reverse else '+', r.reference_start, r.reference_end, r.get_blocks(), r.get_tag('AS'), not r.is_secondary))
    return d

def loci_of(recs):
    """merge records into loci: [dict(chrom, strand, s, e, blocks (of its best-AS record), AS, primary)], best AS first."""
    out = []
    for c, st, s, e, bl, a, p in sorted(recs, key=lambda x: -x[5]):
        for L in out:
            if L['chrom'] == c and L['strand'] == st and L['s'] < e and s < L['e']:
                L['s'] = min(L['s'], s); L['e'] = max(L['e'], e); L['primary'] |= p; break
        else: out.append(dict(chrom=c, strand=st, s=s, e=e, blocks=bl, AS=a, primary=p))
    return out

def find(loci, rec):
    return next((i for i, L in enumerate(loci) if L['chrom'] == rec[0] and L['strand'] == rec[1] and L['s'] < rec[3] and rec[2] < L['e']), None)

def hap_hits(fa, mmi, hap):
    """query -> [(hap, chrom, s, e, NM, aligned)] from minimap2 against one haplotype index."""
    p = subprocess.Popen(['minimap2', '-a', '-x', 'splice:hq', '-uf', '--eqx', '-N', '20', '-p', '0.5', '--secondary=yes', '-t', str(TH), mmi, fa], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    out = collections.defaultdict(list)
    for l in p.stdout:
        if l[0] == '@': continue
        f = l.split('\t')
        if int(f[1]) & 4 or int(f[1]) & 0x800: continue
        ops = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
        end = int(f[3]) - 1 + sum(int(n) for n, op in ops if op in 'MDN=X'); al = sum(int(n) for n, op in ops if op in 'M=X')
        nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), 10 ** 6)
        bl = []; pos = int(f[3]) - 1; cur = pos
        for n, op in ops:
            if op in 'MD=X': pos += int(n)
            elif op == 'N':
                bl.append((cur, pos)); pos += int(n); cur = pos
        bl.append((cur, pos))
        out[f[0]].append((hap, f[2], int(f[3]) - 1, end, nm, al, '-' if int(f[1]) & 16 else '+', bl))
    p.wait(); return out

REPAIR = os.environ.get('NT_REPAIR') == '1'  # r1785; shortcut: opt-in while other sessions score with the old test, make it the default on the user's word

def same_place(a, b):
    """do two haplotype hits name one place? old test: span overlap; r1785 rule 1: same strand and a shared exon base"""
    if a[0] != b[0] or a[1] != b[1] or not (a[2] < b[3] and b[2] < a[3]): return False
    return not REPAIR or (a[6] == b[6] and any(x < v and u < y for x, y in a[7] for u, v in b[7]))
assert same_place(('m', 'c', 0, 100, 0, 50, '+', [(0, 100)]), ('m', 'c', 50, 60, 0, 10, '-', [(50, 60)])) == (not REPAIR)  # opposite strand
assert same_place(('m', 'c', 0, 100, 0, 50, '+', [(0, 40), (90, 100)]), ('m', 'c', 50, 60, 0, 10, '+', [(50, 60)])) == (not REPAIR)  # inside an intron
assert same_place(('m', 'c', 0, 100, 0, 50, '+', [(0, 55)]), ('m', 'c', 50, 60, 0, 10, '+', [(50, 60)]))

def score(ref, mat, pat, outdir, mm, others):
    os.makedirs(outdir, exist_ok=True)
    R = records(mm); O = {n: records(p) for n, p in others}
    order = sorted(R, key=lambda q: hashlib.sha256(q.encode()).hexdigest())
    took = collections.Counter(); reads = {}
    for q in order:
        loci = loci_of(R[q]); lead = loci[0]['AS'] - loci[1]['AS'] if len(loci) > 1 else None
        st = 'S0' if lead is None else 'S1' if lead > 10 else 'S2'
        if took[st] >= CAPS[st]: continue
        took[st] += 1
        pi = next(i for i, L in enumerate(loci) if L['primary'])
        picks = {'primary': pi, 'eichler-keep': pi if st == 'S1' and pi == 0 else None, 'eichler-move': 0 if st == 'S1' else None}
        extra = []
        for n, d in O.items():
            pr = next((r for r in d.get(q, []) if r[6]), None)
            if pr is None: picks[n] = None; continue
            i = find(loci, pr)
            if i is None: extra.append(dict(chrom=pr[0], strand=pr[1], s=pr[2], e=pr[3], blocks=pr[4], AS=pr[5], primary=False)); i = -len(extra)  # a locus minimap2 has no record at
            picks[n] = i
        keep = {i for i, L in enumerate(loci[:MAXLOCI]) if L['AS'] >= 0.8 * loci[0]['AS']} | {i for i in picks.values() if i is not None and i >= 0}
        reads[q] = dict(st=st, lead=lead, loci=loci, extra=extra, picks=picks, keep=sorted(keep))
    print('reads per stratum: ' + ', '.join(f'{k} {took[k]} (of {sum(1 for q in R if ((lambda l: "S0" if len(l) < 2 else "S1" if l[0]["AS"] - l[1]["AS"] > 10 else "S2")(loci_of(R[q]))) == k)})' for k in ('S0', 'S1', 'S2')))
    # sequences: the reads themselves come from the fastq behind MM.bam (primary records carry SEQ); placements from REF
    fa = pysam.FastaFile(ref); rfa = f'{outdir}/reads.fa'; cfa = f'{outdir}/placements.fa'
    with open(rfa, 'w') as w:
        for r in pysam.AlignmentFile(mm, check_sq=False):
            if r.is_secondary or r.is_supplementary or r.query_name not in reads or r.query_sequence is None: continue
            s = r.query_sequence
            if r.is_reverse: s = s.translate(COMP)[::-1]
            w.write(f">{r.query_name}\n{s}\n")
    with open(cfa, 'w') as w:
        for q, d in reads.items():
            for i in d['keep'] + [-(j + 1) for j in range(len(d['extra']))]:
                L = d['loci'][i] if i >= 0 else d['extra'][-i - 1]
                s = ''.join(fa.fetch(L['chrom'], a, b) for a, b in L['blocks']).upper()
                if L['strand'] == '-': s = s.translate(COMP)[::-1]
                w.write(f">{q}|{i}\n{s}\n")
    H = collections.defaultdict(list)
    for hap, mmi in (('mat', mat), ('pat', pat)):
        for src in (rfa, cfa):
            for k, v in hap_hits(src, mmi, hap).items(): H[k] += v
    near_perfect = lambda h: h[4] <= max(2, 0.003 * h[5])
    T = collections.Counter(); rows = []
    for q, d in reads.items():
        hs = H.get(q, []); best = min(hs, key=lambda h: (h[4], -h[5])) if hs else None
        if not best or not near_perfect(best): status, origin = 'no near-perfect hit', set()
        else:
            org = [h for h in hs if h[4] == best[4] and h[5] >= 0.98 * best[5]]; origin = set(); unloc = set()
            for i in d['keep'] + [-(j + 1) for j in range(len(d['extra']))]:
                ch = H.get(f'{q}|{i}', []); bh = {}
                for h in ch: bh[h[0]] = min(bh.get(h[0], (10 ** 9, 0)), (h[4], h[5]))
                loc = [h for h in ch if h[4] <= bh[h[0]][0] + 2 and bh[h[0]][0] <= max(2, 0.003 * bh[h[0]][1])]
                if any(same_place(a, b) for a in loc for b in org): origin.add(i)
                if not loc: unloc.add(i)
            status = 'decidable' if len(origin) == 1 else 'tie' if origin else 'outside'
            if REPAIR:  # r1785 rule 2: a near-tied competitor the haplotypes cannot locate makes the read unscorable
                pi = d['picks']['primary']; comp = {i for i, L in enumerate(d['loci'][:MAXLOCI]) if L['AS'] >= d['loci'][0]['AS'] - 10} | {pi}
                if comp & unloc: status = 'competitor unlocatable'
        T[(d['st'], status)] += 1
        for arm, i in d['picks'].items():
            if status == 'decidable': T[(d['st'], arm, 'abstain' if i is None else 'right' if i in origin else 'wrong')] += 1
        if status == 'decidable' and d['st'] == 'S1' and d['picks']['primary'] != 0:  # reading 2: the best-AS locus is not the primary
            T[('moved', 'primary', d['picks']['primary'] in origin)] += 1; T[('moved', 'eichler-move', 0 in origin)] += 1
        rows.append((q, d['st'], d['lead'], status, ';'.join(f"{a}={'NA' if i is None else i}" for a, i in d['picks'].items()), ','.join(map(str, sorted(origin))),
                     ';'.join(f"{i}:{L['chrom']}:{L['s']}-{L['e']}:{L['strand']}:AS{L['AS']}" for i, L in enumerate(d['loci'][:MAXLOCI]))))
    with open(f'{outdir}/reads.tsv', 'w') as w:
        w.write('read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n')
        for r in rows: w.write('\t'.join(map(str, r)) + '\n')
    arms = ['primary', 'eichler-keep', 'eichler-move'] + [n for n, _ in others]
    print('\nstratum\treads\tno near-perfect hit\toutside\ttie\tdecidable\tcompetitor unlocatable')
    for st in ('S0', 'S1', 'S2'): print(f"{st}\t{took[st]}\t" + '\t'.join(str(T[(st, s)]) for s in ('no near-perfect hit', 'outside', 'tie', 'decidable', 'competitor unlocatable')))
    print('\nstratum\tarm\tdecidable\tright\twrong\tabstain\twrong_share_of_decided\tright_share_of_decidable')
    for st in ('S0', 'S1', 'S2'):
        for a in arms:
            r, x, ab = T[(st, a, 'right')], T[(st, a, 'wrong')], T[(st, a, 'abstain')]; n = r + x + ab
            print(f"{st}\t{a}\t{n}\t{r}\t{x}\t{ab}\t{x / (r + x):.4f}\t{r / n:.4f}" if r + x and n else f"{st}\t{a}\t{n}\t{r}\t{x}\t{ab}\t-\t-")
    m = T[('moved', 'primary', True)] + T[('moved', 'primary', False)]
    print(f"\nreading 2: decidable S1 reads whose best-AS locus is not the primary: {m}; primary right {T[('moved', 'primary', True)]}, eichler-move right {T[('moved', 'eichler-move', True)]}")

if __name__ == '__main__':
    if sys.argv[1] == 'select': select(*sys.argv[2:4])
    else: score(*sys.argv[2:7], [a.split('=', 1) for a in sys.argv[7:]])
