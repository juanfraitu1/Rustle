#!/usr/bin/env python3
"""One read as the seed of a locus (register 1987; driver: RUSTLE_ONE_READ_LOCI=1 in stage `multicopy`).

    one_read_loci.py BASE.gtf CAND.gtf BAM AS_TABLE OUT.gtf [THREADS [supported]]   |   one_read_loci.py --selftest
    one_read_loci.py clean BASE.gtf FLOOR1.gtf BAM GENOME.fa OUT.gtf [THREADS [MIN_INTRONS=3 [AS_TABLE]]]   (row 1990, see clean())

CAND.gtf = multicopy_rescue.py candidates (spliced transcripts of the one-read-floor assembly where BASE has no transcript
on that strand). OUT.gtf = BASE + the candidates that
  exactly ONE primary alignment has as its intron chain (the lone read), and
  that read has no other alignment within MS_BAND of its best ms (AS_TABLE columns best_ms, second_ms): a place a read
  reaches only by a tie is mostly a relative's shadow (0.11-0.51 correct against 0.84-0.94 on simulated reads).
Kept transcript lines get one_read_locus "supported" when every intron of the chain is also in ANOTHER primary read (reads
that disagree on the full chain, so no two of them made a transcript) and one_read_locus "1" otherwise; with the word
`supported` only the first kind is written. Prints one summary line (the driver logs it).
Supported loci on a human library not used to choose the rule: all introns seen again in two replicates 0.94 (two-read
chains of the assembly 0.945), all introns annotated 0.78 (0.61), chain exact or a contiguous part of an annotated one
0.71; the others 0.33-0.48 / 0.26-0.32 / 0.25. Most supported loci are FRAGMENTS of a real transcript (exact chain 12%).
What it is worth (row 1987): simulated fly / mouse / zebrafish 0.94 / 0.89 / 0.84 correct, +8 to +17% correct chains;
a real human library 13% annotated, 20% seen again in two replicate libraries. Opt-in for that reason."""
import collections, multiprocessing, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import oneread

MS_BAND = 10
_G = {}

def read_chain(a):
    p = a.reference_start; ch = []
    for op, n in a.cigartuples:
        if op == 3: ch.append((p, p + n))
        if op in (0, 2, 3, 7, 8): p += n
    return tuple(ch)

def _scan(contig):
    import pysam
    want = _G['want'][contig]; out = collections.defaultdict(list); ins = {i for c in want for i in c}; n = collections.Counter()
    for a in pysam.AlignmentFile(_G['bam']).fetch(contig):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or 'N' not in a.cigarstring: continue
        c = read_chain(a); n.update(i for i in c if i in ins)
        if c in want and len(out[c]) < 2: out[c].append(a.query_name)
    return {(contig, c): v for c, v in out.items()}, {(contig, i): k for i, k in n.items()}

def reads_with(bam, keys, threads=4, introns=False):
    """{(contig, chain): names of up to two primary reads with exactly that chain} for the wanted (contig, chain) keys;
    with introns=True also {(contig, intron): primary reads with it} for the introns of those chains"""
    import pysam
    want = collections.defaultdict(set)
    for c, ch in keys: want[c].add(ch)
    _G.update(bam=bam, want=want); out = {}; n = {}
    contigs = [c for c in pysam.AlignmentFile(bam).references if c in want]
    with multiprocessing.get_context('fork').Pool(max(1, min(threads, len(contigs) or 1))) as pool:
        for d, k in pool.imap_unordered(_scan, contigs): out.update(d); n.update(k)
    return (out, n) if introns else out

def rivals(table, names):
    """names of `names` with a second alignment within MS_BAND of the best ms"""
    out = set()
    for l in open(table):
        if l[0] == '#': continue
        f = l.split('\t', 8)
        if f[0] in names and f[6].lstrip('-').isdigit() and f[7].lstrip('-').isdigit() and 0 <= int(f[7]) >= int(f[6]) - MS_BAND: out.add(f[0])
    return out

def keep(cand, reads, tied):
    """cand = oneread.models(CAND); -> transcript ids to keep"""
    out = []
    for tid, (c, s, ex, _) in cand.items():
        r = reads.get((c, oneread.chain(ex)), ())
        if len(r) == 1 and r[0] not in tied: out.append(tid)
    return out

def supported(c, ex, n):
    """every intron of the chain is in the lone read and at least one other primary read"""
    return all(n.get((c, i), 0) >= 2 for i in oneread.chain(ex))

def main(base, cand, bam, table, out, threads=4, only=''):
    assert only in ('', 'supported'), only
    C = oneread.models(cand); R, N = reads_with(bam, {(c, oneread.chain(ex)) for c, s, ex, _ in C.values()}, int(threads), introns=True)
    tied = rivals(table, {v[0] for v in R.values() if len(v) == 1}); K = keep(C, R, tied)
    S = {t for t in K if supported(C[t][0], C[t][2], N)}
    if only: K = [t for t in K if t in S]
    with open(out + '.tmp', 'w') as o:
        o.write(open(base).read())
        for t in K:
            for l in C[t][3]: o.write(l.rstrip('\n').rstrip(';') + f'; one_read_locus "{"supported" if t in S else 1}";\n' if l.split('\t')[2] == 'transcript' else l)
    os.replace(out + '.tmp', out)
    print(f'{len(K)} of {len(C)} one-read candidates kept as loci, {len(S)} with every intron in another read (lone read {sum(len(v) == 1 for v in R.values())}, of them with a rival within {MS_BAND} ms {len(tied)})')

def clean_ok(C, others, base):
    """others = {read chain: reads} sharing an intron with C, base = assembled chains sharing one. True when exactly one read
    has C, every other read's chain is a contiguous part of C (no conflicting read, no longer read) and no assembled chain holds C"""
    return others.get(C, 0) == 1 and all(oneread.run(r, C) for r in others) and not any(oneread.run(C, b) for b in base)

def clean(base, f1, bam, fasta, out, threads=4, min_introns=3, table=''):
    """Row 1990, rule A: OUT.gtf = BASE + the spliced chains of FLOOR1 (the one-read-floor assembly) that BASE lacks and that
      one primary read has, and no read contradicts: every read sharing an intron has the same chain or a contiguous part of it;
      no read's chain and no transcript of BASE contains the chain;
      have a longest ORF >= 900 nt that does not stop > 50 nt before the last junction (gtf_ptc), and >= MIN_INTRONS introns;
      with AS_TABLE (the driver gives it): the lone read has no other alignment within MS_BAND of its best ms. Row 1991: in
      planted families at two reads per copy this removed every clean chain on a silent copy (4 of 4) for 16 of 96 correct ones.
    Ids get the prefix R1_, transcript lines one_read_chain "clean". Held-out libraries (human, 3 mouse; 2 gorillas as a check):
    0.36-0.48 exact (gorillas 0.24-0.30) against 0.12-0.35 for two-read chains; cost <= 0.003 of exact-chain precision."""
    import pysam, re, gtf_ptc
    have = collections.defaultdict(set); bhas = collections.defaultdict(set)
    for c, s, ex, _ in oneread.models(base).values():
        if len(ex) > 1:
            ch = oneread.chain(ex); have[c].add(ch)
            for i in ch: bhas[(c, i)].add(ch)
    fa = pysam.FastaFile(fasta); cand = collections.OrderedDict(); seen = 0
    for c, s, ex, lines in oneread.models(f1).values():
        ch = oneread.chain(ex) if len(ex) > 1 else ()
        if len(ch) < int(min_introns) or s not in '+-' or ch in have[c] or (c, ch) in cand: continue
        seen += 1; o, p = gtf_ptc.call(ex, s, lambda a, b: fa.fetch(c, a, b))
        if not p and o >= oneread.MIN_ORF: cand[(c, ch)] = lines
    want = {(c, i) for c, ch in cand for i in ch}; prim = collections.Counter()
    contigs = sorted({c for c, _ in want} & set(pysam.AlignmentFile(bam).references)); oneread._SCAN.update(bam=bam, want=want)
    if contigs:
        with multiprocessing.get_context('fork').Pool(max(1, min(int(threads), len(contigs)))) as pool:
            for part in pool.imap_unordered(oneread._scan, contigs, chunksize=1): prim.update(part)
    has = collections.defaultdict(set)
    for c, ch in prim:
        for i in ch: has[(c, i)].add(ch)
    ok = [k for k in cand if clean_ok(k[1], {r: prim[(k[0], r)] for i in k[1] for r in has[(k[0], i)]}, {b for i in k[1] for b in bhas.get((k[0], i), ())})]
    nr = 0
    if table and ok:
        R = reads_with(bam, set(ok), int(threads)); tied = rivals(table, {v[0] for v in R.values()})
        nr = len(ok); ok = [k for k in ok if R.get(k) and R[k][0] not in tied]; nr -= len(ok)
    n = 0
    with open(out + '.tmp', 'w') as o:
        o.write(open(base).read())
        for k in ok:
            lines = cand[k]; n += 1
            for l in lines:
                l = re.sub(r'(gene_id|transcript_id) "', r'\1 "R1_', l)
                o.write(l.rstrip('\n').rstrip(';') + '; one_read_chain "clean";\n' if l.split('\t')[2] == 'transcript' else l)
    os.replace(out + '.tmp', out)
    print(f'{n} one-read chains added (no conflicting read, complete reading frame, >= {min_introns} introns) of {seen} the assembly lacks, {len(cand)} with the frame; {nr} left out for a rival placement of the read')

def selftest():
    A = ((10, 20), (30, 40), (50, 60)); part = ((30, 40), (50, 60))
    assert clean_ok(A, {A: 1, part: 3}, set()) and not clean_ok(A, {A: 2}, set()) and not clean_ok(A, {A: 1, ((10, 20), (30, 45)): 1}, set())
    assert not clean_ok(A, {A: 1, ((1, 5),) + A: 1}, set()) and not clean_ok(A, {A: 1}, {A + ((70, 80),)}) and clean_ok(A, {A: 1}, {part})
    C = {'a': ('c', '+', [(0, 10), (20, 30)], []), 'b': ('c', '+', [(100, 110), (120, 130)], []), 'd': ('c', '+', [(200, 210), (220, 230)], []), 'e': ('c', '+', [(300, 310), (320, 330)], [])}
    R = {('c', ((10, 20),)): ['r1'], ('c', ((110, 120),)): ['r2'], ('c', ((210, 220),)): ['r3', 'r4']}
    assert keep(C, R, {'r2'}) == ['a']          # b: the read has a rival; d: two reads; e: no primary read
    assert supported('c', C['a'][2], {('c', (10, 20)): 2}) and not supported('c', C['a'][2], {('c', (10, 20)): 1})
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.tsv') as t:
        t.write('#h\nr1\t9\t-1\t1\tc\t9\t9\t-1\t9\nr2\t9\t9\t2\tc\t9\t50\t41\t50\nr3\t9\t9\t2\tc\t9\t50\t39\t50\n'); t.flush()
        assert rivals(t.name, {'r1', 'r2', 'r3'}) == {'r2'}
    print('one_read_loci selftest ok')

if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else clean(*sys.argv[2:]) if sys.argv[1] == 'clean' else main(*sys.argv[1:])
