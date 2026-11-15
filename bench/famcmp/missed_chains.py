#!/usr/bin/env python3
"""RefSeq intron chains StringTie or FLAIR output and ours does not: read support and what ours has there instead.

  missed_chains.py REF_GTF BAM OUT_TSV OURS_GTF STRINGTIE_GTF FLAIR_GTF
Support = records whose whole intron chain equals the chain (primary / secondary counted apart). Descriptive, no bar.
"""
import sys, re, gzip, collections, pysam
def chains(path):
    ex = collections.defaultdict(list)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    out = set()
    for (c, s, _), e in ex.items():
        e.sort()
        if len(e) > 1: out.add((c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))))  # intron = 0-based half-open
    return out
ref_gtf, bam_path, out_tsv, ours_p, st_p, fl_p = sys.argv[1:7]
ref = chains(ref_gtf); ours = chains(ours_p); st = chains(st_p); fl = chains(fl_p)
missed = (ref & (st | fl)) - ours; found = ref & ours
by_intron = collections.defaultdict(set)  # (chrom, intron) -> our chains, strand ignored (a read has none)
for c, s, ch in ours:
    for i in ch: by_intron[(c, i)].add(ch)
sub = lambda a, b: any(b[k:k + len(a)] == a for k in range(len(b) - len(a) + 1))  # a is a contiguous run of b
want = collections.defaultdict(set)
for c, s, ch in missed | found: want[c].add(ch)
prim = collections.Counter(); sec = collections.Counter(); part = collections.Counter()
near = collections.Counter(); dpart = collections.Counter(); isup = collections.Counter()  # near-exact reads, partial reads no chain of ours explains, primary records per intron
FUZZ = 10; fuzzy = collections.defaultdict(list)  # (chrom, n introns, first donor // 64) -> missed chains
for c, s, ch in missed: fuzzy[(c, len(ch), ch[0][0] // 64)].append(ch)
ours_other = {(c, ch) for c, s, ch in ours}  # strand-blind: identical chain of ours on any strand
first = collections.defaultdict(list)  # (chrom, first intron) -> missed chains, for the partial-read count
for c, s, ch in missed:
    for i in ch: first[(c, i)].append(ch)
bam = pysam.AlignmentFile(bam_path)
for a in bam:
    if a.is_unmapped or a.is_supplementary or 'N' not in a.cigarstring: continue
    c = a.reference_name; ch = []; p = a.reference_start
    for op, n in a.cigartuples:
        if op == 3: ch.append((p, p + n))
        if op in (0, 2, 3, 7, 8): p += n
    ch = tuple(ch)
    if not a.is_secondary:
        for i in ch: isup[(c, i)] += 1
    if ch in want[c]: (sec if a.is_secondary else prim)[(c, ch)] += 1
    elif not a.is_secondary:
        explained = None
        for m in {m for i in ch for m in first.get((c, i), ())}:
            if len(ch) < len(m) and sub(ch, m):
                part[(c, m)] += 1
                if explained is None: explained = any(sub(ch, o) for o in {o for i in ch for o in by_intron.get((c, i), ())})
                if not explained: dpart[(c, m)] += 1
        for k in (-1, 0, 1):  # same intron count, every junction within FUZZ bp, not exact
            for m in fuzzy.get((c, len(ch), ch[0][0] // 64 + k), ()):
                if all(abs(x[0] - y[0]) <= FUZZ and abs(x[1] - y[1]) <= FUZZ for x, y in zip(ch, m)): near[(c, m)] += 1
def rel(c, ch):
    cand = {o for i in ch for o in by_intron.get((c, i), ())}
    if not cand: return 'no model shares an intron'
    if any(len(o) > len(ch) and sub(ch, o) for o in cand): return 'ours has a longer chain containing it'
    if any(len(o) < len(ch) and sub(o, ch) for o in cand): return 'ours has a shorter run of it'
    return 'ours has another isoform there'
bin_ = lambda n: '0' if n == 0 else '1' if n == 1 else '2' if n == 2 else '3-4' if n < 5 else '5+'
tab = collections.Counter(); who_tab = collections.Counter()
with open(out_tsv, 'w') as w:
    w.write('chrom\tstrand\tn_introns\tfound_by\tprimary_full\tsecondary_full\tprimary_partial\tours\tours_other_strand\tnear_exact\tmin_intron_reads\tpartial_unexplained\tfirst_intron\n')
    for c, s, ch in sorted(missed):
        who = 'both' if (c, s, ch) in st and (c, s, ch) in fl else 'StringTie' if (c, s, ch) in st else 'FLAIR'
        r = rel(c, ch); n = prim[(c, ch)]; tab[(bin_(n), r)] += 1; who_tab[(who, bin_(n))] += 1
        w.write(f'{c}\t{s}\t{len(ch)}\t{who}\t{n}\t{sec[(c, ch)]}\t{part[(c, ch)]}\t{r}\t{int((c, ch) in ours_other)}\t{near[(c, ch)]}\t{min(isup[(c, i)] for i in ch)}\t{dpart[(c, ch)]}\t{ch[0][0]}-{ch[0][1]}\n')
print(f'ref chains {len(ref)}; ours has {len(ours)} chains, finds {len(found)} (precision {len(found) / len(ours):.3f}); missed but found by StringTie or FLAIR {len(missed)}')
B = ['0', '1', '2', '3-4', '5+']; R = ['no model shares an intron', 'ours has a longer chain containing it', 'ours has a shorter run of it', 'ours has another isoform there']
print('full-chain primary records\t' + '\t'.join(R) + '\ttotal')
for b in B: print(b, *(tab[(b, r)] for r in R), sum(tab[(b, r)] for r in R), sep='\t')
print('found by\t' + '\t'.join(B))
for k in ('StringTie', 'FLAIR', 'both'): print(k, *(who_tab[(k, b)] for b in B), sep='\t')
print('for comparison, chains ours finds, by full-chain primary records:', {b: sum(1 for c, s, ch in found if bin_(prim[(c, ch)]) == b) for b in B})
print('missed with 0 full records: have partial primary records', sum(1 for c, s, ch in missed if prim[(c, ch)] == 0 and part[(c, ch)] > 0), '; only secondary full records', sum(1 for c, s, ch in missed if prim[(c, ch)] == 0 and sec[(c, ch)] > 0))
