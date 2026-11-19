"""Label an assembly's multi-exon transcripts TRUE (intron chain equals a truth transcript's) or FALSE, and cross with flags."""
import sys, re, collections, csv
truth, prof = sys.argv[1], sys.argv[2]
def chains(p):
    ex, at = collections.defaultdict(list), {}
    for l in open(p):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append((int(f[3]), int(f[4]))); at[t] = (f[0], f[6])
    return {t: (at[t][0], at[t][1], tuple((a[1], b[0]) for a, b in zip(sorted(e), sorted(e)[1:]))) for t, e in ex.items()}
T = set(chains(truth).values()); A = chains(sys.argv[3])
c = collections.Counter()
for r in csv.DictReader(open(prof), delimiter='\t'):
    if int(r['exons']) < 2: continue
    lab = 'TRUE' if A[r['transcript']] in T else 'FALSE'
    c[(lab, 'full_lt10')] += int(r['full_chain']) < 10
    c[(lab, 'n')] += 1; c[(lab, 'flagged')] += r['flags'] != '-'
    for f in r['flags'].split(','):
        if f != '-': c[(lab, f)] += 1
for lab in ('TRUE', 'FALSE'):
    n = c[(lab, 'n')]
    import sys as _s
    names = _s.argv[4:6] if len(_s.argv) > 5 else ('TRUE', 'FALSE')
    print(f"{dict(zip(('TRUE','FALSE'), names))[lab]}: {n} transcripts, < 10 full-chain reads {c[(lab, 'full_lt10')]}, flagged {c[(lab, 'flagged')]} ({c[(lab, 'flagged')] / max(1, n):.1%}); " + ", ".join(
        f"{f} {c[(lab, f)]}" for f in ('WEAK_JUNCTION', 'SHIFTED_DOMINANT', 'THIN_EXACT', 'NEAR_MISS_DOMINANT', 'UNCOVERED', 'NO_FULL_READ')))
