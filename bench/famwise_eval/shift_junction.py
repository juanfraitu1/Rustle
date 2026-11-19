"""Copy of a GTF's multi-exon transcripts with ONE junction moved: the middle intron's donor shifted by SHIFT bp (exon end and
intron start move together). Ids get a _shift suffix."""
import sys, re, collections
src, out, shift = sys.argv[1], sys.argv[2], int(sys.argv[3])
ex, at = collections.defaultdict(list), {}
for l in open(src):
    f = l.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] == 'exon':
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append([int(f[3]), int(f[4])]); at[t] = (f[0], f[6])
with open(out, 'w') as w:
    for t, e in ex.items():
        if len(e) < 2: continue
        e.sort(); k = (len(e) - 1) // 2
        e[k][1] += shift
        if e[k][1] <= e[k][0] or e[k][1] >= e[k + 1][0]: continue
        for a, b in e:
            w.write(f'{at[t][0]}\tshift\texon\t{a}\t{b}\t.\t{at[t][1]}\t.\tgene_id "{t}_shift"; transcript_id "{t}_shift";\n')
