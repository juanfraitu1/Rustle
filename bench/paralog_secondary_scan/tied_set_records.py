#!/usr/bin/env python3
"""Reproduce the EXACT 40 tied reads (|AS_DAZ1 - AS_DAZ3|<=5, best placement per
copy by AS) and dump FULL records for both placements so positional/clipping
channel analysis can run on them. Definition copied verbatim from
tied_reads_daz.py (copy_of = midpoint test; best record per copy = max AS)."""
import re, collections, json, sys

DAZ1 = (42783133, 42859657)   # - strand gene
DAZ3 = (42879918, 42945552)   # + strand gene (inverted paralog)

def copy_of(pos, endp):
    mid = (pos + endp) // 2
    if DAZ1[0] <= mid <= DAZ1[1]: return 'DAZ1'
    if DAZ3[0] <= mid <= DAZ3[1]: return 'DAZ3'
    return None

def ref_len(cig):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSH=X])', cig) if op in 'MDN=X')

def tags(fields):
    d = {}
    for t in fields[11:]:
        p = t.split(':')
        d[p[0]] = p[-1]
    return d

def clips(cig):
    """leading/trailing soft(S)+hard(H) clip lengths."""
    m = re.findall(r'(\d+)([MIDNSH=X])', cig)
    lead = trail = 0
    # leading
    for n, op in m:
        if op in 'SH': lead += int(n)
        else: break
    for n, op in reversed(m):
        if op in 'SH': trail += int(n)
        else: break
    return lead, trail

reads = collections.defaultdict(dict)   # qname -> copy -> best record (max AS)
for line in open('/tmp/daz_aln.sam'):
    f = line.rstrip('\n').split('\t')
    qn, flag, pos, cig = f[0], int(f[1]), int(f[3]), f[5]
    if cig == '*': continue
    endp = pos + ref_len(cig)
    cp = copy_of(pos, endp)
    if cp is None: continue
    tg = tags(f)
    lead, trail = clips(cig)
    rec = dict(AS=int(tg.get('AS', -1)), ms=int(tg.get('ms', -1)),
               NM=int(tg.get('NM', -1)), de=float(tg.get('de', 'nan')),
               s1=int(tg.get('s1', -1)), s2=int(tg.get('s2', -1)),
               flag=flag, pos=pos, endp=endp, cig=cig,
               strand='-' if flag & 0x10 else '+', alen=ref_len(cig),
               sec=bool(flag & 0x100), supp=bool(flag & 0x800),
               lead_clip=lead, trail_clip=trail,
               mapq=int(f[4]), seqlen=(len(f[9]) if f[9] != '*' else -1))
    cur = reads[qn].get(cp)
    if cur is None or rec['AS'] > cur['AS']:
        reads[qn][cp] = rec

both = {q: c for q, c in reads.items() if 'DAZ1' in c and 'DAZ3' in c}
tied = {q: c for q, c in both.items()
        if abs(c['DAZ1']['AS'] - c['DAZ3']['AS']) <= 5}

print(f"reads in BOTH: {len(both)}; tied (|dAS|<=5): {len(tied)}", file=sys.stderr)
out = {q: c for q, c in tied.items()}
json.dump(out, open('/tmp/daz_tied_records.json', 'w'))
print("\n".join(sorted(tied.keys())))
