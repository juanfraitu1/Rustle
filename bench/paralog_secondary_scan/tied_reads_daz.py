#!/usr/bin/env python3
"""Extract DAZ multimapping reads whose two placements (DAZ1 vs DAZ3) have very
similar / identical alignment scores, and ask: is there ANY signal left that
separates them?  Uses only BAM tags (AS, ms, NM, de, s1/s2) -- all present on
secondaries too -- so no realignment needed.

Input: /tmp/daz_aln.sam  (samtools view GGO.bam NC_073248.2:42783133-42945552)
"""
import re, collections

DAZ1 = (42783133, 42859657)   # - strand
DAZ3 = (42879918, 42945552)   # + strand (inverted paralog)

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

reads = collections.defaultdict(dict)   # qname -> copy -> best record (max AS)
for line in open('/tmp/daz_aln.sam'):
    f = line.rstrip('\n').split('\t')
    qn, flag, pos, cig = f[0], int(f[1]), int(f[3]), f[5]
    if cig == '*': continue
    cp = copy_of(pos, pos + ref_len(cig))
    if cp is None: continue
    tg = tags(f)
    rec = dict(AS=int(tg.get('AS', -1)), ms=int(tg.get('ms', -1)),
               NM=int(tg.get('NM', -1)), de=float(tg.get('de', 'nan')),
               s1=int(tg.get('s1', -1)), s2=int(tg.get('s2', -1)),
               strand='-' if flag & 0x10 else '+', alen=ref_len(cig),
               sec=bool(flag & 0x100))
    cur = reads[qn].get(cp)
    if cur is None or rec['AS'] > cur['AS']:
        reads[qn][cp] = rec

# keep reads aligned to BOTH copies
both = {q: c for q, c in reads.items() if 'DAZ1' in c and 'DAZ3' in c}
print(f"reads with a placement in BOTH DAZ1 and DAZ3: {len(both)}")

def delta(a, b):  # signed: (favoring-DAZ1 means AS1>AS3 ; de1<de3)
    return a - b

rows = []
for q, c in both.items():
    a1, a3 = c['DAZ1'], c['DAZ3']
    dAS = a1['AS'] - a3['AS']
    dNM = a3['NM'] - a1['NM']          # >0 => DAZ1 fewer edits => favors DAZ1
    dde = a3['de'] - a1['de']          # >0 => DAZ1 lower divergence => favors DAZ1
    dms = a1['ms'] - a3['ms']
    rows.append(dict(q=q, dAS=dAS, dNM=dNM, dde=dde, dms=dms,
                     AS1=a1['AS'], AS3=a3['AS'], de1=a1['de'], de3=a3['de'],
                     NM1=a1['NM'], NM3=a3['NM']))

# --- how tied are these reads on AS? ---
def bucket(x):
    x = abs(x)
    return ('0 (identical)' if x == 0 else '1-5' if x <= 5 else
            '6-20' if x <= 20 else '21-100' if x <= 100 else '>100')
asb = collections.Counter(bucket(r['dAS']) for r in rows)
print("\n|AS(DAZ1) - AS(DAZ3)|  distribution:")
for k in ['0 (identical)', '1-5', '6-20', '21-100', '>100']:
    print(f"   {k:>14}: {asb.get(k,0)}")

# --- among the AS-TIED reads (|dAS|<=5), what still separates them? ---
tied = [r for r in rows if abs(r['dAS']) <= 5]
print(f"\nAS-tied reads (|dAS|<=5): {len(tied)}")
sep_de = [r for r in tied if abs(r['dde']) >= 0.001]
sep_nm = [r for r in tied if abs(r['dNM']) >= 2]
sep_any = [r for r in tied if abs(r['dde']) >= 0.001 or abs(r['dNM']) >= 2]
print(f"   still separated by de  (|dde|>=0.001): {len(sep_de)}")
print(f"   still separated by NM  (|dNM|>=2):     {len(sep_nm)}")
print(f"   separated by de OR NM:                 {len(sep_any)}  "
      f"({100*len(sep_any)/max(1,len(tied)):.0f}%)")

# --- the hardest: AS IDENTICAL ---
ident = [r for r in rows if r['dAS'] == 0]
print(f"\nAS-IDENTICAL reads (dAS==0): {len(ident)}")
ide = [r for r in ident if abs(r['dde']) >= 0.001 or abs(r['dNM']) >= 1]
print(f"   of those, de or NM still differs: {len(ide)}")

# --- minimap2's own ambiguity flag: s2 close to s1 on the primary ---
amb = [q for q, c in both.items()
       for cp in c.values() if not cp['sec'] and cp['s2'] > 0 and cp['s1'] > 0
       and cp['s2'] >= 0.9 * cp['s1']]
print(f"\nprimaries where minimap2's 2nd-best chain s2 >= 0.9*s1 (it FLAGS the tie itself): {len(set(amb))}")

# --- show a few concrete tied examples ---
print("\nexamples of AS-tied reads that a finer signal still resolves:")
print(f"   {'read':<14} {'AS1/AS3':>9} {'dAS':>4} {'de1':>7} {'de3':>7} {'NM1/NM3':>9} -> call")
for r in sorted(tied, key=lambda r: abs(r['dAS']))[:8]:
    call = 'DAZ1' if (r['dde'] > 0 or r['dNM'] > 0) else 'DAZ3' if (r['dde'] < 0 or r['dNM'] < 0) else '??'
    print(f"   {r['q'][-12:]:<14} {str(r['AS1'])+'/'+str(r['AS3']):>9} {r['dAS']:>4} "
          f"{r['de1']:>7.4f} {r['de3']:>7.4f} {str(r['NM1'])+'/'+str(r['NM3']):>9} -> {call}")
