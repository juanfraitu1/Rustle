#!/usr/bin/env python3
"""Parse /tmp/daz1_vs_daz3.sam (DAZ3rc aligned to DAZ1 span, cs=long) to:
  - locate diagnostic SNVs (genomic DAZ1 coords)
  - locate the copy-specific structural block(s) (deletions >= 30 bp in DAZ1
    that are absent in DAZ3)
  - build a piecewise DAZ1->DAZ3rc-query coordinate map over the MATCHED
    (shared) backbone, so a DAZ1 junction can be tested for existence in DAZ3.
Writes /tmp/copymap.json
"""
import re, json

DAZ1_OFF = 42783133            # genomic start of the DAZ1 contig (1-based)
DAZ1 = (42783133, 42859657); DAZ3 = (42879918, 42945552)

rec = open('/tmp/daz1_vs_daz3.sam').readline().rstrip('\n').split('\t')
cs = next(t[5:] for t in rec[11:] if t.startswith('cs:Z:'))

refpos = int(rec[3])           # 1-based position in DAZ1 contig (local)
qpos   = 1                     # 1-based position in DAZ3rc query (local)

diag = []          # genomic DAZ1 coords of substitutions
struct = []        # (genomic_start, genomic_end) of DAZ1-only deletions >=30bp
shared_segs = []   # list of (daz1_local_start, daz1_local_end, delta) where
                   # delta = qpos - refpos so daz3_local = daz1_local + delta
                   # over a contiguous matched run (exclusive end, 1-based starts)

def push_match(rlen):
    """record a shared backbone segment of length rlen starting at current
    refpos/qpos; delta constant within it."""
    global refpos, qpos
    delta = qpos - refpos
    shared_segs.append((refpos, refpos + rlen, delta))
    refpos += rlen; qpos += rlen

for op, val in re.findall(r'([=:*+\-~])([A-Za-z0-9]+)', cs):
    if op == '=':                       # matched run (bases listed)
        push_match(len(val))
    elif op == ':':                     # matched run (count)
        push_match(int(val))
    elif op == '*':                     # substitution: 1 ref + 1 qry base
        diag.append(DAZ1_OFF + refpos - 1)
        # treat as a 1bp shared (aligned) column with same delta
        delta = qpos - refpos
        shared_segs.append((refpos, refpos + 1, delta))
        refpos += 1; qpos += 1
    elif op == '-':                     # deletion: ref bases present, qry absent
        n = len(val)
        if n >= 30:
            struct.append((DAZ1_OFF + refpos - 1, DAZ1_OFF + refpos - 1 + n))
        refpos += n                     # qpos unchanged (DAZ1-only)
    elif op == '+':                     # insertion: qry bases, no ref
        qpos += len(val)
    elif op == '~':                     # intron (none expected in genomic-genomic)
        n = int(re.search(r'\d+', val).group())
        refpos += n

diag.sort()
struct.sort()
out = dict(daz1_off=DAZ1_OFF, diag=diag, struct=struct, shared_segs=shared_segs,
          ref_consumed=refpos-int(rec[3]), qry_consumed=qpos-1)
json.dump(out, open('/tmp/copymap.json', 'w'))

print(f"diagnostic SNVs: {len(diag)}")
print(f"  genomic SNV coords: {diag}")
print(f"structural DAZ1-only blocks (>=30bp): {struct}")
for s, e in struct:
    print(f"  block {s}-{e}  len={e-s} bp")
print(f"shared matched segments: {len(shared_segs)}  total matched ref bp = "
      f"{sum(e-s for s,e,_ in shared_segs)}")
print(f"ref bp consumed = {out['ref_consumed']}  (DAZ1 span = {DAZ1[1]-DAZ1[0]})")
