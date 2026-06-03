#!/usr/bin/env python3
"""CHANNEL = splice-graph contribution, over the 40 tied DAZ reads.

Hypothesis: a tied read contributes EQUALLY to both copies' shared splice-graph
backbone -> structure is recoverable regardless of copy assignment; only
abundance is ambiguous.

For each tied read we have its BEST-AS placement in DAZ1 (-) and in DAZ3 (+).
A read's "splice-graph contribution" to a copy = the set of introns (junctions)
its CIGAR carves into that copy + the aligned exon blocks.

We test junction identity ACROSS copies two independent ways:

  (A) DIRECT: does the read's OWN DAZ3 placement carve the SAME junctions
      (after mapping DAZ1->DAZ3 through the copy-copy alignment) as its DAZ1
      placement?  i.e. is each DAZ1 junction matched by a DAZ3 junction of the
      same read.

  (B) COPY-IDENTITY: independent of the read, is each DAZ1 junction a position
      that is colinear/shared in the copy-copy alignment (both donor and
      acceptor land in matched backbone, with the SAME delta) -> the junction
      exists structurally in BOTH copies.

Then we tally:
  (1) total junction-support (read x junction incidences) and % shared/identical.
  (2) total aligned bp and % on shared (non-copy-specific) sequence.
  (3) bp overlapping the copy-specific 10.8 kb block (must be 0).
  (4) one representative read's exon chain + per-exon/junction existence in both.
"""
import re, json, collections

DAZ1 = (42783133, 42859657); DAZ3 = (42879918, 42945552)
cm = json.load(open('/tmp/copymap.json'))
DAZ1_OFF = cm['daz1_off']
struct = [tuple(x) for x in cm['struct']]          # genomic DAZ1-only blocks
shared_segs = cm['shared_segs']                     # (l1_start,l1_end,delta) local DAZ1
# DAZ3rc query is the reverse complement of the DAZ3 plus-strand span.
DAZ3_LEN = DAZ3[1] - DAZ3[0]                         # span length used for faidx

TIED = set(l.strip() for l in open('/dev/stdin')) if False else None

def copy_of(pos, endp):
    mid = (pos + endp) // 2
    if DAZ1[0] <= mid <= DAZ1[1]: return 'DAZ1'
    if DAZ3[0] <= mid <= DAZ3[1]: return 'DAZ3'
    return None

def ref_len(cig):
    return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',cig) if o in 'MDN=X')

def exon_blocks(pos, cig):
    """genomic aligned exon intervals [start,end) ; M/=/X consume both, D small
    deletions stay within an exon, N splits exons."""
    out=[]; cur=pos; seg_start=pos
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',cig):
        n=int(n)
        if o in 'M=X': cur+=n
        elif o=='D':   cur+=n        # small del -> still same exon
        elif o=='N':                  # intron -> close current exon
            out.append((seg_start,cur)); cur+=n; seg_start=cur
        # I/S/H consume no ref
    out.append((seg_start,cur))
    return out

def introns(pos, cig):
    """genomic intron intervals [donor,acceptor) from N ops (the junctions)."""
    out=[]; cur=pos
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',cig):
        n=int(n)
        if o in 'M=XD': cur+=n
        elif o=='N': out.append((cur,cur+n)); cur+=n
    return out

# ----- copy-copy coordinate map (local DAZ1) -> delta lookup -----
import bisect
seg_starts=[s for s,_,_ in shared_segs]
def daz1_local_delta(local):
    """return delta (q-r) if local DAZ1 pos is inside a matched backbone seg,
    else None (it's in the structural block / unaligned)."""
    i=bisect.bisect_right(seg_starts, local)-1
    if i<0: return None
    s,e,d=shared_segs[i]
    return d if s<=local<e else None

def in_struct(g_pos):
    return any(s<=g_pos<e for s,e in struct)

def bp_in_struct(blocks):
    tot=0
    for s,e in blocks:
        for ss,ee in struct:
            lo=max(s,ss); hi=min(e,ee)
            if lo<hi: tot+=hi-lo
    return tot

# ----- load tied set + per-read best placements (BAM tags) -----
TIED=set(l.strip() for l in open('/tmp/tied40.txt'))
reads=collections.defaultdict(dict)
for line in open('/tmp/daz_aln.sam'):
    f=line.rstrip('\n').split('\t'); cig=f[5]
    if cig=='*': continue
    qn=f[0]
    if qn not in TIED: continue
    pos=int(f[3]); cp=copy_of(pos,pos+ref_len(cig))
    if cp is None: continue
    AS=next((int(t[5:]) for t in f[11:] if t.startswith('AS:i:')),-1)
    cur=reads[qn].get(cp)
    if cur is None or AS>cur['AS']:
        reads[qn][cp]=dict(AS=AS,pos=pos,cig=cig,
                           strand='-' if int(f[1])&0x10 else '+')

assert all('DAZ1' in c and 'DAZ3' in c for c in reads.values())
print(f"tied reads loaded with both placements: {len(reads)}")

# ====== (B) JUNCTION COPY-IDENTITY via colinear backbone ======
# A DAZ1 junction (genomic donor g_d, acceptor g_a) is "shared/exists in DAZ3"
# if BOTH endpoints lie in matched backbone with the SAME delta AND no
# structural block lies strictly between them being skipped differently.
def daz1_junction_shared(g_d, g_a):
    ld = g_d - DAZ1_OFF + 1            # local DAZ1 1-based (donor = exon end+1 boundary)
    la = g_a - DAZ1_OFF + 1
    # use the bases just inside each exon next to the junction for delta lookup
    dd = daz1_local_delta(ld-1)        # last base of upstream exon
    da = daz1_local_delta(la)          # first base of downstream exon
    if dd is None or da is None: return False, dd, da
    return (dd==da), dd, da

junc_incidence=0; junc_shared=0
per_read_junc=[]
all_daz1_junctions=collections.Counter()
for qn,c in reads.items():
    j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
    sh=0
    for (gd,ga) in j1:
        all_daz1_junctions[(gd,ga)]+=1
        ok,_,_=daz1_junction_shared(gd,ga)
        junc_incidence+=1; sh+=ok; junc_shared+=ok
    per_read_junc.append((qn,len(j1),sh))

print("\n=== (1) JUNCTION SUPPORT over 40 tied reads (channel B: copy-identity) ===")
print(f"total junction-support (read x junction incidences) = {junc_incidence}")
print(f"  of which on junctions identical/shared between DAZ1 & DAZ3 = {junc_shared}"
      f"  ({100*junc_shared/max(1,junc_incidence):.1f}%)")
print(f"distinct DAZ1 junctions used by tied reads = {len(all_daz1_junctions)}")
ds=sum(1 for j in all_daz1_junctions if daz1_junction_shared(*j)[0])
print(f"  distinct junctions that are shared/colinear in both copies = {ds}"
      f"  ({100*ds/max(1,len(all_daz1_junctions)):.1f}%)")

# ====== (A) DIRECT cross-placement junction match (same read in both copies) ======
# Map each DAZ1 junction endpoint to DAZ3rc-query local, then to DAZ3 genomic
# (plus strand), and check the read's own DAZ3 placement has a junction there.
def daz1_to_daz3_genomic(g):
    """genomic DAZ1 pos -> genomic DAZ3 plus-strand pos via colinear backbone +
    reverse-complement of the DAZ3rc query."""
    local=g-DAZ1_OFF+1
    d=daz1_local_delta(local)
    if d is None: return None
    q_local=local+d                 # 1-based in DAZ3rc query
    # DAZ3rc = revcomp of DAZ3 span; query pos q maps to span pos (LEN - q +1)
    span_pos = DAZ3_LEN - q_local + 1    # 1-based within the DAZ3 faidx span
    return DAZ3[0] + span_pos - 1        # genomic plus-strand
direct_inc=0; direct_match=0
TOL=12
for qn,c in reads.items():
    j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
    j3=set(introns(c['DAZ3']['pos'],c['DAZ3']['cig']))
    j3_ends=[]
    for (d,a) in j3: j3_ends.append((d,a))
    for (gd,ga) in j1:
        direct_inc+=1
        # map both ends to DAZ3 genomic (inverted -> donor/acceptor swap)
        m1=daz1_to_daz3_genomic(gd); m2=daz1_to_daz3_genomic(ga)
        if m1 is None or m2 is None: continue
        lo,hi=sorted((m1,m2))
        # match if read's DAZ3 placement has a junction with both ends within TOL
        for (dd,aa) in j3_ends:
            jl,jh=sorted((dd,aa))
            if abs(jl-lo)<=TOL and abs(jh-hi)<=TOL:
                direct_match+=1; break

print("\n=== (1b) DIRECT cross-placement check (same read's DAZ3 junctions) ===")
print(f"DAZ1 junction incidences whose mapped position has a matching DAZ3 "
      f"junction in the SAME read (+-{TOL}bp): {direct_match}/{direct_inc} "
      f"({100*direct_match/max(1,direct_inc):.1f}%)")

# ====== (2) ALIGNED BP and % on shared sequence; (3) structural block ======
tot_bp=0; shared_bp=0; struct_bp=0
for qn,c in reads.items():
    bl=exon_blocks(c['DAZ1']['pos'],c['DAZ1']['cig'])
    for s,e in bl:
        tot_bp += e-s
        # shared bp = within matched backbone segments (genomic DAZ1)
        for ss,ee,_ in shared_segs:
            gss=DAZ1_OFF+ss-1; gee=DAZ1_OFF+ee-1
            lo=max(s,gss); hi=min(e,gee)
            if lo<hi: shared_bp+=hi-lo
    struct_bp += bp_in_struct(bl)
print("\n=== (2)/(3) ALIGNED BP over 40 tied reads (DAZ1 placement) ===")
print(f"total aligned exon bp = {tot_bp}")
print(f"  on shared (non-copy-specific) backbone = {shared_bp}  "
      f"({100*shared_bp/max(1,tot_bp):.2f}%)")
print(f"  on copy-specific 10.8 kb block = {struct_bp}  "
      f"({100*struct_bp/max(1,tot_bp):.4f}%)")

# ====== (4) one representative read: exon chain + per-feature existence ======
# pick a spliced read with the most junctions (best illustrates the backbone)
rep=max(reads.items(), key=lambda kv: len(introns(kv[1]['DAZ1']['pos'],kv[1]['DAZ1']['cig'])))
qn,c=rep
bl=exon_blocks(c['DAZ1']['pos'],c['DAZ1']['cig'])
j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
print(f"\n=== (4) REPRESENTATIVE read: {qn} ===")
print(f"  DAZ1 placement: pos={c['DAZ1']['pos']} strand={c['DAZ1']['strand']} "
      f"AS={c['DAZ1']['AS']}; DAZ3 placement: pos={c['DAZ3']['pos']} "
      f"strand={c['DAZ3']['strand']} AS={c['DAZ3']['AS']}")
print(f"  exon chain ({len(bl)} exons), each tested for existence in BOTH copies:")
for i,(s,e) in enumerate(bl):
    in_back = all(daz1_local_delta(p-DAZ1_OFF+1) is not None
                  for p in (s, e-1))
    in_blk = bp_in_struct([(s,e)])
    print(f"    exon{i+1}: {s}-{e} len={e-s}  in_shared_backbone={in_back} "
          f"bp_in_copyblock={in_blk}")
print(f"  junctions ({len(j1)}), each tested for existence in BOTH copies:")
for i,(gd,ga) in enumerate(j1):
    ok,dd,da=daz1_junction_shared(gd,ga)
    m1=daz1_to_daz3_genomic(gd); m2=daz1_to_daz3_genomic(ga)
    daz3map = None if m1 is None or m2 is None else tuple(sorted((m1,m2)))
    print(f"    junc{i+1}: DAZ1 {gd}-{ga} (intron {ga-gd}bp)  shared_in_both={ok}"
          f"  -> DAZ3 genomic {daz3map}")
