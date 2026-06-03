#!/usr/bin/env python3
"""Decide whether read 50660416 carries TRUE copy-specific (DAZ1-only) sequence
or whether minimap2 merely parked query bases in the block on the secondary.

Method: walk both CIGARs with query coordinates. For the 3 'in-block' DAZ1 exons
(query intervals q_a..q_b), find where those SAME query bases land in the DAZ3
(primary) placement. If in DAZ3 they map to SHARED backbone (a real DAZ3 exon),
the bases are NOT copy-specific and the DAZ1 in-block placement is an artifact.
If in DAZ3 those query bases are SOFT-CLIPPED or unaligned, the read truly has
DAZ1-only content -> real separating signal.
"""
import re,json,bisect,subprocess
DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
cm=json.load(open('/tmp/copymap.json')); DAZ1_OFF=cm['daz1_off']
struct=[tuple(x) for x in cm['struct']]; S0,S1=struct[0]
QN='m64076_221110_210557/50660416/ccs'
out=subprocess.run(['grep',QN,'/tmp/daz_aln.sam'],capture_output=True,text=True).stdout
def copy_of(p,e):
    m=(p+e)//2
    return 'DAZ1' if DAZ1[0]<=m<=DAZ1[1] else 'DAZ3' if DAZ3[0]<=m<=DAZ3[1] else None
def ref_len(c): return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',c) if o in 'MDN=X')

# pick best per copy + grab flag for strand
best={}
for line in out.splitlines():
    f=line.split('\t');flag=int(f[1]);p=int(f[3]);cig=f[5]
    cp=copy_of(p,p+ref_len(cig))
    if cp is None: continue
    AS=int(next(t[5:] for t in f[11:] if t.startswith('AS:i:')))
    if cp not in best or AS>best[cp][0]: best[cp]=(AS,p,cig,flag)

def qmap(pos,cig,flag,readlen=None):
    """yield (qstart,qend,rstart,rend,op) per CIGAR op in QUERY orientation as
    stored (we use the query coordinate along the FORWARD read). For a reverse
    alignment minimap2 stores CIGAR in ref-forward; soft-clip lengths tell us
    the query offset. We return query coords in the read's own (reverse for -)
    SAM orientation, which is consistent within one record."""
    q=0; r=pos; res=[]
    toks=re.findall(r'(\d+)([MIDNSH=X])',cig)
    for n,o in toks:
        n=int(n)
        if o in 'M=X':
            res.append((q,q+n,r,r+n,'M')); q+=n; r+=n
        elif o=='I': res.append((q,q+n,r,r,'I')); q+=n
        elif o=='S': res.append((q,q+n,r,r,'S')); q+=n
        elif o=='H': res.append((q,q,r,r,'H'))
        elif o=='D': res.append((q,q,r,r+n,'D')); r+=n
        elif o=='N': res.append((q,q,r,r+n,'N')); r+=n
    return res,q  # total query length consumed (incl soft clips)

# The two records may be in opposite SAM orientation (DAZ1 - vs DAZ3 +).
# To compare query bases we need a common coordinate = position along the
# ORIGINAL read. SAM stores query as-aligned; for a reverse record the stored
# query is revcomp, so stored-q maps to original via (Ltot-1 - q). We normalize
# both to "original forward read" coordinate.
recs={}
for cp,(AS,p,cig,flag) in best.items():
    seg,Lstored=qmap(p,cig,flag)
    rev=bool(flag&0x10)
    recs[cp]=dict(p=p,cig=cig,flag=flag,rev=rev,seg=seg,Ltot=Lstored)

# get total read length from a primary/hard-clip-free record if possible
L=max(r['Ltot'] for r in recs.values())
def to_orig(cp,qs,qe):
    r=recs[cp]
    if r['rev']: return (L-qe, L-qs)
    return (qs,qe)

# in-block DAZ1 exons -> their query (original) ranges
print("=== read 50660416: in-block DAZ1 exons, traced to DAZ3 placement ===")
print(f"read length (stored) = {L}")
d1=recs['DAZ1']; d3=recs['DAZ3']
inblk_refs=[(42829829,42829881),(42830441,42830533),(42830807,42830954)]
# build DAZ3 ref-coverage by query (original) interval
d3_cov=[]
for qs,qe,rs,re_,op in d3['seg']:
    if op=='M':
        oqs,oqe=to_orig('DAZ3',qs,qe)
        d3_cov.append((oqs,oqe,rs,re_))
d3_clip=[]
for qs,qe,rs,re_,op in d3['seg']:
    if op in ('S','H'):
        oqs,oqe=to_orig('DAZ3',qs,qe); d3_clip.append((oqs,oqe))

for (rs,re_) in inblk_refs:
    # find query (orig) range of this DAZ1 exon
    qranges=[]
    for qs,qe,a,b,op in d1['seg']:
        if op=='M' and not (b<=rs or a>=re_):
            # clip to exon
            la=max(a,rs); lb=min(b,re_)
            qa=qs+(la-a); qb=qs+(lb-a)
            qranges.append(to_orig('DAZ1',qa,qb))
    print(f"\n  DAZ1 in-block exon {rs}-{re_}  -> orig-read query {qranges}")
    for oqs,oqe in qranges:
        # where do these query bases go in DAZ3?
        mapped=[(max(oqs,a),min(oqe,b),c,d) for a,b,c,d in d3_cov if not(b<=oqs or a>=oqe)]
        clipped=[(max(oqs,a),min(oqe,b)) for a,b in d3_clip if not(b<=oqs or a>=oqe)]
        mbp=sum(e-s for s,e,_,_ in mapped); cbp=sum(e-s for s,e in clipped)
        print(f"    query {oqs}-{oqe} ({oqe-oqs}bp): in DAZ3 -> aligned {mbp}bp "
              f"to ref {[(c,d) for _,_,c,d in mapped]}, soft/hard-clipped {cbp}bp")
        if mapped:
            for s,e,c,d in mapped:
                inb=' (DAZ3 ref also in a copyblock? no DAZ3 has no block)'
                print(f"        -> DAZ3 ref {c}-{d}  (shared backbone exon)")
print("\nINTERPRETATION:")
print("  If those query bases ALIGN in DAZ3 to shared backbone -> the read does")
print("  NOT carry DAZ1-only sequence; the in-block DAZ1 exons are minimap2")
print("  parking the same bases at a paralogous offset = ARTIFACT, not signal.")
print("  If those query bases are CLIPPED/unaligned in DAZ3 -> real DAZ1-only")
print("  content = genuine separating signal.")
