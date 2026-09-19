#!/usr/bin/env python3
"""Drop transcripts whose intron chain is a CONTIGUOUS SUB-CHAIN of another transcript's chain on the same
strand/contig -- i.e. incomplete-splice-match fragments of a longer emitted transcript (§6p4: 5'-truncation).
Single-exon transcripts are dropped only if contained in a multi-exon transcript's span. Keeps the longest."""
import sys,re,collections
inp,out=sys.argv[1],sys.argv[2]
rows=collections.defaultdict(list); attr={}
for l in open(inp):
    if l.startswith('#'): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='exon': continue
    t=re.search(r'transcript_id "([^"]+)"',f[8]).group(1)
    rows[t].append((f[0],f[6],int(f[3])-1,int(f[4])))
chain={}; span={}
for t,ex in rows.items():
    ex.sort(key=lambda x:x[2]); ch=ex[0][0]; st=ex[0][1]
    chain[t]=(ch,st,tuple((ex[i][3],ex[i+1][2]) for i in range(len(ex)-1)))
    span[t]=(ch,st,ex[0][2],ex[-1][3])
by=collections.defaultdict(list)
for t,(ch,st,c) in chain.items(): by[(ch,st)].append(t)
drop=set()
for key,ts in by.items():
    multi=[t for t in ts if chain[t][2]]
    multi.sort(key=lambda t:-len(chain[t][2]))
    for i,a in enumerate(multi):
        if a in drop: continue
        ca=chain[a][2]
        for b in multi:
            if b==a or b in drop: continue
            cb=chain[b][2]
            if len(cb)>=len(ca): continue
            # contiguous sub-chain?
            if any(ca[k:k+len(cb)]==cb for k in range(len(ca)-len(cb)+1)): drop.add(b)
    for t in ts:
        if chain[t][2] or t in drop: continue
        s=span[t]
        if any(span[m][2]<=s[2] and s[3]<=span[m][3] for m in multi if m not in drop): drop.add(t)
kept=0
with open(out,"w") as fo:
    for l in open(inp):
        if l.startswith('#'): fo.write(l); continue
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        m=re.search(r'transcript_id "([^"]+)"',f[8])
        if m and m.group(1) in drop: continue
        fo.write(l)
print(f"{len(rows)} transcripts -> dropped {len(drop)} sub-chain fragments -> {len(rows)-len(drop)} kept")
