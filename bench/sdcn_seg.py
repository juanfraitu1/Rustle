import collections
B='/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed'
iv=collections.defaultdict(list); n=0
for l in open(B):
    if l.startswith('#'): continue
    f=l.split('\t')
    L=min(int(f[2])-int(f[1]),int(f[5])-int(f[4])); x=float(f[20])
    if L<1000 or x<0.90: continue
    n+=1; iv[f[0]].append((int(f[1]),int(f[2]))); iv[f[3]].append((int(f[4]),int(f[5])))
print('pairs kept',n)
tot=0; seg=0; segbp=0; out=open('segments.bed','w')
for c in sorted(iv):
    xs=iv[c]; ev=collections.Counter()
    for s,e in xs: ev[s]+=1; ev[e]-=1
    cov=0; prev=None
    for p in sorted(ev):
        if prev is not None and cov>0 and p-prev>=1000:
            seg+=1; segbp+=p-prev; out.write(f"{c}\t{prev}\t{p}\t{cov}\n")
        cov+=ev[p]; prev=p
    end=-1
    for s,e in sorted(xs):
        if s>end: tot+=e-s; end=e
        elif e>end: tot+=e-end; end=e
print('SD union bp',tot,'elementary segments >=1kb',seg,'bp',segbp)
