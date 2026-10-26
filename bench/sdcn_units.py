import re, random, collections, subprocess
SEX={'NC_073247.2','NC_073248.2','NC_011120.1'}
span={}
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv'):
    if l.startswith('family_id'): continue
    f=l.split('\t'); span[f"{f[0]}|{f[1]}"]=(f[3],int(f[4]),int(f[5]))
ex={}
for l in open('/mnt/linuxdisk/tmp/hapcn/ref.paf'):
    f=l.rstrip('\n').split('\t'); q=f[0].split('|'); key=f"{q[0]}|{q[1]}"
    c,s,e=span[key]
    if f[5]!=c or int(f[7])>=e or int(f[8])<=s: continue
    cg=[x for x in f[12:] if x.startswith('cg:Z:')][0][5:]
    p=int(f[7]); a=p
    for n,o in re.findall(r'(\d+)([MIDN=X])',cg):
        n=int(n)
        if o=='N':
            if p>a: ex.setdefault((c,a,p),q[0])
            p+=n; a=p
        elif o in 'MD=X': p+=n
    if p>a: ex.setdefault((c,a,p),q[0])
ex={k:v for k,v in ex.items() if k[2]-k[1]>=100 and k[0] not in SEX}
with open('exons.bed','w') as o:
    for (c,a,b),fam in sorted(ex.items()): o.write(f"{c}\t{a}\t{b}\tEX|{fam}|{c}:{a}-{b}\n")
print('exons',len(ex),'families',len(set(ex.values())))
# SD windows
w=[]
for l in open('segments.bed'):
    c,s,e,_=l.split('\t'); s,e=int(s),int(e)
    if c in SEX: continue
    x=s
    while x<e:
        y=min(x+5000,e)
        if y-x>=1000: w.append((c,x,y))
        x=y
rng=random.Random(20261025); smp=sorted(rng.sample(w,6000))
with open('windows.bed','w') as o:
    for c,a,b in smp: o.write(f"{c}\t{a}\t{b}\tSD|{c}:{a}-{b}\n")
print('windows total',len(w),'sampled',len(smp))
