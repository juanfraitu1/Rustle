import collections
src={l.split()[0]:l.split()[1] for l in open('src.tsv')}
SEX={'NC_073247.2','NC_073248.2','NC_011120.1'}
auto={r:s for r,s in src.items() if r not in SEX}
hap_of={s:('mat' if '_mat_' in s else 'pat') for s in auto.values()}
hits={g:collections.defaultdict(list) for g in ('ref','mat','pat')}
for g in hits:
    for l in open(f'{g}.paf'):
        f=l.split('\t')
        if (int(f[3])-int(f[2]))/int(f[1])>=0.8 and int(f[9])/int(f[10])>=0.90: hits[g][f[0]].append((f[5],int(f[7]),int(f[8])))
def cnt(xs):
    iv=collections.defaultdict(list)
    for t,s,e in xs: iv[t].append((s,e))
    out=collections.Counter()
    for t,v in iv.items():
        end=-1
        for s,e in sorted(v):
            if s>=end: out[t]+=1
            end=max(end,e)
    return out
units=[l.split('\t')[3].strip() for l in open('units.bed')]
for kind in ('EX','SD'):
    agree=dis=0; byc=collections.Counter(); tot=collections.Counter()
    for u in units:
        if not u.startswith(kind): continue
        r=cnt(hits['ref'][u]); h={g:cnt(hits[g][u]) for g in ('mat','pat')}
        R=sum(r[c] for c in auto); H=sum(h[hap_of[s]][s] for c,s in auto.items())
        if R==0 and H==0: continue
        k=min(R,6); tot[k]+=1
        if R==H: agree+=1
        else: dis+=1; byc[k]+=1
    print(kind,'same-sequence control: agree',agree,'disagree',dis,f'({100*dis/(agree+dis):.1f}%)','by ref copies',{k:f'{byc[k]}/{tot[k]}' for k in sorted(tot)})
