import collections
SEXREF={'NC_073247.2','NC_073248.2','NC_011120.1'}
placed={x for l in open('/mnt/linuxdisk/tmp/p18/chrom_map.tsv') for x in l.rstrip('\n').split('\t')[1:]}
src={l.split()[0]:l.split()[1] for l in open('src.tsv')}
def sex(g,t): return t in SEXREF if g=='ref' else t.split('_')[0] in ('chrX','chrY','chrM')
FILES={'EX':{g:f'ex_sr_{g}.paf' for g in ('ref','mat','pat')},'SD':{g:f'{g}.paf' for g in ('ref','mat','pat')}}
copies=collections.defaultdict(list)
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv'):
    if l.startswith('family_id'): continue
    f=l.split('\t'); copies[f[3]].append((int(f[4]),int(f[5])))
def loci(xs):
    iv=collections.defaultdict(list)
    for t,s,e in xs: iv[t].append((s,e))
    k=collections.Counter()
    for t,v in iv.items():
        end=-1
        for s,e in sorted(v):
            if s>=end: k[t]+=1
            end=max(end,e)
    return k
S=[(C,T) for C in (0.5,0.8,0.95) for T in (0.90,0.95,0.98)]
for kind in ('EX','SD'):
    H={g:collections.defaultdict(list) for g in FILES[kind]}
    for g,p in FILES[kind].items():
        for l in open(p):
            f=l.split('\t')
            if f[0].startswith(kind): H[g][f[0]].append(((int(f[3])-int(f[2]))/int(f[1]), int(f[9])/int(f[10]), f[5], int(f[7]), int(f[8])))
    U=sorted(set(H['ref'])|set(H['mat'])|set(H['pat']))
    U=[u for u in U if not any(len(H[g][u])>=400 for g in H)]
    print(f'\n== {kind}: {len(U)} units')
    for placed_only in (False,True):
        stable=None
        for s in S:
            n={}
            for u in U:
                n[u]=[]
                for g in ('ref','mat','pat'):
                    k=loci([(t,a,b) for cov,idn,t,a,b in H[g][u] if cov>=s[0] and idn>=s[1] and not sex(g,t) and (g=='ref' or not placed_only or t in placed)])
                    n[u].append(sum(k.values()))
                if s==(0.8,0.90):
                    hd=n
            d={u for u in U if n[u][1]!=n[u][2]}; stable=d if stable is None else stable&d
            if s==(0.8,0.90):
                line=f'  {"placed contigs only" if placed_only else "all contigs"} C=0.8 T=0.90: mat!=pat {len(d)}/{len(U)} ({100*len(d)/len(U):.1f}%); ref!=mat {100*sum(n[u][0]!=n[u][1] for u in U)/len(U):.1f}%, ref!=pat {100*sum(n[u][0]!=n[u][2] for u in U)/len(U):.1f}%'
                if kind=='EX':
                    fams={u.split('|')[1] for u in U}; df={u.split('|')[1] for u in d}
                    line+=f'; families with >=1 differing exon {len(df)}/{len(fams)} ({100*len(df)/len(fams):.1f}%)'
                print(line)
                print('    |mat-pat|:',dict(sorted(collections.Counter(min(abs(n[u][1]-n[u][2]),5) for u in d).items())),
                      ' ref copies of differing units:',dict(sorted(collections.Counter(min(n[u][0],6) for u in d).items())),
                      ' all units:',dict(sorted(collections.Counter(min(n[u][0],6) for u in U).items())))
                if kind=='SD':
                    def oncat(u):
                        c,x=u.split('|')[1].split(':'); a,b=map(int,x.split('-'))
                        return any(s0<b and e0>a for s0,e0 in copies[c])
                    g1=[u for u in U if oncat(u)]; g0=[u for u in U if not oncat(u)]
                    print(f'    overlapping a catalog copy {sum(u in d for u in g1)}/{len(g1)}, not {sum(u in d for u in g0)}/{len(g0)};'
                          f' low-copy (ref 1-2) {sum(u in d for u in U if n[u][0]<=2)}/{sum(1 for u in U if n[u][0]<=2)}')
        print(f'    differ under all 9 settings: {len(stable)} ({100*len(stable)/len(U):.1f}%)')
