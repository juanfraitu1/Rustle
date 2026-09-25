"""Score mcl_families clusters against the gorilla referee at the PAIR level, by annotated-mRNA nucleotide identity band.
usage: referee_band_score.py clusters.tsv loci.gff3 label   (paths relative to /mnt/linuxdisk/tmp/gw22/sec)
truth pairs = referee same-family gene pairs, both genes expressed (>=2 primary reads); band = best minimap2 identity
between the two genes' ANNOTATED mRNAs (independent of our output); 'none' = no alignment record at all."""
import sys, re, collections, itertools
S='/mnt/linuxdisk/tmp/gw22/sec/'
clusters, gff3, label = sys.argv[1:4]
fam={}
for ln in open(S+'ref/NC_073244.2.tsv'):
    f=ln.rstrip('\n').split('\t')
    if f[0]!='Gene Name' and len(f)>1: fam[f[0]]=f[1]
expr=set(l.split('\t')[0] for l in open(S+'ref/NC_073244.2.expressed.tsv') if not l.startswith('Gene'))
ident={}
for ln in open(S+'ref/referee_mrna.paf'):
    f=ln.split('\t')
    if f[0]==f[5]: continue
    k=frozenset((f[0],f[5])); idn=int(f[9])/max(1,int(f[10]))
    if k not in ident or idn>ident[k]: ident[k]=idn
BANDS=[(0.90,1.01,'>=90'),(0.80,0.90,'80-90'),(0.70,0.80,'70-80'),(0.60,0.70,'60-70'),(0.0,0.60,'<60')]
def band(k):
    if k not in ident: return 'none'
    for lo,hi,n in BANDS:
        if lo<=ident[k]<hi: return n
truth={}
for fid in set(fam.values()):
    gs=sorted(g for g,v in fam.items() if v==fid and g in expr)
    for p in itertools.combinations(gs,2): truth[frozenset(p)]=band(frozenset(p))
genes={}
for ln in open(S+'ref/NC_073244.2.genes.gff'):
    f=ln.rstrip('\n').split('\t')
    if len(f)<9 or f[0]!='NC_073244.2' or f[2] not in ('gene','pseudogene'): continue
    m=re.search(r'(?:^|;)Name=([^;]+)',f[8])
    if m: genes[m.group(1)]=(int(f[3])-1,int(f[4]))
gl=sorted((s,e,g) for g,(s,e) in genes.items())
def gene_of(s,e):
    best=None
    for gs,ge,g in gl:
        if ge<=s: continue
        if gs>=e: break
        o=min(e,ge)-max(s,gs)
        if o>0 and (best is None or o>best[0]): best=(o,g)
    return best[1] if best else None
cg=collections.defaultdict(set)
for ln in open(S+clusters):
    f=ln.rstrip('\n').split('\t')
    if f[0]=='cluster_id': continue
    g=gene_of(int(f[6]),int(f[7]))
    if g: cg[f[0]].add(g)
pairs={frozenset(p) for gs in cg.values() for p in itertools.combinations(sorted(gs),2)}
largest=max((len(v) for v in cg.values()),default=0)
by=collections.defaultdict(lambda:[0,0])
for k,b in truth.items():
    by[b][1]+=1; by[b][0]+=(k in pairs)
# precision against the COMPLETE referee (expression is irrelevant to whether two genes are one family)
jd={k for k in pairs if all(g in fam for g in k)}; tp=sum(1 for k in jd if len({fam[g] for g in k})==1)
big=max(cg.items(),key=lambda x:len(x[1]))[1]
comp=collections.Counter(fam.get(g,'not-in-referee') for g in big)
order=['>=90','80-90','70-80','60-70','<60','none']
rec=' · '.join(f"{b} {by[b][0]}/{by[b][1]}" for b in order if b in by)
print(f"{label:16s} largest {largest:3d} genes {dict(comp.most_common(3))} | prec {tp}/{len(jd)}={tp/len(jd) if jd else 0:.3f} | recall by annotated-mRNA identity: {rec}")
