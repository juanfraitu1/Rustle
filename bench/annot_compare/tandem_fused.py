#!/usr/bin/env python3
"""Descriptive (2026-11-03, no bar): which r1495 fused loci join paralog / same-stem RefSeq genes, what blocks the F1 cut there,
and the read-level evidence at the crossing junction. usage: tandem_fused.py OUT.tsv JUNCTION_READS.tsv (run from the repo root;
famtruth env vars; r1495 arm A under /anvme/workspace/v138da13-rustle/fusedgw/A). Output of the 2026-11-03 run: ~/fused/tandem_fused.tsv."""
import sys, csv, re, collections, itertools, gzip, pysam
sys.path.insert(0,'bench/denovo_eval'); import famtruth as F
G=F.load_cat(); pairs,univ=F.truth_pairs(G,F.PARALOGS,set(F.LEVELS['primates']))
name_of={g:d['name'] for g,d in G.items()}
ppairs={frozenset((name_of[a],name_of[b])) for a,b in pairs if a in name_of and b in name_of}
idx=F.exon_index(G)
stem=lambda n: re.sub(r'(L?\d+[A-Z]?|P\d*)$','',n)
W='/anvme/workspace/v138da13-rustle/fusedgw/A/'
RS='/home/hpc/v138da/v138da13/work/jfim/rustle_data/analysis/benchmark_collapse/A119b.chm13v2.0_RefSeq.gtf.gz'
A=lambda k,a:(re.search(k+r' "([^"]+)"',a) or [None,''])[1]
gspan={}
for l in gzip.open(RS,'rt'):
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='exon': continue
    n=A('gene_name',f[8]) or A('gene_id',f[8]); s,e=int(f[3]),int(f[4])
    gspan[n]=(f[0],min(gspan[n][1],s),max(gspan[n][2],e)) if n in gspan else (f[0],s,e)
tx={}
for l in open(W+'x.gtf'):
    f=l.rstrip('\n').split('\t')
    if len(f)<9: continue
    g,t=A('gene_id',f[8]),A('transcript_id',f[8]); r=tx.setdefault(t,{'g':g,'ex':[],'reads':0,'strand':f[6],'chrom':f[0]})
    if f[2]=='transcript': r['reads']=int(A('reads',f[8]) or 0)
    elif f[2]=='exon': r['ex'].append((int(f[3]),int(f[4])))
for r in tx.values(): r['ex'].sort()
byg=collections.defaultdict(list)
for t,r in tx.items(): byg[r['g']].append(t)
def comps(ts):
    par={t:t for t in ts}
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for a,b in itertools.combinations(ts,2):
        if tx[a]['strand']==tx[b]['strand'] and any(x[0]<=y[1] and y[0]<=x[1] for x in tx[a]['ex'] for y in tx[b]['ex']): par[find(a)]=find(b)
    c=collections.defaultdict(list)
    for t in ts: c[find(t)].append(t)
    return list(c.values())
bam=pysam.AlignmentFile('/home/hpc/v138da/v138da13/work/jfim/rustle_data/reads/A119b.t2t.bam')
def junction_reads(chrom,s,e):
    names=[]
    for a in bam.fetch(chrom,s-1,e):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        pos=a.reference_start
        for op,ln in a.cigartuples:
            if op==3 and pos==s-1 and pos+ln==e: names.append(a.query_name); break
            if op in (0,2,3,7,8): pos+=ln
    return names
# CAT genes inside a locus span -> truth pairs inside
def cat_genes_in(chrom,s,e):
    out=set()
    for g,d in G.items():
        if d['chrom']==chrom and d['exons'] and d['exons'][0][0]>=s-1 and d['exons'][-1][1]<=e: out.add(g)
    return out
truth_set=set(frozenset(p) for p in pairs)
rows=[r for r in csv.DictReader(open(W+'audit.tsv'),delimiter='\t') if r['kind']=='fused']
out=open(sys.argv[1],'w'); jr=open(sys.argv[2],'w')
blk=collections.Counter(); prize=0; prize_loci=0
print('class locus n_genes n_tx kind blocker crossing_junction reads_TJ n_junction_reads genes',file=out)
for r in rows:
    gs=r['genes'].split(',')
    tp=any(frozenset(p) in ppairs for p in itertools.combinations(gs,2)); sm=any(v>=2 for v in collections.Counter(stem(g) for g in gs if not g.startswith('LOC')).values())
    cl='paralog' if tp else ('same_stem' if sm else None)
    cg=cat_genes_in(r['chrom'],int(r['start']),int(r['end'])); tpin=sum(1 for p in itertools.combinations(cg,2) if frozenset(p) in truth_set)
    if tpin: prize+=tpin; prize_loci+=1
    if cl is None: continue
    T=byg[r['locus']]; spans=[gspan[g] for g in gs if g in gspan]
    if max(len(tx[t]['ex']) for t in T)==1:
        blk[(cl,'mono_blob')]+=1
        names=[a.query_name for a in bam.fetch(r['chrom'],int(r['start']),int(r['end'])) if not(a.is_secondary or a.is_supplementary)]
        for n in names: jr.write(f"{r['locus']}\tblob\t{n}\n")
        print(cl,r['locus'],len(gs),len(T),'mono_blob','no_junction','NA',tx[T[0]]['reads'],len(names),','.join(gs[:2]),file=out); continue
    def ng(t): return sum(1 for c,s,e in spans if tx[t]['ex'][0][0]<=e and tx[t]['ex'][-1][1]>=s)
    spanning=[t for t in T if ng(t)>=2]; private=[t for t in T if ng(t)==1]
    best=None
    for t in spanning:
        ex=tx[t]['ex']
        for (a1,a2),(b1,b2) in zip(ex,ex[1:]):
            s,e=a2+1,b1-1
            if not any((ge>=s and ge<=e) or (gs_>=s and gs_<=e) for c,gs_,ge in spans): continue
            TJ=[u for u in T if any(x[1]+1==s and y[0]-1==e for x,y in zip(tx[u]['ex'],tx[u]['ex'][1:]))]
            R=[u for u in T if u not in TJ and tx[u]['strand']==tx[t]['strand']]
            up=down=straddle=0
            for comp in comps(R):
                L=any(x[0]<s for u in comp for x in tx[u]['ex']); Rt=any(x[1]>e for u in comp for x in tx[u]['ex'])
                if L and Rt: straddle+=1
                elif L: up+=1
                elif Rt: down+=1
            st='STRUCTURAL' if (straddle==0 and up>=1 and down>=1) else ('straddle' if straddle else 'no_side')
            cand=(st,s,e,sum(tx[u]['reads'] for u in TJ))
            if best is None or (st=='STRUCTURAL' and best[0]!='STRUCTURAL'): best=cand
    kind='span_only' if not private else 'with_private'
    if best is None: blk[(cl,kind,'no_crossing')]+=1; print(cl,r['locus'],len(gs),len(T),kind,'no_crossing','NA','NA','NA',','.join(gs[:2]),file=out); continue
    names=junction_reads(r['chrom'],best[1],best[2])
    for n in names: jr.write(f"{r['locus']}\t{best[1]}-{best[2]}\t{n}\n")
    blk[(cl,kind,best[0])]+=1
    print(cl,r['locus'],len(gs),len(T),kind,best[0],f"{r['chrom']}:{best[1]}-{best[2]}",best[3],len(names),','.join(gs[:2]),file=out)
print('BLOCKERS',sorted(blk.items()))
print('PRIZE: truth pairs (primates, coding, CAT) with both genes inside one fused locus:',prize,'in',prize_loci,'loci of',len(rows))
