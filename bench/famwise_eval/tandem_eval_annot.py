"""Member size on REAL arrays against the curated annotation, protein-coding members only (their annotated full-length
transcripts are real gene models). A119b testis, every array of row 1820 with >= 2 protein-coding units. Per member with
>= 5 sense reads: the rule (tandem.member_sizes), StringTie and Rustle primaries (the model of the member = the same-strand model
with the largest exonic overlap with the member's annotated span), each against the member's longest annotated transcript span."""
import sys, csv, collections, statistics
sys.path.insert(0,'/storage/home/jxi21/Rustle/bench/famwise_eval')
import pysam
from tandem import member_sizes
from undecidable import read_gtf_tx, load_annot
W='/storage/group/kdm16/default/jxi21/fritz/ws'
tx,genes=load_annot(f'{W}/r1797_famsim_mimic/annot/chm13v2.0_CAT_Liftoff.slim.gff3.gz')
span=collections.defaultdict(int)
for t,v in tx.items():
    span[v['gene']]=max(span[v['gene']], v['exons'][-1][1]-v['exons'][0][0])
arr=collections.defaultdict(list)
for r in csv.DictReader(open(f'{W}/r1820_undecidable_copies/arrays/units.tsv'),delimiter='\t'):
    if r['group']=='protein_coding':
        arr[r['array']].append((int(r['start']),int(r['end']),r['strand'],r['unit'],r['gene']))
bam=pysam.AlignmentFile(f'{W}/r1820_undecidable_copies/arrays/P.bam')
G={'stringtie':read_gtf_tx('st_P.gtf'),'rustle_P':read_gtf_tx(f'{W}/r1820_undecidable_copies/arrays/asm_P.gtf')}
TANDEM={'CT47A','CT45A','GAGE','MAGEA','TSPY','AMY','SSX','PRAMEF','SPATA31D','DEFA','XAGE','SPANX','CT45','RBMY','DAZ','TBC1D'}
detail=[]
res=collections.defaultdict(list); fam=collections.defaultdict(lambda: collections.defaultdict(list))
for an,us in arr.items():
    us.sort()
    if len(us)<2: continue
    chrom=an.split('@')[1].split(':')[0]; st=collections.Counter(u[2] for u in us).most_common(1)[0][0]
    us=[u for u in us if u[2]==st]
    if len(us)<2: continue
    units=[(u[0],u[1]) for u in us]
    P=statistics.median(b[0]-a[0] for a,b in zip(units,units[1:]))
    s0,e0=units[0][0]-P,units[-1][1]+P
    reads=[(r.reference_start,r.reference_end,'-' if r.is_reverse else '+',any(o==3 for o,_ in r.cigartuples)) for r in bam.fetch(chrom,max(0,int(s0)),int(e0)) if not (r.is_secondary or r.is_supplementary or r.is_unmapped)]
    est,_=member_sizes(reads,units,P,st)
    for (a,b,_,gid,name),m in zip(us,est):
        if m['source']=='none': continue
        truth=span[gid] or (b-a)
        row={'rule':m['size']}
        for k,g in G.items():
            best=None
            for v in g.values():
                if v['chrom']!=chrom or v['strand']!=st: continue
                ov=sum(max(0,min(y,b)-max(x,a)) for x,y in v['exons'])
                if ov>0 and (best is None or ov>best[0]): best=(ov,v['exons'][-1][1]-v['exons'][0][0])
            row[k]=best[1] if best else 0
        f=an.split('@')[0]
        if f in TANDEM: detail.append((f,name,truth,row['rule'],row['stringtie'],row['rustle_P'],m['reads']))
        for k,v in row.items():
            res[k].append(abs(v-truth)/truth); fam[f][k].append(abs(v-truth)/truth)
n=len(res['rule'])
print(f"protein-coding array members with >= 5 sense reads: {n} in {len(fam)} families")
for k in ('rule','stringtie','rustle_P'):
    e=res[k]; print(f"  {k}: within 10% {sum(x<=0.10 for x in e)}/{n}, within 20% {sum(x<=0.20 for x in e)}/{n}, median error {statistics.median(e):.3f}, > 2x off {sum(x>1 for x in e)}")
print('per family (members; rule / stringtie / rustle_P within 20%):')
for f,d in sorted(fam.items(), key=lambda kv: -len(kv[1]['rule'])):
    m=len(d['rule']); print(f"  {f}: {m}; " + " / ".join(str(sum(x<=0.2 for x in d[k])) for k in ('rule','stringtie','rustle_P')))

print('tandem-family members: family gene | annotated | rule | stringtie | rustle_P | reads | rule vs annotated')
sg=collections.Counter()
for f,g,t,r,st,rp,n in sorted(detail):
    d='within 20%' if abs(r-t)<=0.2*t else 'shorter' if r<t else 'longer'; sg[d]+=1
    print(f'  {f} {g} | {t} | {r} | {st} | {rp} | {n} | {d}')
print('rule vs annotation in tandem families:', dict(sg))
for k,i in (('rule',3),('stringtie',4),('rustle_P',5)):
    print(f"  {k}: tandem-family members within 20% of annotation {sum(abs(d[i]-d[2])<=0.2*d[2] for d in detail)}/{len(detail)}")
