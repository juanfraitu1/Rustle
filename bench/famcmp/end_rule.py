'''PREREG r1660 (docs/PREREG_tes_mode_2026-11-12.md): end-choice rules on full-chain primary reads vs RefSeq ends.
  end_rule.py BAM REFSEQ OURS_GTF ISOQUANT_GTF x'''
import sys,re,gzip,collections,statistics,multiprocessing,pickle
sys.path.insert(0,'/home/hpc/v138da/v138da13/Rustle/bench/famcmp')
from flair_like import best_sites
import pysam
BAM,REF,OURS,IQ,OUT=sys.argv[1:6]
def load(p):
    ex=collections.defaultdict(list)
    for l in (gzip.open if p.endswith('.gz') else open)(p,'rt'):
        f=l.split('\t')
        if len(f)>8 and f[2]=='exon': ex[(f[0],f[6],re.search(r'transcript_id "([^"]*)"',f[8]).group(1))].append((int(f[3])-1,int(f[4])))
    out=collections.defaultdict(list)
    for (c,s,_),e in ex.items():
        e.sort()
        if len(e)>1: out[(c,s,tuple((a[1],b[0]) for a,b in zip(e,e[1:])))].append((e[0][0],e[-1][1]))
    return out
ref,ours,iq=load(REF),load(OURS),load(IQ)
common=set(ref)&set(ours)  # correct chains of ours; IsoQuant scored where it has the chain
want=collections.defaultdict(dict); want=collections.defaultdict(dict)
for c,s,ch in common: want[c][ch]=s
def work(c):
    g=collections.defaultdict(list)
    for r in pysam.AlignmentFile(BAM).fetch(c):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        pos=r.reference_start; ch=[]
        for op,ln in r.cigartuples:
            if op==3: ch.append((pos,pos+ln))
            if op in(0,2,3,7,8): pos+=ln
        ch=tuple(ch)
        if ch in want[c]: g[ch].append((r.reference_start,pos))
    return c,dict(g)
if __name__=='__main__':
    with multiprocessing.Pool(32) as p: reads=dict(p.imap_unordered(work,list(want)))
    rules={'mode (ties to longer)':lambda v,lo:(lambda c:min(c,key=lambda x:(-c[x],x if lo else -x)))(collections.Counter(v)),
           'median':lambda v,lo:int(statistics.median(v)),'most extended':lambda v,lo:min(v) if lo else max(v),
           'FLAIR window peak':lambda v,lo:best_sites(collections.Counter(v),len(v),lo)[0][0],
           '10th pct inward':lambda v,lo:sorted(v)[int(0.1*len(v))] if lo else sorted(v)[-1-int(0.1*len(v))]}
    res=collections.defaultdict(lambda:{"5'":[],"3'":[]}); n=0
    for (c,s,ch) in common:
        v=reads.get(c,{}).get(ch)
        if not v: continue
        n+=1
        for name,G in (('ours (emitted)',ours),('IsoQuant (emitted)',iq)):
            if (c,s,ch) not in G: continue
            for which in ("5'","3'"):
                lo=(s=='+')==(which=="5'")
                res[name][which].append(min(abs((a if lo else b)-(ra if lo else rb)) for a,b in G[(c,s,ch)] for ra,rb in ref[(c,s,ch)]))
        for name,f in rules.items():
            for which in ("5'","3'"):
                lo=(s=='+')==(which=="5'"); x=f([a if lo else b for a,b in v],lo)
                res[name][which].append(min(abs(x-(ra if lo else rb)) for ra,rb in ref[(c,s,ch)]))
    print('correct chains of ours with full-chain reads:',n)
    print("rule\tend\tmedian |d|\twithin 10\twithin 50\twithin 100")
    for name,d in res.items():
        for which,ab in d.items(): print(name,which,statistics.median(ab),'%.3f'%(sum(x<=10 for x in ab)/len(ab)),'%.3f'%(sum(x<=50 for x in ab)/len(ab)),'%.3f'%(sum(x<=100 for x in ab)/len(ab)),sep='\t')
