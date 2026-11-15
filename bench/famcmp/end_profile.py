'''Why do end rules differ between samples? (descriptive, 2026-11-12, follows r1660)
  end_profile.py BAM REFSEQ OURS_GTF OUT_TSV
Per RefSeq-correct chain of OURS with full-chain primary reads, and per end (5', 3'): reads, signed error of the most
extended end and of the mode (+ = beyond RefSeq, - = short), share of reads within 5 bp of the mode, nearest RefSeq id.'''
import sys,re,gzip,collections,multiprocessing
import pysam
BAM,REF,OURS,OUT=sys.argv[1:5]
def load(p,ids=False):
    ex=collections.defaultdict(list)
    for l in (gzip.open if p.endswith('.gz') else open)(p,'rt'):
        f=l.split('\t')
        if len(f)>8 and f[2]=='exon': ex[(f[0],f[6],re.search(r'transcript_id "([^"]*)"',f[8]).group(1))].append((int(f[3])-1,int(f[4])))
    out=collections.defaultdict(list)
    for (c,s,t),e in ex.items():
        e.sort()
        if len(e)>1: out[(c,s,tuple((a[1],b[0]) for a,b in zip(e,e[1:])))].append((e[0][0],e[-1][1],t))
    return out
ref,ours=load(REF),load(OURS); common=set(ref)&set(ours); want=collections.defaultdict(dict)
for c,s,ch in common: want[c][ch]=s
def work(c):
    g=collections.defaultdict(list)
    for r in pysam.AlignmentFile(BAM).fetch(c):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        pos=r.reference_start; ch=[]; cg=r.cigartuples
        for op,ln in cg:
            if op==3: ch.append((pos,pos+ln))
            if op in(0,2,3,7,8): pos+=ln
        ch=tuple(ch)
        if ch in want[c]: g[ch].append((r.reference_start,pos,cg[0][1] if cg[0][0]==4 else 0,cg[-1][1] if cg[-1][0]==4 else 0))
    return c,dict(g)
if __name__=='__main__':
    with multiprocessing.Pool(32) as p, open(OUT,'w') as w:
        w.write("contig\tstrand\tintrons\treads\tend\text_err\tmode_err\tmode_share5\tclip_median\trefseq\tmode_pos\n")
        for c,g in p.imap_unordered(work,list(want)):
            for ch,v in g.items():
                s=want[c][ch]
                for which in ("5'","3'"):
                    lo=(s=='+')==(which=="5'"); xs=[a if lo else b for a,b,_,_ in v]; cnt=collections.Counter(xs)
                    mode=min(cnt,key=lambda x:(-cnt[x],x if lo else -x)); ext=min(xs) if lo else max(xs)
                    sgn=lambda x,r:(r-x) if lo else (x-r)   # + = beyond RefSeq
                    ra,rb,rid=min(ref[(c,s,ch)],key=lambda t:abs(mode-(t[0] if lo else t[1]))); r=ra if lo else rb
                    cl=sorted(x[2] if lo else x[3] for x in v)[len(v)//2]
                    w.write(f"{c}\t{s}\t{len(ch)}\t{len(v)}\t{which}\t{sgn(ext,r)}\t{sgn(mode,r)}\t{sum(abs(x-mode)<=5 for x in xs)/len(xs):.3f}\t{cl}\t{rid}\t{mode}\n")
