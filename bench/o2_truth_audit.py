"""S1: audit the junction-anchored truth. For each anchored read (anchors from ORIG family's model junctions),
realign the READ (minimap2 -x splice) to every copy's core hull; a read that aligns to a hull other than its
anchor with the same intron count and NM <= NM(anchor)+2 has an anchor that is an annotation gap, not a truth."""
import csv,collections,pysam,subprocess,sys,re
orig,hulls_fam=sys.argv[1],sys.argv[2]; BAMP='/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam'; BAM=pysam.AlignmentFile(BAMP)
FA=pysam.FastaFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa')
# anchors (as in score_fixed_truth)
junc=collections.defaultdict(set); cp={}
for r in csv.DictReader(open(f'{orig}/copies.tsv'),delimiter='\t'):
    cp[r['copy_idx']]=(r['chrom'],int(r['start']),int(r['end'])); bl=[tuple(map(int,x.split('-'))) for x in r['exons'].split(',')]
    for (a,b),(c,e) in zip(bl,bl[1:]): junc[r['copy_idx']].add((r['chrom'],b,c))
owner=collections.defaultdict(set)
for i,js in junc.items():
    for j in js: owner[j].add(i)
uniqj={j for j,o in owner.items() if len(o)==1}
A0={r['read_name']:r for r in csv.DictReader(open(f'{orig}/A.assignments.tsv'),delimiter='\t')}
c0=cp['0'][0]; s0=min(v[1] for v in cp.values()); e0=max(v[2] for v in cp.values()); anchor={}; mapq={}; seq={}; injs={}
for a in BAM.fetch(c0,s0-50000,e0+50000):
    if a.flag&2308 or a.query_name not in A0: continue
    pos=a.reference_start; hits=set(); ajs=[]
    for op,l in a.cigartuples:
        if op in (0,7,8,2): pos+=l
        elif op==3:
            j=(a.reference_name,pos,pos+l); ajs.append(j)
            if j in uniqj: hits|=owner[j]
            pos+=l
    if len(hits)==1:
        anchor[a.query_name]=next(iter(hits)); mapq[a.query_name]=a.mapping_quality; seq[a.query_name]=a.query_sequence
        aj=[j for j in ajs if j in uniqj]; injs[a.query_name]=aj
# hull sequences of the HULLS family (copy_idx aligned with orig by tid)
H=list(csv.DictReader(open(f'{hulls_fam}/copies.tsv'),delimiter='\t')); tid2idx={r['tid']:r['copy_idx'] for r in csv.DictReader(open(f'{orig}/copies.tsv'),delimiter='\t')}
hull={}
for r in H:
    h=r.get('core_hull','NA')
    if h!='NA': hull[tid2idx.get(r['tid'],'?')]=tuple(map(int,h.split('-')))
with open('audit_hulls.fa','w') as f:
    for i,(c,s_,e_) in cp.items(): f.write(f">{i}\n{FA.fetch(c,s_,e_)}\n")
with open('audit_reads.fa','w') as f:
    for n,s in seq.items(): f.write(f">{n}\n{s}\n")
paf=subprocess.run(['minimap2','-x','splice','-N','60','-p','0.1','--secondary=yes','-c','-t','4','audit_hulls.fa','audit_reads.fa'],capture_output=True,text=True).stdout
hits=collections.defaultdict(list)
for l in paf.split('\n'):
    f=l.split('\t')
    if len(f)<12: continue
    nm=next((int(x[5:]) for x in f[12:] if x.startswith('NM:i:')),-1); cs=next((x for x in f[12:] if x.startswith('cg:Z:')),'')
    introns=cs.count('N') if cs else 0; qcov=(int(f[3])-int(f[2]))/int(f[1])
    hits[f[0]].append((f[5],nm,introns,round(qcov,2)))
verdict=collections.Counter(); rows=[]
for n in anchor:
    hs=hits.get(n,[]); an=[h for h in hs if h[0]==anchor[n]]
    if not an: verdict['anchor-unalignable']+=1; rows.append((n,mapq[n],anchor[n],'anchor-unalignable',None)); continue
    a=min(an,key=lambda h:h[1]); alts=[h for h in hs if h[0]!=anchor[n] and h[3]>=a[3]-0.05 and h[2]==a[2] and h[1]<=a[1]+2]
    inhull=all(hull.get(anchor[n],(0,0))[0]<=j[1] and j[2]<=hull.get(anchor[n],(0,10**12))[1] for j in injs[n]) if anchor[n] in hull else None
    v=('ambiguous-anchor' if alts else 'valid-anchor')+('' if inhull in (True,None) else '/junction-outside-core'); verdict[v]+=1; rows.append((n,mapq[n],anchor[n],v,(a,alts[:3])))
print(f"anchored reads {len(anchor)} | verdicts: {dict(verdict)}")
for band,sel in (('MAPQ60',[r for r in rows if r[1]>=60]),('MAPQ<60',[r for r in rows if r[1]<60])):
    print(f"  {band}: n={len(sel)} {dict(collections.Counter(r[3] for r in sel))}")
# the disagreements seen so far
disc={'SRR27178663.848890','SRR27438212.7897311','SRR27438212.1595710','SRR27438212.2780452','SRR27438212.983160','SRR27178663.1690333','SRR27438212.8049629','SRR27438212.2432386','SRR27178663.742461','SRR27438212.9782780','SRR27438212.2714427','SRR27438212.2621222','SRR27438212.2824742'}
print("disagreement reads (read, mapq, anchor, verdict, anchor-hit, alternatives):")
for r in rows:
    if r[0] in disc: print('  ',r)
with open('truth_audit.tsv','w') as o:
    o.write('read\tmapq\tanchor\tverdict\tanchor_hit\talternatives\n')
    for r in rows: o.write('\t'.join(map(str,r))+'\n')
