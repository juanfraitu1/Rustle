"""O2-2: junction-anchored truth. A read whose intron (donor,acceptor) matches a junction of exactly ONE copy is anchored to it."""
import csv,collections,pysam,sys
d=sys.argv[1]; BAM=pysam.AlignmentFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
cp={}; junc=collections.defaultdict(set)
for r in csv.DictReader(open(f'{d}/copies.tsv'),delimiter='\t'):
    cp[r['copy_idx']]=(r['chrom'],int(r['start']),int(r['end']))
    bl=[tuple(map(int,x.split('-'))) for x in r['exons'].split(',')]
    for (a,b),(c,e) in zip(bl,bl[1:]): junc[r['copy_idx']].add((r['chrom'],b,c))
owner=collections.defaultdict(set)
for i,js in junc.items():
    for j in js: owner[j].add(i)
uniqj={j for j,o in owner.items() if len(o)==1}
import os
arms={arm:{r['read_name']:r for r in csv.DictReader(open(f'{d}/{arm}.assignments.tsv'),delimiter='\t')} for arm in ('A','B','C') if os.path.exists(f'{d}/{arm}.assignments.tsv')}
names=set(arms['A'])
anchor={}; mapq={}
c0,s0,e0=cp['0'][0],min(v[1] for v in cp.values()),max(v[2] for v in cp.values())
for a in BAM.fetch(c0,s0-50000,e0+50000):
    if a.flag&2308 or a.query_name not in names: continue
    pos=a.reference_start; hits=set()
    for op,l in a.cigartuples:
        if op in (0,7,8,2): pos+=l
        elif op==3:
            j=(a.reference_name,pos,pos+l)
            if j in uniqj: hits|=owner[j]
            pos+=l
    if len(hits)==1: anchor[a.query_name]=next(iter(hits)); mapq[a.query_name]=a.mapping_quality
print(f"{d}: copies {len(cp)}, copy-specific junctions {len(uniqj)}/{len(owner)}; junction-anchored reads {len(anchor)} of {len(names)} molecules")
def tab(label,sel):
    for arm in arms:
        st=collections.Counter(); ok=tot=0
        for n in sel:
            r=arms[arm][n]; st[r['status']]+=1
            if r['status']=='assigned': tot+=1; ok+=(r['assigned_copy']==anchor[n])
        print(f"   {label:34s} arm {arm}: n={len(sel):4d} assigned {st['assigned']:4d} tied {st['tied']:4d} amb {st['ambiguous']:3d} | junction agreement {ok}/{tot}")
tab('all anchored',list(anchor))
tab('anchored, MAPQ<60',[n for n in anchor if mapq[n]<60])
changed=[n for n in anchor if 'B' in arms and (arms['A'][n]['status']!=arms['B'][n]['status'] or arms['A'][n]['assigned_copy']!=arms['B'][n]['assigned_copy'])]
tab('anchored, changed A->B (re-threaded)',changed)
