"""score assignments in NEW family dir against the FIXED junction-anchored truth built from the ORIGINAL family's model junctions"""
import csv,collections,pysam,sys
orig,new=sys.argv[1],sys.argv[2]; BAM=pysam.AlignmentFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
junc=collections.defaultdict(set); cp={}
for r in csv.DictReader(open(f'{orig}/copies.tsv'),delimiter='\t'):
    cp[r['copy_idx']]=(r['chrom'],int(r['start']),int(r['end'])); bl=[tuple(map(int,x.split('-'))) for x in r['exons'].split(',')]
    for (a,b),(c,e) in zip(bl,bl[1:]): junc[r['copy_idx']].add((r['chrom'],b,c))
owner=collections.defaultdict(set)
for i,js in junc.items():
    for j in js: owner[j].add(i)
uniqj={j for j,o in owner.items() if len(o)==1}
A={r['read_name']:r for r in csv.DictReader(open(f'{new}/A.assignments.tsv'),delimiter='\t')}
A0={r['read_name']:r for r in csv.DictReader(open(f'{orig}/A.assignments.tsv'),delimiter='\t')}
c0=cp['0'][0]; s0=min(v[1] for v in cp.values()); e0=max(v[2] for v in cp.values()); anchor={}; mapq={}
for a in BAM.fetch(c0,s0-50000,e0+50000):
    if a.flag&2308 or a.query_name not in A0: continue
    pos=a.reference_start; hits=set()
    for op,l in a.cigartuples:
        if op in (0,7,8,2): pos+=l
        elif op==3:
            j=(a.reference_name,pos,pos+l)
            if j in uniqj: hits|=owner[j]
            pos+=l
    if len(hits)==1: anchor[a.query_name]=next(iter(hits)); mapq[a.query_name]=a.mapping_quality
print(f"{new} scored on {orig}'s fixed truth: anchored reads {len(anchor)}")
for lab,sel in (('all anchored',list(anchor)),('MAPQ<60',[n for n in anchor if mapq[n]<60])):
    for tag,AA in (('orig',A0),('new ',A)):
        st=collections.Counter(); ok=tot=0; bad=[]
        for n in sel:
            r=AA.get(n)
            if r is None: st['absent']+=1; continue
            st[r['status']]+=1
            if r['status']=='assigned': tot+=1; ok+=(r['assigned_copy']==anchor[n]); bad+= [] if r['assigned_copy']==anchor[n] else [(n,anchor[n],r['assigned_copy'],r['min_p_value'])]
        print(f"   {lab:12s} {tag}: n={len(sel):4d} assigned {st['assigned']:4d} tied {st['tied']:4d} amb {st['ambiguous']:3d} absent {st['absent']:3d} | agreement {ok}/{tot} {bad[:4]}")
