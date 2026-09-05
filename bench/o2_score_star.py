import csv,collections,pysam,sys
M='/mnt/linuxdisk/home/juanfraitu/mcl_ann'; D=f'{M}/adj/d3'; BAM=pysam.AlignmentFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
def score(cpf,asg,label):
    cp={r['copy_idx']:(r['chrom'],int(r['start']),int(r['end'])) for r in csv.DictReader(open(cpf),delimiter='\t')}
    A={r['read_name']:r for r in csv.DictReader(open(asg),delimiter='\t')}
    truth={};mapq={}
    for i,(c,s,e) in cp.items():
        for a in BAM.fetch(c,s,e):
            if a.flag&2308 or a.query_name not in A: continue
            o=sum(min(b1,e)-max(b0,s) for b0,b1 in a.get_blocks() if b1>s and b0<e)
            if o>0 and (a.query_name not in truth or o>truth[a.query_name][1]): truth[a.query_name]=(i,o,a.query_alignment_length); mapq[a.query_name]=a.mapping_quality
    st=collections.Counter(A[n]['status'] for n in truth); q60=[n for n in truth if mapq[n]>=60]; low=[n for n in truth if mapq[n]<60]
    asg60=[n for n in q60 if A[n]['status']=='assigned']; agree=sum(1 for n in asg60 if A[n]['catalog_copy_idx']==truth[n][0])
    inside=[n for n in asg60 if truth[n][1]>=0.5*truth[n][2]]; agree_in=sum(1 for n in inside if A[n]['catalog_copy_idx']==truth[n][0])
    lowst=collections.Counter(A[n]['status'] for n in low)
    print(f'{label}: unit reads {len(truth)} | assigned {st["assigned"]} tied {st["tied"]} amb {st["ambiguous"]} | q60 agree {agree}/{len(asg60)} (reads >=50% inside their copy: {agree_in}/{len(inside)}) | MAPQ<60 {len(low)}: {dict(lowst)}')
for fam in sys.argv[1:]:
    score(f'{M}/sweep_v7/{fam}/copies.tsv',f'{M}/sweep_v7/{fam}/A.assignments.tsv',f'{fam} OFF')
    score(f'{M}/sweep_v7/{fam}/copies.tsv',f'{D}/{fam}_star.assignments.tsv',f'{fam} read-star')
if 'fam_MCL1_073242' in sys.argv[1:]:
    core={r['copy_idx']:(r['chrom'],int(r['start']),int(r['end'])) for r in csv.DictReader(open(f'{M}/o2scale/fam_NPIPcore_073242/copies.tsv'),delimiter='\t')}
    v7={r['copy_idx']:(r['chrom'],int(r['start']),int(r['end'])) for r in csv.DictReader(open(f'{M}/sweep_v7/fam_MCL1_073242/copies.tsv'),delimiter='\t')}
    ovl=lambda a,b: a[0]==b[0] and a[1]<b[2] and b[1]<a[2]
    m={}
    for ci,cc in core.items():
        hits=[vi for vi,vv in v7.items() if ovl(cc,vv)]; m[ci]=hits[0] if len(hits)==1 else None
    truth=[r for r in csv.DictReader(open(f'{M}/o2scale/truth_valid.tsv'),delimiter='\t')]
    for label,asg in [('OFF',f'{M}/sweep_v7/fam_MCL1_073242/A.assignments.tsv'),('read-star',f'{D}/fam_MCL1_073242_star.assignments.tsv')]:
        A={r['read_name']:r for r in csv.DictReader(open(asg),delimiter='\t')}
        st=collections.Counter(); wrong=[]; right=0
        for t in truth:
            r=A.get(t['read'])
            if r is None: st['absent']+=1; continue
            st[r['status']]+=1
            if r['status']=='assigned':
                exp=m.get(t['anchor_copy'])
                if exp is None: st['unmappable']+=1
                elif r['catalog_copy_idx']==exp: right+=1
                else: wrong.append((t['read'][-8:],t['mapq'],exp,r['catalog_copy_idx'],r['n_decisive'],r['min_p_value']))
        print(f'NPIP 62 valid anchors, {label}: {dict(st)} | right {right} wrong {len(wrong)} {wrong[:6]}')
