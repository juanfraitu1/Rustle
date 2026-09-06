import csv,collections,glob,os,pysam
M='/mnt/linuxdisk/home/juanfraitu/mcl_ann'; BAM=pysam.AlignmentFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
def score(d):
    cp={r['copy_idx']:(r['chrom'],int(r['start']),int(r['end']),r['tid'].startswith('ADM_')) for r in csv.DictReader(open(f'{d}/copies.tsv'),delimiter='\t')}
    A={r['read_name']:r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'),delimiter='\t')}
    truth={};mapq={}
    for i,(c,s,e,adm) in cp.items():
        for a in BAM.fetch(c,s,e):
            if a.flag&2308 or a.query_name not in A: continue
            o=sum(min(b1,e)-max(b0,s) for b0,b1 in a.get_blocks() if b1>s and b0<e)
            if o>0 and (a.query_name not in truth or o>truth[a.query_name][1]): truth[a.query_name]=(i,o,a.query_alignment_length); mapq[a.query_name]=a.mapping_quality
    base=[n for n in truth if not cp[truth[n][0]][3]]  # reads whose placement is an ORIGINAL unit (paired set)
    st=collections.Counter(A[n]['status'] for n in base); rej=sum(1 for n in base if A[n].get('origin_rejected')=='1')
    q60=[n for n in base if mapq[n]>=60 and A[n]['status']=='assigned']; agree=sum(1 for n in q60 if A[n]['catalog_copy_idx']==truth[n][0])
    adm=[n for n in truth if cp[truth[n][0]][3]]; adm_asg=[n for n in adm if A[n]['status']=='assigned']; adm_ok=sum(1 for n in adm_asg if A[n]['catalog_copy_idx']==truth[n][0])
    return dict(base=len(base),assigned=st['assigned'],tied=st['tied'],amb=st['ambiguous'],rej=rej,q60=len(q60),agree=agree,adm=len(adm),adm_asg=len(adm_asg),adm_ok=adm_ok)
T=collections.Counter(); rows=[]
for d in sorted(glob.glob(f'{M}/sweep_v9_adm/fam_*')):
    if not os.path.exists(f'{d}/A.done'): continue
    f=os.path.basename(d); a=score(f'{M}/sweep_v9/{f}'); b=score(d)
    rows.append((f,a,b))
    for k,v in a.items(): T['a'+k]+=v
    for k,v in b.items(): T['b'+k]+=v
print(f'families rerun: {len(rows)}')
for tag,name in [('a','v9 (before admission)'),('b','v9 + admitted copies')]:
    t={k[1:]:v for k,v in T.items() if k.startswith(tag)}
    print(f"  {name}: paired unit reads {t['base']} | assigned {t['assigned']} tied {t['tied']} amb {t['amb']} | origin-rejected {t['rej']} | MAPQ-60 agreement {t['agree']}/{t['q60']} = {t['agree']/max(1,t['q60']):.4f}" + (f" | reads at admitted loci {t['adm']}: assigned {t['adm_asg']}, to the admitted locus {t['adm_ok']}" if tag=='b' else ''))
for f,a,b in sorted(rows,key=lambda r:-r[1]['rej'])[:8]: print(f"  {f}: rejected {a['rej']} -> {b['rej']} | assigned {a['assigned']} -> {b['assigned']} | admitted-locus reads {b['adm']} assigned-to-it {b['adm_ok']}/{b['adm_asg']}")
