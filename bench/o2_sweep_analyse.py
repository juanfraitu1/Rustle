import csv,collections,glob,os,statistics as st,math,pysam,re
BAM=pysam.AlignmentFile(os.environ.get("O2_BAM", "/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam"))
fam_rows=[]; unit_rows=[]
for d in sorted(glob.glob('fam_*')):
    if not os.path.exists(f'{d}/A.done') or not os.path.exists(f'{d}/A.assignments.tsv'): continue
    done=open(f'{d}/A.done').read().strip(); err=open(f'{d}/A.err').read()
    t=re.search(r'elapsed=([\d.]+) rss_kb=(\d+)',err); el,rss=(float(t.group(1)),int(t.group(2))//1024) if t else (None,None)
    cp={r['copy_idx']:(r['chrom'],int(r['start']),int(r['end'])) for r in csv.DictReader(open(f'{d}/copies.tsv'),delimiter='\t')}
    fc={r['copy_idx']:r for r in csv.DictReader(open(f'{d}/forecast.tsv'),delimiter='\t')}

    A={}
    for r in csv.DictReader(open(f'{d}/A.assignments.tsv'),delimiter='\t'):
        r['assigned_copy']=r['catalog_copy_idx']; A[r['read_name']]=r   # native column (§6et)
    # hygiene: reads with an aligned block inside a unit; truth-by-placement = unit with max block overlap
    truth={}; mapq={}
    for i,(c,s,e) in cp.items():
        for a in BAM.fetch(c,s,e):
            if a.flag&2308 or a.query_name not in A: continue
            o=sum(min(b1,e)-max(b0,s) for b0,b1 in a.get_blocks() if b1>s and b0<e)
            if o>0 and (a.query_name not in truth or o>truth[a.query_name][1]): truth[a.query_name]=(i,o); mapq[a.query_name]=a.mapping_quality
    st_all=collections.Counter(A[n]['status'] for n in truth); q60=collections.Counter(A[n]['status'] for n in truth if mapq[n]>=60); qlow=collections.Counter(A[n]['status'] for n in truth if mapq[n]<60)
    agree=sum(1 for n in truth if mapq[n]>=60 and A[n]['status']=='assigned' and A[n]['assigned_copy']==truth[n][0]); asg60=q60['assigned']
    per_unit=collections.defaultdict(collections.Counter)
    for n,(i,o) in truth.items(): per_unit[i][A[n]['status']]+=1
    for i in cp:
        c=per_unit[i]; tot=sum(c.values())
        if tot>=5:
            ni=fc[i]['nearest_ident']; unit_rows.append((d,i,float(ni) if ni!='NA' else float('nan'),c['assigned']/tot,tot,float(fc[i]['rep_frac']) if fc[i]['rep_frac']!='NA' else float('nan'),fc[i]['source']))
    fam_rows.append(dict(fam=d,n=len(cp),exit=done.split()[0],reads_in_units=len(truth),assigned=st_all['assigned'],tied=st_all['tied'],amb=st_all['ambiguous'],
        q60=sum(q60.values()),q60_assigned=asg60,q60_agree=agree,qlow=sum(qlow.values()),qlow_assigned=qlow['assigned'],elapsed=el,rss_mb=rss,
        forecast_unassignable=sum(1 for i in cp if fc[i]['nearest_ident']!='NA' and float(fc[i]['nearest_ident'])>=0.99)))
ok=[r for r in fam_rows if r['exit']=='exit=0']
print(f"families with results: {len(fam_rows)}, completed: {len(ok)}, aborted: {[r['fam'] for r in fam_rows if r['exit']!='exit=0']}")
T=lambda k: sum(r[k] for r in ok)
print(f"reads with a block inside a unit: {T('reads_in_units')} | assigned {T('assigned')} tied {T('tied')} ambiguous {T('amb')}")
print(f"MAPQ-60: {T('q60')} reads, assigned {T('q60_assigned')}, placement agreement {T('q60_agree')}/{T('q60_assigned')} | MAPQ<60: {T('qlow')} reads, assigned {T('qlow_assigned')}")
print(f"wall: total {sum(r['elapsed'] for r in ok if r['elapsed']):.0f}s, max {max(r['elapsed'] for r in ok if r['elapsed']):.0f}s; RSS max {max(r['rss_mb'] for r in ok if r['rss_mb'])} MB")
# forecast vs observed at unit level
def spearman(x,y):
    def rk(v):
        s=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v); i=0
        while i<len(s):
            j=i
            while j+1<len(s) and v[s[j+1]]==v[s[i]]: j+=1
            for k in range(i,j+1): r[s[k]]=(i+j)/2+1
            i=j+1
        return r
    rx,ry=rk(x),rk(y); mx,my=st.mean(rx),st.mean(ry); num=sum((a-mx)*(b-my) for a,b in zip(rx,ry)); den=math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry)); return num/den if den else float('nan')
u=[r for r in unit_rows if not math.isnan(r[2])]
hi=[r for r in u if r[2]>=0.99]; lo=[r for r in u if r[2]<0.99]
print(f"units with >=5 reads: {len(u)}; forecast-unassignable (nearest_ident>=0.99): {len(hi)} median assigned share {st.median([r[3] for r in hi]) if hi else 'na'}; others {len(lo)} median {st.median([r[3] for r in lo]) if lo else 'na'}; Spearman(nearest_ident, assigned share) = {spearman([r[2] for r in u],[r[3] for r in u]):.3f}")
rep=[r for r in u if not math.isnan(r[5]) and r[5]>=0.8]; print(f"units with rep_frac>=0.8 (repeat-derived exon chains) among scored: {len(rep)}")
with open('sweep_families.tsv','w') as o:
    o.write('\t'.join(fam_rows[0].keys())+'\n')
    for r in fam_rows: o.write('\t'.join(str(v) for v in r.values())+'\n')
with open('sweep_units.tsv','w') as o:
    o.write('fam\tcopy_idx\tnearest_ident\tassigned_share\treads\trep_frac\tsource\n')
    for r in unit_rows: o.write('\t'.join(str(v) for v in r)+'\n')
