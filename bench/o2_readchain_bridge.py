import csv,collections,pysam,os,sys
src,dst=sys.argv[1],sys.argv[2]; os.makedirs(dst,exist_ok=True)
FA=pysam.FastaFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa'); BAM=pysam.AlignmentFile('/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
GIANT=50000; MINSUP=3
rows=list(csv.DictReader(open(f'{src}/copies.tsv'),delimiter='\t')); fid=rows[0]['family_id']+'_rc'
RC=str.maketrans('ACGTacgtNn','TGCAtgcaNn'); out=[]; note=collections.Counter()
for r in rows:
    c,s,e=r['chrom'],int(r['start']),int(r['end'])
    reads=[]
    for a in BAM.fetch(c,s,e):
        if a.flag&2308: continue
        blocks=[]; introns=[]; pos=a.reference_start
        for op,l in a.cigartuples:
            if op in (0,7,8): blocks.append((pos,pos+l)); pos+=l
            elif op==2: pos+=l
            elif op==3: introns.append((pos,pos+l)); pos+=l
        mb=[]
        for b in blocks:
            if mb and b[0]<=mb[-1][1] and not any(i[0]==mb[-1][1] for i in introns): mb[-1]=(mb[-1][0],b[1])
            else: mb.append(b)
        if not any(b[1]>s and b[0]<e for b in mb): continue
        ts=a.get_tag('ts') if a.has_tag('ts') else None
        strand=('-' if a.is_reverse else '+') if ts is None else (ts if not a.is_reverse else {'+':'-','-':'+'}[ts])
        reads.append((mb,introns,strand))
    sup=collections.Counter(i for _,ins,_ in reads for i in ins)
    # mis-chain cut: drop blocks beyond a giant unsupported intron (keep the piece overlapping the copy)
    cov=collections.Counter(); strands=collections.Counter()
    for mb,ins,st in reads:
        keep=[]; cur=[mb[0]]
        for i,intr in enumerate(ins):
            if intr[1]-intr[0]>GIANT and sup[intr]<MINSUP: keep.append(cur); cur=[]
            if i+1<len(mb): cur.append(mb[i+1])
        keep.append(cur)
        seg=next((sg for sg in keep if any(b[1]>s and b[0]<e for b in sg)),[])
        for b0,b1 in seg:
            for x in range(max(b0,s),min(b1,e)): cov[x]+=1
        strands[st]+=1
    blocks=[]; 
    for x in sorted(cov):
        if cov[x]>=MINSUP:
            if blocks and x==blocks[-1][1]: blocks[-1][1]=x+1
            else: blocks.append([x,x+1])
    blocks=[tuple(b) for b in blocks]
    if len(reads)<MINSUP or not blocks:
        blocks=[tuple(map(int,x.split('-'))) for x in r['exons'].split(',')]; strand=r['strand']; note['gff_fallback']+=1
    else:
        strand=strands.most_common(1)[0][0]; note['read_chain']+=1
    seq=''.join(FA.fetch(c,a,b) for a,b in blocks)
    if strand=='-': seq=seq.translate(RC)[::-1]
    out.append((r,blocks,strand,seq,len(reads)))
with open(f'{dst}/copies.tsv','w') as t, open(f'{dst}/copies.fa','w') as f:
    t.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n')
    for i,(r,blocks,strand,seq,n) in enumerate(out):
        t.write(f"{fid}\t{i}\t{r['tid']}\t{r['chrom']}\t{blocks[0][0]}\t{blocks[-1][1]}\t{len(blocks)}\t{strand}\t{max(n,1)}\t{','.join(f'{a}-{b}' for a,b in blocks)}\n")
        f.write(f">{fid}|{i}|{r['chrom']}:{blocks[0][0]}-{blocks[-1][1]}|{strand}|nexon={len(blocks)}\n{seq}\n")
open(f'{dst}/regions','w').write(open(f'{src}/regions').read())
print(dst,dict(note),'copies',len(out),'; exonic length GFF vs read-chain (first 6):',[(sum(int(y.split('-')[1])-int(y.split('-')[0]) for y in r['exons'].split(',')),sum(b-a for a,b in bl)) for r,bl,_,_,_ in out[:6]])
