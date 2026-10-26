import subprocess, sys, pysam, collections
REF='/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta'
HAP={'mat':'/mnt/linuxdisk/tmp/p18/mat.fa','pat':'/mnt/linuxdisk/tmp/p18/pat.fa'}
cmap={}
for l in open('/mnt/linuxdisk/tmp/p18/chrom_map.tsv'):
    r,*h=l.rstrip('\n').split('\t'); cmap[r]={('mat' if '_mat_' in x else 'pat'):x for x in h}
rf=pysam.FastaFile(REF)
loci=[tuple(x.split(':')) for x in sys.argv[1:]]
for c,s,e in loci:
    s,e=int(s),int(e)
    for g in ('mat','pat'):
        hc=cmap[c][g]; hf=pysam.FastaFile(HAP[g])
        L=rf.fetch(c,s-8000,s-200); R=rf.fetch(c,e+200,e+8000); X=rf.fetch(c,s,e)
        open('q.fa','w').write(f'>L\n{L}\n>R\n{R}\n')
        tgt=f't_{g}_{hc}.fa'
        import os
        if not os.path.exists(tgt): open(tgt,'w').write(f'>{hc}\n{hf.fetch(hc)}\n')
        out=subprocess.run(['minimap2','-c','-x','asm5','--secondary=no',tgt,'q.fa'],capture_output=True,text=True).stdout
        hits={}
        for l in out.splitlines():
            f=l.split('\t'); 
            if int(f[11])>=30 and f[0] not in hits: hits[f[0]]=(int(f[7]),int(f[8]),f[4])
        if 'L' not in hits or 'R' not in hits: print(c,s,e,g,hc,'flank not placed',hits); continue
        a=min(hits['L'][1],hits['R'][1]); b=max(hits['L'][0],hits['R'][0])
        gap=hf.fetch(hc,min(a,b),max(a,b)) if a!=b else ''
        open('x.fa','w').write(f'>X\n{X}\n'); open('g.fa','w').write(f'>G\n{gap or "N"*50}\n')
        o2=subprocess.run(['minimap2','-c','-x','sr','g.fa','x.fa'],capture_output=True,text=True).stdout.splitlines()
        best=max((int(l.split('\t')[9])/int(l.split('\t')[1]) for l in o2), default=0)
        print(f'{c}:{s}-{e} in {g} ({hc}): flanks at {hits["L"][:2]} / {hits["R"][:2]}, between {abs(b-a)-400} bp; exon present: {"yes" if best>0.8 else "NO"} (match {best:.2f}); ref locus length {e-s}')
