#!/usr/bin/env bash
# Genome-wide member map-back (fixes the -p0.5 self-suppression): find each member's genome paralogs,
# not just listed siblings. Answers: do the 7 benchmark-singletons link to a real genome paralog?
set -uo pipefail
SAM=/home/juanfra/miniforge3/bin/samtools
MM2=/home/juanfra/miniforge3/bin/minimap2
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
[ -f "$FA" ] || FA=/mnt/wsl/PHYSICALDRIVE0p2/home/juanfraitu/winloci_data/chm13v2.0.fa
O=/home/juanfra/winloci_scratch/dna_complement
BED=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.chr.bed
ATTR=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/member_attribution.final.tsv

echo "map members -> whole genome, -p 0.02 -N 100 (keep partial+divergent paralogs the self-hit suppressed)"
$MM2 -c --eqx -N 100 -p 0.02 "$FA" "$O/members.fa" 2>/dev/null > "$O/gmap.paf"
echo "PAF lines: $(wc -l < "$O/gmap.paf")"

python3 - "$BED" "$ATTR" "$O/gmap.paf" "$O" <<'PY'
import sys
from collections import defaultdict
BED,ATTR,PAF,O=sys.argv[1:5]
members=[]; hdr2m={}; length={}; famsize=defaultdict(int)
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t"); s=int(s);e=int(e)
    fam=name.split("|")[1]; gene=name.split("|")[0]; m=(c,s,e,fam,gene)
    members.append(m); hdr2m[f"{c}:{s+1}-{e}"]=m; length[(c,s,e)]=e-s; famsize[fam]+=1
rna={}
for i,ln in enumerate(open(ATTR)):
    if i==0: continue
    p=ln.rstrip("\n").split("\t"); rna[(p[0],int(p[3]))]=(p[12],p[12].startswith("FOUND"))
# per member: distinct non-self genome paralog loci with a >=1kb (or >=50% short) block, keyed by identity
par=defaultdict(dict)   # member -> {locuskey: bestid}
for ln in open(PAF):
    f=ln.rstrip("\n").split("\t")
    if len(f)<12: continue
    q=f[0]; qlen=int(f[1]); t=f[5]; ts=int(f[7]); te=int(f[8]); matches=int(f[9]); alnlen=int(f[10])
    qm=hdr2m.get(q)
    if not qm: continue
    if matches < min(1000, 0.5*qlen): continue
    qc,qS,qE=qm[0],qm[1],qm[2]
    if t==qc and ts< qE and te> qS: continue   # self locus
    ident=matches/alnlen if alnlen else 0
    lk=(t, ts//5000)  # 5kb-binned locus
    k=(qc,qS,qE)
    par[k][lk]=max(par[k].get(lk,0), ident)
def best(m):
    v=par.get((m[0],m[1],m[2]),{})
    return max(v.values()) if v else 0.0
def nloci(m,th):
    return sum(1 for id_ in par.get((m[0],m[1],m[2]),{}).values() if id_>=th)
N=len(members)
print("\n=== GENOME-WIDE DNA recall (>=1 non-self paralog locus, >=1kb block) ===")
for th in (0.98,0.95,0.90,0.85,0.80,0.0001):
    d=sum(1 for m in members if best(m)>=th)
    print(f"  id>={th:<7}: {d}/{N} = {100*d/N:.1f}%")
# compare to the 10 within-subset 'no sibling block' members
print("\n=== the 10 benchmark-singletons/complex — now linked GENOME-WIDE? ===")
watch=["DUX4L50","BOLA2B","AC243829.6","AL590399.2","AC119751.3","CR381670.1","DEFB104B","ANKRD36B","GOLGA8A","AC124944.2"]
rows=[]
for m in members:
    if m[4] in watch:
        b=best(m); n90=nloci(m,0.90); n95=nloci(m,0.95)
        rows.append((m,b,n90,n95))
for (m,b,n90,n95) in rows:
    print(f"  {m[4]:14s} fam_size={famsize[m[3]]} len={length[(m[0],m[1],m[2])]:>7}  best_par_id={b:.3f}  n_par(>=0.95)={n95} n_par(>=0.90)={n90}")
# dump per-member genome table for verification
with open(f"{O}/member_genome_paralogs.tsv","w") as fh:
    fh.write("family\tgene\tlocus\tlen\tfam_size\trna_found\tbest_par_id\tn_par_ge95\tn_par_ge90\tbest_par_locus\n")
    for m in members:
        v=par.get((m[0],m[1],m[2]),{})
        bl=max(v.items(),key=lambda x:x[1])[0] if v else ("-",0)
        rf=rna.get((m[3],m[1]),("NA",False))[1]
        fh.write(f"{m[3]}\t{m[4]}\t{m[0]}:{m[1]}-{m[2]}\t{length[(m[0],m[1],m[2])]}\t{famsize[m[3]]}\t{rf}\t{best(m):.3f}\t{nloci(m,0.95)}\t{nloci(m,0.90)}\t{bl[0]}:{bl[1]*5000}\n")
print(f"\nwrote {O}/member_genome_paralogs.tsv")
print("DONE_GMAP")
PY
