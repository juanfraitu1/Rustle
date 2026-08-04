#!/usr/bin/env bash
# Per-member DNA recoverability (Soto's map-back step) vs RNA recovery -> complementarity cross-tab.
set -uo pipefail
SAM=/home/juanfra/miniforge3/bin/samtools
MM2=/home/juanfra/miniforge3/bin/minimap2
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
[ -f "$FA" ] || FA=/mnt/wsl/PHYSICALDRIVE0p2/home/juanfraitu/winloci_data/chm13v2.0.fa
BED=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.chr.bed
ATTR=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/member_attribution.final.tsv
O=/home/juanfra/winloci_scratch/dna_complement; rm -rf "$O"; mkdir -p "$O"

echo "[1] member sequences (reference) -> multi-fasta"
awk -F'\t' '{split($4,a,"|"); printf "%s:%d-%d\n", $1, $2+1, $3}' "$BED" > "$O/regions.txt"
$SAM faidx -r "$O/regions.txt" "$FA" > "$O/members.fa" 2>/dev/null
echo "    members: $(grep -c '^>' "$O/members.fa")"

echo "[2] map members back to CHM13 (Soto's map-back params)"
$MM2 -c --end-bonus 5 --eqx -N 50 -p 0.5 "$FA" "$O/members.fa" 2>/dev/null > "$O/members.paf"
echo "    PAF lines: $(wc -l < "$O/members.paf")"

echo "[3] cross-tab RNA vs DNA per member"
python3 - "$BED" "$ATTR" "$O/members.paf" "$O" <<'PY'
import sys
from collections import defaultdict
BED, ATTR, PAF, O = sys.argv[1:5]

# members: key=(chrom, start0, end) -> (fam, gene)
members=[]; by_fam=defaultdict(list)
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t")
    fam=name.split("|")[1]; gene=name.split("|")[0]
    m=(c,int(s),int(e),fam,gene); members.append(m); by_fam[fam].append(m)

# RNA verdict per (fam,start)
rna={}
for i,ln in enumerate(open(ATTR)):
    if i==0: continue
    p=ln.rstrip("\n").split("\t")
    rna[(p[0],int(p[3]))]=(p[12], p[12].startswith("FOUND"))

def overlaps(c,s,e, c2,s2,e2): return c==c2 and s<e2 and e>s2

# PAF: query = chrom:start1-end (1-based). Parse hits.
# For each member query, find hits (id>=0.98, qcov>=0.9) that land on a SIBLING member locus (same fam, diff locus).
hits_sib=defaultdict(bool); best_sib=defaultdict(float); n_par=defaultdict(int); par_seen=defaultdict(set)
for ln in open(PAF):
    f=ln.split("\t")
    if len(f)<12: continue
    q=f[0]; qlen=int(f[1]); qs=int(f[2]); qe=int(f[3])
    tname=f[5]; ts=int(f[7]); te=int(f[8]); matches=int(f[9]); alnlen=int(f[10])
    ident=matches/alnlen if alnlen else 0
    qcov=(qe-qs)/qlen if qlen else 0
    if ident<0.98 or qcov<0.9: continue
    # query member coords
    qc, rng = q.split(":"); qS=int(rng.split("-")[0])-1; qE=int(rng.split("-")[1])
    # find this query's member record
    qm=next((m for m in members if m[0]==qc and m[1]==qS and m[2]==qE), None)
    if qm is None: continue
    fam=qm[3]
    # is the target a non-self genomic locus?
    if overlaps(qc,qS,qE, tname,ts,te):  # self hit
        continue
    par_seen[(qc,qS,qE)].add((tname, ts//1000))  # distinct paralog loci (kb-binned)
    # sibling member of same fam overlapping the target?
    for sm in by_fam[fam]:
        if (sm[0],sm[1],sm[2])==(qc,qS,qE): continue
        if overlaps(sm[0],sm[1],sm[2], tname,ts,te):
            hits_sib[(qc,qS,qE)]=True
            best_sib[(qc,qS,qE)]=max(best_sib[(qc,qS,qE)], ident)

# write per-member table + cross-tab
rows=[]
for (c,s,e,fam,gene) in members:
    rv=rna.get((fam,s),("NA",False)); bucket=rv[0]; rf=rv[1]
    dsib=hits_sib[(c,s,e)]; bid=best_sib[(c,s,e)]; npar=len(par_seen[(c,s,e)])
    rows.append((fam,gene,f"{c}:{s}-{e}",bucket,rf,dsib,round(bid,3),npar))

with open(f"{O}/member_dna_vs_rna.tsv","w") as fh:
    fh.write("family\tgene\tlocus\trna_bucket\trna_found\tdna_sibling_hit\tdna_best_sib_id\tn_paralog_loci\n")
    for r in rows: fh.write("\t".join(str(x) for x in r)+"\n")

N=len(rows)
rf=sum(1 for r in rows if r[4]); df=sum(1 for r in rows if r[5])
both=sum(1 for r in rows if r[4] and r[5])
rna_only=sum(1 for r in rows if r[4] and not r[5])
dna_only=sum(1 for r in rows if r[5] and not r[4])
neither=sum(1 for r in rows if not r[4] and not r[5])
union=both+rna_only+dna_only
print(f"\n=== COMPLEMENTARITY (n={N} members) ===")
print(f"RNA-recovered:        {rf}/{N} = {100*rf/N:.1f}%")
print(f"DNA-recovered (sib):  {df}/{N} = {100*df/N:.1f}%")
print(f"RNA ∪ DNA:            {union}/{N} = {100*union/N:.1f}%")
print(f"  both:      {both}")
print(f"  RNA-only:  {rna_only}   (DNA misses: too divergent for >=98% sib link)")
print(f"  DNA-only:  {dna_only}   (RNA misses: not-expressed / K=0 / mischain)")
print(f"  neither:   {neither}")
print(f"\n=== DNA-only recoveries by RNA-miss bucket (what DNA fills in) ===")
b=defaultdict(int)
for r in rows:
    if r[5] and not r[4]: b[r[3]]+=1
for k,v in sorted(b.items(),key=lambda x:-x[1]): print(f"  {v:3d}  {k}")
print(f"\n=== RNA-only (DNA-missed) members — which does DNA fail on? ===")
for r in rows:
    if r[4] and not r[5]:
        print(f"  {r[0]:8s} {r[1]:14s} {r[2]:24s} best_sib_id={r[6]} n_par={r[7]}")
print(f"\n=== neither (hard core) ===")
for r in rows:
    if not r[4] and not r[5]:
        print(f"  {r[0]:8s} {r[1]:14s} {r[2]:24s} bucket={r[3]} best_sib_id={r[6]} n_par={r[7]}")
PY
echo "DONE_COMPLEMENT"
