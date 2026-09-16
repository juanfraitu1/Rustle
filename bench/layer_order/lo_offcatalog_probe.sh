#!/usr/bin/env bash
# DNA probe for universe genes outside D's catalogs (not a D build): gene spans vs the 52 genes of D MCL0 + MCL4.
# Writes integrate_slim/probe/{q,t}.regions, q_vs_t.paf (asm20 defaults), short_probe_opts.out (short genes, sensitive seeds).
set -euo pipefail
P=/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim/probe
REF=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0.fa
mkdir -p $P && cd $P
python3 - <<'PY'
import csv
U={r["gene_id"]:r for r in csv.DictReader(open("../../light/universe.corrected.tsv"),delimiter="\t")}
q=["USP6NL","LOC100420289","LOC124905656","TBC1D3P6","NPIPB1P"]
with open("q.regions","w") as fq, open("t.regions","w") as ft, open("names.tsv","w") as fn:
    for g,r in U.items():
        reg=f'{r["chrom"]}:{int(r["start"])+1}-{r["end"]}'
        fn.write(f'{reg}\t{r["name"]}\t{r["D_group"]}\t{r["family_side"]}\n')
        if r["name"] in q: fq.write(reg+"\n")
        elif r["D_group"] in ("D|c15_17_22|MCL4","D|c16_19_20|MCL0"): ft.write(reg+"\n")
PY
samtools faidx $REF -r q.regions > q.fa
samtools faidx $REF -r t.regions > t.fa
minimap2 -c -x asm20 --eqx -N 50 -t 4 t.fa q.fa 2> mm2.err > q_vs_t.paf
printf "chr1:15436125-15436587\nchr4:158717685-158718287\n" > short.regions
samtools faidx $REF -r short.regions > short.fa
for opt in "-k13 -w5" "-x sr" "-x asm20 -k13 -w5 -m 50 -s 50"; do
  minimap2 -c $opt --eqx -N 50 -t 4 t.fa short.fa 2>/dev/null > s.paf
  echo "opts [$opt]: $(wc -l < s.paf) records"
  awk -F'\t' '{printf "  %s -> %s id %.3f qcov %.3f blen %d\n",$1,$6,$10/$11,($4-$3)/$2,$11}' s.paf | sort -k6,6r | head -6
done > short_probe_opts.out
# probe.out summary (records >= 300 bp and >= 0.70 identity, pooled identity, union coverage of the query body)
python3 - <<'PY' > probe.out
import collections
names={}
for l in open("names.tsv"):
    reg,n,d,f=l.rstrip("\n").split("\t"); names[reg]=(n,d)
acc=collections.defaultdict(lambda:[0,0,[]]); allq=collections.defaultdict(int); qlen={}
for l in open("q_vs_t.paf"):
    f=l.split("\t"); nm,bl=int(f[9]),int(f[10]); qlen[f[0]]=int(f[1]); allq[f[0]]+=1
    if bl<300 or nm/bl<0.70: continue
    a=acc[(f[0],f[5])]; a[0]+=nm; a[1]+=bl; a[2].append((int(f[2]),int(f[3])))
def cov(iv,L):
    iv=sorted(iv); t=0; cs,ce=iv[0]
    for s,e in iv[1:]:
        if s<=ce: ce=max(ce,e)
        else: t+=ce-cs; cs,ce=s,e
    return (t+ce-cs)/L
for q in open("q.regions").read().split():
    hits=[(names[t][0],names[t][1],a[0]/a[1],cov(a[2],qlen[q])) for (qq,t),a in acc.items() if qq==q]
    hits.sort(key=lambda x:-x[2]*x[3])
    grp=collections.Counter(h[1] for h in hits)
    print(f"{names[q][0]} ({q}): PAF records {allq.get(q,0)}; targets with >=300bp/>=0.70 records {len(hits)} by D group {dict(grp)}; "
          f"top: " + "; ".join(f"{n} id {i:.3f} qcov {c:.3f}" for n,_,i,c in hits[:3]))
PY
