#!/usr/bin/env bash
# Reproduce Soto's SHARED-EXON criterion at exon resolution.
# For each Soto member locus: collect overlapping annotated exons (full GFF), extract their sequences,
# map genome-wide, and test whether each exon RECURS (>=98% id, >=99% of exon covered) inside a SIBLING
# member's locus (same family). A member "shares an exon" if >=1 of its exons recurs at a sibling.
set -uo pipefail
SAM=/home/juanfra/miniforge3/bin/samtools
MM2=/home/juanfra/miniforge3/bin/minimap2
FA=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
GFF=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz
ATTR=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/member_attribution.final.tsv
O=/home/juanfra/winloci_scratch/shared_exon; rm -rf "$O"; mkdir -p "$O"

echo "[1] per-member overlapping exons -> regions (dedup), only exons >=80bp (Soto uses >=99% cov; tiny exons unreliable)"
python3 - "$ATTR" "$GFF" "$O" <<'PY'
import gzip,sys
from collections import defaultdict
ATTR,GFF,O=sys.argv[1:4]
mem=[]
for i,ln in enumerate(open(ATTR)):
    if i==0: continue
    p=ln.rstrip("\n").split("\t")
    mem.append((p[0],p[1],p[2],int(p[3]),int(p[4])))  # fam,gene,chrom,start,end
# exons per chrom
exons=defaultdict(list)
for ln in gzip.open(GFF,'rt'):
    if ln[0]=='#': continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="exon": continue
    exons[f[0]].append((int(f[3]),int(f[4])))
for c in exons: exons[c]=sorted(set(exons[c]))
# assign exons to member by coordinate overlap; write regions labeled memberidx:exonidx
regf=open(f"{O}/exon_regions.txt","w")
mapf=open(f"{O}/exon_member.tsv","w")  # exonID  fam  memberidx  gene  chrom  estart  eend
cnt=0; memex=defaultdict(int)
for mi,(fam,gene,c,s,e) in enumerate(mem):
    seen=set()
    for (a,b) in exons.get(c,[]):
        if b< s: continue
        if a> e: break
        # clip exon to member span, require >=80bp
        cs,ce=max(a,s),min(b,e)
        if ce-cs+1 < 80: continue
        key=(a,b)
        if key in seen: continue
        seen.add(key)
        eid=f"m{mi}_e{cnt}"
        regf.write(f"{c}:{a}-{b}\n")
        mapf.write(f"{eid}\t{fam}\t{mi}\t{gene}\t{c}\t{a}\t{b}\n")
        cnt+=1; memex[mi]+=1
regf.close(); mapf.close()
# member table
with open(f"{O}/members.tsv","w") as fh:
    for mi,(fam,gene,c,s,e) in enumerate(mem):
        fh.write(f"{mi}\t{fam}\t{gene}\t{c}\t{s}\t{e}\t{memex[mi]}\n")
print(f"  members={len(mem)}  exon-regions(>=80bp)={cnt}  members-with-any-exon={sum(1 for v in memex.values() if v)}")
PY

echo "[2] extract exon sequences (labeled)"
# faidx by region, then relabel headers to exon IDs in order
$SAM faidx -r "$O/exon_regions.txt" "$FA" > "$O/exon_raw.fa" 2>/dev/null
python3 - "$O" <<'PY'
import sys
O=sys.argv[1]
ids=[l.split("\t")[0] for l in open(f"{O}/exon_member.tsv")]
out=open(f"{O}/exons.fa","w"); i=-1; buf=[]
def flush():
    if i>=0 and i<len(ids) and buf:
        out.write(f">{ids[i]}\n"+"".join(buf)+"\n")
for ln in open(f"{O}/exon_raw.fa"):
    if ln[0]=='>':
        flush(); i+=1; buf=[]
    else: buf.append(ln.strip())
flush(); out.close()
print(f"  wrote {i+1} exon seqs")
PY

echo "[3] map exons genome-wide (sensitive)"
$MM2 -c --eqx -N 100 -p 0.02 "$FA" "$O/exons.fa" 2>/dev/null > "$O/exons.paf"
echo "  PAF lines: $(wc -l < "$O/exons.paf")"

echo "[4] shared-exon test: exon recurs >=98% id, >=99% of exon covered, inside a SIBLING member locus"
python3 - "$O" <<'PY'
import sys
from collections import defaultdict
O=sys.argv[1]
# member loci
M={}  # mi -> (fam,gene,chrom,start,end,nexon)
famloci=defaultdict(list)
for ln in open(f"{O}/members.tsv"):
    mi,fam,gene,c,s,e,nx=ln.rstrip("\n").split("\t")
    mi=int(mi); M[mi]=(fam,gene,c,int(s),int(e),int(nx)); famloci[fam].append(mi)
# exon -> member
ex2mi={}; exlen={}
for ln in open(f"{O}/exon_member.tsv"):
    eid,fam,mi,gene,c,a,b=ln.rstrip("\n").split("\t")
    ex2mi[eid]=int(mi); exlen[eid]=int(b)-int(a)+1
# parse PAF: does exon recur at a sibling locus?
mem_shared=defaultdict(bool)   # mi -> has a shared exon
exon_shared=defaultdict(bool)
for ln in open(f"{O}/exons.paf"):
    f=ln.split("\t")
    if len(f)<12: continue
    q=f[0]; qlen=int(f[1]); qs=int(f[2]); qe=int(f[3]); t=f[5]; ts=int(f[7]); te=int(f[8]); matches=int(f[9]); alnlen=int(f[10])
    if q not in ex2mi: continue
    mi=ex2mi[q]; fam=M[mi][0]
    ident=matches/alnlen if alnlen else 0
    cov=(qe-qs)/qlen if qlen else 0
    if ident<0.98 or cov<0.99: continue
    # does target land inside a SIBLING member locus (same fam, different member)?
    for smi in famloci[fam]:
        if smi==mi: continue
        _,_,sc,ss,se,_=M[smi]
        if t==sc and ts< se and te> ss:
            mem_shared[mi]=True; exon_shared[q]=True; break
# outcomes
import csv
attr={}
for i,ln in enumerate(open("/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/member_attribution.final.tsv")):
    if i==0: continue
    p=ln.rstrip("\n").split("\t"); attr[(p[0],p[1])]=(p[12],p[12].startswith("FOUND"),p[13]=="True")
rows=[]
for mi,(fam,gene,c,s,e,nx) in M.items():
    bucket,found,topup=attr.get((fam,gene),("NA",False,False))
    rows.append(dict(mi=mi,fam=fam,gene=gene,nexon=nx,shared=mem_shared[mi],found=found,topup=topup,bucket=bucket))
with open(f"{O}/member_sharedexon.tsv","w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["fam","gene","n_exon_ge80","has_shared_exon","rna_found","topup","bucket"])
    for r in rows: w.writerow([r['fam'],r['gene'],r['nexon'],r['shared'],r['found'],r['topup'],r['bucket']])

N=len(rows)
have_exon=[r for r in rows if r['nexon']>0]
print(f"\n=== SHARED-EXON reconstruction (Soto criterion: exon recurs >=98%id/>=99%cov at a sibling member) ===")
print(f"members with >=1 exon (>=80bp): {len(have_exon)}/{N}")
print(f"members with a SHARED exon (Soto-groupable): {sum(1 for r in rows if r['shared'])}/{N} = {100*sum(1 for r in rows if r['shared'])/N:.1f}%")
print(f"  (among members that HAVE an exon: {sum(1 for r in have_exon if r['shared'])}/{len(have_exon)} = {100*sum(1 for r in have_exon if r['shared'])/max(len(have_exon),1):.1f}%)")
# cross-tab shared-exon vs RNA
def grp(r): return "FOUND" if r['found'] else ("topup" if r['topup'] else r['bucket'].replace("MISS:",""))
by=defaultdict(lambda:[0,0])
for r in rows:
    g=grp(r); by[g][0]+= 1; by[g][1]+= 1 if r['shared'] else 0
print("\n=== has-shared-exon by RNA outcome ===")
order=["FOUND","topup","not-expressed","mischain","collapse-K0","genuine-miss","seeding-gap","coverage-limited","thin-single-exon","resolved-elsewhere"]
for g in order:
    if g in by: n,s=by[g]; print(f"  {g:20s} shared-exon {s}/{n}")
# the interesting cell: RNA-missed members that DO have a shared exon (Soto groups, RNA doesn't)
missed_shared=[r for r in rows if not r['found'] and not r['topup'] and r['shared']]
missed_noshare=[r for r in rows if not r['found'] and not r['topup'] and not r['shared']]
print(f"\n=== RNA-missed WITH a shared exon (Soto would group them; RNA-specific miss): {len(missed_shared)} ===")
for r in sorted(missed_shared,key=lambda x:x['bucket']): print(f"   {r['fam']:8s} {r['gene']:14s} {r['bucket']:22s} n_exon={r['nexon']}")
print(f"\n=== RNA-missed WITHOUT any shared exon (hard for Soto's method too): {len(missed_noshare)} ===")
for r in sorted(missed_noshare,key=lambda x:x['bucket']): print(f"   {r['fam']:8s} {r['gene']:14s} {r['bucket']:22s} n_exon={r['nexon']}")
PY
echo "DONE_SHAREDEXON"
