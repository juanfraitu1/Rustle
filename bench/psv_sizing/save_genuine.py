#!/usr/bin/env python3
"""Save the strand-aware genuine multimapper-recoverable copies (the validation target)
to /tmp/cre_guided/genuine_strand.tsv: tx, family, chrom, start, end, strand, own, wrong."""
import csv, os, collections, re
import pysam
BAM="/home/juanfra/winloci_scratch/GGO.bam"
EX="/tmp/cre_guided/exons.tsv"; UNI="/tmp/cre_guided/universe.tsv"; REC="/tmp/cre_guided/gd_recovery.tsv"
GFF="/home/juanfra/winloci_scratch/GGO_tx.gff"; K=2
exon_map={}
for row in csv.DictReader(open(EX),delimiter="\t"):
    s=[int(x) for x in row["exon_starts"].split(",") if x]; e=[int(x) for x in row["exon_ends"].split(",") if x]
    exon_map[row["transcript_id"]]=(row["chrom"],list(zip(s,e)))
fam_members=collections.defaultdict(list); fam_of={}
for row in csv.DictReader(open(UNI),delimiter="\t"):
    fam_members[row["family_id"]].append(row["transcript_id"]); fam_of[row["transcript_id"]]=row["family_id"]
recovered=set()
for row in csv.DictReader(open(REC),delimiter="\t"):
    if row.get("fsm")=="true" or row.get("ism")=="true": recovered.add(row["ref_transcript"])
strand={}; want=set(exon_map)
for l in open(GFF):
    f=l.split("\t")
    if len(f)<8 or f[2] not in ("mRNA","transcript","lnc_RNA","ncRNA","primary_transcript","transcript_region"): continue
    m=re.search(r'Name=([^;]+)', f[8])
    if m and m.group(1) in want: strand[m.group(1)]=f[6]
missed=[t for t in fam_of if t in exon_map and t not in recovered]
bam=pysam.AlignmentFile(BAM,"rb")
def best(chrom,exons):
    b={}
    for s,e in exons:
        try: it=bam.fetch(chrom,max(0,s),e)
        except: continue
        for r in it:
            if r.is_unmapped: continue
            try: a=r.get_tag("AS")
            except KeyError: continue
            c=b.get(r.query_name)
            if c is None or a>c[0]: b[r.query_name]=(a,r.is_reverse)
    return b
rows=[]
for tx in missed:
    if tx not in strand: continue
    chrom,exons=exon_map[tx]; fam=fam_of[tx]
    sisters=[exon_map[s] for s in fam_members[fam] if s!=tx and s in exon_map]
    if not sisters: continue
    cas=best(chrom,exons); sas={}
    for sch,sex in sisters:
        for rd,(a,_) in best(sch,sex).items():
            if rd not in sas or a>sas[rd]: sas[rd]=a
    want_rev=(strand[tx]=="-")
    own=wrong=0
    for rd,(a,isrev) in cas.items():
        s=sas.get(rd)
        if s is None or a>s:
            if isrev==want_rev: own+=1
            else: wrong+=1
    if own>=K:
        st_,en_=exons[0][0],exons[-1][1]
        rows.append((tx,fam,chrom,st_,en_,strand[tx],own,wrong))
rows.sort(key=lambda r:-r[6])
with open("/tmp/cre_guided/genuine_strand.tsv","w") as fh:
    fh.write("tx\tfamily\tchrom\tstart\tend\tstrand\town\twrong\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
print(f"wrote {len(rows)} genuine copies to /tmp/cre_guided/genuine_strand.tsv")
PY=0
