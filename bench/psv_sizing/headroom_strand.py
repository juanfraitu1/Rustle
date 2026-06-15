import csv, os, collections
import pysam
BAM="/home/juanfra/winloci_scratch/GGO.bam"
EX="/tmp/cre_guided/exons.tsv"; UNI="/tmp/cre_guided/universe.tsv"; REC="/tmp/cre_guided/gd_recovery.tsv"
REF="/tmp/cre_guided/ref.gtf"; K=2
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
# copy strand from ref.gtf (one pass)
strand={}
import re
want=set(exon_map)
for l in open(REF):
    if "\ttranscript\t" not in l: continue
    m=re.search(r'transcript_id "([^"]+)"', l)
    if m and m.group(1) in want:
        strand[m.group(1)]=l.split("\t")[6]
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
buck=collections.Counter(); genuine=[]
for tx in missed:
    if tx not in strand: continue
    chrom,exons=exon_map[tx]; fam=fam_of[tx]
    sisters=[exon_map[s] for s in fam_members[fam] if s!=tx and s in exon_map]
    if not sisters: continue
    cas=best(chrom,exons); sas={}
    for sch,sex in sisters:
        for rd,(a,_) in best(sch,sex).items():
            if rd not in sas or a>sas[rd]: sas[rd]=a
    want_rev = (strand[tx]=="-")
    dec_own=dec_wrong=0
    for rd,(a,isrev) in cas.items():
        s=sas.get(rd)
        if s is None or a>s:   # AS-decisive
            if isrev==want_rev: dec_own+=1
            else: dec_wrong+=1
    if not cas: buck["unexpressed"]+=1
    elif dec_own>=K: buck["GENUINE (own-strand decisive)"]+=1; genuine.append((fam,tx,strand[tx],dec_own,dec_wrong))
    elif dec_wrong>=K: buck["strand-confounded (wrong-strand only)"]+=1
    else: buck["thin/tied"]+=1
print(f"missed copies analyzed: {sum(buck.values())}")
for k,v in buck.most_common(): print(f"  {k:<38}{v}")
print(f"\n>>> STRAND-AWARE genuine multimapper-recoverable: {buck['GENUINE (own-strand decisive)']} (was 23 strand-blind)")
print("\nGenuine (family, tx, strand, own-strand-decisive, wrong-strand):")
for r in sorted(genuine,key=lambda x:-x[3]): print(f"  {r[0]:<18} {r[1]:<18} {r[2]} own={r[3]:<4} wrong={r[4]}")
