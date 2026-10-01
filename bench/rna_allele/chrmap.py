# pri chromosome -> identical haplotype chromosome (by length, then verified by uppercase-sequence md5) and its B (other-haplotype) partner
import csv, hashlib, subprocess, sys
D="/mnt/linuxdisk/tmp/rna_allele"; H="/mnt/linuxdisk/home/juanfraitu/gorilla_haps"
PRI="/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
hap={}
for h in ("mat","pat"):
    for name,num,ln in csv.reader(open(f"{D}/{h}.len.tsv"),delimiter="\t"): hap[(h,num)]=(name,int(ln))
pri=[l.split("\t")[:2] for l in open(PRI+".fai")]
def md5(fa,name):
    p=subprocess.Popen(["samtools","faidx",fa,name],stdout=subprocess.PIPE); m=hashlib.md5(); next(p.stdout)
    for line in p.stdout: m.update(line.strip().upper())
    p.wait(); return m.hexdigest()
out=[]
for name,ln in pri:
    ln=int(ln); hits=[(h,num) for (h,num),(n,l) in hap.items() if l==ln]
    if not hits: continue
    h,num=hits[0]; src=hap[(h,num)][0]
    ok = md5(PRI,name)==md5(f"{H}/{h}.fa",src)
    other="pat" if h=="mat" else "mat"
    B=hap.get((other,num),("",0))[0] if num not in ("X","Y") else ""
    out.append((name,num,h,src,"identical" if ok else "DIFFERS",other if B else "",B))
    print(*out[-1],sep="\t",flush=True)
with open(f"{D}/chrmap.tsv","w") as f:
    f.write("pri\tchrom\tsame_hap\tsame_name\tseq_check\tB_hap\tB_name\n")
    for r in out: f.write("\t".join(r)+"\n")
