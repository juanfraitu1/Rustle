#!/usr/bin/env python3
"""SEMI-SYNTHETIC top-up for NPIPB12: does clean, locus-internal read evidence recover it?

§24b showed NPIPB12 is missed despite 273 reads because 86% of its spliced reads are mis-chained out of the
locus across a 104 kb intron. This tests the counterfactual directly: simulate reads from NPIPB12's OWN
annotated transcript (not a paralog's -- that would be circular), map them exactly as the real data is
mapped, merge with the real chr16 BAM, and re-run the catalogue.

⚠ Detection alone would be a trivial result -- reads simulated from a reference transcript always align back
cleanly. The INFORMATIVE question is whether NPIPB12 then JOINS the NPIP family, because its homology edges
cap at 0.25 coverage (§24b). Two possible outcomes, both worth knowing:
   copy + family  -> mis-chaining was the whole story
   copy, no family -> mis-chain AND the coverage floor both bind
"""
import gzip, re, subprocess, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sim_reads import simulate_reads

HFA="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
HGFF="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
SAM="/home/juanfra/miniforge3/bin/samtools"; MM2="/home/juanfra/miniforge3/bin/minimap2"
O="/home/juanfra/winloci_scratch/npipb12_topup"; os.makedirs(O, exist_ok=True)
GENE="NPIPB12"; CHROM="chr16"; S,E = 29765424, 29788033
N_READS=int(sys.argv[1]) if len(sys.argv)>1 else 60

# --- NPIPB12's own transcripts, longest first
tx=dict()
with gzip.open(HGFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon" or f[0]!=CHROM: continue
        g=re.search(r'(?:^|;)gene=([^;]+)',f[8]); t=re.search(r'(?:^|;)transcript_id=([^;]+)',f[8])
        if g and t and g.group(1).upper()==GENE:
            tx.setdefault(t.group(1),[]).append((int(f[3])-1,int(f[4])))
if not tx:
    sys.exit(f"no annotated transcript for {GENE}")
best=max(tx.items(), key=lambda kv: sum(b-a for a,b in kv[1]))
tid, exons = best[0], sorted(best[1])
print(f"{GENE}: {len(tx)} annotated transcripts; using {tid} with {len(exons)} exons, "
      f"{sum(b-a for a,b in exons)} bp spliced, span {exons[0][0]}-{exons[-1][1]}")

seq=[]
for a,b in exons:
    r=subprocess.run([SAM,"faidx",HFA,f"{CHROM}:{a+1}-{b}"],capture_output=True,text=True).stdout
    seq.append("".join(r.split("\n")[1:]))
mrna="".join(seq)

# IsoSeq-like: low error, some 5'/3' degradation -- matched to the real library's character
reads=simulate_reads(mrna, N_READS, err=0.003, indel=0.001, seed=12, trunc_frac=0.15)
with open(f"{O}/sim.fastq","w") as fh:
    for i,(s,q) in enumerate(reads):
        fh.write(f"@SIMTOPUP|{GENE}|{i}\n{s}\n+\n{q}\n")
print(f"simulated {len(reads)} reads from {GENE}'s own transcript -> {O}/sim.fastq")
