#!/usr/bin/env python3
"""A KNOWN-GENE-FAMILY benchmark, independent of Soto, on the gorilla substrate.

Why not Soto: its members are segmental-duplication BLOCKS. Only 126/362 are protein_coding; 121 are
pseudogenes and 110 carry no named annotation at all. So "did we recover the family" there is partly a
question about duplicated sequence, not about genes.

Why gorilla: the only FULL-GENOME long-RNA alignment available is GGO_mm.bam (mGorGor1, 10.7M reads). Both
human IsoSeq BAMs are restricted to Soto regions -- they have literally 0 reads at GSTM, HBB, PCDHB, KRT,
TUBB and S100A -- so no non-Soto human family can be tested with them. Gorilla is also the thesis substrate.

Ground truth by ORTHOLOGY, since the gorilla annotation carries only LOC ids:
  1. define families in HUMAN by curated gene symbol (CHM13 RefSeq has symbols)
  2. extract each member's exon-sum from chm13v2.0.fa
  3. map those to GGO.fasta -- the gorilla loci they hit are the family's gorilla members

That ground truth is independent of the method under test (which uses RNA reads + self-alignment), so it is
not circular. ⚠ But it is an INFERENCE: it inherits cross-species mapping ambiguity in duplicated regions,
which is the same difficulty one level up. Treat member counts as approximate and report them as such.
"""
import gzip, re, subprocess, sys, json
from collections import defaultdict

HFA="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
HGFF="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
GFA="/home/juanfra/winloci_scratch/GGO.fasta"
GBAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
SAM="/home/juanfra/miniforge3/bin/samtools"; MM2="/home/juanfra/miniforge3/bin/minimap2"
OUT="/home/juanfra/winloci_scratch/knownfam"
subprocess.run(["mkdir","-p",OUT],check=False)

# Curated multi-copy families spanning several duplication mechanisms. Deliberately EXCLUDES families that
# sit inside Soto regions (AMY, DEFB, NPIP, GOLGA, NBPF, TBC1D3) so this is a genuinely independent set.
FAMILIES = {
 "GSTM":   r"^GSTM\d+$",          "PCDHB":  r"^PCDHB\d+$",       "PCDHGA": r"^PCDHGA\d+$",
 "KRTAP":  r"^KRTAP\d+-\d+$",     "S100A":  r"^S100A\d+$",       "TUBB":   r"^TUBB\d[A-Z]?$",
 "HBB":    r"^HB[BDGE]\d?$",      "UGT1A":  r"^UGT1A\d+$",       "CYP2C":  r"^CYP2C\d+$",
 "SERPINA":r"^SERPINA\d+$",       "APOL":   r"^APOL\d+$",        "CCL":    r"^CCL\d+$",
 "MMP":    r"^MMP\d+$",           "H2BC":   r"^H2BC\d+$",        "H4C":    r"^H4C\d+$",
 "RPS4":   r"^RPS4[XY]\d?$",      "TUBA":   r"^TUBA\d[A-Z]?$",   "SIGLEC": r"^SIGLEC\d+$",
}

# --- human family members with symbols
genes=defaultdict(list)
with gzip.open(HGFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="gene": continue
        m=re.search(r'(?:^|;)Name=([^;]+)',f[8]); b=re.search(r'gene_biotype=([^;]+)',f[8])
        if not m or not b or b.group(1)!="protein_coding": continue
        for root,rx in FAMILIES.items():
            if re.match(rx,m.group(1)):
                genes[root].append((m.group(1),f[0],int(f[3])-1,int(f[4]))); break

# --- exon-sum per human member (mRNA-like query for cross-species mapping)
exons=defaultdict(list)
with gzip.open(HGFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon": continue
        m=re.search(r'(?:^|;)gene=([^;]+)',f[8])
        if m: exons[(f[0],m.group(1))].append((int(f[3])-1,int(f[4])))

def merge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o

q=open(f"{OUT}/queries.fa","w"); nq=0
for root,mem in sorted(genes.items()):
    if len(mem)<2: continue
    for sym,c,s,e in mem:
        ex=merge(exons.get((c,sym),[]))
        if not ex: continue
        seq=[]
        for a,b in ex:
            r=subprocess.run([SAM,"faidx",HFA,f"{c}:{a+1}-{b}"],capture_output=True,text=True).stdout
            seq.append("".join(r.split("\n")[1:]))
        s2="".join(seq)
        if len(s2)>=300:
            q.write(f">{root}|{sym}\n{s2}\n"); nq+=1
q.close()
print(f"human family members with an exon-sum >=300bp: {nq} across {sum(1 for r,m in genes.items() if len(m)>=2)} families")
for root,mem in sorted(genes.items()):
    if len(mem)>=2: print(f"  {root:9s} {len(mem):>3} human members")
print(f"\nwrote {OUT}/queries.fa")
