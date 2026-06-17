#!/usr/bin/env python3
"""Two-tier capstone: overlay the RNA layer on the DNA/protein-defined families.

DNA tier = protein clusters (mmseqs2 on translated CDS) -> formal multi-copy families, INCLUDING the
ancient/diverged ones the RNA-similarity definition missed (and currently-silent coding copies).
RNA overlay = per DNA-family copy: transcribed? at what level? -> 'what actually transcribes'.

Produces, per family, the 3-number summary (DNA copies -> transcribed -> well-expressed) and the
genome-vs-transcriptome aggregate, plus the ancient-family gain over the RNA-only definition.
"""
import csv
import os
from collections import defaultdict

BASE = os.path.dirname(__file__)
CLUST = "/home/juanfra/winloci_scratch/protfam_cluster.tsv"
COV = "/home/juanfra/winloci_scratch/gene_cov.tsv"
RNA_FAM = os.path.join(BASE, "families.tsv")           # RNA-tier families (POA)
OUT_MD = os.path.join(BASE, "dna_rna_overlay.md")
OUT_TSV = os.path.join(BASE, "dna_families_overlay.tsv")
T_EXPR = 5       # transcribed
T_WELL = 40      # well-expressed
CURATED = {  # textbook families; ANCIENT ones (TAS2R/DEFB/SIGLEC/CRYBG) were missed by RNA tier
    "RABL2": ["RABL2A", "RABL2B"], "DAZ": ["DAZ1", "DAZL"],
    "APOBEC3": ["APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"],
    "RFPL": ["RFPL1", "RFPL2", "RFPL3", "RFPL4A"],
    "CRYBG": ["CRYBG1", "CRYBG2", "CRYBG3"],
}
PREFIX = ["TAS2R", "DEFB", "MAGEA", "PRAMEF", "SIGLEC"]
ANCIENT = {"TAS2R*", "DEFB*", "SIGLEC*", "CRYBG"}      # the RNA-tier misses


def main():
    # DNA families from protein clusters
    clust = defaultdict(set)
    with open(CLUST) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            clust[rep].add(mem); clust[rep].add(rep)
    dna_fams = [c for c in clust.values() if len(c) >= 2]
    g2dna = {g: i for i, c in enumerate(dna_fams) for g in c}
    allgenes = set(g for c in clust.values() for g in c)

    cov = {}
    for line in open(COV):
        if line.startswith("gene\t"):
            continue
        f = line.rstrip("\n").split("\t"); cov[f[0]] = int(f[5])

    # RNA-tier families (for the gain comparison)
    g2rna = {}
    if os.path.exists(RNA_FAM):
        with open(RNA_FAM) as fh:
            fh.readline()
            for line in fh:
                fid, n, genes = line.rstrip("\n").split("\t")
                for g in genes.split(","):
                    g2rna[g] = fid

    # curated sets (resolve prefixes against the protein-clustered gene universe)
    curated = {f: [m for m in ms if m in allgenes] for f, ms in CURATED.items()}
    for pre in PREFIX:
        curated[pre + "*"] = sorted(g for g in allgenes if g.startswith(pre))
    curated = {f: ms for f, ms in curated.items() if len(ms) >= 2}

    def recovered_in(g2fam, members):
        ch = defaultdict(int)
        for m in members:
            if m in g2fam:
                ch[g2fam[m]] += 1
        return (max(ch.values()) if ch else 0) >= 2

    # overlay per DNA family
    rows = []
    tot_copies = tot_tx = tot_well = 0
    for i, c in enumerate(sorted(dna_fams, key=lambda c: -len(c))):
        members = sorted(c)
        ntx = sum(1 for g in members if cov.get(g, 0) >= T_EXPR)
        nwell = sum(1 for g in members if cov.get(g, 0) >= T_WELL)
        tot_copies += len(members); tot_tx += ntx; tot_well += nwell
        rows.append((f"DFAM{i}", len(members), ntx, nwell, members))

    # ancient-family gain: curated families recovered by DNA vs RNA
    dna_rec = {f: recovered_in(g2dna, ms) for f, ms in curated.items()}
    rna_rec = {f: recovered_in(g2rna, ms) for f, ms in curated.items()}
    gained = [f for f in curated if dna_rec[f] and not rna_rec[f]]

    L = []
    def P(s=""): L.append(s)
    P("# Two-tier overlay: RNA layer on DNA/protein-defined multi-copy families")
    P()
    P("**DNA tier** = protein clusters (mmseqs2 easy-cluster, ≥30% id / 50% cov, on translated CDS) →")
    P("formal multi-copy families, INCLUDING ancient/diverged families the RNA-similarity definition")
    P("missed. **RNA overlay** = per copy, is it transcribed (real GGO IsoSeq), and how well.")
    P()
    P(f"- DNA multi-copy families (protein clusters ≥2): **{len(dna_fams):,}** "
      f"({tot_copies:,} gene copies)")
    P(f"- RNA-tier families (POA, for comparison): {len({v for v in g2rna.values()}):,}")
    P()
    P("## What actually transcribes (genome vs transcriptome)")
    P(f"- of {tot_copies:,} copies in DNA-defined families: **transcribed (≥{T_EXPR} reads): "
      f"{tot_tx:,} ({100*tot_tx/tot_copies:.0f}%)** ; well-expressed (≥{T_WELL}): {tot_well:,} "
      f"({100*tot_well/tot_copies:.0f}%) ; **silent: {tot_copies-tot_tx:,} "
      f"({100*(tot_copies-tot_tx)/tot_copies:.0f}%)**")
    P("- i.e. the DNA tier enumerates every copy; the RNA layer shows which are live — many copies are")
    P("  present in the genome but transcriptionally silent.")
    P()
    P("## Ancient-family gain (the whole point of the DNA tier)")
    P(f"- curated families recovered by the DNA tier that the RNA tier MISSED: "
      f"**{', '.join(gained) if gained else 'none'}**")
    P("| family | DNA tier | RNA tier | DNA copies | transcribed |")
    P("|---|---|---|---|---|")
    for f, ms in sorted(curated.items()):
        ch = defaultdict(int)
        for m in ms:
            ch[g2dna.get(m, -1)] += 1
        comp = max((k for k in ch if k >= 0), key=lambda k: ch[k], default=None)
        fammem = [m for m in ms if comp is not None and g2dna.get(m) == comp]
        ntx = sum(1 for m in fammem if cov.get(m, 0) >= T_EXPR)
        tag = " (ANCIENT)" if f in ANCIENT else ""
        P(f"| {f}{tag} | {'YES' if dna_rec[f] else 'no'} | {'YES' if rna_rec[f] else 'no'} | "
          f"{len(fammem) if fammem else 0} | {ntx} |")
    P()
    P("## Per-family 3-number summary (sample: largest + curated)")
    P("| DNA family | copies | transcribed | well-expressed | example members |")
    P("|---|---|---|---|---|")
    for fid, n, ntx, nwell, members in rows[:12]:
        P(f"| {fid} | {n} | {ntx} | {nwell} | {', '.join(members[:5])}{'…' if n>5 else ''} |")
    P()
    P("## Honest scope")
    P("- DNA tier = protein clustering of annotated CODING genes → catches ancient coding families +")
    P("  currently-silent coding copies. Non-coding/pseudogene + UNANNOTATED copies need a genomic")
    P("  self-alignment pass (next extension).")
    P("- 'transcribed' uses the REAL IsoSeq (not the ideal-coverage synthetic) — so silent = silent in")
    P("  this sample; ideal-coverage GGO would test detectability, not biology.")
    P("- copy-resolvability (which transcribed copies are distinguishable per-read) follows from the")
    P("  identifiability theorem: dispersed copies resolve by locus; co-located need PSVs (rare here).")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(OUT_TSV, "w") as fh:
        fh.write("dna_family\tn_copies\tn_transcribed\tn_well_expressed\tgenes\n")
        for fid, n, ntx, nwell, members in rows:
            fh.write(f"{fid}\t{n}\t{ntx}\t{nwell}\t{','.join(members)}\n")
    print("\n".join(L))
    print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")


if __name__ == "__main__":
    main()
