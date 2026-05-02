#!/usr/bin/env python3
"""Per-paralog presence/absence breakdown across multiple assemblers.

For each family (GOLGA6/GOLGA8/AMY/TBC1D3), shows which reference paralog
loci are recovered by which tool (StringTie / rustle-em / rustle-em+snp).
This tells us *which* copies are missed and which are caught by SNP-aware EM.

Usage: per_paralog_breakdown.py <st.gtf> <em.gtf> <em_snp.gtf>

Format of gffcompare's .loci: tab-separated 4 columns
   XLOC_id  chrom[strand]start-end  ref_ids  query_ids
Where ref_ids = comma-list of "gene-XXX|rna-YYY" pairs (or "-" if no ref).
"""
import os, re, sys, subprocess, tempfile
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit("usage: per_paralog_breakdown.py <st.gtf> <em.gtf> <em_snp.gtf>")

GFFCOMPARE = "/storage/work/jxi21/.pixi_cache/pkgs/gffcompare-0.12.10-h9948957_0/bin/gffcompare"
queries = {"ST": sys.argv[1], "rustle-EM": sys.argv[2], "rustle-EM+SNP": sys.argv[3]}
families = ["GOLGA6", "GOLGA8", "AMY", "TBC1D3"]

def parse_ref_loci(reffile):
    loci = {}
    with open(reffile) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip().split("\t")
            if len(p) < 9 or p[2] != "gene": continue
            attrs = p[8]
            m_id = re.search(r"ID=([^;]+)", attrs)
            m_name = re.search(r"Name=([^;]+)", attrs)
            if not m_id: continue
            loci[m_id.group(1)] = {
                "chrom": p[0], "start": int(p[3]), "end": int(p[4]),
                "name": m_name.group(1) if m_name else m_id.group(1),
            }
    return loci

def parse_loci_file_for_caught_genes(loci_path):
    """Return set of ref gene_ids that are caught (i.e. an XLOC has both
    ref and query columns populated and the gene_id appears there)."""
    caught = set()
    with open(loci_path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 4: continue
            ref_col, q_col = p[2], p[3]
            if ref_col == "-" or q_col == "-": continue
            # ref_col looks like: gene-LOC101150678|rna-XM_055363285.2,gene-LOC101150678|rna-XM_055363287.2
            for entry in ref_col.split(","):
                entry = entry.strip()
                if not entry: continue
                gene = entry.split("|", 1)[0]   # "gene-LOC101150678"
                caught.add(gene)
    return caught

# (label, family) -> set of caught ref gene IDs
caught = {}

with tempfile.TemporaryDirectory() as tmp:
    for fam in families:
        ref = f"ref_{fam}.gff"
        for label, q in queries.items():
            prefix = os.path.join(tmp, f"cmp_{fam}_{label.replace('+','_')}")
            cmd = [GFFCOMPARE, "-r", ref, "-o", prefix, "--no-merge", q]
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            loci_path = prefix + ".loci"
            caught[(label, fam)] = parse_loci_file_for_caught_genes(loci_path)

# Build per-family per-paralog table.
for fam in families:
    print(f"\n=== {fam} ===")
    ref_loci = parse_ref_loci(f"ref_{fam}.gff")
    header = f"{'paralog (gene)':30} {'chrom':14} {'range':25} {'ST':>4} {'EM':>4} {'EM+SNP':>8}"
    print(header); print("-"*len(header))
    for gid, info in sorted(ref_loci.items(), key=lambda kv: (kv[1]['chrom'], kv[1]['start'])):
        st_hit  = 1 if gid in caught[("ST", fam)] else 0
        em_hit  = 1 if gid in caught[("rustle-EM", fam)] else 0
        snp_hit = 1 if gid in caught[("rustle-EM+SNP", fam)] else 0
        rng = f"{info['start']}-{info['end']}"
        nm = info['name'][:28]
        print(f"{nm:30} {info['chrom']:14} {rng:25} {st_hit:>4} {em_hit:>4} {snp_hit:>8}")
    print("-"*len(header))
    n_loci = len(ref_loci)
    st_n  = sum(1 for gid in ref_loci if gid in caught[("ST", fam)])
    em_n  = sum(1 for gid in ref_loci if gid in caught[("rustle-EM", fam)])
    snp_n = sum(1 for gid in ref_loci if gid in caught[("rustle-EM+SNP", fam)])
    print(f"{'TOTAL':30} {'':14} {f'{n_loci} ref paralogs':>25} {st_n:>4} {em_n:>4} {snp_n:>8}")
