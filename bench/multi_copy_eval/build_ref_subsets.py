#!/usr/bin/env python3
"""Build per-family reference GFF subsets for GOLGA6/GOLGA8/AMY/TBC1D3."""
import re
import sys

GFF = "/scratch/jxi21/Assembler/GGO_genomic.gff"

FAMILIES = {
    # Family name -> regex matching the gene `description` attribute
    "GOLGA6": re.compile(r"golgin subfamily A member 6|golgin A6 family", re.IGNORECASE),
    "GOLGA8": re.compile(r"golgin subfamily A member 8", re.IGNORECASE),
    "AMY":    re.compile(r"alpha-amylase|amylase alpha [12]", re.IGNORECASE),
    "TBC1D3": re.compile(r"TBC1 domain family member 3[A-Z]?-?(like)?$|TBC1 domain family member 3[A-Z]?\b(?!.*23)", re.IGNORECASE),
}

# Pass 1: collect gene IDs per family.
gene_ids_per_fam = {f: set() for f in FAMILIES}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene": continue
        attrs = parts[8]
        m_id = re.search(r"ID=([^;]+)", attrs)
        m_desc = re.search(r"description=([^;]+)", attrs)
        if not m_id or not m_desc: continue
        gid, desc = m_id.group(1), m_desc.group(1)
        for fam, rx in FAMILIES.items():
            if rx.search(desc):
                gene_ids_per_fam[fam].add(gid)

for fam, ids in gene_ids_per_fam.items():
    print(f"[{fam}] {len(ids)} matching genes", file=sys.stderr)

# Pass 2: collect mRNA IDs whose Parent is one of those gene IDs.
mrna_ids_per_fam = {f: set() for f in FAMILIES}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9: continue
        if parts[2] not in ("mRNA", "transcript", "ncRNA", "lnc_RNA"): continue
        attrs = parts[8]
        m_id = re.search(r"ID=([^;]+)", attrs)
        m_par = re.search(r"Parent=([^;]+)", attrs)
        if not m_id or not m_par: continue
        rid, par = m_id.group(1), m_par.group(1)
        for fam, ids in gene_ids_per_fam.items():
            if par in ids:
                mrna_ids_per_fam[fam].add(rid)

for fam, ids in mrna_ids_per_fam.items():
    print(f"[{fam}] {len(ids)} matching mRNAs", file=sys.stderr)

# Pass 3: write per-family GFFs containing gene + mRNA + exon/CDS records
# whose ID or Parent matches.
out_files = {f: open(f"ref_{f}.gff", "w") for f in FAMILIES}
counts = {f: {"gene":0, "mRNA":0, "exon":0, "CDS":0, "other":0} for f in FAMILIES}

with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            for f in out_files.values(): f.write(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9: continue
        attrs = parts[8]
        m_id = re.search(r"ID=([^;]+)", attrs)
        m_par = re.search(r"Parent=([^;]+)", attrs)
        rid = m_id.group(1) if m_id else None
        par = m_par.group(1) if m_par else None
        ftype = parts[2]
        for fam, gids in gene_ids_per_fam.items():
            mids = mrna_ids_per_fam[fam]
            keep = False
            if rid and rid in gids: keep = True
            elif par and par in gids: keep = True
            elif rid and rid in mids: keep = True
            elif par and par in mids: keep = True
            if keep:
                out_files[fam].write(line)
                key = ftype if ftype in counts[fam] else "other"
                counts[fam][key] += 1

for f in out_files.values(): f.close()

print()
print("=== Final reference subset counts ===", file=sys.stderr)
for fam, cnt in counts.items():
    print(f"  ref_{fam}.gff: gene={cnt['gene']} mRNA={cnt['mRNA']} exon={cnt['exon']} CDS={cnt['CDS']} other={cnt['other']}", file=sys.stderr)
