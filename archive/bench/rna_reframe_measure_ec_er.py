#!/usr/bin/env python3
r"""Per-copy O2 read-resolvability proxy (NOT the O1 de-tie drop): MAPQ within transcript-homology (E_r)
families. Classifies an E_r family EC_DROPPED iff <2 of its loci carry >=5% MAPQ0 primary reads -- a per-locus
HARD-ambiguity proxy, which OVER-states the O1 family-drop ~3x (96.6% here) vs the faithful de-tie
edge-formation measure in rna_reframe_measure_detie.py (29.8%). Use detie for the O1-drop; use this only for
the per-read/per-locus MAPQ mass (0.51% within E_r loci)."""
import csv, sys, os, pysam
from collections import defaultdict

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
# Output TSV: repo-relative bench/ by default; override with RNA_REFRAME_OUT.
OUT = os.environ.get("RNA_REFRAME_OUT", "bench/rna_reframe_ec_er_rows.tsv")
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FAM = "bench/denovo_families.tsv"

# 1. coords for each DN locus
coord = {}
for r in csv.DictReader(open(META), delimiter='\t'):
    coord[r['id']] = (r['chrom'], int(r['start']), int(r['end']), int(r['n_reads']))

# 2. E_r families (multi-copy)
fams = []
for r in csv.DictReader(open(FAM), delimiter='\t'):
    mem = r['members'].split(',')
    fams.append((r['family_id'], int(r['n_copies']), mem))

bam = pysam.AlignmentFile(BAM, "rb")
refs = set(bam.references)

CAP = 2500          # cap primary reads per locus
MAXWIN = 400_000    # cap window length (avoid giant-span artifacts dominating)
MAPQ0_FRAC = 0.05   # a locus "carries ambiguity" if >=5% of its primary reads are MAPQ0

def locus_mapq(chrom, s, e):
    """Return (n_primary, n_mapq0) over a locus window (primary alignments only)."""
    if chrom not in refs:
        return (0, 0)
    e = min(e, s + MAXWIN)
    n = m0 = 0
    try:
        for a in bam.fetch(chrom, max(0, s), e):
            if a.is_secondary or a.is_supplementary or a.is_unmapped:
                continue
            n += 1
            if a.mapping_quality == 0:
                m0 += 1
            if n >= CAP:
                break
    except ValueError:
        return (0, 0)
    return (n, m0)

# per-family classification
rows = []
size_cap = int(sys.argv[1]) if len(sys.argv) > 1 else 40   # skip giant repeat blobs above this
gw_n = gw_m0 = 0
for fid, ncop, mem in fams:
    if ncop > size_cap:
        rows.append((fid, ncop, "SKIP_BLOB", 0, 0, 0, 0))
        continue
    loci_amb = 0        # loci carrying >=MAPQ0_FRAC ambiguity
    loci_reads = 0
    tot_n = tot_m0 = 0
    covered = 0
    for dn in mem:
        if dn not in coord:
            continue
        c, s, e, nr = coord[dn]
        n, m0 = locus_mapq(c, s, e)
        if n == 0:
            continue
        covered += 1
        tot_n += n; tot_m0 += m0
        if m0 / n >= MAPQ0_FRAC:
            loci_amb += 1
    gw_n += tot_n; gw_m0 += tot_m0
    # classify
    if covered < 2:
        cls = "UNCOVERED"      # can't evaluate (loci not expressed / off-BAM)
    elif loci_amb >= 2:
        cls = "EC_VISIBLE"     # >=2 loci carry >=5% MAPQ0 (hard-ambiguity proxy)
    else:
        cls = "EC_DROPPED"     # PROXY ONLY (per-locus hard-ambiguity); overstates the O1-drop ~3x
                               # -- the true O1 de-tie drop (29.8%) is rna_reframe_measure_detie.py
    frac = tot_m0 / tot_n if tot_n else 0.0
    rows.append((fid, ncop, cls, covered, loci_amb, tot_n, round(frac, 4)))

# summarize
from collections import Counter
cls_fam = Counter(r[2] for r in rows)
cls_cop = Counter()
for fid, ncop, cls, cov, amb, tn, fr in rows:
    cls_cop[cls] += ncop
print("=== per-family classification (size_cap=%d) ===" % size_cap)
print("families by class:", dict(cls_fam))
print("copies by class  :", dict(cls_cop))
evaluable = cls_fam['EC_VISIBLE'] + cls_fam['EC_DROPPED']
if evaluable:
    print("EC_DROPPED / evaluable families: %d/%d = %.1f%%" %
          (cls_fam['EC_DROPPED'], evaluable, 100*cls_fam['EC_DROPPED']/evaluable))
    evc = cls_cop['EC_VISIBLE'] + cls_cop['EC_DROPPED']
    print("EC_DROPPED / evaluable copies  : %d/%d = %.1f%%" %
          (cls_cop['EC_DROPPED'], evc, 100*cls_cop['EC_DROPPED']/evc))
print("genome-wide (within E_r loci) primary reads: %d, MAPQ0 %d = %.2f%% ambiguous (=> %.2f%% MAPQ>0 resolvable)" %
      (gw_n, gw_m0, 100*gw_m0/max(1,gw_n), 100*(1-gw_m0/max(1,gw_n))))
# dump for reuse
with open(OUT,"w") as fh:
    fh.write("family_id\tn_copies\tclass\tcovered\tloci_amb\ttot_reads\tmapq0_frac\n")
    for r in rows:
        fh.write("\t".join(map(str,r))+"\n")
