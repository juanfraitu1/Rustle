#!/usr/bin/env python3
"""make_dna_family_manifest.py — DNA-family manifest for rustle --family-manifest.
Builds the family-manifest TSV (family_id, gene_id, chrom, start, end, strand) from the
genome cDNA all-vs-all paralog families. The one-side-intronless RETROCOPY filter is applied
by DEFAULT (drops 1931 processed-retrocopy edges; keeps real intronless tandems; mega 389->238)
-- see bench/family_def_retrocopy_filter.py. Opt out (over-merged truth) with
FAMILY_DEF_NO_RETROCOPY_FILTER=1. Domain-bridge over-merge still partly remains (max comp 238);
the junction-POSITION co-threading metric is the next de-bridging lever.
Run: python bench/make_dna_family_manifest.py > /home/juanfra/winloci_scratch/dna_family_manifest.tsv
"""
import sys
from collections import Counter
import pysam
sys.path.insert(0, "bench")
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph
from family_def_artifact_filter import partition_families, DISABLE as ART_OFF

coord = {}
for line in open("/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 4:
        coord[p[3]] = (p[0], int(p[1]), int(p[2]))
Hd, _ = dna_homology()
DG, prov, _ = build_family_graph(Hd)   # retrocopy filter on by default (env opt-out)
sys.stderr.write(f"[manifest] retrocopy filter={prov['filter_retrocopy']}: "
                 f"{prov['retrocopy_dropped']} retrocopy edges dropped, "
                 f"{prov['both_spliced_kept']+prov['both_intronless_kept']} kept\n")

# VG-topological artifact filter (default ON; opt out FAMILY_DEF_NO_ARTIFACT_FILTER=1 to KEEP+tag the
# repeat/over-merge/lncRNA families so their graphs can be inspected). Genomic seq (capped) feeds the
# k-mer self-recurrence = cycle signal; the graph itself gives cliqueness + copy-count.
fa = pysam.FastaFile("/home/juanfra/winloci_scratch/GGO.fasta")
def seq_provider(g):
    if g not in coord: return None
    c, s, e = coord[g]
    try: return fa.fetch(c, max(0, s), min(e, s + 20000))   # cap: a tandem repeat shows in the first ~20 kb
    except Exception: return None

fams = partition_families(DG, seq_provider)
side = open("/home/juanfra/winloci_scratch/dna_family_artifact_class.tsv", "w")
side.write("family_id\tclass\tn\tdensity\trecur\tkept\n")
tally = Counter()
print("family_id\tgene_id\tchrom\tstart\tend\tstrand")
for i, (members, klass, metrics, kept) in enumerate(fams):
    fid = f"DNAFAM_{i}"
    tally[klass] += 1
    side.write(f"{fid}\t{klass}\t{metrics['n']}\t{metrics['density']}\t{metrics['recur']}\t{int(kept)}\n")
    if not kept:
        continue
    rows = [(g, coord[g]) for g in members if g in coord]
    if len(rows) >= 2:
        for g, (c, s, e) in rows:
            print(f"{fid}\t{g}\t{c}\t{s}\t{e}\t.")
side.close()
n_drop = sum(1 for _, _, _, k in fams if not k)
sys.stderr.write(f"[manifest] artifact filter={'OFF (kept+tagged)' if ART_OFF else 'ON'}: "
                 f"classes={dict(tally)}; {n_drop} artifact families "
                 f"{'KEPT (disabled)' if ART_OFF else 'DROPPED'}; sidecar -> dna_family_artifact_class.tsv\n")
