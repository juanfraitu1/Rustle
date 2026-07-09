#!/usr/bin/env python3
"""AIRTIGHT proof: homology-primary finds families with NO read-conflict — copies diverged enough that the
aligner maps every read as a PRIMARY (no de-tie edge can form). Plants one 3-copy family at 20% divergence
(identity ~0.80), spaced 50 kb apart (distinct loci), plus 2 single-copy controls, on one synthetic chrom.
Read names carry TRUE family/copy => non-circular. Then: map -> confirm all-primary -> run BOTH catalogs.
Expectation: --cross-chrom (read-conflict) forms NO 'divergent' family; --homology-primary forms it (1 fam, 3
copies) from sequence homology alone; single-copy controls form nothing in either.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_genome import make_gene, copy_layout, simulate_reads, Chrom  # reuse the validated planting machinery

OUT = "/home/juanfra/winloci_scratch"
DIVERGENCE = float(sys.argv[1]) if len(sys.argv) > 1 else 0.20   # per-copy divergence (CLI-overridable)
SPACER = 50000            # 50 kb between copies -> distinct loci, no co-location

chrom = Chrom("simC")
chrom.add_bg(30000)
reads = []
truth = []

# the divergent 3-copy family: ONE base gene, three copies each independently diverged by DIVERGENCE
exons, introns = make_gene(6)
for ci in range(3):
    g, spliced, _intr = copy_layout(exons, introns, DIVERGENCE)
    start = chrom.add_copy(g)
    chrom.add_bg(SPACER)
    truth.append(("divergent", ci, start, start + len(g)))
    for i, rq in enumerate(simulate_reads(spliced, 40, err=0.003, indel=0.0008, seed=hash(("div", ci)) & 0xffff)):
        reads.append((f"SIM|divergent|{ci}|{i}", rq))

# single-copy controls: must NOT form a family in either catalog
for s in range(2):
    e2, i2 = make_gene(6)
    g, spliced, _ = copy_layout(e2, i2, 0.0)
    start = chrom.add_copy(g)
    chrom.add_bg(SPACER)
    truth.append((f"single{s}", 0, start, start + len(g)))
    for i, rq in enumerate(simulate_reads(spliced, 30, err=0.003, indel=0.0008, seed=hash(("sing", s)) & 0xffff)):
        reads.append((f"SIM|single{s}|0|{i}", rq))

with open(f"{OUT}/simap.fasta", "w") as f:
    f.write(f">simC\n{chrom.seq()}\n")
with open(f"{OUT}/simap_reads.fastq", "w") as f:
    for name, (r, q) in reads:
        f.write(f"@{name}\n{r}\n+\n{q}\n")

print(f"planted: 1 divergent family (3 copies @ {DIVERGENCE:.0%} divergence, {SPACER//1000}kb apart) + 2 single-copy controls; {len(reads)} reads")
for fam, ci, s, e in truth:
    print(f"  {fam} copy{ci}: simC:{s}-{e}")
