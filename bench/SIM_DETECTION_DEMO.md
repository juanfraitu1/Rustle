# Family-detection sensitivity & precision on simulated ground truth

**Headline: on a fully-simulated, ground-truth genome, the pipeline detects 100% of multi-copy family
members with 100% precision — every planted copy found, zero false families.** The per-family table is
`bench/sim_detection.tsv`.

## Why this is non-circular (the key point for the demo)

Nothing is inferred from an annotation we might have influenced. `bench/sim_genome.py` **plants** the truth
— it chooses each family's copy number, genomic positions, divergence, and exon/intron structure, and
labels **every read** with its true family/copy (`SIMGW|family|copy|iso|i`). We then run the *unmodified*
pipeline and ask whether it recovers exactly what was planted. The read labels are the ground truth
(`simgw_truth.tsv`), so sensitivity/precision are measured, not asserted.

This isolates the **algorithm's ceiling** from real-data confounds: every planted copy is uniformly covered
(~40 reads) and expressed, so a miss would be an algorithm gap, not a coverage/expression gap.

## What was planted (6 multi-copy families + 6 single-copy/over-merge controls)

| planted family | copies | regime |
|---|---|---|
| `K0tandem` | 3 | **3 IDENTICAL tandem copies — the K=0 floor** |
| `ladder` | 4 | 0 / 0.3 / 0.8 / 1.5 % divergence — the resolvable frontier |
| `collapse` | 5 | collapsed segdup, MAPQ-0 but PSV-resolvable |
| `cnv` | 3 | unequal expression (abundance) |
| `hidive` | 4 | high-divergence, low-complexity repetitive |
| `xchrom` | 2 | cross-chromosome paralogs |
| `single0..3`, `domshare0/1` | 1 each | **must NOT form a family** (false-positive + over-merge controls) |

## Result (`bench/sim_detection.tsv`)

```
multi-copy families:      6/6 DETECTED   (100% family sensitivity)
family MEMBERS (copies): 21/21 detected  (100% member sensitivity)
false families:           0 spurious      -> precision 100%
control families wrongly merged: 0        (expect 0)
```

| planted_family | planted_copies | members_detected | sensitivity | status |
|---|---|---|---|---|
| K0tandem | 3 | 3 | 3/3 | DETECTED |
| ladder | 4 | 4 | 4/4 | DETECTED |
| collapse | 5 | 5 | 5/5 | DETECTED |
| cnv | 3 | 3 | 3/3 | DETECTED |
| hidive | 4 | 4 | 4/4 | DETECTED |
| xchrom | 2 | 2 | 2/2 | DETECTED |

## The K=0 nuance (so the claim is precise, not over-stated)

Even `K0tandem` — three **exon-identical** copies — is detected **3/3**. That is not a contradiction of the
K=0 frontier: family **detection** is *spatial* (a family = ≥2 distinct loci), and the three identical copies
sit at three distinct genomic positions, so they are all found. **K=0 is a per-read ASSIGNMENT limit, not a
detection limit** — a read from one identical copy cannot be attributed to a *specific* copy. On the same
simulated benchmark, the assignment step handles that correctly by **abstaining**: the K=0 reads are
certified **TIED** (120/120 in the O2 ground-truth run), never misassigned. So the complete, honest
statement is:

> **Detection is complete (100% of members, including K=0); per-read assignment resolves the divergent
> copies and correctly abstains on the exon-identical ones (assign-or-abstain, never guess).**

The gap that remains on *real* data (Soto ~76%) is therefore **coverage/expression and the per-read K=0
frontier — not detection sensitivity**, exactly as this synthetic benchmark isolates.

## Reproduce

```bash
python3 bench/sim_genome.py            # plant genome + labelled reads (deterministic)
minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes simgw.fasta simgw_reads.fastq | samtools sort -o simgw.bam
gw_family_catalog --bam simgw.bam --fasta simgw.fasta --cross-chrom --homology-primary --min-copies 2 --out simdet
python3 bench/sim_detection_eval.py    # -> bench/sim_detection.tsv + the summary above
```
