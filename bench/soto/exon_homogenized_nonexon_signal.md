# Exon-homogenized copies carry recoverable non-exon signal — PPIAL4 (ID_431)

**Question (advisor):** for the "exon-homogenized K=0" floor, have we checked whether retained
introns / UTRs / soft-clipped flank distinguish copies that are exon-identical?

**Answer: no, not fully — and when we do, ~1/3 of the reads become distinguishable.**

## The family
PPIAL4 = Soto ID_431, 5 co-family copies on chr1 (dispersed 120.9–148.7 Mb); 4 classified
`expressed-K=0: exon-homogenized (reads tie in exon, genome <99.9%)`.

## The genomic divergence is abundant, not absent
PPIAL4A vs PPIAL4E genomic sequence (core + 1.5 kb flank each side), `minimap2 --eqx asm20`:
**94.8% identity — ~191 distinguishing positions over 3.7 kb** (de = 0.0145 gap-compressed).
The exons are homogenized (reads tie there); the divergence lives in intron/UTR/flank.

## A third of the reads DISTINGUISH when non-exon sequence is used
942 unique family reads (all 5 loci, deduped) aligned to a two-copy full-genomic reference
(copyA + copyE), `minimap2 map-hifi --eqx`:

| outcome | reads | % |
|---|---|---|
| **DISTINGUISH** (uniquely assigned, MAPQ > 0) | **302** | **32%** |
| tie (MAPQ 0, exon-confined) | 640 | 67% |

Of the 302 distinguishing reads, **290 (96%) carry soft-clip flank** — the distinguishing base
is in NON-EXON sequence (UTR / intron / flank), carried by full-length IsoSeq reads that extend
past the identical exon. Distinguishing reads align at de ≈ 5–6% (they span the divergent bases).

## Why the pipeline currently calls this K=0
`copy_assign_pipeline.rs`: PSVs + the per-read feature vector are built in **spliced/exon space**
(exon_map/gen2off; introns mapped into spliced coords). The distinguishing-column search looks
only at the exon body — which is exactly where gene conversion homogenized the copies. The
non-exon divergence (191 positions) is never searched.

## Prior context
`k0_flank_experiment`: exon-confined reads resolve 0/120; flank-bearing (same) copies resolve
120/120 — but only ~1.1% of primaries carry ≥20 bp soft-clip. Here 32% distinguish, because
full-length IsoSeq carries UTR (dense, on every read) in addition to sparse soft-clip flank.

## The rescue lever (concrete, testable)
Extend the distinguishing-column search **beyond the exon body** — UTR positions (dense in
full-length IsoSeq; gene-conversion-resistant), retained-intron bases, and aggregated soft-clip
flank. On PPIAL4 that reaches ~1/3 of reads. Wiring this into the PSV/assignment path could move
some of the 31 exon-homogenized members off the DNA-only floor.

**Caveats (honest):** 2-copy test (adding all 5 copies returns some ambiguity); MAPQ>0 from
genomic alignment demonstrates the information exists — wiring it into the copy-assignment PSV
search is the follow-up; 32% is PPIAL4-specific (other families will differ).

Reproduce: `bench/soto/` region extraction + `minimap2` as logged in the session.
