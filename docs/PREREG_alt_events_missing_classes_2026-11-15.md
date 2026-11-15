# PREREG row 1744 (2026-11-15): the event classes row 1714 did not measure, and whether alignment loses short exons

User question: are all alternative splicing events modelled as the annotation expects, and would a graph (vg) help with
the tricky ones. Row 1714 covered shifted exon ends, alternative first / last exons, nested starts / ends and a cassette
control. Not measured: retained introns, mutually exclusive exons, microexons. Measurement only; written and committed
before the job ran. No default changes.

## Classes added to `bench/famcmp/alt_events.py`

Pairs of annotated RefSeq intron chains, as in row 1714 (both expressed = >= 2 full-chain primary reads each; read at
minor >= 5 reads; share with both chains emitted):

- **RI**: one chain equals the other without one internal intron (a terminal one is the nested class of row 1714).
- **MXE**: the chains differ in one internal exon replaced by a disjoint exon between the same outer splice sites.
- **MICRO**: cassette pairs whose skipped exon is <= 30 bp (a subset of CASSETTE).

**Skip reads** (the alignment question). The pair measure cannot see an exon the aligner misses: reads that skip it do
not spell the chain, so the pair is never "expressed". Hence a second table: for each annotated chain and each internal
exon that no annotated chain skips, the full-chain reads and the reads spelling the same chain with that exon skipped,
units with >= 5 such reads, by exon length (<= 15, 16-30, 31-60, > 60 bp). Unit = chain x exon (an exon shared by
isoforms counts once per chain). A skip read is a missed alignment or a real unannotated skip; the > 60 bp bin is the
control for the latter.

## Inputs

Row 1714's: OR6737 (`r1628/orrun/or.gtf`, `.final.gtf`), A119b (`r1617/default/a.gtf`, `.final.gtf`), with the benchmark
StringTie and FLAIR GTFs; KB3781 (`r1622/kb/kb.gtf`) reported without tools. Workspace `r1744_altevents2`.

## Reading rules (minor >= 5 reads; OR6737 and A119b decide, KB3781 reported)

1. **Per pair class (RI, MXE, MICRO):** fewer than 30 expressed pairs on a sample = too few to read there. Otherwise
   both-emitted share within 0.05 of that sample's CASSETTE share = modelled like ordinary splicing; lower = a gap,
   reported with the split major-only / minor-only / neither.
2. **Alignment (skip reads):** an alignment problem a graph could address is flagged when the skip-read share of the
   <= 30 bp exons (both short bins pooled) exceeds the > 60 bp control by >= 0.05 absolute on both deciding samples.
   Flagged = a graph test is worth registering, with that share as the number to beat. Not flagged = no alignment loss
   on short exons for a graph to recover; vg is not pursued for splicing events on this evidence.

## Limits

RefSeq pairs only; an RI chain's reads can be unspliced precursor RNA; a skipped short exon can also be absorbed as
an insertion or mismatches beside the junction, which gives the same skip chain and is counted; exons skipped by reads
that also differ elsewhere are not counted. The assemblies are the row 1714 files (default assembly rules have not
changed since).
