# PREREG — what counts as a FOUND copy: spliced support, not subregion overlap (written before any run, 2026-10-04)

## The problem (user, 2026-10-04 00:14, on the NPIP Read Pools page F3gJty4egn598SCZ9RBiM1)

The page's "own node" (and its headline 23 / 21 / 24 copies with a node of their own) calls a copy found when a de novo locus of the
NPIP clusters has a representative whose exons OVERLAP the copy's exons on the same strand. Figure 6d / Figure 8 call a Liftoff locus
"read-supported" when >= 2 primary reads have an aligned block on its exon union (`figures/_liftoff.py` `SUPPORT_READS`). Neither asks
whether any read is a spliced transcript of the copy: a read inside one exon, a read spanning an intron unspliced, or a 6.6-kb unspliced
stub (NPIPB4's P-arm "own node", `docs/NPIP_READ_SUPPORT_2026-10-03.md`) all count. The tools benchmark (`docs/COPY_RECOVERY_TOOLS_*.md`)
already uses the strict notion (E2 = copies with >= 2 exact-chain reads). This prereg fixes ONE standard for "found" and applies it first to
the page, then to the figure denominators.

## Definitions (fixed now)

- **Reads at a copy:** primary records (`-F 2308`) of the sample's BAM whose alignment overlaps the copy's annotated exon union on the
  copy's strand (the read's alignment strand; the libraries are stranded FLNC).
- **Junction:** an `N` CIGAR operation of >= 50 bp between two aligned blocks, keyed by its exact (donor, acceptor) coordinates. Annotated
  introns of < 50 bp are not counted either (CAT NPIPB4 carries 13-61 bp lift artefacts, RefSeq NPIPB4 1-2 bp frameshift gaps; the same
  floor on both sides keeps the comparison fair).
- **Supported junction (annotation-free):** a junction inside the copy's span carried by >= 3 reads at the copy (the 10-03 note's rule).
- **Structural-support read:** a read at the copy carrying >= k supported junctions, with k = min(2, the copy's annotated intron count
  over its transcripts); for an intronless copy (k = 0: a retrocopy or a processed pseudogene is unspliced by nature) a read whose aligned
  blocks cover >= 50% of the exon union.
- **Spliced-expressed copy:** >= 2 structural-support reads (the house floor of 2, register row 1073).
- **Found copy (by a locus set):** spliced-expressed AND some locus of the set on the copy's strand has a representative whose junctions
  include >= k of the copy's supported junctions (k = 0: the representative's exons cover >= 50% of the exon union). "Own node" / "read-
  supported" as previously defined is reported beside as the OLD rule.
- Reported beside, not decided on: reads matching >= 2 ANNOTATED introns (the annotation-bound variant); exact-chain reads (E2, the tools
  benchmark's notion: the read's junction set equals a transcript's intron set).

## Substrates and what is run

1. **Human NPIP, the page's 25 CAT/Liftoff chr16 copies** (`copy_recovery_tools_cat/ann/copies.hsa.tsv` family NPIP, transcripts from its
   `truth.hsa.gtf`), BAM `winloci_data/A119b.t2t.bam`, the three arms' loci `/mnt/linuxdisk/tmp/readpool_npip/{P,GOOD,ALL}.gff3` and the
   page's NPIP-cluster node sets (`pagedata.json`).
2. **Gorilla NPIP, 25 copies** (`copy_recovery_tools/ann/copies.ggo.tsv`, `truth.ggo.gtf`; OR6737 testis `winloci_data/GGO_mm.bam`):
   spliced-expressed counts only (no arm loci there; the tools benchmark found no exact-chain read at any gorilla NPIP copy).
3. Scorer `bench/copy_support.py` (per-copy table + summary JSON). Light: indexed region fetches.

## Pre-registered statements

- **P1 (the user's claim):** under the strict rule the number of found NPIP copies is LOWER than the page's own-node count in at least one
  arm, and the copies that drop are the ones the 10-03 note flagged as fragment-dominated (NPIPB4, NPIPB13, NPIPB10P, LOC124907807,
  LOC128966608) or others — the list is reported, whatever it is.
- **P2 (reporting, no bar):** per copy and arm: reads at copy, unspliced, 1-junction, structural-support reads, annotated-intron-matched
  reads, exact-chain reads, old own node, strict found. The page's headline is re-stated under the strict rule.
- **Decision:** the strict rule becomes THE "found copy" definition for every copy-recovery claim in the repo (the page, Figure 6d / 8's
  denominators, future O1 recovery tables); the old overlap counts stay only as "overlap" columns. Re-running the figure data under it is a
  follow-up with its own cost (noted in `docs/PENDING_2026-10-04.md`), not part of this prereg's run.

## Not changed

Node admission in the assembler (single-exon loci stay admitted: retrocopies are real members), the tools benchmark's E2, O2's certificate.
