# PREREG — what counts as a FOUND copy: spliced support, not subregion overlap (written before any run, 2026-10-04)

## The problem (user, 2026-10-04 00:14, on the NPIP Read Pools page F3gJty4egn598SCZ9RBiM1)

The page's "own node" (and its headline 23 / 21 / 24 copies with a node of their own) calls a copy found when a de novo locus of the
NPIP clusters has a representative whose exons OVERLAP the copy's exons on the same strand. Figure 6d / Figure 8 call a Liftoff locus
"read-supported" when >= 2 primary reads have an aligned block on its exon union (`figures/_liftoff.py` `SUPPORT_READS`). Neither asks
whether any read is a spliced transcript of the copy: a read inside one exon, a read spanning an intron unspliced, or a 6.6-kb unspliced
stub (NPIPB4's P-arm "own node", `docs/archive/2026-10/NPIP_READ_SUPPORT_2026-10-03.md`) all count. The tools benchmark (`docs/COPY_RECOVERY_TOOLS_*.md`)
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

## Amendment A (2026-10-04 11:45, user correction, written before the re-score): support must be the COPY's own introns

The rule above lets the reads define the junction set (any junction carried by >= 3 reads at the copy). That admits a set of fragments
that is internally consistent but shares nothing with the copy's annotated transcript: at NPIPB4, 212 reads carry >= 2 "supported"
junctions while only 144 carry >= 2 of CAT's introns, and the loci built from such reads were counted as found. The user's standard:
**the introns and exons must align with the transcript the read is said to support; otherwise it is a fragment that is there by
happenstance.** Re-registered definitions, replacing "supported junction" wherever the FOUND verdict is concerned:

- **Annotated intron set A(c):** the introns (>= 50 bp, exact donor/acceptor) of every annotated transcript of the copy (human: CAT/Liftoff
  v2.0; gorilla: RefSeq). Where RefSeq also annotates a human copy (the NPIP page), RefSeq's set is reported beside — CAT and RefSeq NPIPB4
  share 2 of 20 introns, so the choice of annotation is part of the answer.
- **Support read:** a read at the copy with >= k junctions in A(c), k = min(2, |A(c)|); k = 0 (intronless copy): aligned blocks cover
  >= 50% of the exon union.
- **Spliced-expressed copy:** >= 2 support reads.
- **FOUND copy:** spliced-expressed AND a same-strand locus whose representative's junctions include >= k introns of A(c) (k = 0: the
  representative's exons cover >= 50% of the exon union). Locus level (any transcript of the locus) reported beside.
- The annotation-free reading of the original rule is kept as a REPORTED column (it answers "is there a consistent spliced structure here
  at all"), never as the FOUND verdict. Exact-chain reads (E2) stay reported.

Applied to: the NPIP page (3 arms, 25 copies; the page is restated once more), and the representative-rule run's H3 on all 7 contigs
(`/mnt/linuxdisk/tmp/rep_rule/`, both arms) — the latter is a re-score of saved outputs, no families re-run; `docs/LOCUS_REPRESENTATIVE_
RULE_2026-10-04.md` gains the annotation-anchored H3 table and its decision clause (b) is re-evaluated under it (R_J is opt-in already;
the clause can only confirm or add a second violation). Scorer: `bench/copy_support.py` gains the `ann_*` columns; nothing else changes.

## Amendment B (2026-10-04 12:00, user correction, written before the re-score): the exon-intron CHAIN must align, and both annotations count

Amendment A counts a read for a copy when >= 2 of its junction coordinates are introns of an annotated transcript — anywhere in the read,
in any order, with anything in between. The user's standard is stricter and is the right one: **a read supports a transcript only where its
exon-intron chain aligns with the transcript's** — a matched intron must be followed (or preceded) by the transcript's next intron with the
annotated exon between them spanned exactly, i.e. a read that happens to hit the same junctions with other structure in between is not a read
of that transcript. Re-registered:

- **Chain support:** a read's junction chain (its junctions in order, inside the copy's span) is a **contiguous sub-chain of the intron chain of
  some annotated transcript of the copy** (an incomplete splice match in SQANTI's sense), with >= k junctions, k = min(2, the transcript's
  intron count); every junction of the read inside the span must belong to that sub-chain (no extra junction, no skipped intron); k = 0
  (intronless copy): coverage >= 50% of the exon union as before. Junction coordinates are exact; introns < 50 bp are dropped on both sides.
- **Both annotations:** for the human copies the transcript models are the union of the CAT/Liftoff v2.0 models and the RefSeq models of the
  same gene (the RefSeq-era 26-copy table `copy_recovery_tools/ann/copies.hsa.tsv` + `truth.hsa.gtf`, matched by gene name); a read supports
  the copy if it chain-matches a model of EITHER annotation; the two are also reported separately. Where only one annotation exists
  (gorilla RefSeq; the held-out human genes until a RefSeq mapping is wired) the single set is used and said so.
- **Spliced-expressed:** >= 2 chain-support reads. **FOUND:** spliced-expressed AND a same-strand locus whose representative's junction chain
  inside the copy's span is a contiguous sub-chain (>= k junctions) of a model of the copy. Locus level (any transcript of the locus) beside.
- Amendment A's coordinate-match columns and the original read-defined columns stay reported; neither is the verdict.

Applied as Amendment A was: the NPIP page (both annotations) and the representative-rule H3 on the 7 contigs (single annotation).

## Amendment C (2026-10-04 12:25, user: "ensure that only reads that truly support the intron chains are counted"; written before the run): expressed chains, annotation-free

`docs/archive/2026-10/NPIP_CHAIN_COMPARISON_2026-10-04.md` showed that at NPIP the annotated chains are not what is expressed (27% of multi-junction reads are an
annotated chain; the dominant chain is unannotated at 13 of 25 copies, by uniquely placed reads at reference divergence). Amendment B's verdict
therefore measures the annotation as much as the nodes. Amendment C keeps B's chain logic and replaces the annotated models by the reads'
own expressed chains:

- **Read chain:** the read's junctions inside the copy's span (exon union of the CAT ∪ RefSeq models; `N` >= 50 bp, exact coordinates), in
  order. **Unique read:** MAPQ > 0 (the aligner's own ambiguity call; a tied read is counted only where O2 assigns it — reported beside with
  tied reads included).
- **Expressed chain of a copy:** a chain of >= 2 junctions carried identically by >= 3 unique reads (the floor every junction rule of this
  prereg uses; counts at >= 2 and >= 5 reported beside, fixed now, not chosen after).
- **A read TRULY SUPPORTS the copy** iff its chain (>= 2 junctions) equals an expressed chain of the copy or is a contiguous sub-chain of one
  (a 5'-truncated read of the same isoform). Any other read — unspliced, one junction, a chain not nested in an expressed chain — does not.
  **Spliced-expressed copy:** >= 1 expressed chain (hence >= 3 support reads).
- **FOUND:** spliced-expressed AND a same-strand locus whose representative's in-span chain (>= 2 junctions) equals or is a contiguous
  sub-chain of an expressed chain of the copy; locus level (any transcript of the locus) beside. Intronless copies: the coverage rule as before.
- **Reported, no rule:** each expressed chain's class against the CAT ∪ RefSeq models (FSM / ISM / NIC / NNC); the number of expressed chains per
  copy and the reads on the dominant one; Amendment B's annotated-chain verdicts beside.
- Applied to the NPIP page (3 arms; the page is restated under C with B beside) and to the representative-rule H3 on the 7 contigs (saved arms,
  no families re-run; models of the contig's single annotation for the class column only). Scorer: `bench/copy_support.py` (`xc_*` columns,
  `--chain-floor 3`).

## Amendment D (2026-10-04 12:40, user: "a read that starts at the TSS, has up to 3 junctions and introns in common with the intron chain it is supporting"; written before the run)

Amendments B and C accept any contiguous sub-chain, so a 3' fragment of the right isoform counts as support and a 3' fragment representative
counts as found. The user's standard anchors support at the transcript start:

- **TSS-anchored support:** a read supports an annotated transcript t iff (i) its 5' end (strand-aware: the alignment start on `+`, the
  alignment end on `-`) lies within **±150 bp** of t's TSS (the library's measured 5' dispersion, register §6w3; ±50 bp and ±300 bp reported
  beside, fixed now) and (ii) its first m junctions, read 5' -> 3', equal t's first m introns, m = min(3, |introns(t)|) — so the first m exons
  after the TSS are spanned exactly; a read with fewer than m junctions does not support t. Models: CAT ∪ RefSeq (each reported alone).
  Intronless transcripts (m = 0): 5' end within the tolerance and the read covers >= 50% of the exon.
- **Spliced-expressed copy:** >= 2 TSS-anchored support reads for some model of the copy.
- **FOUND:** spliced-expressed AND a same-strand locus whose representative satisfies (i) and (ii) itself against a model of the copy; locus
  level (any transcript of the locus) beside.
- Reported beside, no rule: the same test against the reads' own chains (the dominant expressed chain of Amendment C with its carrying reads'
  modal 5' end as the TSS), and Amendments A-C's counts.
- Applied to the NPIP page (restated under D) and to the representative-rule H3 on the 7 contigs (saved arms; single annotation).

## Amendment E (2026-10-04 13:58, written before the run): the transcription start is where the CAPPED reads begin

Amendment D' took a chain's modal 5' end as its start; at 7 of 24 NPIP copies that start carries no cap signal while other reads of the copy
do (`docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, "The cap signal"). Re-registered:

- **Cap read:** a primary whose RNA 5' end carries a 1-3 bp untemplated G (leading soft clip of G on `+`, trailing soft clip of C on `-`;
  `--polish-tss`'s CAP signal). **Capped start (TSS_cap):** a 20-bp bin of cap reads' 5' ends holding >= 3 cap reads; its position = the
  median 5' end of those reads; a locus may have several.
- **Anchored expressed chain:** among the reads whose 5' end lies within ±150 bp of a TSS_cap, an identical >= 2-junction chain carried by >= 3
  uniquely placed reads; its first m introns (m = min(3, length)) are the chain's anchor.
- **A read supports the copy** iff its 5' end lies within ±150 bp of a TSS_cap and its first m junctions equal an anchored expressed chain's
  first m introns at that start. **Spliced-expressed:** >= 2 such reads. **FOUND:** a same-strand locus whose representative's 5' end lies within
  ±150 bp of a TSS_cap and whose first m junctions equal such an anchor; locus level beside. Reported beside: Amendment D (annotated TSS) and
  D' (modal start); the number of capped starts per copy and their positions.
- Applied to the NPIP page and to the representative-rule H3 on the 7 contigs (saved arms). Scorer `tc_*` columns.
