# What counts as a FOUND copy — spliced support at the 25 human NPIP copies (and gorilla's 25): the page's "own node" overstates recovery about 2x, and the loss is the REPRESENTATIVE, not the reads — 2026-10-04

Prereg: `docs/PREREG_spliced_copy_support_2026-10-04.md` (c5508e2a, before the run). Scorer `bench/copy_support.py`; per-copy tables
`docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv`, `docs/SPLICED_COPY_SUPPORT_ggo_npip.tsv`; work dir `/mnt/linuxdisk/tmp/readpool_npip/`
(`support_hsa.*`, `support_ggo.*`). Prompted by the user on the NPIP Read Pools page (F3gJty4egn598SCZ9RBiM1): a read somewhere between
a copy's exons and introns is not a spliced transcript of the copy.

## The rule (registered)

Reads at a copy = same-strand primaries overlapping its annotated exon union. Junction = `N` >= 50 bp, exact coordinates. Supported
junction = carried by >= 3 reads at the copy (annotation-free). Structural-support read = >= 2 supported junctions (k = min(2, annotated
introns); k = 0 -> covers >= 50% of the exon union). Spliced-expressed copy = >= 2 such reads. FOUND = spliced-expressed AND a same-strand
locus whose REPRESENTATIVE carries >= 2 of the copy's supported junctions. Reported beside: the locus-level reading (any transcript of the
locus), reads matching >= 2 annotated introns, exact-chain reads (E2's notion), the page's own-node rule (exon overlap).

## Human NPIP (A119b, 25 CAT/Liftoff chr16 copies; `copy_support.py` summary)

| | P (primaries) | GOOD (+ tied secondaries, the default) | ALL |
|---|---|---|---|
| page: copies with an own node in the NPIP clusters (exon overlap) | 23 | 21 | 24 |
| **FOUND, registered: own node whose representative carries >= 2 supported junctions** | **13** | **9** | **7** |
| locus level: some transcript of the locus carries >= 2 supported junctions | 20 | 19 | 16 |
| any same-strand locus overlapping (chr16-wide, not only NPIP clusters) | 24 | 22 | 25 |

- **Reads are not the problem: 24 of 25 copies are spliced-expressed**, most with hundreds of structural-support reads (NPIPA9 875,
  NPIPB14P 1,047, NPIPA1 641); the exception is **NPIPB13** (131 reads, 70 unspliced, 8 supported junctions but no read carries two of them).
  NPIPB4 is fragment-dominated (433 of 897 unspliced) yet still spliced-expressed (212 support reads).
- **P1 HOLDS: the strict count is lower than the own-node count in every arm** (13 < 23, 9 < 21, 7 < 24). The copies that drop in P are
  NPIPA7, NPIPA9, NPIPB4, NPIPB5, NPIPB9, NPIPB10P, NPIPB12, NPIPB13, NPIPB14P, LOC124907808 (and NPIPB2 had no node). Only two of the
  10-03 note's fragment-dominated five (NPIPB4, NPIPB13) are among them; the rest are highly expressed, well-spliced copies.
- **Where the structure is lost: the representative.** At NPIPA9 the P locus holds 42 transcripts with 38 supported junctions, and its
  representative (`DN_chr16_18372354_2`, the most-read transcript) is a 2-exon fragment with 1 junction; NPIPB14P 1 vs 15, NPIPA7 1 vs 12,
  NPIPB9 1 vs 10, LOC124907834 2 vs 17, NPIPB4 1 vs 2 (P) / 2 vs 16-18 (GOOD, ALL). The locus-level reading recovers 20 / 19 / 16: the
  loci contain the spliced copies; the one exon chain handed to the family step and to O2's copy table does not. The 5'-truncated
  library makes the most-read transcript a 3' fragment.
- **ALL is the worst arm by structure (7), not the best (24):** its extra nodes are 2-exon pieces; even at the locus level it finds 16.
  The default GOOD arm loses NPIPA6 and NPIPB6 (merged into neighbours, as the page said) and its reps are fragments at 12 copies.
- Exact-chain reads (the tools benchmark's E2): 12 copies have >= 2 (NPIPB2 193, NPIPB6 107, NPIPB15 74, LOC128966608 37, NPIPB4 27,
  NPIPB7 21, NPIPB14P 20, LOC124907807 14, NPIPA9 13, NPIPA1 10, PKD1P6-NPIPP1 5, LOC124907808 5); 13 copies have 0-1: the CAT models
  are not what the reads splice at half the copies, which is why the registered rule is annotation-free.

## Gorilla NPIP (OR6737 testis, 25 copies)

8 of 25 copies spliced-expressed (NPIPB4 9 support reads, NPIPB13 9, NPIPB14P 7, NPIPB15 7, NPIPB8 5, NPIPB2 4, NPIPA5 3, LOC124907807 2);
0-66 reads per copy, 0 exact-chain reads at any copy (as `docs/COPY_RECOVERY_TOOLS_2026-09-29.md` found). The testis library is too thin
for a locus-level test there; the fibroblast library has no NPIP.

## Decision (as registered) and what follows

- "Found copy" now means the strict rule; the page's headline is restated (its three answer cards keep their numbers as "own node"
  counts, with the strict counts beside them); overlap counts stay only as overlap.
- The lever is the representative rule, not the read pool: a locus's representative should be the transcript carrying the most supported
  junctions (or the exon union of its >= 2-read transcripts), not the most-read one. That is a change to the copy table O1 and O2 consume
  (`docs/THESIS_OBJECTIVES.md` O1 row: "positional exon sum ... of the locus's most-read transcript") and to `project_locus_representatives_study`'s
  parked verdict "rep fine for families/O2" — it is not fine at NPIP. To be pre-registered (families and figures re-scored under it; held-out
  chromosomes), listed in `docs/PENDING_2026-10-04.md`.
- Figure 6d / 8's "read-supported" denominators (>= 2 overlapping reads) and Figure 7's recovery are to be re-scored under spliced
  expression; same pending entry.

Register rows 1232-1233.
