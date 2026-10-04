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

## Amendment A re-score (2026-10-04, user correction): support = the copy's OWN annotated introns

The user's objection to the first rule: a read spliced somewhere inside NPIPB4 is not a read of NPIPB4 unless its introns are NPIPB4's;
the first rule let the reads define the junction set. Amendment A (prereg, committed 8735ceeb before this re-score) anchors support and
FOUND on the annotated introns (exact splice sites, >= 50 bp; k = min(2, annotated introns); intronless copies by coverage). Scorer columns
`ann_support_reads`, `ann_expressed`, `<arm>_ann_found`, `<arm>_locus_ann_found` (`bench/copy_support.py`); the read-defined columns stay as
the annotation-free reading. Tables: `docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv` (CAT-anchored, the page's 25 copies, with the page's node
sets) and `docs/SPLICED_COPY_SUPPORT_hsa_npip_refseq.tsv` (the RefSeq-era 26-copy table, RefSeq models, no node sets).

**NPIP, CAT-anchored: FOUND = 14 / 10 / 8 of 25 (P / GOOD / ALL), within the page's NPIP-cluster nodes 13 / 9 / 8** (read-defined rule:
14 / 10 / 7 and 13 / 9 / 7); 24 of 25 copies spliced-expressed under both rules. RefSeq-anchored (26 copies): 15 / 10 / 9, all 26
spliced-expressed. So the correction changes the per-copy verdicts, not the headline: own node 23 / 21 / 24 vs found 13 / 9 / 8.

| copy | reads | support reads CAT / RefSeq | exact chains CAT / RefSeq | P found / own node | GOOD | ALL |
|---|---|---|---|---|---|---|
| NPIPB2 | 432 | 348 / 347 | 193 / 7 | 0 / 0 | 0 / 0 | 0 / 1 |
| NPIPA2 | 348 | 132 / 287 | 0 / 140 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA1 | 806 | 631 / 576 | 10 / 0 | 1 / 1 | 1 / 1 | 0 / 1 |
| PKD1P6-NPIPP1 | 253 | 139 / 361 | 5 / 0 | 1 / 0 | 1 / 0 | 1 / 1 |
| NPIPA5 | 147 | 110 / 110 | 0 / 30 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA6 | 204 | 142 / 119 | 1 / 3 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPA7 | 285 | 184 / 184 | 0 / 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA8 | 192 | 140 / 140 | 0 / 0 | 1 / 1 | 0 / 1 | 0 / 1 |
| NPIPA9 | 998 | 865 / 790 | 13 / 8 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC128966608 | 1104 | 585 / 234 | 37 / 12 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB4 | 897 | 199 / 161 | 27 / 2 | 0 / 1 | 1 / 1 | 1 / 1 |
| NPIPB5 | 785 | 317 / 46 | 1 / 2 | 0 / 1 | 0 / 1 | 1 / 1 |
| NPIPB6 | 684 | 563 / 568 | 107 / 127 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPB7 | 164 | 105 / 105 | 21 / 7 | 1 / 1 | 0 / 1 | 0 / 1 |
| NPIPB8 | 174 | 128 / 27 | 0 / 5 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB9 | 568 | 433 / 455 | 0 / 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB10P | 70 | 28 / 27 | 0 / 0 | 0 / 1 | 1 / 1 | 1 / 1 |
| NPIPB11 | 150 | 82 / 80 | 0 / 1 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB12 | 54 | 9 / 9 | 1 / 1 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907834 | 634 | 302 / 296 | 0 / 20 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB13 | 131 | 0 / 52 | 0 / 1 | 0 / 1 | 0 / 1 | 0 / 0 |
| NPIPB14P | 1255 | 1023 / 1188 | 20 / 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB15 | 220 | 170 / 171 | 74 / 3 | 1 / 1 | 0 / 1 | 0 / 1 |
| LOC124907808 | 73 | 46 / 46 | 5 / 4 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907807 | 96 | 41 / 41 | 14 / 13 | 1 / 1 | 0 / 1 | 0 / 1 |

- **NPIPB2:** its reads ARE NPIPB2's transcript (348 of 432 carry >= 2 CAT introns; 193 exact CAT chains) — the failure is the nodes: none in P
  and GOOD, and ALL's six loci are 2-exon fragments carrying none of its introns (found 0 in every arm, as the user read the page).
- **NPIPB4:** the P-arm own node is the unspliced 6.6-kb stub (0 annotated introns in its representative -> not found); the GOOD and ALL
  arms hold loci whose representatives carry >= 2 CAT introns of NPIPB4 (found). 199 of 897 reads carry >= 2 CAT introns (RefSeq: 161).
- **The annotation is part of the answer:** NPIPA2 has 0 exact CAT chains and 140 exact RefSeq chains (CAT support 132 reads, RefSeq 287);
  NPIPB5 the reverse (317 vs 46); NPIPB13 has 0 CAT-anchored support reads and 52 RefSeq-anchored ones. The per-copy verdicts should be read
  against both models until the human annotation question is settled (CAT is the default since 2026-10-01).
- The register row 1232 headline (23/21/24 -> 13/9/7) is superseded by 13/9/8 (annotation-anchored); row 1237 records the amendment.

## Amendment B re-score (2026-10-04 12:00, user correction): the exon-intron CHAIN must align; CAT and RefSeq models both count

Amendment A credited a read when >= 2 of its junction coordinates were annotated introns, in any order with anything in between. Amendment B
(prereg 1e2801b0, before this re-score) requires the read's junction chain inside the copy's span to be a **contiguous sub-chain of an
annotated transcript's intron chain** (an incomplete/full splice match: every junction of the read in the span belongs to it, the exon between
two matched introns is spanned exactly), >= 2 junctions, against the union of the CAT/Liftoff v2.0 and RefSeq models of the gene. FOUND = >= 2
such reads AND a node whose representative's in-span junction chain is such a sub-chain. Scorer columns `chain_support_reads` (+ `_ann1` CAT,
`_ann2` RefSeq), `chain_expressed`, `<arm>_chain_found`, `<arm>_locus_chain_found`. The page's per-copy table is `docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv`
(all three rules' columns).

**NPIP: FOUND = 10 / 7 / 7 of 25 (P / GOOD / ALL), within the page's NPIP-cluster nodes 9 / 6 / 6** (Amendment A: 13 / 9 / 8; the first rule:
13 / 9 / 7; own node 23 / 21 / 24). Every copy has >= 2 chain-support reads, but at most copies they are a small minority of the reads:

| copy | reads | chain-support reads: both / CAT / RefSeq | coordinate rule (A) | exact chains | P found / own | GOOD | ALL |
|---|---|---|---|---|---|---|---|
| NPIPB2 | 432 | 216 / 199 / 39 | 348 | 193 | 0 / 0 | 0 / 0 | 0 / 1 |
| NPIPA2 | 348 | 201 / 1 / 201 | 132 | 0 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA1 | 806 | 155 / 155 / 38 | 631 | 10 | 0 / 1 | 0 / 1 | 0 / 1 |
| PKD1P6-NPIPP1 | 253 | 37 / 36 / 32 | 139 | 5 | 1 / 0 | 1 / 0 | 0 / 1 |
| NPIPA5 | 147 | 83 / 9 / 83 | 110 | 0 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA6 | 204 | 18 / 18 / 7 | 142 | 1 | 0 / 1 | 0 / 0 | 0 / 1 |
| NPIPA7 | 285 | 36 / 36 / 36 | 184 | 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA8 | 192 | 2 / 2 / 2 | 140 | 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA9 | 998 | 55 / 55 / 39 | 865 | 13 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC128966608 | 1104 | 338 / 338 / 66 | 585 | 37 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB4 | 897 | 49 / 49 / 31 | 199 | 27 | 0 / 1 | 1 / 1 | 1 / 1 |
| NPIPB5 | 785 | 133 / 132 / 11 | 317 | 1 | 0 / 1 | 0 / 1 | 1 / 1 |
| NPIPB6 | 684 | 206 / 164 / 206 | 563 | 107 | 0 / 1 | 0 / 0 | 0 / 1 |
| NPIPB7 | 164 | 38 / 28 / 15 | 105 | 21 | 1 / 1 | 0 / 1 | 0 / 1 |
| NPIPB8 | 174 | 74 / 74 / 6 | 128 | 0 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB9 | 568 | 102 / 6 / 102 | 433 | 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB10P | 70 | 2 / 2 / 2 | 28 | 0 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB11 | 150 | 6 / 5 / 6 | 82 | 0 | 0 / 1 | 0 / 1 | 1 / 1 |
| NPIPB12 | 54 | 3 / 2 / 3 | 9 | 1 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907834 | 634 | 131 / 72 / 127 | 302 | 0 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB13 | 131 | 23 / 0 / 23 | 0 | 0 | 1 / 1 | 0 / 1 | 1 / 0 |
| NPIPB14P | 1255 | 43 / 43 / 1 | 1023 | 20 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB15 | 220 | 141 / 117 / 141 | 170 | 74 | 1 / 1 | 0 / 1 | 0 / 1 |
| LOC124907808 | 73 | 33 / 19 / 33 | 46 | 5 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907807 | 96 | 33 / 31 / 33 | 41 | 14 | 1 / 1 | 0 / 1 | 0 / 1 |

- **The user's concern is confirmed at the read level.** At NPIPA8, 2 of 192 reads chain-match any annotated model (Amendment A counted 140:
  their junction coordinates are annotated introns, their chains are not an annotated chain); NPIPA9 55 of 998 (A: 865); NPIPB14P 43 of 1,255
  (A: 1,023); NPIPB4 49 of 897; NPIPA6 18 of 204; NPIPB10P 2 of 70; NPIPB12 3 of 54. The reads at these copies splice through annotated
  splice sites in chains the annotation does not have — alternative or novel isoforms, 5'-variable structures, or chains of a sibling copy.
  Only NPIPB2 (216), NPIPB6 (206), NPIPA2 (201, all RefSeq), LOC128966608 (338), NPIPB15 (141), NPIPB5 (133), LOC124907834 (131), NPIPB9 (102,
  RefSeq) have >= 100 reads that are an annotated chain.
- **Which annotation matters, copy by copy:** NPIPA2 1 CAT / 201 RefSeq; NPIPB9 6 / 102; NPIPA5 9 / 83; NPIPB13 0 / 23; NPIPB2 199 / 39;
  NPIPB5 132 / 11; NPIPB14P 43 / 1. Neither annotation alone describes what is expressed at NPIP; the union is used, as registered.
- **Nodes:** the primaries-only arm keeps 9 chain-found copies among its 23 own nodes; "+ all" 6 of 24. NPIPB2 again 0 in every arm (ALL's
  loci are fragments); NPIPB4 0 / 1 / 1; NPIPA8, NPIPA9, NPIPB14P, NPIPB9, NPIPB10P, NPIPB11, NPIPB12, LOC124907808 are own nodes in every arm
  and chain-found in none — the representatives at those copies are not an annotated chain (fragments, or chains with a junction the
  annotation lacks).
- Locus level (any transcript of the locus is an annotated sub-chain): 17 / 16 / 16. The loci still carry annotated chains at two thirds of the
  copies; the representative does not.

Register row 1238. The headline the page now carries: own node 23 / 21 / 24 vs found 9 / 6 / 6.

## Amendment C re-score (2026-10-04 12:25): only reads that TRULY support an expressed intron chain are counted

Amendment C (prereg dd710b37, before the run): the copy's expressed chains are its reads' own identical >= 2-junction chains carried by >= 3
UNIQUELY placed reads (MAPQ > 0); a read truly supports the copy iff its chain equals or is a contiguous sub-chain of an expressed chain; FOUND
= the node's representative is such a sub-chain. No annotation enters the verdict; each expressed chain is classed against CAT ∪ RefSeq for
the record. Columns `xc_*` in `docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv` (the C run carries every rule's columns).

**NPIP: 24 of 25 copies have an expressed chain (NPIPB12 none: 54 reads, no three identical unique chains). FOUND = 14 / 10 / 7 of 25 by any
locus, 13 / 9 / 6 within the page's own nodes (P / GOOD / ALL); locus level 21 / 19 / 17.** Amendment B (annotated chains): 9 / 6 / 6.

| copy | reads | expressed chains | truly-supporting reads (unique / incl. tied) | dominant expressed chain: reads, junctions, class vs CAT ∪ RefSeq | B support | P found / own node | GOOD | ALL |
|---|---|---|---|---|---|---|---|---|
| NPIPB2 | 432 | 19 | 301 / 301 | 167, 5j, FSM | 216 | 0 / 0 | 0 / 0 | 0 / 1 |
| NPIPA2 | 348 | 16 | 236 / 236 | 79, 6j, FSM | 201 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA1 | 806 | 52 | 402 / 402 | 67, 6j, NIC | 155 | 1 / 1 | 1 / 1 | 0 / 1 |
| PKD1P6-NPIPP1 | 253 | 13 | 57 / 57 | 6, 3j, NNC | 37 | 1 / 0 | 1 / 0 | 0 / 1 |
| NPIPA5 | 147 | 6 | 86 / 86 | 38, 6j, ISM | 83 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA6 | 204 | 7 | 47 / 61 | 10, 7j, NIC | 18 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPA7 | 285 | 5 | 43 / 141 | 10, 5j, ISM | 36 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA8 | 192 | 1 | 4 / 100 | 3, 4j, NNC | 2 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA9 | 998 | 48 | 447 / 469 | 30, 6j, NIC | 55 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC128966608 | 1104 | 23 | 320 / 412 | 135, 11j, ISM | 338 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB4 | 897 | 13 | 96 / 101 | 21, 6j, NIC | 49 | 0 / 1 | 1 / 1 | 1 / 1 |
| NPIPB5 | 785 | 19 | 197 / 219 | 108, 2j, ISM | 133 | 0 / 1 | 0 / 1 | 1 / 1 |
| NPIPB6 | 684 | 29 | 511 / 513 | 184, 7j, NNC | 206 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPB7 | 164 | 8 | 69 / 72 | 21, 6j, FSM | 38 | 1 / 1 | 0 / 1 | 0 / 1 |
| NPIPB8 | 174 | 6 | 86 / 92 | 27, 15j, ISM | 74 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB9 | 568 | 22 | 396 / 401 | 194, 7j, NNC | 102 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB10P | 70 | 3 | 24 / 25 | 10, 6j, NNC | 2 | 0 / 1 | 1 / 1 | 0 / 1 |
| NPIPB11 | 150 | 6 | 49 / 49 | 11, 6j, NIC | 6 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB12 | 54 | 0 | 0 / 0 | 1, 8j, NNC | 3 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907834 | 634 | 17 | 256 / 257 | 68, 2j, ISM | 131 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB13 | 131 | 5 | 42 / 42 | 9, 6j, NNC | 23 | 1 / 1 | 0 / 1 | 1 / 0 |
| NPIPB14P | 1255 | 46 | 918 / 918 | 400, 6j, NIC | 43 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB15 | 220 | 4 | 116 / 144 | 70, 6j, FSM | 141 | 1 / 1 | 0 / 1 | 0 / 1 |
| LOC124907808 | 73 | 1 | 5 / 33 | 5, 5j, ISM | 33 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907807 | 96 | 2 | 22 / 31 | 14, 6j, FSM | 33 | 1 / 1 | 0 / 1 | 0 / 1 |

- **The reads do carry consistent intron chains at almost every copy** — NPIPB14P 46 expressed chains and 918 truly-supporting unique reads
  (dominant chain 400 reads, 6 junctions, NIC: annotated splice sites in an unannotated combination), NPIPB6 511 (dominant 184 reads, NNC),
  NPIPA9 447, NPIPA1 402, NPIPB9 396, NPIPB2 301 — and the nodes do not: at NPIPB14P, NPIPB9, NPIPA9, NPIPB2 no arm's representative is a
  sub-chain of any expressed chain. The copy is expressed and structured; what we hand on is a fragment.
- **Where the reads are shared with siblings the unique-read floor exposes it:** NPIPA8 4 unique truly-supporting reads (100 with tied reads),
  NPIPA7 43 / 141, LOC124907808 5 / 33, NPIPA6 47 / 61, NPIPB15 116 / 144 — these copies' chains are also their siblings' chains; the tied
  reads count only through O2.
- Amendment B vs C per copy: B's annotated-chain support is lower than C's everywhere except where the annotation happens to be the expressed
  chain (NPIPA2, NPIPA5, NPIPB15, LOC124907807); at NPIPB14P 43 vs 918, NPIPA9 55 vs 447, NPIPA1 155 vs 402.

Register row 1240. The page's headline is now the Amendment C number (13 / 9 / 6 within own nodes 23 / 21 / 24), with B, A and the first rule
beside.
