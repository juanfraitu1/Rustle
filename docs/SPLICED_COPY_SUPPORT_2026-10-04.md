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

## Amendment D re-score (2026-10-04 12:40, user: "a read that starts at the TSS, has up to 3 junctions and introns in common with the intron chain it is supporting")

Amendment D (prereg 8affcfb2, before the run): a read supports a transcript iff its 5' end lies within ±150 bp of the transcript's TSS
(strand-aware) AND its first m junctions equal the transcript's first m introns, m = min(3, introns); found = the node's representative does
the same. Registered against the annotated models (CAT ∪ RefSeq); the same test against the reads' own expressed chains (Amendment C's chains,
TSS = the modal 5' end of their unique carriers, 20-bp bins) was registered as the reading beside (D'). Columns `td_*` (annotated) and `te_*`
(expressed) in `docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv`.

**Against the annotated TSS (D): 5 of 25 copies have >= 2 supporting reads; FOUND 1 / 1 / 1 (NPIPB8). Against the expressed TSS (D'): 24 of 25
copies supported, FOUND 10 / 6 / 2 within the own nodes 23 / 21 / 24 (any locus 11 / 7 / 2); locus level 20 / 19 / 16.**

| copy | reads | D: from the ANNOTATED TSS (±150; CAT / RefSeq; ±300) | D': from the EXPRESSED TSS (unique) | expressed TSS(s) | P found / own | GOOD | ALL |
|---|---|---|---|---|---|---|---|
| NPIPB2 | 432 | 34 (27 / 34; 35) | 305 (305) | 12012752, 11977729, 11973618 | 0 / 0 | 0 / 0 | 0 / 1 |
| NPIPA2 | 348 | 0 (0 / 0; 0) | 227 (227) | 14741044, 14746189, 14749609 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA1 | 806 | 0 (0 / 0; 0) | 463 (463) | 14938049, 14925111, 14925612 | 1 / 1 | 1 / 1 | 0 / 1 |
| PKD1P6-NPIPP1 | 253 | 0 (0 / 0; 0) | 68 (68) | 15129275, 15128694, 15130861 | 1 / 0 | 1 / 0 | 0 / 1 |
| NPIPA5 | 147 | 1 (1 / 0; 1) | 82 (82) | 15382704, 15382704, 15400736 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPA6 | 204 | 0 (0 / 0; 0) | 71 (70) | 16337769, 16339536, 16339733 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPA7 | 285 | 0 (0 / 0; 0) | 130 (28) | 16342074, 16395878, 16395769 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPA8 | 192 | 0 (0 / 0; 0) | 121 (10) | 18339521 | 1 / 1 | 0 / 1 | 0 / 1 |
| NPIPA9 | 998 | 0 (0 / 0; 2) | 557 (541) | 18386676, 18391557, 18382592 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC128966608 | 1104 | 22 (15 / 7; 24) | 284 (281) | 21685230, 21689577, 21685160 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB4 | 897 | 0 (0 / 0; 0) | 97 (97) | 22326971, 22359378, 22365631 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB5 | 785 | 0 (0 / 0; 0) | 72 (72) | 22759094, 22696417, 22761474 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB6 | 684 | 1 (1 / 0; 1) | 516 (514) | 28637415, 28637415, 28637415 | 1 / 1 | 0 / 0 | 0 / 1 |
| NPIPB7 | 164 | 0 (0 / 0; 1) | 48 (48) | 28771997, 28751524, 28771999 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB8 | 174 | 52 (52 / 0; 73) | 86 (86) | 28903704, 28903501, 28903794 | 1 / 1 | 1 / 1 | 1 / 1 |
| NPIPB9 | 568 | 0 (0 / 0; 0) | 360 (360) | 28932688, 28969226, 28992124 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB10P | 70 | 0 (0 / 0; 0) | 19 (19) | 29037610, 29319354, 29319354 | 0 / 1 | 1 / 1 | 1 / 1 |
| NPIPB11 | 150 | 0 (0 / 0; 0) | 44 (44) | 29679826, 29666321, 29679826 | 1 / 1 | 1 / 1 | 0 / 1 |
| NPIPB12 | 54 | 1 (0 / 1; 1) | 0 (0) |  | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907834 | 634 | 4 (0 / 4; 4) | 170 (170) | 30523997, 30519740, 30511567 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB13 | 131 | 0 (0 / 0; 0) | 25 (25) | 30626433, 30625966, 30621808 | 0 / 1 | 0 / 1 | 0 / 0 |
| NPIPB14P | 1255 | 0 (0 / 0; 0) | 729 (729) | 75875797, 75799850, 75843724 | 0 / 1 | 0 / 1 | 0 / 1 |
| NPIPB15 | 220 | 0 (0 / 0; 0) | 138 (129) | 80195303, 80195303, 80195303 | 1 / 1 | 0 / 1 | 0 / 1 |
| LOC124907808 | 73 | 2 (2 / 0; 2) | 6 (4) | 80308465 | 0 / 1 | 0 / 1 | 0 / 1 |
| LOC124907807 | 96 | 0 (0 / 0; 0) | 24 (21) | 80424051, 80432162 | 1 / 1 | 0 / 1 | 0 / 1 |

- **Why the annotated TSS fails: the models' first exons.** The CAT NPIP models have first "exons" of 12-53 kb (NPIPB2 47,479 bp; LOC128966608
  52,989; NPIPB4 41,574; NPIPA1 25,831; NPIPA9 25,056; NPIPB5 33,609), RefSeq's at several copies too (NPIPB9 34,034; NPIPB7 27,579; NPIPB6 19,184);
  the models whose intron chain the full-length reads carry (their FSM match) place the TSS 3.7-51.5 kb upstream of where every read starts —
  while the reads' 5' ends cluster to the base (NPIPB2: 196 FSM reads start at one position 35,027 bp downstream of that model's TSS, IQR
  35,027-35,032; NPIPA2 8,676 [8,676, 8,676]; NPIPB6 6,464; NPIPB14P 5,781). The reads start at a consistent transcription start; the annotation's
  first exon absorbs the upstream region. Post hoc; the per-copy first-exon lengths are in the session record and the `te_tss` column lists
  the expressed starts (e.g. NPIPB2 12,012,752; NPIPB6 28,637,415; NPIPB15 80,195,303).
- **Under the user's standard the nodes fall further:** the primaries-only arm's representative starts at the expressed TSS with the first
  three introns at 10 of its 23 own nodes, "+ good" 6 of 21, "+ all" 2 of 24. The loci hold such a transcript at 20 / 19 / 16.
- Sibling-shared copies under the unique-read floor: NPIPA7 130 supporting reads but 28 unique, NPIPA8 121 / 10, LOC124907808 6 / 4.

Register row 1241. The page's headline is D' (10 / 6 / 2), with D (1 / 1 / 1), C, B, A and the first rule beside.

## Do we know where NPIP's 5' end is? CAT vs RefSeq vs the reads (2026-10-04 13:50, user question)

Per copy: the CAT and RefSeq model spans, CAT's 5' extension beyond RefSeq (strand-aware; positive = CAT reaches further upstream), and the
distance from the reads' modal transcription start (Amendment D', the dominant expressed chain's unique carriers) to the nearest CAT and
RefSeq TSS. Session table (also the basis of the page's tooltips):

| copy | CAT span | RefSeq span | CAT 5' extension | reads' TSS -> CAT (bp) | -> RefSeq (bp) |
|---|---|---|---|---|---|
| NPIPB2 | 49,450 | 49,446 | 0 | 4 | 4 |
| NPIPA2 | 17,354 | 22,936 | -5,533 | 5,422 | 111 |
| NPIPA1 | 31,910 | 14,597 | +17,309 | 16,849 | 460 |
| PKD1P6-NPIPP1 | 33,210 | 35,998 | +17,605 | 30,136 | 12,531 |
| NPIPA5 | 17,420 | 18,025 | -601 | 3,139 | 3,740 |
| NPIPA6 | 36,006 | 18,742 | +17,282 | 14,771 | 2,511 |
| NPIPA7 | 14,894 | 14,827 | +85 | 49,204 | 49,289 |
| NPIPA8 | 18,795 | 18,818 | 0 | 4,458 | 4,458 |
| NPIPA9 | 31,193 | 18,745 | +12,448 | 16,888 | 4,440 |
| LOC128966608 | 57,390 | 22,731 | +34,678 | 52,864 | 18,186 |
| NPIPB4 | 46,003 | 22,541 | +23,449 | 9,401 | 32,850 |
| NPIPB5 | 43,879 | 32,704 | +13,262 | 9,250 | 22,512 |
| NPIPB6 | 20,809 | 22,119 | -1,272 | 6,464 | 7,736 |
| NPIPB7 | 20,987 | 41,603 | -20,573 | 13,959 | 6,614 |
| NPIPB8 | 35,746 | 10,633 | +25,174 | 54 | 25,120 |
| NPIPB9 | 21,030 | 37,372 | -16,285 | 99,677 | 83,392 |
| NPIPB10P | 14,083 | 14,485 | -821 | 282,178 | 281,357 |
| NPIPB11 | 22,428 | 25,154 | -2,682 | 5,996 | 8,678 |
| NPIPB12 | 22,608 | 23,125 | -409 | – (no expressed chain) | – |
| LOC124907834 | 28,891 | 22,360 | +8,932 | 14,728 | 5,796 |
| NPIPB13 | – (no CAT gene) | 25,274 | – | – | 8,237 |
| NPIPB14P | 19,826 | 89,796 | -69,970 | 70,166 | 196 |
| NPIPB15 | 14,206 | 15,802 | -1,561 | 341 | 1,220 |
| LOC124907808 | 14,188 | 15,984 | -1,761 | 1,572 | 189 |
| LOC124907807 | 14,169 | 17,110 | -2,906 | 336 | 2,570 |

- **CAT does not lengthen NPIP systematically.** It reaches further 5' than RefSeq at 9 of 24 copies, by 9-35 kb (NPIPA1, PKD1P6-NPIPP1, NPIPA6,
  NPIPA9, LOC128966608, NPIPB4, NPIPB5, NPIPB8, LOC124907834 — the copies with 12-53 kb first exons); RefSeq reaches further at 6 (NPIPB14P by
  70 kb, NPIPB7 21 kb, NPIPB9 16 kb, NPIPA2 5.5 kb, NPIPB11, LOC124907807); the median difference is -204 bp. The two annotations disagree
  with each other in both directions.
- **Neither annotation's TSS is where the reads start:** the reads' modal start is within 150 bp of a CAT TSS at 2 of 23 copies and of a
  RefSeq TSS at 2 of 24 (within 1 kb: 4 and 5).
- **The library does reach transcription starts elsewhere:** of the 660 chr16 CAT genes with >= 10 reads and >= 3 introns, 586 (89%) have
  >= 2 reads starting within 150 bp of the CAT TSS and carrying its first three introns (Amendment D as registered); NPIP 5 of 25. The NPIP
  discrepancy is annotation, not 5' truncation.
- **So: the annotated 5' ends of NPIP are not known, and the reads' own starts are the only defensible estimate** — tight to the base at 24
  of 25 copies (NPIPB2 196 full-chain reads at one position; NPIPA2, NPIPB6, NPIPB14P likewise). Caveat on three copies: the dominant
  expressed chain starts 49-282 kb away from the copy (NPIPA7 49 kb, NPIPB9 100 kb, NPIPB10P 282 kb) — transcription units that run through
  a neighbouring copy or gene (NPIPA6 -> NPIPA7 is the one the read-pool page saw as a merged locus); their TSS belongs to the neighbour, and
  "the copy's 5' end" is not a well-posed question for them without the family context.

Register row 1242.

## The cap signal marks where NPIP transcription really starts (2026-10-04 13:55)

The A119b library carries the template-switching cap signature (`--polish-tss`'s CAP signal, `docs/ASSEMBLY_POLISH.md` addendum 2): a 1-3 bp
untemplated G at the RNA 5' end, present on molecules that were reverse-transcribed to the cap. Per copy, the fraction of same-strand primaries
carrying it, at the reads' modal start (±20 bp, Amendment D') and >200 bp away from it:

| locus | all reads | at the modal start | > 200 bp away |
|---|---|---|---|
| VPS4A (control) | 407/688 (59%) | 317/420 (75%) | 7/134 (5%) |
| NPIPB2 | 206/436 (47%) | 189/306 (62%) | 1/100 (1%) |
| NPIPA2 | 218/366 (60%) | 192/216 (89%) | 5/118 (4%) |
| NPIPA5 | 85/153 (56%) | 77/86 (90%) | 2/60 (3%) |
| NPIPA8 | 125/192 (65%) | 114/134 (85%) | 0/46 (0%) |
| NPIPA9 | 102/999 (10%) | 76/112 (68%) | 23/809 (3%) |
| NPIPB4 | 80/901 (9%) | 59/71 (83%) | 21/817 (3%) |
| NPIPB6 | 422/693 (61%) | 410/501 (82%) | 7/185 (4%) |
| NPIPB7 | 48/168 (29%) | 29/45 (64%) | 11/113 (10%) |
| NPIPB8 | 99/198 (50%) | 35/45 (78%) | 28/108 (26%) |
| NPIPB10P | 28/92 (30%) | 25/26 (96%) | 2/64 (3%) |
| NPIPB11 | 58/181 (32%) | 50/56 (89%) | 5/117 (4%) |
| NPIPB13 | 20/255 (8%) | 13/17 (76%) | 7/230 (3%) |
| NPIPB14P | 602/1264 (48%) | 346/425 (81%) | 160/664 (24%) |
| NPIPB15 | 133/231 (58%) | 127/147 (86%) | 5/83 (6%) |
| LOC124907834 | 122/644 (19%) | 106/130 (82%) | 15/467 (3%) |
| LOC124907807 | 31/106 (29%) | 24/25 (96%) | 6/77 (8%) |
| NPIPA1 | 177/807 (22%) | 0/10 (0%) | 177/774 (23%) |
| NPIPA6 | 22/204 (11%) | 0/22 (0%) | 22/139 (16%) |
| NPIPA7 | 109/285 (38%) | 0/10 (0%) | 109/265 (41%) |
| LOC128966608 | 116/1111 (10%) | 0/14 (0%) | 114/1074 (11%) |
| NPIPB5 | 46/786 (6%) | 3/19 (16%) | 37/732 (5%) |
| NPIPB9 | 321/624 (51%) | 6/8 (75%) | 307/605 (51%) |
| LOC124907808 | 22/83 (27%) | 1/8 (12%) | 20/74 (27%) |
| PKD1P6-NPIPP1 | 9/255 (4%) | 0/6 (0%) | 8/215 (4%) |

- **At 17 of 24 copies the reads' modal start is a capped transcription start** (62-96% of the reads there carry the cap G, against 0-10%
  of the reads starting elsewhere; the control gene 75% vs 5%). The 5' ends of NPIP are knowable from the data: they are where the capped
  molecules begin, 3.7-51 kb downstream of the CAT first exons.
- **At 7 copies the modal start picked by Amendment D' is not the capped start** (NPIPA1, NPIPA6, NPIPA7, LOC128966608, NPIPB5, LOC124907808
  with 0-16% capped there while 11-41% of the other reads are capped; NPIPB9 has 51% capped reads away from the one start chosen): D' took the
  first expressed chain's modal 5' end, which at these copies is a downstream fragment start or the start of a chain running in from a
  neighbour. PKD1P6-NPIPP1 has almost no capped reads at all (9/255).
- **So the fix is to define the start from the capped reads, not from the modal 5' end:** cluster the 5' ends of cap-clipped reads per locus,
  call each cluster of >= 3 a transcription start, and anchor support, the representative and the evaluation there (Amendment E below).

Register row 1243.
