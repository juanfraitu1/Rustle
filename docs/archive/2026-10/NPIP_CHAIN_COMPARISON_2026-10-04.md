# Is the annotation what is expressed at NPIP? Chain-level comparison of the reads with the CAT and RefSeq models at the 25 chr16 copies — 2026-10-04

Prompted by the user after Amendment B of `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` (12:15: "this is a big issue, let's do the chain-level
comparison and ensure that only reads that truly support the intron chain are counted"). Descriptive, no decision rule; what is reported was
fixed in `bench/npip_chains.py`'s docstring before the run. Data: A119b primaries (`-F 2308`, same strand) on the exon union of the CAT ∪
RefSeq models of each copy; junction = `N` >= 50 bp, exact coordinates; a read's chain = its junctions inside the span, in order. Classes
against a model set: **FSM** = equals a model's intron chain; **ISM** = contiguous sub-chain of one (FSM + ISM = Amendment B's support);
**NIC** = every splice site annotated, chain not a sub-chain of any model; **NNC** = at least one unannotated splice site; 1-junction and
unspliced reads apart. Models: CAT/Liftoff v2.0 (`copy_recovery_tools_cat/ann/truth.hsa.gtf`) and RefSeq (`copy_recovery_tools/ann/truth.hsa.gtf`),
matched by gene name; the "union" class uses both. Tables: `docs/NPIP_CHAIN_COMPARISON_copies.tsv`, `docs/NPIP_CHAIN_COMPARISON_chains.tsv`
(every chain with >= 2 reads: class, novel junctions, recurrence). Work dir `/mnt/linuxdisk/tmp/readpool_npip/chains_hsa.*`.

## Headline

- **14,115 reads at the 25 copies: 693 FSM (4.9%), 1,413 ISM (10.0%), 1,475 NIC (10.5%), 4,218 NNC (29.9%), 1,331 one-junction (9.4%),
  4,985 unspliced (35.3%).** Of the 7,799 reads with >= 2 junctions, **27% are an annotated chain, 19% recombine annotated splice sites, 54%
  use at least one splice site neither annotation has.**
- **The dominant read chain is an annotated chain at 8 of 25 copies** (FSM/ISM vs CAT or RefSeq), NIC at 6, NNC at 7, a single junction at 4.
  Copies whose most frequent multi-junction chain is NNC with recurrent novel junctions (>= 3 reads each): PKD1P6-NPIPP1 (64 reads), NPIPA7
  (50), NPIPA8 (45), NPIPB6 (184), NPIPB7 (13), NPIPB9 (194), NPIPB13 (8).
- **The NNC reads are, at most copies, uniquely placed and reference-like:** pooled MAPQ-0 fraction 11.8% (FSM 2.1%, ISM 15.8%, NIC 2.0%),
  median divergence `de` 0.0039 (FSM 0.0034). At NPIPB2, NPIPA2, NPIPA1, PKD1P6-NPIPP1, NPIPB6, NPIPB11, NPIPB14P the NNC reads are 0-1% MAPQ 0:
  they are this locus's reads, spliced in chains the annotation does not have — unannotated isoforms (or annotation errors), not coin-toss
  placements. The exceptions are the copies whose siblings share the structure: NPIPA7 (NNC 77% MAPQ 0), NPIPA8 (81%), LOC124907808 (63%),
  NPIPB15 (45%), NPIPA6 (40%) — there a chain cannot say which copy a read came from, only PSVs can (O2).
- **Which annotation fits is copy-specific:** NPIPB2's dominant chain (167 reads) is a CAT FSM and a RefSeq NNC; NPIPA2's (79) a RefSeq FSM
  and a CAT NNC; NPIPA5's an ISM of RefSeq only; LOC124907807's an FSM of both. CAT carries 1-21 models per copy, RefSeq 1-14; neither covers
  what is expressed.

## Per copy

| copy | reads | FSM | ISM | NIC | NNC | 1-junction | unspliced | annotated-chain share of multi-junction reads | dominant chain: reads / class vs CAT / vs RefSeq / novel junctions recurrent | NNC reads MAPQ 0 |
|---|---|---|---|---|---|---|---|---|---|---|
| NPIPB2 | 432 | 174 | 42 | 40 | 101 | 35 | 40 | 61% | 167 / FSM / NNC / 0 | 0.0 |
| NPIPA2 | 351 | 140 | 38 | 41 | 79 | 16 | 37 | 60% | 79 / NNC / FSM / 0 | 0.0 |
| NPIPA1 | 806 | 10 | 145 | 197 | 297 | 52 | 105 | 24% | 67 / NIC / NNC / 0 | 0.003 |
| PKD1P6-NPIPP1 | 504 | 4 | 62 | 58 | 258 | 52 | 70 | 17% | 64 / NNC / NNC / 1 | 0.0 |
| NPIPA5 | 148 | 30 | 53 | 12 | 18 | 14 | 21 | 73% | 38 / NNC / ISM / 0 | 0.0 |
| NPIPA6 | 204 | 1 | 17 | 70 | 58 | 14 | 44 | 12% | 10 / NIC / NNC / 0 | 0.397 |
| NPIPA7 | 285 | 0 | 36 | 6 | 152 | 31 | 60 | 19% | 50 / NNC / NNC / 1 | 0.774 |
| NPIPA8 | 192 | 0 | 2 | 1 | 146 | 9 | 34 | 1% | 45 / NNC / NNC / 1 | 0.806 |
| NPIPA9 | 998 | 13 | 42 | 389 | 438 | 50 | 66 | 6% | 33 / 1J:ISM / 1J:ISM / 0 | 0.11 |
| LOC128966608 | 1104 | 37 | 301 | 56 | 234 | 126 | 350 | 54% | 155 / ISM / NNC / 0 | 0.22 |
| NPIPB4 | 897 | 22 | 27 | 35 | 166 | 214 | 433 | 20% | 21 / NIC / NNC / 0 | 0.146 |
| NPIPB5 | 800 | 2 | 128 | 17 | 241 | 221 | 191 | 34% | 121 / ISM / NNC / 0 | 0.136 |
| NPIPB6 | 690 | 131 | 75 | 12 | 351 | 49 | 72 | 36% | 184 / NNC / NNC / 1 | 0.008 |
| NPIPB7 | 234 | 11 | 9 | 0 | 99 | 46 | 69 | 17% | 13 / NNC / NNC / 1 | 0.021 |
| NPIPB8 | 174 | 0 | 74 | 31 | 36 | 18 | 15 | 52% | 27 / ISM / NNC / 0 | 0.119 |
| NPIPB9 | 3328 | 0 | 95 | 12 | 390 | 71 | 2760 | 19% | 194 / NNC / NNC / 1 | 0.027 |
| NPIPB10P | 70 | 0 | 2 | 0 | 36 | 15 | 17 | 5% | 10 / 1J:ISM / 1J:ISM / 0 | 0.051 |
| NPIPB11 | 156 | 1 | 5 | 16 | 79 | 16 | 39 | 6% | 11 / NIC / NNC / 0 | 0.0 |
| NPIPB12 | 54 | 2 | 1 | 2 | 8 | 24 | 17 | 23% | 7 / 1J:ISM / 1J:NNC / 0 | 0.28 |
| LOC124907834 | 701 | 4 | 106 | 85 | 133 | 104 | 269 | 34% | 52 / NNC / NNC / 0 | 0.107 |
| NPIPB13 | 173 | 1 | 16 | 0 | 45 | 34 | 77 | 27% | 8 / noann / NNC / 1 | 0.067 |
| NPIPB14P | 1411 | 17 | 52 | 388 | 770 | 68 | 116 | 6% | 72 / NIC / NNC / 0 | 0.002 |
| NPIPB15 | 224 | 74 | 60 | 7 | 31 | 24 | 28 | 78% | 71 / FSM / ISM / 0 | 0.447 |
| LOC124907808 | 78 | 6 | 13 | 0 | 28 | 11 | 20 | 40% | 9 / 1J:ISM / 1J:ISM / 0 | 0.633 |
| LOC124907807 | 101 | 13 | 12 | 0 | 24 | 17 | 35 | 51% | 12 / FSM / FSM / 0 | 0.281 |


## What this means for the counting rule and for O1

1. **The counting rule stands as amended (Amendment B):** a read supports a transcript only when its chain is a sub-chain of that transcript's
   chain; a node is found only when its representative is such a chain. `bench/copy_support.py` implements it; it is the FOUND verdict on the
   NPIP page (9 / 6 / 6) and in the representative-rule tables.
2. **But at NPIP the annotation is not the transcript to support.** Against CAT ∪ RefSeq only 27% of the multi-junction reads are annotated
   chains, and at 13 of 25 copies the dominant expressed chain is NIC or NNC with reference-level divergence and unique placement. Scoring
   "found copies" against these models measures agreement with an incomplete annotation as much as it measures our nodes.
3. **The reads' own dominant chains are the expressed transcripts** (the first rule's instinct, at chain level instead of junction level): an
   expressed chain = a read chain carried by >= N reads (N = 3, the junction floor; the exact count is a parameter to pre-register), its
   sub-chains inherit support. A node "truly supports" copy X when its representative is a sub-chain of an expressed chain of X; a read
   "truly supports" X when its chain is. That is annotation-free, chain-level (no fragment that merely overlaps counts), and it is what the
   de novo assembler already enforces per transcript (floor 2 per exact chain, strict junctions). Pre-registering it as Amendment C, with the
   annotated-chain reading reported beside, is the next step; the representative question is then asked against expressed chains, not models.
4. **Where siblings share structure (NPIPA6-A9, LOC124907808, NPIPB15), chains do not separate copies.** Every read there is a chain of
   several copies at once; only PSV evidence (O2's certificate) places it. The chain rule should therefore count a tied read for a copy only
   when O2 assigns it there, or report it as ambiguous support.

Register row 1239.
