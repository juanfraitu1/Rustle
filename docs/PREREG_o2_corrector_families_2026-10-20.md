# Pre-registration: does the sample-aware corrector recover a share of minimap2's errors across many families? (2026-10-20)

**User (2026-10-20):** after the vault survey found the sample-aware corrector to be the only O2 method with held-out gains:
"ok can we try to prove that it works in other families, maybe not DAZ but as long as it can recover a small percent.."

## The rule (frozen; nothing changes)

`bench/o3_graph/o2correct.py donor`, version 2, unchanged since bfa15f12, default flags (no `--v1`, no `--guard`):
- **Domain:** a read whose alignments within 10 AS of its best (Eichler) lie on two or more copies of one family and on no
  locus off it, minimap2's primary among them.
- **The sample's copies:** per copy, its confident reads (primary on it, every alignment on another copy of its family more
  than 9 AS lower) give its events (substitutions, left-shifted indels): PSVs as this sample carries them. No phasing, no
  haplotypes, no assembly.
- **Decision:** the tied copy whose alignment has the fewest events its sample reads do not explain takes the read, and
  overrules minimap2 only when it dominates; otherwise minimap2's pick stays.

## Substrates, none ever seen by this rule

- **D (primary): a never-used donor, simulated from its own assembly.** HG03516 (HPRC year 1, Esan, Nigeria; hifiasm
  phased `year1_f1_assembly_v2_genbank`, maternal and paternal), never mentioned in this project. The recipe is register
  1346's, unchanged:
  - every CHM13 copy of the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`), ± 20 kb, aligned to each
    haplotype (`minimap2 -x asm20 -c --secondary=yes -N 10`);
  - `donorsim.py orthologs`, then `donorsim.py build ... A119b_testis_v4_ref.json --per-copy 150 --seed 3516`: reads
    simulated from the donor's own transcripts, so its alleles, indels and missing copies are real;
  - aligned to CHM13 v2.0 (`minimap2 -a -x splice:hq -uf --secondary=yes -N 100 -p 0.5`).
  - Truth: the CHM13 copy whose donor ortholog produced the read.
- **A (replication, real reads): a third A119b fragment cut.** `o2frag.py sources` with `--seed 3`, excluding every
  source read of seeds 1 and 2; fragments `--per-read 3 --min 300 --max 1200 --seed 3`; truth = the source read's copy,
  kept when clear by 30 AS (register 1343's design and known bias). The corrector learns from the whole library
  (`--teach A119b.t2t.bam`), each fragment's source read left out.
- **N (no harm, one genome):** simulation round 27's six CROSS worlds (`fast/heldout27/*/CROSS`; INTACT and SILENT repeat
  their reads), never run through the corrector.

## Scoring

`o2correct.py score OUT.tsv TRUTH.tsv --among-candidates --families qualify_nonpar.json` (D also `--donor`): a read is
scored when its truth is among its tied copies. Per substrate: scored reads, minimap2 right, corrector right, fixes
(corrector right where minimap2 was wrong), breaks (the reverse), and per family whether it fixes more than it breaks.
N: `o2correct.py score-worlds` on the six CROSS worlds, PAR families excluded.

## Bars (D and A each)

- **P (power):** at least 300 scored reads from at least 11 families, minimap2 wrong on at least 50.
- **R (recovers a share):** fixes − breaks >= 2% of minimap2's wrong calls, and >= 20 reads.
- **S (safety):** fixes >= 5 × breaks.
- **F (families):** at least 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up)
  where it breaks more.

**N (no harm):** on the one-genome worlds, breaks <= fixes + 2.

**Verdict:** **D PASS** = P, R, S, F on D, and N. FAIL = P holds and a bar fails (or N fails). UNDERPOWERED = P fails.
A gets its own verdict on P, R, S, F; a claim across a real library needs A to PASS as well.

**What a PASS would claim:** for reads tied within 10 AS between near-identical copies, reading each copy's PSVs from the
sample's own unambiguous reads recovers a share of minimap2's wrong placements, many times more than it breaks, in more
than 10 families of a donor it never saw. **It would not claim:** reads identical in the donor too (they stay minimap2's
coin toss: no fixes, no breaks), copies the reference lacks (O3), the assembly step, or anything past 10 AS.

## Rules

- Each substrate is built and scored once; bugs → amendment + rerun.
- Heavy jobs (the donor download and alignments, the CHM13 splice alignments, the library pass) in the foreground, one at
  a time, each command under 10 minutes; the corrector's library pass resumes from its cache.
- The metric-auditor checks the verdicts on a frozen copy before they are reported. Human numbers only.
- Push only when asked.

## Amendments

### Amendment 1 (2026-10-20, before any data): gorilla KB3781 only — replaces the substrates, bars D / A / N

**User:** "wait no, I think this is one of the over-reaching problems my advisor complains, can we stick with gorilla
fibroblasts and its T2T genome?"

- **Withdrawn, unrun:** D (HG03516), A (A119b seed 3), N (round 27). The interrupted command had started HG03516's
  maternal download; the 84 MB partial file was deleted unread. Nothing else was built.
- **The substrate (G2):** KB3781's fibroblast Iso-Seq on the GGO primary (`GCA_029281585.2_flnc_mm.bam`), in **the 383
  GWFAM families never scored** (`rna_allele/fibro_cat.copies.tsv` less the 283 of `personal/qual.json`, used by registers
  1347-1350). Every alignment of every read touching their copies is cut from the library, as `personal/sub.bam` was.
  DAZ-like reads (identical on two copies in his own genome too) are undetermined and not scored.
- **The rule:** the frozen corrector, as above, with that family table.
- **The answer key (evaluation only; the corrector never sees it):** each domain read aligned to his maternal and paternal
  assemblies; its copy is the one whose alignment has the best intron-free score (`ms`), at least one substitution (5
  points) ahead of any other copy or locus (`o2window.py truth-ms`), through each copy's ALT ortholog (`o2personal.py
  queries` / `orthologs`, unchanged). **refcons** = his reference-haplotype copy is among the best (his copy equals the
  reference); **altonly** = only the other haplotype's.
- **Scoring:** the frozen `o2correct.py score --among-candidates` on the determined reads (truth as `read, copy`); then the
  same on the refcons reads alone.
- **Bars (G2):** P (>= 300 scored reads from >= 11 families, minimap2 wrong on >= 50), R (fixes − breaks >= 2% of
  minimap2's wrong calls and >= 20), S (fixes >= 5 × breaks), F (>= 11 families fixing more than breaking, <= a fifth as
  many rounded up the reverse), **H (no harm where his copy equals the reference: on refcons reads, breaks <= fixes + 2)**.
  PASS = all five; FAIL = P holds and another fails; UNDERPOWERED = P fails.
- **Known caveats, stated now:** on refcons reads the key is the reference's own ordering (near self-scoring; H is a harm
  check, not a gain); in the 283 families tested before, the corrector's gorilla gains were concentrated (r1347: 82 fixes,
  76 in SM392).

### Amendment 2 (2026-10-20, before any corrector run on these families): the thesis's mainstay families only

**User:** "I don't need the GWFAM, my advisor is very picky and he does not believe any of the method works, we need to stick
to the YAGs or mainstays that we have tried earlier, unless those families don't have nearly AS tied multi-mapping reads";
on the census below and the proposed test: "Yes".

**Census (substrate only; no rule run):** KB3781's whole fibroblast library, reads with a primary on a copy / reads whose
alignments on two or more of the family's copies lie within 10 AS:
- YAGs on his Y (RBMY 13, TSPY 6, HSFY 5, BPY2 2, DAZ 2, CDY 1 copies): **0 reads** (testis-only genes): untestable here.
- TBC1D3 (14 copies) 5 / 3; SMN (RefSeq, 3) 1,478 / 0; GOLGA6 (2) 0 / 0; FAM90A (1) 4 / 0.
- **NPIP** (the curated 25-copy set, `fibro_cat` family `NPIP`): 1,040 / **182** (44 exact).
- **LRPAP1** (its 22-copy family, `fibro_cat` `SM22`; the 10-04 case-study copies): 2,555 / **38**.
- SM21 (204 tied) is left out: genes riding in the NPIP duplication (PLA2G10, PARN, SRRM2-like), not a mainstay.

**Replaces Amendment 1's family set and bars.** Families: NPIP and LRPAP1 (47 copies; neither among the 283 scored in
registers 1347-1350). Everything else in Amendment 1 holds (his library, the frozen corrector, his two assemblies as the
answer key by the intron-free score, refcons / altonly).

**Bars (pooled over the two families; each family reported):**
- **P:** >= 50 scored reads and minimap2 wrong on >= 10; else UNDERPOWERED.
- **R:** fixes > breaks, and fixes − breaks >= 2% of minimap2's wrong calls.
- **S:** fixes >= 5 × breaks.
- **H:** on refcons reads, breaks <= fixes + 2.

PASS = P, R, S, H; FAIL = P holds and another fails.

## Results (code 67383d28; run once; `o3_graph/mainstay/`): **UNDERPOWERED**

KB3781's fibroblast reads touching NPIP's 25 and LRPAP1's 22 copies: 5,821 reads, 77,090 alignments. The corrector's
domain (tied within 10 AS): **185 reads** (NPIP 164, LRPAP1 21); it learns from 3,037 confident reads on 33 copies and
overrules minimap2 on 4.

**What his genome says about the 185** (each copy's counterpart on his other haplotype lifted for 20 of the 47 copies):

| family | tied | determined, his copy = the reference | determined, other parent's copy | identical on two copies in his genome | no traceable copy |
|---|---|---|---|---|---|
| NPIP | 164 | 75 (1 outside the tie) | 6 | 48 | 35 |
| LRPAP1 | 21 | 0 | 1 | 2 | 18 |

**Scored (truth among the tied copies): 81 reads; minimap2 right 79, corrector right 79; fixes 0, breaks 0.**
- P fails (minimap2 wrong on 2, bar 10): **UNDERPOWERED.** R and S cannot be met with 0 fixes; H holds (refcons 74: 0
  breaks).
- minimap2's 2 errors are NPIP reads (one his reference copy NPIP.14 placed on NPIP.16; one his other parent's NPIP.18
  placed on NPIP.17); the corrector kept minimap2's pick on both.
- The corrector's 4 overrules are all LRPAP1 reads whose best place in his genome has no traceable reference copy, so
  the key cannot score them. Descriptive look (not a registered result): two are reads whose best match is his maternal
  copy of SM22.1, minimap2's pick, by one substitution (the overrule would be a break); two leave a copy his genome
  scores 5 and 30 points worse for one of two equally best copies.

**Reading.** In the matched sample, for the thesis's mainstay families, minimap2 is already right on 79 of the 81 tied
reads his own genome can place. The rest are reads identical on two copies in his genome too (50: no method can place
them) or in duplications his assemblies cannot trace (53). There is almost nothing left for O2 to recover here. YAGs are
not expressed in fibroblasts and cannot be tested in this sample. Audit pending.

## Audit (metric-auditor, 2026-10-21): UNDERPOWERED holds with caveats; the reading above is withdrawn

Every number reproduces (0 status mismatches over 185 reads). UNDERPOWERED survives every orthology-respecting repair of
the lift once a key artefact is removed (minimap2 errors 4-9 against the bar of 10). Corrections:

- **Self-scored majority:** 74 of the 81 scored reads are refcons, keyed by the reference's own best alignment (73 / 74).
  The only informative reads are the **7 from his other parent's copy: minimap2 right on 6.**
- **Key artefact:** one of minimap2's two "errors" (SRR27178663.1441099) and the "outside the tie" read are truths built on
  stitched alignments (272-418 kb gaps; `truth_ms` never applied `stitched()`). **Genuine scored errors: 1**
  (SRR27438212.7310158, NPIP.17 vs NPIP.18).
- **The 53 off-copy reads are outcome-linked:** all have minimap2's copy >= 5 ms below their best placement (median 262).
  - **29 NPIP reads come from one maternal-only locus, mat chr18 (CM054600.2) 6.72-6.82 Mb**: 0-3 substitutions there,
    18-167 at their best paternal placement. A copy the reference lacks (O3's case).
  - 24 are traceable but the lift refused their copies: 11 NPIP copies tripped the 50-bp exon guard (MAPQ 31-60,
    99.77-99.97% identity), 11 SM22 dropped as ambiguous, 6 SM22 refused as interchromosomal. Repaired lifts give 105-113
    scored reads and 4-9 genuine minimap2 errors.
- **The corrector could not act on 66 of 185 reads** (minimap2's copy has no confident reads, so dominance never passes),
  including the one genuine error. NPIP.3's sample carries 1.35% events, consistent with a missing copy's reads teaching
  their nearest copy foreign alleles (unverified).
- **The 4 overrules (LRPAP1):** two are exact reference ties where the maternal locus is 99.83% / 99.94% identical to the
  two candidates (undecidable); in the other two his genome scores the corrector's copy 4-5 ms below minimap2's. The
  descriptive "two likely breaks" and "one of two equally best" were wrong.
- "50 identical" is 29 exact ties, 5 within 1-4 ms, 16 tied with an unlifted locus. "YAGs 0 reads" is 0 primary
  alignments (107 secondary, from 51 reads at 2.7-16.7% divergence). The census counts were not saved and reproduce
  only approximately (NPIP 175-180 tied; LRPAP1 22-61).
- Wording not to use: "minimap2 is right on 79 of 81"; "almost nothing left to recover"; "53 his assemblies cannot trace";
  "no method can place the 50"; "H holds" as safety; "minimap2 wrong on 2" without the artefact.

**What stands:** too few minimap2 errors among NPIP / LRPAP1 tied reads that his genome places on a reference copy to test
the corrector (1 genuine; 4-9 under repaired lifts). **What the data does show:** the largest group of minimap2's
mis-placements here is 29 NPIP reads from a maternal NPIP copy absent from the reference — O3, not a choice among
reference copies.

### Amendment 3 (2026-10-21, before rescoring): the repaired answer key

**User:** "lets try 1" (fix the answer key and rescore the frozen corrector output; no new alignments).

**Already seen, stated plainly:** the r1360 audit's exploratory bracket. Under several lift repairs the scored set was
105-113 reads with 4-9 genuine minimap2 errors (UNDERPOWERED in each). This rescore applies one convention, chosen by
principle and fixed here; it is a confirmation, not a blind test. **Expected: UNDERPOWERED.**

**The repaired key** (`--repaired` on `o2personal.py orthologs` and `o2window.py truth-ms`; without it both reproduce
r1360's files byte for byte):
1. **Key alignments:** an alignment of a read to his assemblies is no placement when it is stitched across two copies'
   exons (`stitched()`, >= 50 bp on each) or bridges a gap of >= 100 kb (`no_placement`, the audit's recipe).
2. **Counterparts (identity convention):** a copy's counterpart on his other haplotype is the primary alignment (MAPQ
   >= 20) of its ± 20 kb query **on any contig that is not the reference's own** — the homologous chromosome or another;
   its exons are lifted without the 50-bp growth guard.
3. **Overlaps:** when two copies' counterparts overlap, the copy whose query aligned with the lower divergence (`de`)
   keeps it; the other gets none (instead of both dropped).

Unchanged: the frozen corrector output (`arm.tsv`), the scorer, the bars (P: minimap2 wrong on >= 10 of >= 50 scored; R,
S, H) and the verdict rules.
