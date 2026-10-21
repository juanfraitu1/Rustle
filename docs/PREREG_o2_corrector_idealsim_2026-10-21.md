# Pre-registration: the frozen corrector on ideal-expression reads simulated from both of KB3781's haplotypes (2026-10-21)

**User (2026-10-21):** after register 1361 (his real fibroblast reads: UNDERPOWERED, minimap2 wrong on 3 of 107 tied reads),
asked whether reads simulated at optimal expression from every known gene would work better; on the proposal to simulate
from his non-reference haplotype: "ok lets do that".

**Question.** When every known copy is expressed from both of his haplotypes, does the frozen sample-aware corrector fix a
share of minimap2's wrong placements among tied reads, many times more than it breaks? The answer is a ceiling under ideal
expression, from a simulation, not a real-read result.

## The rule (frozen; nothing changes)

`bench/o3_graph/o2correct.py donor`, version 2, unchanged since bfa15f12, default flags, run as for register 1360
(`o2correct.py donor BAM COPIES.tsv GGO.fasta OUT.tsv --cache DIR --budget 500`). Domain, sample and decision are those of
`PREREG_o2_corrector_families_2026-10-20.md`. The corrector sees only the reads aligned to the GGO reference and the copies
table, never his haplotype assemblies.

## Families and copies

- **NPIP** (25 copies, catalog family `NPIP`), **LRPAP1** (22, `SM22`), **TBC1D3** (14, `TBC1D3`): 61 copies, each with its
  catalog exon chain (`rna_allele/fibro_cat.copies.tsv`, one model per copy).
- **SMN is dropped.** The census's "SMN (RefSeq, 3)" (Amendment 2 of PREREG_o2_corrector_families) matched the words "motor
  neuron": MNX1, SMNDC1 and SMN2. The gorilla annotation has one SMN gene (SMN2), in no catalog family, so SMN is not a
  multi-copy family here. This corrects that census line.
- **The YAGs stay out.** His Y has no second haplotype, so in his own genome every YAG tie is a true tie.

## Simulation (reads only)

- **Copy instances.**
  - **Reference instance:** each copy on the reference haplotype. Per chromosome the GGO reference is identical to one of his
    haplotype contigs (`rna_allele/chrmap.tsv`).
  - **Alternate instance:** the copy's counterpart on his other haplotype, by Amendment 4's positional key of
    PREREG_o2_corrector_families: ± 500 kb, `o2personal.py queries / orthologs --repaired --flank 500000`; the longest primary
    with MAPQ >= 20 on a non-reference contig, no exon-growth guard, overlaps to the lower divergence. It is recomputed for the
    61 copies; the 47 NPIP / LRPAP1 copies must reproduce `mainstay/positional/orthologs.tsv`.
- **Contested counterparts.** A counterpart that overlapped another copy's counterpart before the lower-divergence rule resolved
  it, winner and losers alike (for example maternal chr12 25.8-25.9 Mb, register 1361). Its reads are simulated, because his
  genome has the copy, but they are **never scored**. A copy that lost the overlap has a reference instance only.
- **Ideal expression.** Every instance gets the same **100 reads** over its one model: no expression spread, no silent copies
  (profile `A119b_testis_v4_ref` with `expr_log_sd` 0 and without the paralog-expression keys).
- **Reads.** `donorsim.py build`, the simulator of register 1346: the instance's own sequence (the reference, or his
  haplotype assembly through the lifted exons), with the profile's retained introns, truncation, end jitter and HiFi errors.
  - The read model is A119b testis Iso-Seq, the audited profile of register 1346 and the one-genome rounds. It is not fitted
    to his fibroblast library: `calibrate.py` is bound to CHM13 and its annotation.
  - His alleles, indels and conversions are real; only expression and sequencing are simulated.
- **Not simulated.** Copies his haplotypes carry but the reference lacks, such as the maternal-only NPIP locus of register
  1360. They are not known genes. Leaving them out favours the corrector: no foreign alleles teach its sample.
- **Two replicates:** seeds 3781 and 3782.
- **Alignment:** the real library's command, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, against the
  GGO reference (prebuilt `winloci_data/GGO.splice.mmi`), with minimap2 2.30 (the library used 2.31).

## Truth and scoring

- **Truth:** the reference copy whose instance produced the read. It is exact; no answer key is built from alignments.
- **Scoring:** `o2correct.py score ARM.tsv TRUTH.tsv --donor --among-candidates`. A read is scored when its truth is among its
  tied copies; reads from contested instances are left out.
- **Strata:**
  - **ref:** reads from reference instances, one genome. This is the no-harm control: minimap2's errors there should be exact
    sequence ties.
  - **alt:** reads from his other haplotype, the informative ones.
  - **all:** both together.

## Bars (each replicate, pooled over the three families)

- **P (power):** >= 50 scored reads and minimap2 wrong on >= 10; else UNDERPOWERED.
- **R (recovers a share):** fixes > breaks, and fixes − breaks >= 2% of minimap2's wrong calls.
- **S (safety):** fixes >= 5 × breaks.
- **H (no harm):** on ref reads, breaks <= fixes + 2.
- **F (families):** fixes > breaks in at least 2 of the 3 families, and breaks > fixes in none.

**Verdict:**
- **PASS:** P, R, S, H and F hold in both replicates.
- **FAIL:** a replicate where P holds and another bar fails.
- **UNDERPOWERED:** anything else.

## Expectations, written before any read is built

- **P:** expected to hold. The real library had 24 informative reads; ideal expression gives 100 reads to every alternate
  instance.
- **R, S, F:** no prediction. The corrector's earlier donor simulation (HG00733, register 1346) gave 408 fixes / 40 breaks;
  his haplotypes differ from the reference far less than a stranger's differ from CHM13.

## What a PASS would claim

With every known copy expressed, the sample's own unambiguous reads carry enough of his copy-specific variants for the frozen
corrector to recover a share of minimap2's wrong placements among tied reads, in more than one family, many times more than
it breaks. That is a ceiling for NPIP, LRPAP1 and TBC1D3 in KB3781.

**It would not claim:**
- real-read performance (register 1361: UNDERPOWERED);
- reads from copies the reference lacks (O3);
- contested loci;
- the YAGs;
- other read models or expression levels.

## Rules

- Each replicate is built and scored once; bugs → amendment + rerun.
- Heavy commands in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Gorilla only; never pooled with human numbers. Push only when asked.
- Code added for this run, behind new options with defaults unchanged:
  - `o2personal.py orthologs` writes the contested counterparts beside its output;
  - `donorsim.py build --per-instance` gives each haplotype instance its own read count;
  - `bench/o3_graph/o2idealsim.py` builds the copies, the simulation's instances, the truth strata and the bars.

## Results (code 9acb0233; run once per replicate; `/mnt/linuxdisk/tmp/o3_graph/idealsim/`): **FAIL**

- **Counterparts:** 53 of 61 copies (TBC1D3 13 of 14). The 47 NPIP / LRPAP1 copies reproduce Amendment 4's
  `orthologs.tsv` exactly. 14 LRPAP1 counterparts are contested: 7 simulated but unscored, 7 without an alternate instance.
- **Instances:** 61 reference, 46 alternate scored, 7 alternate contested; 100 reads each = 11,400 reads per replicate, all
  aligned (11,400 primaries, 0 unmapped).
- **Corrector:** 1,544 / 1,595 reads tied within 10 AS; its sample has 7,958 / 7,742 confident reads on 57 / 56 copies.

| replicate | scored | minimap2 wrong | fixes | breaks | families better / worse | P | R | S | H | F | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3781 | 1,319 | 471 | 8 | 7 | 1 / 1 | ✓ | ✗ | ✗ | ✓ | ✗ | FAIL |
| 3782 | 1,365 | 503 | 10 | 4 | 3 / 0 | ✓ | ✗ | ✗ | ✓ | ✓ | FAIL |

The corrector changes ~1% of the tied calls and recovers a net 1 and 6 of minimap2's 471 and 503 errors (0.2%, 1.2%).

**Where minimap2's errors are (descriptive; 3781 / 3782):**
- **Reference instances, one genome:** 147 / 142 errors, all but one exact AS ties (identical copies over the read; the
  coin toss the set answers). Near ties: minimap2 right on 386 / 386 and 392 / 393.
- **His other haplotype, exact ties:** 188 / 192.
- **His other haplotype, near ties:** 136 / 169. Here the corrector's own evidence points to minimap2's pick: the true copy
  has more unexplained events than minimap2's copy on 116 / 153 reads, equal on 11 / 8, fewer on 9 / 8.
  - They concentrate in four instances (120 of 136, 129 of 169), each sending its reads to one neighbour:
    - TBC1D3.1 → .0: 65 / 73;
    - TBC1D3.7 → .6: 20 / 28;
    - TBC1D3.0 → .1: 16 / 16;
    - NPIP.1 → .2 (and .13): 19 / 12.
  - TBC1D3.0 and .1 are an inverted pair 105 kb apart. The reference segment spanning both (31.88-32.03 Mb) aligns reversed
    on his maternal chr4 (32.28-32.42 Mb), between collinear flanks: an inversion between his haplotypes with its ends in
    the two copies.
  - TBC1D3.6 and .7 lie inside a region whose whole ± 500 kb aligns reversed on his maternal chr4.
  - In these pairs his copy at a copy's position is closer in sequence to the neighbouring reference copy. Position and
    sequence disagree on provenance, and no read-level evidence follows position there.

Audit: below (FAIL holds with caveats; the TBC1D3 reading above is withdrawn).

## Audit (metric-auditor, 2026-10-21): FAIL holds with caveats; the TBC1D3 diagnosis above is withdrawn

Every number reproduces: scored, minimap2 wrong, fixes, breaks, families and bars, plus the diagnosis counts (116 / 11 / 9
and 153 / 8 / 8). Checked and found sound:
- **Order:** design 09:04, code 09:07, first output 09:07; the design text is unchanged.
- **Simulation:** 114 instances × 100 reads; the ideal profile changes only the expression keys.
- **Blindness:** no truth leaks through read names (a seeded permutation); the corrector reads only the BAM, the copies table
  and GGO.
- **Read sources:** 13 spot-checked reads match their own instance best.
- **Robustness:** no reasonable reading turns FAIL into PASS. Adding the contested reads gives 11 / 8 and 10 / 5. The
  31 / 36 reads whose truth lies outside the tied copies are minimap2 errors either way.

Corrections:
- **The power is positional only.** On just 3 / 471 and 9 / 503 of minimap2's errors does the true copy have strictly fewer
  edits (NM + clipped bases) than minimap2's pick. Under a fewest-edit key the test would be UNDERPOWERED. Of the
  corrector's scored overrules (19 / 16), one per replicate moves a read to the copy with fewer edits. **FAIL means the
  reads do not carry the positional provenance minimap2 misses; no read-based rule can recover it.**
- **TBC1D3: "closer to the neighbouring copy" is withdrawn.** His maternal copies are mosaics of the pair's paralog
  differences, equidistant from the two reference copies over the copy models.
  - Edits to its own copy vs the neighbour: maternal .1 4 / 4, .0 5 / 5, .7 2 / 2, .6 21 / 21.
  - His maternal .0 and .1 differ at one base (exon 14); reads that miss the 5′ exons look like .0.
  - These two copies give 86 / 471 and 97 / 503 of minimap2's errors.
  - The .0 / .1 inversion is supported: a 144.6-kb reversed piece between collinear pieces of 526 and 411 kb. Its ends are
    not pinned to the copies.
  - The .6 / .7 whole-region reversal keeps the positional correspondence and explains nothing.
- **"No read-level evidence follows position" overreaches.** On the alternate-haplotype near-tie errors, NM is higher on the
  true copy for 96 / 129 reads, equal for 39 / 36, lower for 1 / 4.
- **F fails in 3781 on LRPAP1.** All 3 SM22 breaks are reads of his maternal SM22.1 moved to SM22.8. That copy is closer to
  SM22.8 (25-mer containment 0.865 vs 0.792), and 72 / 74 of its reads sit confidently on SM22.8, where they enter the
  corrector's sample (register 1343's mechanism).
- **Small corrections:**
  - The corrector overrules 37 / 1,544 and 18 / 1,595 tied calls (2.4% / 1.1%), not "~1%".
  - In 3782 the fourth instance is TBC1D3.5 → .3 (13 reads), not NPIP.1 (12).
- **Outside O2's domain:** 526 / 518 untied reads from his other haplotype sit confidently on the wrong copy by the
  positional key. Lifted NPIP exons change by up to −5.7 kb, because `--repaired` has no growth guard.
- **Not checked:** whether the near-identity of his maternal TBC1D3.0 / .1 is an assembly collapse. Checking it needs his
  HiFi reads at the exon-14 base, a heavy run.

**What stands:** FAIL. With every known NPIP, LRPAP1 and TBC1D3 copy expressed from both of his haplotypes, minimap2 is
wrong by position on ~500 tied reads per replicate. On all but 3-9 of them, though, the read is at least as close to
minimap2's pick as to its own copy, and the frozen corrector recovers a net 1 and 6. In his own genome, tied reads carry
almost no sequence evidence for the copy minimap2 misses: ideal expression adds power, not information.
