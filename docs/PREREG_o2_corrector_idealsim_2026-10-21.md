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
