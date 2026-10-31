# Multi-copy assembly mode: family assembly, then PSV placement (design, 2026-10-31)

> **Outcome (2026-10-31):** implemented as `bench/stgaps/famasm.py`; held-out HG002 **FAIL** (M1 recall +1.4 PASS, M2 precision 14.0 -> 8.5 FAIL, M3 FAIL, M4 PASS); A119b dev FAIL on the same bars. Not adopted. Record: `docs/LEDGER_O1.md` `### FAMASM RESULT`, register 1479.

User (2026-10-31, advisor's request): a multi-copy assembly mode that borrows information across family members, both reads and
structure, running on top of the normal single-copy assembly. Approach A of the brainstorm was approved.

## Why (what the ledger already knows)
- The single- vs multi-copy recall gap is in what the reads show at the own copy: only 19% (gorilla KB3781) / 27% (human A119b) of
  multi-copy reference transcripts have any primary read carrying their whole chain at their own locus, against 38% / 33% for
  single-copy genes; once a chain is observed twice at its copy we recover it as well as a single-copy one (87.3 vs 87.0 human)
  (`bench/stgaps/RESULTS.md`, register 1428). "Own locus" there is the PRIMARY placement: reads at a paralog count as unobserved.
- Borrowing STRUCTURE alone (family layer FW, PREREG_family_layer_2026-10-26: a single-read chain admitted when the lifted family
  chain has 2 molecules anywhere) FAILED its held-out: +0.1 recall point. Admitting every single-read chain gains 2 points for
  3.5-7.2 precision points. Single-read chains are not the gap; mis-placed reads are.
- Nothing may be emitted at a copy without reads that belong to it (registers 647, 1285; r1269: secondary-seeded transcripts are
  paralog shadows). The per-read copy certificate exists: O2 / `copy_phasing` (PSV certificate, correlation clustering default).

## The mode (FAMASM)
Input: BASE = the shipped `assemble` GTF (unchanged, every line kept); the O1 families (`<p>.fam.copies.tsv`, copy loci); the BAM;
the genome; the per-read copy calls of `copy_phasing run` on the families (Assigned / Undecided / Tied per read, with the copy).
Per family with >= 2 copies:
1. **Pool.** Every read with a primary or good secondary record (AS >= 0.98 x the molecule's genome-wide best, the shipped pool rule)
   on any copy locus of the family. One molecule once.
2. **Lift.** The family's copy map (fw_layer.py step 2: each copy locus aligned to the most-read copy with `minimap2 -x asm20 -c --cs`;
   a junction lifts when both its boundary bases fall on aligned columns). Each pooled read's intron chain, taken from the record on
   the copy O2 certified it to (or its primary when uncertified), is rewritten in family junctions. A read with an unliftable
   junction is kept with that junction as a copy-private junction (it can only support chains at that copy).
3. **Family chains.** The unchanged pass-1 rule over the pooled, lifted reads: exact family-chain groups, floor 2 molecules; the
   `--direct-transfrags`-style compatible-read support is NOT used for chain discovery (r1264 / r1265 precision), only in step 4.
4. **Per-copy admission (must-observe, certified).** A family chain is emitted at copy c when (a) >= 2 reads CERTIFIED to c carry
   the whole chain, or (b) reads certified to c that are contiguous, end-compatible sub-chains of it jointly cover every junction
   with >= 2 reads per junction (the FW/compat idea, now on certified reads only). Undecided and Tied reads count for no copy.
   A chain admitted at c is written in c's coordinates by inverting the lift (every family junction maps back to c; a chain with a
   copy-private junction of another copy cannot be admitted at c).
5. **Output.** FAMASM = BASE + admitted chains absent from BASE at that copy (same contig, exact intron chain), each tagged
   `famasm "1"; family "F"; certified_reads "n"; via "whole|compat"`. BASE lines are never altered or removed. Single-copy loci and
   one-copy families are untouched by construction (byte-identical, checked).
Reads are borrowed (pooled across copies and placed by certificate), structure is borrowed (chains are family-wide), evidence is
never borrowed (admission needs reads certified to the copy).

## Measurement (bars written here; numbers after)
Strata from `bench/stgaps/stratify.py` (annotation-only copy status). Scorer: gffcompare `=` on reference multi-exon transcripts,
multi-copy stratum recall and precision; `observed.py`-style recall on observed transcripts where "observed at the copy" now admits
reads certified to the copy; the mainstay and pair bars of the family thread must not move (the mode changes the GTF the families
are built from, so `mcl_families` is rerun on FAMASM and compared to BASE's families).
- M1 gain: multi-copy recall FAMASM - BASE >= 1.0 point (FW's bar was 0.5 and it gained 0.1; this mode moves reads, so it must do
  more than the structure-only layer to be worth its complexity).
- M2 precision: multi-copy precision FAMASM >= BASE - 0.5 point.
- M3 single-copy untouched: every BASE line kept; single-copy recall / precision |FAMASM - BASE| <= 0.1 point (expected 0.0).
- M4 families: mainstays recalled >= BASE and unrelated foreign <= BASE; pair precision >= BASE - 0.01 and recall >= BASE - 0.01.
- Descriptive: admitted chains by `via`, by family size, certified reads per admitted chain; admitted chains with no reference
  transcript (candidate unannotated isoforms, not gain); runtime.
Dev: human A119b (reported). **Held-out: HG002 Kinnex on HG002 T2T Q100 paternal** (`~/data/hg002/`, Ensembl 2022_07 truth): the one
substrate where PSV certificates can be checked against haplotypes; it was used today for family gates (r1473, r1475) and that is
disclosed; this is a new question (assembly at family copies). ADOPT as opt-in (`--multicopy` on the assemble stage) iff M1-M4 on
HG002. Otherwise a register row and the prototype stays in `bench/`.

## Shape of the implementation
Prototype in Python first (`bench/stgaps/famasm.py`, reusing `fw_layer.py`'s copy map and chain parsing, `copy_phasing run`'s per-read
calls, `observed.py` / `stratify.py` for scoring), one command per step so each product can be inspected: `pool`, `lift`, `chains`,
`admit`, `score`. A Rust port follows only after the held-out, as the assemble stage's `--multicopy` layer (NEXT TARGETS: one binary).
Tests before the run: single-copy byte-identity on a two-family fixture; a chain admitted via `whole` at one copy and refused at the
other for lack of certified reads; a `compat` admission needing every junction covered twice; an unliftable junction stays private.

## Out of scope
Changing BASE's chains or polish; quantification; moving reads between families; a family consensus transcript (approach C);
reads with no certificate counting toward any copy.
