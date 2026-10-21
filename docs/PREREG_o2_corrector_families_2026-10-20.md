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

(none yet)
