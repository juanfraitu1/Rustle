# Figure 9: The closed loop — tied reads given to their assigned copy and re-assembled

**Status.** Pre-registered in `docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md` before any pass-2 assembly or real-read
union-test count existed. The pipeline exists (`tools/rustle_reassemble.sh`, `bench/loop_home.py`, the assembler's
opt-in `RUSTLE_READ_HOME_TABLE`); no verdict sample has been run. Every panel says "not built" until its table exists.
The development numbers below (human chr16 simulation) are disclosed in the pre-registration's amendment 1 and are not
part of any verdict. `python3 figures/fig_loop.py summary` prints every number this caption quotes.

**Claim (to be filled from the tables; the bars were fixed before they existed).**
- LOOP1: on at least 3 of the 4 held-out samples (human testis, gorilla KB3781, chimpanzee PTR, orangutan PPY), giving
  each assigned read to its copy and re-assembling gains reference intron chains where reads tie (tie fraction > 0.5,
  gains in at least 2 read-sharing groups) without lowering intron-chain precision or the fraction of de novo loci at
  a Liftoff locus by more than 0.2 points.
- LOOP2: on each held-out sample's simulation, placing every read at its true copy (the best any assignment could do)
  recovers at least 1.0 point, or 5 copies, more exact intron chains than the default; otherwise the placement ceiling
  is closed.
- Predicted before looking: no effect or a small one (the union test assigned 0 reads on both development
  simulations; tied reads are 1-3% of reads).

## Caption

**Figure 9 | Re-assembling after copy assignment (annotation-free; Rustle-internal arms, not a tool comparison).**
Pass 1 is Rustle's default: loci assembled from each read's primary alignment plus secondary alignments scoring
≥ 98% of the read's best alignment score genome-wide, then the default de novo families (locus representatives
compared by their exon sums, MCL). A read whose best alignment score is reached at several copies (equal best score)
is then tested once over every candidate copy, across families and outside them (the union test). Only a read the test
**assigns** moves: in pass 2 it is taken only at its assigned copy, and every other alignment of it, the primary
included, is dropped; reads left unassigned keep their pass-1 alignments. The families are not recomputed from pass 2.
Arms: Rustle, primary alignments only (P); Rustle default; the loop on the default (assigned reads at their copy,
everything else as the default); the loop on primaries (assigned reads at their copy, everything else primaries only);
and, on simulations only, every read at its true copy (the ceiling of any assignment).
**a**, Alignment records the loop moves per sample: removed away from the assigned copy, and admitted at the assigned
copy beyond the default's read pool. Below 20 records a sample cannot change more than 19 read groups and is not scored
(pre-registered stop rule). **b**, Simulated reads (known source copy): multi-exon copies whose exact intron chain an
arm assembles. **c**, Intron-chain sensitivity and precision against the sample's RefSeq annotation (gffcompare 0.12.10,
as Fig. 1). **d**, Reference chains matched by the loop minus those matched by the default, by the tie fraction of
their reads (Fig. 3's bins: a second alignment ≥ 98% of the best score), with the number of read-sharing groups that
carry them. **e**, Read-supported Liftoff loci (Fig. 8) found by each arm's de novo loci, on two fixed universes:
loci with ≥ 2 reads whose primary alignment lies on the locus's exons (as Fig. 8), and loci with ≥ 2 reads with any
candidate alignment there (primary, or secondary ≥ 98% of the read's best score), which counts the copies whose reads
are tied.

**Samples and exposure.** Verdict: human testis (T2T-CHM13 v2.0), gorilla KB3781 fibroblast cell line (mGorGor1,
GCF_029281585.2; used before for missing-copy work only), chimpanzee PTR (mPanTro3), orangutan PPY (mPonPyg2): never
used to develop the assembler, the seeding rule or copy assignment. Development: human A119b (copy-assignment rules,
chr16) and gorilla OR6737 testis (seeding rule, chr20 (NC_073244.2)). Species and samples are never pooled. Units of
independence: reference chain (c), read-sharing group or gene (d), Liftoff source record (e), source copy (b).

## Methods

1. Pass 1: `tools/rustle_pipeline.sh assemble` and `families` (with `--bam`, so the families' per-locus units are
   written in the copy-assignment contract).
2. `tools/rustle_reassemble.sh union`: `copy_assign --families UNITS --union-certificate` over every contig (exact
   alignment-score ties; ties across contigs are invisible to it and left as they are). Only the union test's verdicts
   are read (`union_certificate.tsv`: `assigned_family` → the unit `family_id:copy_idx`; `assigned_outside` → the
   outside locus); `tied`, `ambiguous` and `no_result` reads never move.
3. `home`: `bench/loop_home.py union` writes one home span per assigned read and leaves out a read with no alignment
   record whose aligned block (M/=/X) overlaps its home (counted; such a read would need a lifted record, not built).
4. `pass2`: the assemble stage again with `RUSTLE_READ_HOME_TABLE`: a listed read's record enters the assembler's read
   pool iff an aligned block overlaps its home span (any flag, any score); its other records are dropped. The copy
   assignment's own read input is never filtered. Unset, every output is byte-identical to the pre-loop binary
   (checked with `cmp` on the chr16 simulation and slice; an empty table is inert).
5. `g0`: the filter's own counters (records kept at home, admitted at home beyond the default pool, removed away from
   home), before the assembler's coordinate de-duplication.
6. Scoring: panel b, `bench/loop_home.py score-sim` (gffcompare against the simulated copies' chains; a single-exon
   copy counts when an arm locus covers ≥ 50% of its bases); panels c-e with the Fig. 1, Fig. 3 and Fig. 8 code on
   each arm's GTF.

## Development results (human chr16 simulation, legacy copy catalog; not a verdict)

`${work}/o2sim/human` (28,453 simulated reads from 1,398 catalog copies, mapped genome-wide); arms run with the
current binaries; `rustle_figures_dev/loop/d1_sim_human/m5.tsv`:

| arm | transcripts | multi-exon copies, exact chain (of 217) | single-exon copies covered (of 1,181) | transcripts matching a copy / at a copy |
|---|---|---|---|---|
| Rustle, primaries only | 437 | 191 | 226 | 401 / 435 |
| Rustle default | 528 | 188 | 247 | 393 / 493 |
| loop on the default | 528 (identical file) | 188 | 247 | 393 / 493 |
| loop on primaries | 437 (identical file) | 191 | 226 | 401 / 435 |
| every read at its true copy | 409 | 191 | 209 | 386 / 409 |

- The union test assigns none of this simulation's reads (991 left unassigned, 2 ambiguous; 1,255 of its 1,263
  MAPQ-0 reads have an NM-identical twin, Fig. 4), so both loop arms reproduce their base arm byte for byte.
- Placing every read at its true copy moves 8,554 records out of the default's pool and admits 1,518 at home: +3
  multi-exon copies (+1.4 points) over the default, precision 94.4% vs 79.7% of transcripts at a copy, but 38 fewer
  single-exon copies. Mechanism of the loss, read from the logs: the polish's single-exon read floor is data-adaptive
  (28 reads on chr16 here); a copy's own simulated reads (10-23) fall below it once its echo reads from sibling copies
  are gone, whereas in the default the echoes lift it over the floor.

## Caveats

- The loop moves only reads the union test assigns. On both development simulations it assigns none; its real-read
  rate is unmeasured until the verdict runs.
- Placements the reads do not have are not created (no lifted records).
- The single-exon floor interaction above applies to any arm that removes reads; it is reported, not tuned.
