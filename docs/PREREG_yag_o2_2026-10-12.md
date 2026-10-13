# PREREG — does O2 work on the Y ampliconic genes (YAGs)? (2026-10-12)

**Why.** The user: "YAGs are definitely an important test" for O2. The Y is HAPLOID, so every difference between two Y
copies is a PSV, never a heterozygous allele. Known before this prereg (10-06 audit, `rustle_eval_2026-10-06/03_yag.md`):
real human Y — O2 assigns 12 of 3,641 contested molecules; 10-06 simulation (human chrY roster, 0.13% error) — O2
assigns 0 of 2,425 MAPQ-0 reads, including the 117 reads whose copies are only 98-99% identical (that last fact is
unexplained: such copies should usually be distinguishable).

## Arm A — simulation with truth (human chrY; reuses the 10-06 products, nothing re-simulated)

Reads `o2_default_roster_2026-10-06/chrY/sim.bam` (10,224 reads from the 458-copy chrY roster, names = `FAM|copy|i`,
genome-mapped), families `asg.fam.copies.tsv` / `.fa`, MCL0 (the Yq12 satellite) skipped — the exact 10-06 command,
with the CURRENT binary and `--dump-star`. Per molecule O2 decided on (any status), from `star_reads.tsv`:
- **distinguishable** = the truth copy t is a candidate and every other candidate c has >= 2 columns where t's and c's
  aligned bases differ and the read covers both (O2's own m >= 2 step; register: "m = 1 never assigns").
- **A-P1 precision:** 0 wrong among assigned (upper 95% bound reported; UNDERPOWERED if < 30 assigned).
- **A-P2 recall:** of distinguishable molecules, >= 0.80 assigned to t.
- **A-P3 calibration:** of NOT distinguishable molecules, <= 0.05 assigned.
- Reported: end-to-end placement (O2's copy where assigned, else the aligner's primary) vs truth; by identity band.

## Arm B — real testis reads, split-half concordance (human A119b, CHM13 chrY)

`copy_assign --dump-star` on `A119b.t2t.bam` restricted to the A119b catalog's chrY families (the default pipeline's
own copies; satellite/array families with > 100 copies skipped as in Arm A). The read-star certificate uses PSV columns
only here (`read_star_junctions false`), with per-column error e/3 = 0.001: for the best copy b and each rival c,
p = P(Binomial(n_bc, 0.001) >= k_bc) (n = columns where b and c differ and the read covers both, k = those matching b);
Assigned iff every rival has p < alpha / (#rivals) and LLR > 0, Tied if some rival has n = 0, else Ambiguous; b = the
maximin pairwise-LLR candidate. Python reimplementation `bench/yag_o2/split_half.py`.
- **B-P0 reproduction:** on the full columns it reproduces O2's own status and copy for >= 0.99 of molecules whose final
  status was not changed by a later gate (origin rejection, tie outside the catalog). Fail -> stop, nothing else read.
- Halves: columns at read positions < L/2 vs >= L/2 (L = read length from the BAM) — disjoint evidence from disjoint
  parts of the same molecule; on a haploid Y no allele can make them disagree.
- **B-P1 split-half concordance:** among molecules where BOTH halves are Assigned, the same copy in >= 0.98 (>= 30 such
  molecules, else UNDERPOWERED).
- **B-P2:** where O2 (full read) and a half both assign, the same copy in >= 0.98.
- Reported: per family; the fraction of contested molecules where at least one half assigns.

## Side check (no bar, report only)

SUN calls vs an external SUN set (human: the QuicK-mer2 k-mer set reproduced in register row 1171), if its files are
on disk.

Gorilla (OR6737 testis) needs a families stage first; added as an amendment with the same bars only if run.
Code `bench/yag_o2/`; outputs `/mnt/linuxdisk/tmp/yag_o2/`.

## Amendments made while running (2026-10-12)

- **Dump bug (found, fixed, TDD):** `--dump-star` kept ONE proof per read name, so a molecule evaluated in two families
  got family B's columns next to family A's verdict. First-pass numbers built on it are RETRACTED (53 "distinguishable"
  simulated reads, a "17 wrong" re-certification). Fix: every registration kept with its candidates' tids; the writer
  picks each row's own family proof, once per (molecule, family). Assignments byte-identical to the unfixed binary.
- **Scope (user, mid-run): "we only need this for the YAGs, so DAZ, PRY, etc."** Families = catalog families with a
  copy overlapping a RefSeq CHM13 chrY gene DAZ / RBMY / TSPY / VCY / XKRY / CDY / HSFY / BPY2 / PRY / PRORY
  (`bench/yag_o2/prep_yag.py`): A119b catalog 11 families / 86 copies (TSPY 46, RBMY 6, HSFY 6+3+2, VCY 6, DAZ 4+2,
  PRY 4, CDY 4, BPY2 3); simulation roster 11 families / 88 copies. (A first Y-only run is not reported as a result.)

## Results

**Arm A (simulation, YAG families):** of 1,460 contested molecules whose true copy is a candidate, 1,416 have a rival
IDENTICAL over every column the read covers and 42 differ at one column; 2 are distinguishable.
- A-P1 precision: 0 assigned -> UNDERPOWERED. A-P2 recall: 0 / 2 -> FAIL as registered (n = 2). **A-P3 calibration:
  0 / 1,458 undecidable assigned -> PASS.** Reproduction of O2 by the Python certificate: 2,404 / 2,404.

**Arm B (A119b testis reads, YAG families):** 7,844 molecules; O2 assigns 4. **B-P0 reproduction 3,445 / 3,445 = 1.0
PASS.** **B-P1 split-half concordance 31 / 57 = 0.544 FAIL** (HSFY 11/23, BPY2 16/22, PRY 2/9, other HSFY 2/3; DAZ and
TSPY: no half assigns). B-P2 3 / 3 -> UNDERPOWERED.

**Reading (post hoc; not tested here).** The haploid argument assumed the reads come from the reference man's Y. They
do not: CHM13's Y is HG002's, the testis RNA is A119b's. Real reads carry MORE distinguishing columns than reads
simulated from the reference copies (386 molecules with an assignable half vs 37), as expected if A119b's YAG copies
differ from HG002's (gene conversion between palindrome arms differs between men), so a read's halves can point to
different REFERENCE copies with no error by O2. O2's full-read certificate abstains on all of those. The split-half
test therefore cannot separate O2 error from donor-vs-reference mosaicism without a matched donor: a T2T Y and testis
Iso-Seq from the SAME man (YAGs are testis-specific; the T2T Ys come from cell lines). The defensible statement for
the YAGs is calibration, not accuracy.
