# Held-out family test — result

Run 2026-09-20 against `docs/PREREG_heldout_families_2026-09-20.md` (pre-registered and committed at
`aa4e7b1d`, truth-regex amendment at `42c298a9`, Test 2 at `ac9718ff` — all **before** the scores below).
Answers `docs/ADVISOR_QUESTIONS.md` §1.3 *"everything you have is one family"* and Q6.

## The headline

**The pre-registered bar is met.** On chromosomes the shipped rule's threshold was never selected on, it
scores **within 0.017** of the chromosomes it was selected on — and it touches **every** truth family.

| arm | chromosomes | Soto families | touched | **pooled F** | exact (F = 1.000) |
|---|---|---|---|---|---|
| ⭐ **held out** | chr2, chr8, chr10 | **20** | **20/20** | **0.7016** | 3 |
| development (LORO set) | chr5, chr7, chr21 | 11 | 11/11 | 0.7183 | 4 |
| | | | | **Δ = −0.017** | |

Pre-registered bar: Δ ≥ −0.10 **HOLDS** · −0.25…−0.10 PARTIAL · < −0.25 FAILS → ⭐ **HOLDS**.
Secondary bar (≥ 1 exact recovery on a zero-exposure chromosome) — met: chr2 has 2, chr8 has 1.

## Per chromosome, worst included

| chromosome | exposure (ledger/register) | families | touched | F | sens | prec | exact |
|---|---|---|---|---|---|---|---|
| **chr2** | **0 / 0** | 7 | 7 | **0.7295** | 0.770 | 0.792 | 2 |
| **chr8** | **0 / 0** | 5 | 5 | 0.6377 | **0.950** | 0.583 | 1 |
| chr10 | 1 / 0 | 8 | 8 | 0.7172 | 0.752 | 0.747 | 0 |
| chr7 | *development (LORO)* | 9 | 9 | **0.7825** | 0.778 | 0.830 | 4 |
| chr5 | *development (LORO)* | 2 | 2 | **0.4296** | 1.000 | 0.281 | 0 |
| chr21 | *development (LORO)* | 0 | — | not scoreable (no Soto family ≥ 3 members) | | | |
| chr16 | *development (NPIP)* | 8 | 7 | 0.6698 | 0.731 | 0.691 | 0 |
| chr17 | *development (TBC1D3)* | 14 | 13 | 0.7328 | 0.758 | 0.733 | 4 |

⭐ **The single worst chromosome in the whole panel is chr5 — a development chromosome**, one of the three
the 0.60 threshold was selected on (F 0.430, precision 0.281 on 2 families). If the rule were fitted to
its selection set, that is the last place you would expect to find the worst score. chr2, never looked
at, is the second best.

⭐ **n is now 53 families across 8 chromosomes**, against the dossier's standing concession of human
n = 3 (NPIP, TBC1D3, AMY). The "everything you have is one family" objection should be updated, not
merely conceded.

## Test 1 was declared VOID before Test 2 was designed

The first pre-registered truth — families by shared gene-symbol root — **fails on the development set
itself** (chr16 F 0.215, chr17 F 0.220), below its own FAIL bar. The reason is visible in what it counts
as a family: chr16's largest "missed families" are **MIR (46), RNA5S (31), ZNF (29), C16orf (12),
SNORA, SNORD, RN7SK, CDH, PRSS** — microRNA and rRNA gene classes, ancient protein-coding families, and
`C16orf`, which is a naming convention ("open reading frame on chromosome 16"), not a family. The method
targets *recent, sequence-similar* duplications and is right to leave those alone; what it does recover
under that truth are real recent duplications (NOMO F = 1.000, HERC2, SMG1 7/7, HBA, CLEC18).

Filtering by GFF biotype removes only **9 of 36** misses, so the truth could not be rescued by scoping.
Per §6 of the pre-registration it is reported void rather than swapped for a better-scoring variant.

## O2 — reproduces exactly, and the honest verdict is unchanged

Re-run today on the Y ampliconic genes (`--families` catalog, 8 families / 30 copies, A119b chrY):

| | recorded §6r4 | today |
|---|---|---|
| contested (O2's subject) | 3,641 | **3,641** |
| assigned | 12 (0.3%) | **12 (0.3%)** |
| tied (abstain) | 3,423 (94.0%) | **3,423 (94.0%)** |
| ambiguous | 206 (5.7%) | **206 (5.7%)** |
| families yielding any assignment | DAZ 9, RBMY 3 | **DAZ 9, RBMY 3** |

⭐ **Q8 is answered by construction**: every row carries a discrete status — `assigned`, `tied` or
`ambiguous` — and there is no fractional weight anywhere in the output. It never splits 1/k.
⚠ **The concession stands**: as a copy-resolution method on the hardest real case it resolves 0.3%.
The defensible claim is that it identifies the ambiguous population, refuses to guess, and says which
families are resolvable at all.

## What this does not show

- **Soto is the truth-source the 0.60 threshold was tuned against** (LORO on chr5/7/21). What is held
  out here is the **chromosome**, not the truth-source. This is not an independent-truth test.
- **chr6 could not be scored at all** — it has no Soto family with ≥ 3 members — so one of the two
  zero-exposure chromosomes contributes nothing to the headline. Its clusters (87 families, 225 members)
  are reported but unscored.
- **chr5's 2 families and chr21's 0** make the development arm thin; Δ = −0.017 rests on 11 families
  against 20, and should not be quoted as a precise effect size.
- ⚠ **Q9 reappeared as a measurement.** Under the amended single-pass root, **NPIPA (7 members) scores as
  wholly missed**: the method puts all 21 NPIP genes in ONE cluster, so bipartite matching assigns that
  cluster to NPIPB and leaves NPIPA unmatched. The method does not separate NPIPA from NPIPB.
