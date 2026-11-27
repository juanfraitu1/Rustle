# PREREG: own-evidence seeding, a precision rule for the assembly with secondaries (row 1984, 2026-11-26)

The user's goal: increase the assembly's precision through empirical rules that are not overfit for apes. Where the
earlier precision work stands (rows 1769-1792, 1876-1879, memory `assembly`): the chain score, its feature families, the
isoform fraction, the sub-chain collapse, 5' pieces, near-miss junctions and read floors are fitted, traded or closed. None
of them looks at WHERE a seed comes from.

## What is new this week and general
`all` assembles with the `good` pool: every secondary alignment at 0.98 or more of the read's best AS seeds a
transcript, wherever it lands. Rows 1952, 1961, 1971 measured what that costs: a transcript no primary alignment
stands behind is right in 0 of 361 cases at 60 planted reads per isoform; a locus only secondaries support is on an
expressed copy in 0 of 80; and in the band below, the seeds that land where the sample has reads of its own fit a true
isoform 84% of the time against 34%. The mechanism is not about apes: a place the sample puts no read of its own at,
or an intron none of its own reads shows, is evidence about a relative, in any genome with duplicated genes.

## Rule (`bench/o2_families/pool_bam.py --rule own:as:0.98:1`)
A secondary alignment seeds only when (a) its AS is at least 0.98 of the read's best, as today, AND (b) at least 50 of
its bases lie where a primary alignment with mapping quality 1 or more is placed, AND (c) each of its introns is an
intron of such a primary alignment. No threshold is new; (b) and (c) use one read, the smallest amount of own evidence.

## Test order: non-ape data first, apes last and only as a check
Arms, `PREFIX.gtf` and `PREFIX.final.gtf` of stage `assemble`, multi-exon chains by exact intron chain
(`r1831_mono_end_rules/score.py`): P = primaries only; G = today's `good` pool; O = the rule.
1. Complete truth: fly, mouse and zebrafish simulated reads.
2. Real, non-ape: mouse ES, heart, adrenal, cortex against the mouse reference annotation.
3. Families with exact truth (human, planted): the half-silent sets and the thin-copy controls of rows 1952 / 1980.
4. Apes (OR6737, KB3781): only after 1-3, as a check that nothing breaks.

## Claims (O against G)
- **V1, complete truth:** on each of the 3 simulated sets wrong chains are not more than G's and correct chains are at
  least 0.998 of G's; on at least 2 of 3 wrong chains are fewer.
- **V2, mouse libraries:** on at least 3 of 4, precision (correct / chains) is not below G's, with correct chains at
  least 0.998 of G's on all 4.
- **V3, families:** in the 3 half-silent planted sets O models fewer silent copies than G in all 3; at the thin copies of
  row 1980's controls O's isoform recall is within 0.03 of G's.
A default change needs V1-V3 and is the user's decision.
## Predictions (mine)
V1 and V2 hold with a small effect: outside duplicated genes few secondaries reach 0.98, so G and O differ by tens of
chains. V3: fewer silent copies holds; the recall bound is at risk (a thin copy whose reads all tie has no placed read,
so nothing is seeded there). If the effect on 1 and 2 is nil, the rule is a rule for the family mode and is reported so.

## Result, steps 1-3 (jobs 4512734-4512741). `PREFIX.gtf`, multi-exon chains: chains / correct (precision).
| dataset | P, primaries only | G, today's pool | O, the rule | O correct / G correct |
|---|---|---|---|---|
| fly simulated | 7,203 / 6,791 (0.943) | 7,851 / 6,796 (0.866) | 7,209 / 6,791 (0.942) | 0.9993 |
| mouse simulated | 27,371 / 25,461 (0.930) | 27,613 / 25,495 (0.923) | 27,403 / 25,470 (0.929) | 0.9990 |
| zebrafish simulated | 14,595 / 12,287 (0.842) | 15,057 / 12,312 (0.818) | 14,656 / 12,301 (0.839) | 0.9991 |
| mouse ES | 44,475 / 17,304 (0.389) | 46,117 / 17,352 (0.376) | 44,862 / 17,313 (0.386) | 0.99775 |
| mouse heart | 22,576 / 13,099 (0.580) | 22,790 / 13,123 (0.576) | 22,606 / 13,104 (0.580) | 0.9986 |
| mouse adrenal | 24,901 / 14,371 (0.577) | 25,174 / 14,402 (0.572) | 24,962 / 14,378 (0.576) | 0.9983 |
| mouse cortex | 31,630 / 15,758 (0.498) | 31,921 / 15,777 (0.494) | 31,685 / 15,765 (0.498) | 0.9992 |
Wrong chains, G -> O: fly 1,055 -> 418; mouse simulated 2,118 -> 1,933; zebrafish 2,745 -> 2,355. Secondary alignments
that seed: 2.7% of all against today's 11.8% (fly), 2.2% / 3.6%, 1.1% / 7.8%, mouse libraries 3.0-7.0% / 5.3-8.9%.
`PREFIX.final.gtf` shows the same order with smaller gaps (fly 0.946 -> 0.974; mouse ES 0.405 -> 0.407).

Families, planted (human; exact truth):
| | P | G | O |
|---|---|---|---|
| all expressed, 60 reads: isoform precision (exact isoforms), seeds 1967-1969 | 0.83 / 0.83 / 0.81 | 0.654 (197) / 0.654 (197) / 0.659 (197) | 0.798 (197) / 0.786 (198) / 0.775 (196) |
| half silent: isoform precision | 0.65 / 0.62 / 0.62 | 0.452 / 0.417 / 0.488 | 0.560 / 0.578 / 0.602 |
| half silent: silent copies with a transcript | 23 / 25 / 24 | 45 / 48 / 40 | 33 / 32 / 27 |
| thin copies (3 reads, row 1980's controls): isoform recall | 0.364 / 0.307 / 0.329 | 0.711 / 0.728 / 0.724 | 0.561 / 0.583 / 0.605 |
| thin copies: isoforms with one own read found | 7% / 6% / 7% | 56% / 60% / 61% | 32% / 40% / 42% |

### Claims
- **V1 HOLDS:** fewer wrong chains on 3 of 3 simulated sets (fly -60%, mouse -9%, zebrafish -14%), correct chains 0.9990-0.9993 of G's.
- **V2 FAILS by four chains:** precision is above G's on 4 of 4 mouse libraries, and correct chains are at least 0.998 of
  G's on three; on mouse ES the ratio is 0.99775 (17,313 against the 17,318 needed).
- **V3 HOLDS in its first half, FAILS in its second:** silent copies modelled 33 / 32 / 27 against 45 / 48 / 40; isoform
  recall at thin copies 0.56-0.61 against 0.71-0.73 (bound 0.03).

### Reading
- **Seeding from every secondary costs precision everywhere, not only in ape families.** Against primaries only, today's
  pool adds 5, 34 and 25 correct chains on the three simulated sets and 643, 208 and 437 wrong ones; on mouse ES 48
  correct for 1,594 more chains. This is the price of `all`'s `good` pool on the whole assembly, measured on non-ape data.
- **The rule takes back about nine tenths of that price and keeps a third to a half of the gain:** its precision is
  within 0.001-0.003 of primaries-only on all seven non-ape sets, with a few more correct chains than primaries-only on six.
- **In families it keeps what is well supported and drops what is not:** at 60 reads the exact isoforms are unchanged and
  precision rises by 12-15 points; silent copies fall by a third. At thin copies it keeps half of the borrowing gain
  (0.33 -> 0.58, today 0.72): a copy whose few reads all tie has no placed read, so nothing may be seeded there.
- So it is a third setting between the two that exist (`primary`, `good`), not a free gain: more precise than `good`
  everywhere, more sensitive than `primary` in families, less sensitive than `good` at the thinnest copies.

## Step 4, the ape check (jobs 4512793 / 4512794; gorilla, RefSeq)
| | P | G | O | O correct / G correct |
|---|---|---|---|---|
| OR6737 | 68,536 / 24,411 (0.356) | 69,647 / 24,563 (0.353) | 69,051 / 24,443 (0.354) | 0.9951 |
| KB3781 | 60,003 / 23,495 (0.392) | 60,650 / 23,554 (0.388) | 60,243 / 23,512 (0.390) | 0.9982 |
The same order as on the non-ape sets (P most precise, G least, O between and close to P). OR6737, a testis sample with
many expressed families, loses 120 of G's correct chains (0.5%) for 596 fewer chains; KB3781 42 for 407. Against
primaries only today's pool buys 152 and 59 correct chains with 1,111 and 647 chains. Nothing breaks; the size is as elsewhere.

## Summary of row 1984
- The rule is not an ape fit: it was specified from the mechanism, tested first on fly, mouse and zebrafish, and gives
  its largest precision gain there (fly 0.866 -> 0.942).
- It is a setting between `primary` and `good`: on nine real or simulated datasets its precision is within 0.001-0.004
  of primaries-only and 0.001-0.076 above today's pool, with 99.5-99.9% of today's correct chains.
- It does not meet the bars for a default as registered (mouse ES by four chains; recall at three-read copies).
- It exists in `bench/o2_families/pool_bam.py` only (`--rule own:as:0.98:1`, then `--seed-pool all`): one extra pass over
  the BAM in Python. Making it a seeding option of the driver, or of the assembler itself, is not done.

## Variant registered after steps 1-4, before any run of it: the anchor also counts tie-break primaries (O0)
The rule missed two bars: mouse ES by four correct chains, and the recall at three-read copies (0.58 against 0.72). Both
losses have one cause that can be named: a copy whose reads ALL tie has no primary with mapping quality 1 or more, so the
rule seeds nothing there. Variant O0 = the same rule with the anchor and the intron check taken from primaries of ANY
mapping quality (`own:as:0.98:1:0`): a place where the aligner put at least one of the sample's reads, by a clear lead
or by the tie-break. Expected price: a silent twin that received tie-break primaries is anchored again.
Same data, same bars V1-V3 as above, O0 against G; plus one more:
- **V4:** O0's precision is at least O's - 0.005 on the seven non-ape sets (it must not give the gain back).
Prediction (mine): thin-copy recall comes back to within 0.03 of G's; silent copies rise to near G's (V3's first half
at risk); V1 and V4 hold; mouse ES passes.

## Result of the variant O0 (jobs 4512903-4512910): the anchor also counts tie-break primaries
`PREFIX.gtf`, chains / correct (precision); O0's correct chains against G's:
| dataset | G | O | O0 | O0 correct / G correct |
|---|---|---|---|---|
| fly simulated | 7,851 / 6,796 (0.866) | 7,209 / 6,791 (0.942) | 7,219 / 6,792 (0.941) | 0.9994 |
| mouse simulated | 27,613 / 25,495 (0.923) | 27,403 / 25,470 (0.930) | 27,441 / 25,491 (0.929) | 0.9998 |
| zebrafish simulated | 15,057 / 12,312 (0.818) | 14,656 / 12,301 (0.839) | 14,730 / 12,308 (0.836) | 0.9997 |
| mouse ES | 46,117 / 17,352 (0.376) | 44,862 / 17,313 (0.386) | 44,929 / 17,321 (0.386) | 0.9982 |
| mouse heart | 22,790 / 13,123 (0.576) | 22,606 / 13,104 (0.580) | 22,644 / 13,113 (0.579) | 0.9992 |
| mouse adrenal | 25,174 / 14,402 (0.572) | 24,962 / 14,378 (0.576) | 24,996 / 14,386 (0.576) | 0.9989 |
| mouse cortex | 31,921 / 15,777 (0.494) | 31,685 / 15,765 (0.498) | 31,726 / 15,772 (0.497) | 0.9997 |
Families (planted): silent copies with a transcript 33 / 32 / 27 (as O; G 45 / 48 / 40); isoform precision at 60 reads
0.785 / 0.776 / 0.767; isoform recall at three-read copies 0.649 / 0.640 / 0.675 (O 0.56-0.61, G 0.71-0.73, P 0.31-0.36).
- **V1 HOLDS** (wrong chains 427, 1,950, 2,422 against G's 1,055, 2,118, 2,745; correct 0.9994-0.9998).
- **V2 HOLDS** (precision above G on 4 of 4; correct at least 0.998 on all 4, mouse ES 0.9982).
- **V3:** fewer silent copies HOLDS (I predicted they would rise: the intron condition keeps them out); thin-copy recall
  still FAILS its bound (0.05-0.09 below G; the bound was 0.03).
- **V4 HOLDS** (precision at most 0.004 below O).
O0 is the better setting: every bar O met, plus mouse ES, plus half of the thin-copy recall O gave up.

## In the pipeline (user's goal; built 2026-11-27)
`--seed-pool own` (also `RUSTLE_SEED_POOL=own`): stage `assemble` first writes the secondaries that pass the rule into
`PREFIX.seedpool.d/pool.bam` (`bench/o2_families/pool_bam.py --rule own:as:0.98:1:Q`), assembles from it with every
remaining secondary admitted, and removes the folder (`RUSTLE_SEEDPOOL_KEEP=1` keeps it). Q = `RUSTLE_SEED_OWN_MAPQ`,
default 0 (variant O0); 1 = arm O. The other stages read the original BAM. Opt-in: `assemble` alone still defaults to
`primary` and `all` to `good`; which pool `all` uses is the user's decision.
Checks (jobs 4512911, 4512978; the second after the option's default became the variant: default = arm O0 exactly, `RUSTLE_SEED_OWN_MAPQ=1` = arm O exactly): on the fly reads the option gives exactly the transcripts of the benchmark arm, with a
score table and without one (the stage builds it); the folder is removed; a width with this pool and a bad
`RUSTLE_SEED_OWN_MAPQ` are refused with exit 2 before anything runs (two new cases in `tests/pipeline_seed_pool.rs`);
whole `cargo test` 1,591 passed, 0 failed.

## Configuration H, registered before any run of it: `own` for the main assembly, `good` only for the mode's rescue
The one bar still missed is the isoform recall at three-read family copies. In the pipeline those copies are served by
the `multicopy` stage, which adds one-read transcripts at family copies the main assembly leaves empty, from a run of
its own with a one-read floor. H = `assemble --seed-pool own`, `families`, then `multicopy --seed-pool good` on the
same prefix: the main assembly `PREFIX.gtf` is the precise one, the family mode's rescue still sees every near tie.
Same planted family sets; `PREFIX.multicopy.gtf` of H against today's.
- **H1:** isoform recall at the three-read copies is within 0.03 of today's in all 3 seeds.
- **H2:** silent copies with a transcript are fewer than today's in all 3 half-silent sets.
- **H3:** at 60 reads per isoform the exact isoforms are at least today's - 1 and isoform precision is above today's, in all 6 sets.
`PREFIX.gtf` of H is arm O0 by construction (V1, V2, V4 above). If H1-H3 hold, every registered bar is met by H.
Prediction (mine): H1 holds; H2 at risk (the rescue may bring the silent copies back); H3 holds.

## Result of configuration H (job output `r1984_own_evidence_seeding/h.*.out`)
| the mode's `PREFIX.multicopy.gtf` | today | `own` for everything | H |
|---|---|---|---|
| isoform recall at three-read copies, seeds 1981-1983 | 0.711 / 0.728 / 0.724 | 0.649 / 0.640 / 0.675 | 0.684 / 0.689 / 0.693 |
| silent copies with a transcript, seeds 1967-1969 | 45 / 48 / 40 | 33 / 32 / 27 | 45 / 48 / 40 |
| 60 reads, all expressed: isoform precision (exact isoforms) | 0.654 (197) / 0.654 (197) / 0.659 (197) | 0.785 / 0.776 / 0.767 (197 each) | 0.785 / 0.776 / 0.767 (197 each) |
| 60 reads, half silent: isoform precision | 0.452 / 0.417 / 0.488 | 0.554 / 0.544 / 0.601 | 0.466 / 0.465 / 0.533 |
H's main assembly `PREFIX.gtf` (half-silent sets): silent copies 32 / 32 / 26, isoform precision 0.560 / 0.544 / 0.605.
- **H1 FAILS narrowly:** 0.027, 0.039 and 0.031 below today's (bound 0.03; one seed of three inside it).
- **H2 FAILS:** the rescue brings every silent copy of today's mode back (45 / 48 / 40); I had marked this as the risk.
- **H3 HOLDS:** exact isoforms unchanged, precision above today's in all six sets.
No configuration meets every bar. What each buys against today's `all`:
- `own` for everything: the main assembly AND the mode more precise, a third fewer silent copies, three-read recall 0.05-0.09 lower.
- H: the main assembly more precise (as O0 on every non-ape set), the mode's product as today's in silent copies and
  0.03-0.04 lower in three-read recall, more precise at 60 reads.
The thin-copy reach and the silent copies are the same reads: a tied read at a place with no evidence of its own is
either a thin copy's only support or a silent copy's shadow, and nothing in it says which (rows 1971, 1980).
A default is the user's decision; by the rule of 2026-11-21 it also needs data this choice was not made on.

## The user's decision (2026-11-27) and its confirmation on unused data, registered before the run
Asked which pool `all` should use, the user chose the MIXED configuration (H): the main assembly from `own`, the other
stages from `good`. Built: `all` with no pool named seeds `PREFIX.gtf` with `own` (`RP_ASM_OWN`); `all --seed-pool good`
is the default before; a filter that cannot run (no pysam / numpy) is a warning and the stage seeds with `good`.
Confirmation on data this choice was not made on, the whole stage chain run both ways (`all` against `all --seed-pool good`):
- fresh planted family sets, seeds 1990 and 1991 (E and S at 60 reads per isoform; E at 3 reads);
- WTC11 replicates 1 and 2 (human cell line; not used for this rule), `PREFIX.gtf` against the human reference annotation.
No unused non-ape reads exist in the workspace (the second zebrafish BAM is the same read set); the non-ape evidence
is V1 / V2 above.
- **U1:** on both WTC11 replicates `PREFIX.gtf` has a higher precision than under `good` and at least 0.998 of its correct chains.
- **U2:** on the fresh planted sets, three-read isoform recall of `PREFIX.multicopy.gtf` is within 0.04 of the old default's
  in both seeds, silent copies with a transcript are not more than the old default's in both S sets, and isoform
  precision at 60 reads is above it in all four sets.
- **U3:** the labels still run and name the same number of expressed copies within 2 (the labels read the original BAM).
Prediction (mine): all hold.

## Result of the confirmation (job 4512998; the first submission, 4512995, is void)
Void first run, said plainly: `all` starts each stage as a child process and my switch lived in the parent only, so the
assemble stage never saw it; "new" and "old" came out identical on every set and the log had no line for the new pool.
Fixed (the parent marks its assemble call with `RP_ASM_MAIN`, the stage clears the mark before it starts any assembly of
its own) and run again; every figure below is from the rerun, where the log shows the pool in use.

| unused data | `all --seed-pool good` (before) | `all` (new default) |
|---|---|---|
| WTC11 r1, `PREFIX.gtf`: chains / correct (precision) | 45,046 / 25,166 (0.5587) | 44,699 / 25,140 (0.5624); 0.9990 of the correct |
| WTC11 r2 | 34,170 / 20,126 (0.5890) | 33,878 / 20,095 (0.5932); 0.9985 |
| planted 1990 / 1991, all expressed, 60 reads: isoform precision (exact) | 0.656 (198) / 0.629 (198) | 0.767 (198) / 0.750 (198) |
| half silent: the mode's isoform precision (exact) | 0.466 (103) / 0.462 (103) | 0.500 (103) / 0.490 (102) |
| half silent: silent copies with a transcript, mode / main assembly | 39 / 38; 48 / 47 | 39 / 28; 48 / 31 |
| three-read copies: isoform recall of the mode | 0.737 / 0.689 | 0.680 / 0.632 |
| labels: loci called expressed | 82, 41 / 83, 43 | 82, 41 / 83, 43 |
- **U1 HOLDS** (both libraries: precision up, 0.9990 and 0.9985 of the correct chains).
- **U2 FAILS in one of its three parts:** three-read recall is 0.057 lower in both seeds (bound 0.04); silent copies are
  not more (equal in the mode, 10 and 16 fewer in the main assembly); precision at 60 reads is higher in all four sets.
- **U3 HOLDS** (the same loci called expressed).
The cost at the thinnest copies is larger on the fresh seeds than on the sets the configuration was chosen on (0.057
against 0.03-0.04). Everything else confirms. The default as built stands in the working tree, uncommitted; whether it
stays is put back to the user with this number.

## Decision (user, 2026-11-27, asked again with the confirmation's numbers): KEEP the mixed default
`all` with no pool named: the main assembly `PREFIX.gtf` is seeded with `own` (`RUSTLE_SEED_OWN_MAPQ` 0), every other
stage with `good`. `all --seed-pool good` = the default before 2026-11-27. `assemble` named alone is unchanged
(`primary`). Known cost, accepted: isoform recall of the family mode at three-read copies 0.057 lower on fresh seeds.
Whole `cargo test` on the final driver: 1,591 passed, 0 failed. Every figure, base run and page made with `all` before
this date used `good` for the main assembly: rerun before comparing.
