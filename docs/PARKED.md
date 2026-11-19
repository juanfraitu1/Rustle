# Parked items (to pick up later, on any cluster)

Each item says where it stands, what is in the repository, what is not (workspace files of the Fritz cluster do not travel), and
the first steps. Newest first. Remove an item when it is closed in `NEGATIVE_RESULTS_REGISTER.md`.

## 5. Tandem arrays: member size (2026-11-19, goal cleared by the user, row 1826)

**Stands at.** `TANDEM_ARRAYS_2026-11-19.md` (sections 1-8), `PREDICTION_tandem_heldout_2026-11-19.md`. Established: StringTie joins
neighbouring units through reads longer than one period (removing them: 50 -> 4 multi-unit models in A119b's regular arrays,
protocadherins unchanged); reads split across units are rare (5 of 15,115); the read classifier (antisense and > 1-period reads set
apart) keeps 330 / 330 real reads on simulated arrays and sets apart every antisense and read-through read. Not established: a
member-size rule for real arrays. The read-end rule (`tandem.py member_sizes`) is exact on simulations and right on near-identical
arrays (CT47A, TSPY in both species) but failed both held-out tests on OR6737 (0.60 vs StringTie 0.80 / Rustle 0.75; 0.55 vs 0.55 /
0.66).
**In the repository.** `bench/famwise_eval/{tandem,tandem_heldout,tandem_eval_annot,array_inject}.py`.
**Not in the repository.** PSU workspace `goal_multicopy_dev/{tandem,arrays_sim,heldout}` (linked from `r1826_tandem_arrays`).
**To do.**
1. A truth better than the longest annotated isoform: a testis Iso-Seq sample of an individual with assembled haplotypes, or
   full-length (capped) reads per member.
2. OR6737 and A119b are spent for this rule; any new variant needs another sample. Do not re-tune on them.
3. Optional, cheap: drop reads longer than one period inside regular arrays before StringTie (removes its unit joins; member sizes
   on A119b 56 -> 59 of 131 within 20%).

## 4. Correlation clustering for copies the reference lacks (2026-11-19, step 0 done, row 1823)

**Stands at.** `CORRELATION_CLUSTERING_2026-11-19.md`: copy_phasing (the Strainy port, correlation clustering) against HG002's
own-haplotype truth: no gain for assignment; for O3 its certified flag hits 13 / 13 families with >= 6 `outside` molecules, with
false flags rising with depth. Development only (all HG002 quarters spent).
**In the repository.** `bench/famwise_eval/cp_truth.py` (inputs from the truth families; O2 / O3 scorer); `src/copy_phasing/`.
**Not in the repository.** PSU workspace `goal_multicopy_dev/cc`, `cc_all` (linked from `r1823_correlation_clustering`).
**To do.**
1. Split `outside` molecules into a copy CHM13 lacks vs an HG002 allele too far from CHM13's copy (the truth records the HG002 place).
2. A judged O3 row with a held-out own-haplotype truth (KB3781 on its own reference; or HG002 reads against `hg002.pat.fa`, where
   maternal-only copies are the positives), bars on recall and on false flags WITHIN depth bins.
3. A depth-aware certificate (cluster size scaled to the pile) if the false flags hold up.
4. Untried variants: Paraphase-style collapse onto one archetype copy; SDA's PSV-level graph. For O3, not for O2.

## 3. Multi-copy families after rows 1820-1822 (parked 2026-11-19, user cleared the goal)

**Stands at.** `MULTICOPY_FAMILIES_2026-11-19.md` (summary), `GOAL_multicopy_families_2026-11-19.md` (the goal and its status),
`MULTICOPY_DEV_LEDGER_2026-11-19.md` (every development number). Rows 1820 (C1 fails), 1821 (R1 fails narrowly, R2 / R3 pass),
1822 (K1, B1 fail on KB3781).
**In the repository.** `RUSTLE_SEED_SPLICED=support` (opt-in), `bench/famwise_eval/{copyset_score,resolution_call,borrow,
family_isoform,member_completeness,member_bed,support_filter,array_score,sim_truth,undecidable}.py`.
**Not in the repository.** PSU workspace folders `r1820_undecidable_copies`, `r1821_resolution_call`, `r1822_borrow_isoformk`,
`goal_multicopy_dev` (rsync line at the end of the summary).
**To do.**
1. The member gap: Rustle drops low-support family members that StringTie keeps (19 of 24 on the HG002 truth quarter get no
   transcript; row 1536's class). Judge a lever on a sample with room (A119b has 835 expressed members, but no truth).
2. Shared isoforms (`family_isoform.py`) as a promotion rule for one-read chains: 17 of 19 right on KB3781, too few. Needs
   bars; KB3781 with `DENOVO_GGO_ENS_MAP` on Fritz gives a larger family truth.
3. Tandem arrays on real data: needs a testis Iso-Seq sample of an individual with assembled haplotypes.
4. `RUSTLE_SEED_SPLICED=support` as a default: the user's call (row 1821: wrong-copy new chains 60 -> 15).
**Data.** No unseen HG002 quarter remains.

## 1. Junction replacement by support ratio, inside the assembler (parked 2026-11-19, user)

**Stands at.** Row 1819, `PREREG_junction_replace_2026-11-19.md`: as a post-step on `PREFIX.gtf`, replacing intron K by the
best-supported intron J within 10 bp at both ends when support(J) >= 20 x support(K) passes its bar on seven datasets (at least
0.995 of the correct chains kept, precision higher on all; wrong chains -3.3% to -7.8% on fly, mouse simulated, mouse ES, heart,
adrenal; -1.2% / -1.4% on A119b / chimpanzee testis; SIRV 53 / 0.589 -> 54 / 0.675). It was only analysed: no GTF is written and
nothing in the assembler changed.
**In the repository.** `bench/famcmp/junction_replace.py` (the analysis; needs `bench/famcmp/curve.py` only). The existing snap:
`--polish-junction-snap off|equiv|reads` in `src/copy_assign_bin.rs` (doc comment at the option; driver
`RUSTLE_POLISH_JUNCTION_SNAP`, default `reads` since row 1800). `reads` moves K onto J only when the reads cannot tell them apart
and uses no support ratio.
**Not in the repository.** The assemblies and tables (Fritz workspace `r1819_junction_replace`, `r1817_transfer_rerun`). On
another cluster: run `tools/rustle_pipeline.sh assemble` at the defaults on the datasets at hand, then the script.
**To do.**
1. Add the ratio condition to the snap as a new opt-in value or option (K is also moved when support(J) >= R x support(K), R
   default 20; support as the snap already pools it), with a test next to the existing snap tests.
2. Build into a separate target directory. The shared Fritz binary (`~/rustle_target/release`, built 2026-11-18 04:20) is older
   than the source: a rebuild also brings in `--polish-spliced-share`, `--polish-mono-shadow-ratio` (row 1801) and the family
   catalog / copy phasing changes of commit 2dfa590d. First confirm that a default `assemble` is unchanged with the new
   condition off (md5 of `PREFIX.gtf` and `PREFIX.final.gtf` against the old binary on one small dataset).
3. Pre-register, then run end to end on the datasets of row 1819 with `PREFIX.gtf` and `PREFIX.final.gtf`; the post-step numbers
   above are the expectation, but inside the assembler the read floor and the later rules see the merged counts.
**Known cost.** Real splice sites a few bases apart: 11 to 160 correct chains per dataset moved onto another correct chain.

## 2. Open after the default changes of 2026-11-18 / 19 (rows 1800-1816)

Today's `assemble` defaults differ from 2026-11-17 in five ways (junction snap on, sub-chain ratio 0.3, extension off, `assemble`
alone seeds with primary alignments, filter = chain score r1813 + 5' truncation rule; `OVERFIT_LEVERS_2026-11-18.md`). Not yet done:
- **Family benchmark at these defaults** (headline pair precision / recall): every family figure predates them.
- **Advisor page and figures:** show the state before row 1793.
- **Row 1801, SIRV rules** (`PREREG_sirv_rules_2026-11-18.md`): draft, bars not accepted; two of its arms need the newer binary.
- **New evidence for the wrong chains that are left** (row 1818): splice-site sequence and alignment quality at a near-miss
  junction; 3'-end support for pieces of longer chains. The 23 read features of `chain_rules_feat.py` are used up
  (rows 1814, 1815, 1818).
- **Promotion of one-read chains** (row 1808) and `all` end to end were not run through the driver at these defaults.
