# PREREG — recovering unannotated family members by projecting the family's core (2026-09-07)

**Written before the experiment was run.** md5 in `mcl_ann/adj/projection/PREREG.md5`.
Companion to `PREREG_annotation_ablation_2026-09-07.md` (md5 65efc5b2), whose arms A0/A0r/A1/A2 have run and
whose A3/A4/A5 are deferred by the user in favour of this experiment.

## What this tests, and what it does NOT
O1's node set is the annotation, so a copy with no annotation record cannot be proposed, clustered or scored
(the open gap G5; the failed detection prediction P3 of `O1_DEFINITION_SWITCH.md` §2). This experiment tests
ONE route out: **for a family that already exists, project its core into the genome and admit what comes back
under the unchanged membership rule.** The definition does not change — a family is still a set of loci each
sharing a duplicated core with at least half the family. Only the PROPOSAL gains a second source.

⚠ **Scope.** This can recover unannotated MEMBERS of a family O1 already found. It cannot discover a family
none of whose members is annotated. That case belongs to the read-proposal route (G5) and to O3.

## Design — leave-one-out over the NPIP family
Substrate: the three gorilla contigs. Family: MCL1 of the ablation control `adj/ablation/A0.units.tsv`
(25 members, reproduces the shipped catalog exactly, P0 held). Truth: the 26 LCR16a loci
(`adj/size/lcr16a.bed`). SEDEF (`GGO_sedef_final.bed`) is the annotation-independent pair instrument.

For each member m of the 25 (each fold treats m as if it had never been annotated):
1. **Survivors** S = the other 24 members.
2. **Recompute the cores over S alone** by the shipped majority rule (a position is core when duplication
   pairs link it to ≥ half of the other survivors). m contributes nothing to any survivor's core.
3. **Project**: the survivors' core-hull sequences are aligned to the three contigs
   (`minimap2 -x asm20 -c -N 50 -p 0.1` against `npip3_contigs.asm20.mmi`). Hits are merged per contig.
4. **Discard** merged intervals overlapping any survivor's annotated span (those loci are already members).
5. **Admit** a remaining interval by the SAME membership rule: its duplication depth to the survivors must
   reach half of them, and its core must cover half its span or half the survivors' median core.
6. **Score** against m: recovered when an admitted interval and m's truth locus overlap by ≥ 50 % of either.

## Metrics
- **R** recovery: folds in which m is recovered, out of 25.
- **B** boundary quality of recovered intervals: median coverage of the truth locus, and the fraction with a
  size ratio in 0.5–2× (the pre-registered in-band rule).
- **F** false admissions: admitted intervals per fold overlapping none of the 26 truth loci.
- **K** the core test's keep: candidate intervals per fold before step 5 versus admitted after it.
- **H** the six known hitchhikers (ABCC1, the two ABCC1-like models, the uncharacterised 16.4 Mb block,
  EIF3C, PLA2G10): the fraction of folds in which each is admitted.

## Predictions
| # | prediction |
|---|---|
| **P1** | R ≥ 20/25 |
| **P2** | B: median truth coverage ≥ 0.80, and ≥ 0.90 of recovered intervals in band 0.5–2× |
| **P3** | F: median ≤ 3 false admissions per fold |
| **P4** | H: each hitchhiker admitted in ≤ 20 % of folds — a record the core rule rejects when annotated must also be rejected when projected |
| **P5** | K: candidates before the core test ≥ 3× the admitted count, i.e. the test does the work, not the alignment |

## Interpretation fixed in advance
- P1 + P4 holding ⟹ the annotation is not needed to FIND a member, only to propose one, and the same
  membership rule filters projections as it filters annotations. That is the answer to "over-reliant on the
  annotation" for members of known families.
- P1 holding while P2 fails ⟹ projection finds copies but does not delineate them; the boundary would then
  have to come from the reads, and that must be reported rather than smoothed over.
- P4 failing ⟹ projection admits what the annotated pipeline rejects, i.e. the route is NOT rule-preserving
  and must not ship.

---

## AMENDMENT 1 (2026-09-07, after the run, before any number was quoted)

⛔ **The pre-registered truth is circular for this experiment and P1–P4 as written are VOID.**
`adj/size/lcr16a.bed` is byte-identical, on its interval columns, to `adj/o1loci/npip_loci.tsv`, whose
`unit` column reads `MCL1:<copy>:<status>`: **the 26 LCR16a intervals ARE this family's own core hulls
projected back to the genome.** Scoring a core projection against them is the tautology the register already
carries three times ("prediction ⊆ its own truth"). The run's headline numbers against that file
(R 25/25, coverage median 1.000, F 0) must NOT be quoted, and are not quoted.

**Substituted reference, fixed here before the re-score was read out:** the held-out member's **own annotated
gene span**, which its fold never saw (its annotation is exactly what the fold removes) and which is produced
by an instrument independent of our cores. Metrics: coverage of that gene span by the projected interval,
and the size ratio, with the same in-band rule. The recovery claim itself (does the projection land on the
held-out locus at all) survives the amendment, because the search set excludes the held-out member's core;
what does not survive is any claim about **boundary accuracy** measured on `lcr16a.bed`.

**P5 fails on its own terms**, independently of the circularity: candidates before the membership test and
admissions after it are both 1 per fold. On NPIP the alignment's coverage floor already isolates the family
and the membership test does no additional work. That is reported, not repaired.
