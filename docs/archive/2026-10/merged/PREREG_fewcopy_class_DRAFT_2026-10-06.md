# DRAFT pre-registration: class-level few-copy O1 test on a held-out substrate (2026-10-06)

**STATUS: DRAFT v0. NOT REGISTERED. Nothing in it has been run.** Written by the adversarial reviewer of wf_dbfd1e66-dd2 before any held-out product was opened; it becomes a pre-registration only when the decisions below are made, its numbers are frozen and it is committed as `docs/PREREG_fewcopy_class_<date>.md`. Amendments are allowed only before outputs exist.
Background and the dev results: `docs/FEWCOPY_IDEAL_CASES_2026-10-06.md`. Frozen input-side criteria (C1-C5): `bench/fewcopy/criteria.md`, sha256 `89363117a885718695d70c6d71a78b75f88c31a62421904c505aeba34c7c1571`.

## Decisions needed from the user before registering

1. **Substrate.** Gorilla OR6737 testis (primary choice of the draft; truth = Ensembl Compara gorilla paralogs, else human Compara 2-4 families projected by orthology as in the 07-29 bench), or a never-run human library (HG002: availability unverified), or human chr13 + chr18 (3 families: UNDERPOWERED by rule 6). Gorilla KB3781 is reported separately, never as independent of OR6737.
2. **Bars.** E1 uses 0.9 (as `PREREG_ideal_expression_2026-10-06.md`) and a NO threshold of 0.7; 0.7 is the one free parameter.
3. **Whether a runner must be written and reviewed first.** The dev instrument scripts live in `/mnt/linuxdisk/tmp/fewcopy_2026-10-06/` and are not committed; a held-out run needs a committed, reviewed `bench/fewcopy/` runner.

## Draft text

0. **Status of dev data.** The 11 human A119b families that pass the criteria and all 377 outcomes are SEEN = DEV / regression. No A119b contig is unexposed to the current default, so none is called held-out.
1. **Held-out (never pooled with human).** Primary: gorilla OR6737 testis (GGO RefSeq GFF); truth = Ensembl Compara gorilla paralogs frozen (sha1) first, else human Compara 2-4 families projected by orthology. Fallbacks as in decision 1.
2. **Class K** = families passing C1-C5 of `criteria.md`, computed on the held-out annotation and BAM from inputs only; no threshold edit or relaxation after counts. Near-miss families are listed, never put into K. |K| per size is printed before any outcome is read.
3. **Controls fixed before outcomes.** (G) annotation-as-loci arm: the same `mcl_families` flags with loci = annotation and no reads; a family is REACHABLE iff G recovers it exactly in S1, S2 and S3. (Neg) per family, 5 matched negative pairs (protein-coding, in no multi-copy Compara family, no >= 100 bp exonic overlap, same contig, gap >= 1 kb, reads >= the family's minimum), seeded.
4. **Arm N** = de novo (`assemble` with f1v2, then `mcl_families` defaults) at frozen binary sha1s, run once per library. Scorer: `family_score` (scipy tie policy named) plus pairwise. "Exact" = S1, S2 and S3 AND every member has its own distinct locus (a locus shared by two members is a miss). Chain fidelity is descriptive only.
5. **Endpoints.** E1 (primary): x = reachable families exact in N, n = reachable families; YES >= ceil(0.9 n); NO < ceil(0.7 n); else PARTLY (exact binomial 90% CI printed). E2: false merges on Neg: N must not exceed G. E3: A-class counterexamples = reachable families whose members all have >= 100 P2 reads and tie <= .07 and that N does not recover exactly; >= 2 means NO. Reported without a bar: sens / prec / bipartite F and pairwise over K and over ALL truth families (always beside K), per size and per depth tier (minimum P2 3-29 / 30-99 / >= 100), G-exact over all K (the definition ceiling), and the gap G minus N (guided against de novo).
6. n_reachable < 8 means UNDERPOWERED: no verdict, and it is not a pass.
7. **Failure** = E1 NO, or E2 violated, or E3 triggered. No rescue by changing instrument, tiers, thresholds, substrate or by re-running; deviations are logged as amendments.
8. **Reading a YES.** It means a necessary condition holds on clean, expressed few-copy families in that species and library. It is not sufficient (ties / O2, large families, other species, unlisted relatives). Dev evidence for planning: all-view exact 4 / 11 on the criteria-passing set, 42 / 292 over all scorable families; the annotation arm equals de novo on the 5 pairs with a 09-25 guided row.
