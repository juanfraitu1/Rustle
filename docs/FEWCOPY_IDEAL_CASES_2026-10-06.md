# Few-copy ideal cases for O1, and are NPIP / TBC1D3 "unreachable"? (2026-10-06)

**Status: EXPLORATORY, DEV, human A119b only, stored products only (no pipeline run). This is not a pre-registered test.** The held-out class-level test is drafted in `docs/PREREG_fewcopy_class_DRAFT_2026-10-06.md` and is NOT registered.
User question (2026-10-06): "do we know of any ideal cases of multi-copy gene families (just 1 or 2 with few copies) that would validate our O1? So far it looks like NPIP and TBC1D3 are a bit unreachable right?"
Method: workflow wf_dbfd1e66-dd2 (four read-only sweeps run in parallel, then one adversarial reviewer that recomputed the decision-bearing numbers). Outputs: `/mnt/linuxdisk/tmp/fewcopy_2026-10-06/` (tables, scripts, the five sweep results; the instrument scripts are NOT committed). Frozen input-side criteria: `bench/fewcopy/criteria.md`, sha256 `89363117a885718695d70c6d71a78b75f88c31a62421904c505aeba34c7c1571`.

## Answer

1. **NPIP and TBC1D3 are not unreachable.** On ideal reads TBC1D3 passes its registered bar (12 of 13 reachable copies, both replicates) and NPIP misses by one copy (12 of 14, bar 13); perfect loci (the annotation through `mcl_families`) reach 14 / 14 and 13 / 13, so both bars are reachable by construction (`docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`, register 1254, recomputed by the reviewer from `g4_s*.copies.tsv` and `gates.json`).
   What is real: 10 of 25 NPIP and 2 of 16 TBC1D3 copies share >= 100 exonic bp with another gene (no per-copy-locus bar is demanded for them, yet 4 and 2 of them were found); real-read recovery is data limited (annotation-anchored FOUND, quoted: NPIP 12 real against 23 ideal of 25, TBC1D3 12 against 16); both verdicts sit on a one-copy margin; the simulation is circular and dev-only.
   They are hard stress tests of fused loci and representative choice, and poor clean validators of O1.
2. **Two clean few-copy cases, chosen from inputs only: GSTA1 / GSTA2 (chr6) and PAGE2 / PAGE2B / PAGE5 (chrX).** Both are exact in all three scoring views under the current default. They are positive controls for a NECESSARY condition on DEV data, not validation (details below).
3. **"Ideal by inputs" does not guarantee recovery, and the class is much harder than the two cases.** Of 377 Compara families with 2-4 genes, 11 pass the five input criteria; only 4 of the 11 are exact in all three views. Over all 292 scorable whole 2-4-gene families on 19 A119b contigs: sens .297, prec .935, F .450.

## Input-side criteria (frozen before any count)

Universe: Compara Primates families with 2, 3 or 4 listed genes (N0 = 377: 292 / 60 / 25), members mapped RefSeq to CAT/Liftoff v2.0 by shared exonic bp (never by name). A family passes only if EVERY member or pair passes.

| step | criterion | cumulative families |
|---|---|---|
| C1 | every member maps strongly to one protein-coding CAT gene; reference transcript >= 3 exons, every intron >= 50 bp and canonical | 129 |
| C2 | not entangled: < 100 exonic bp shared with any other annotated gene | 62 |
| C3 | separate loci: spans >= 1 kb apart | 61 |
| C4 | homologous, not identical: spliced reference transcripts, 0.90 <= identity < 0.999, coverage >= 0.50 (minimap2 asm20) | 35 |
| C5 | observable: >= 3 reads per member carrying its exact annotated chain (E0 rule, primary + good secondaries, region queries of `A119b.t2t.bam`) | **11** |

One code amendment after the first run (spans = exon-union extent; 3 flags changed, the 11 and their order did not); no threshold edited; the criteria hash was re-checked at the end and by me (identical). 1 of 864 members could not be mapped (FAM86B1). 57 families fail exactly one criterion (C5 24, C2 16, C4 11, C1 6; `near_miss.tsv`).

## The 11 families that pass, and what the current default did

Outcome side: stored current-default products (arm D, f1v2 bridge regroup + `--min-cov-shorter 0.70` + GOOD secondaries), scored with `family_score` (S1 = its own view; S2 = any other annotated gene in the cluster counts as a non-member; S3 = membership by >= 1 exonic bp). The outcome sweep was kept blind to the input sweep.

| family | members | contig | identity (min) | min reads | tie fraction | S1 | S2 | S3 | all three |
|---|---|---|---|---|---|---|---|---|---|
| CF150 | HBA1, HBA2 | chr16 | .967 | 369 | .004 | exact | exact | exact | **yes** |
| CF316 | GSTA1, GSTA2 | chr6 | .954 | 111 | .000 | exact | exact | exact | **yes** |
| CF334 | GTF2IRD2, GTF2IRD2B | chr7 | .993 | 18 | .209 | exact | exact | exact | **yes** |
| CF400 | PAGE2, PAGE2B, PAGE5 | chrX | .918 | 12 | .000 | exact | exact | exact | **yes** |
| CF63 | ZNF33A, ZNF33B | chr10 | .921 | 11 | .000 | exact | exact | partial + non-member | no |
| CF9 | FCGR3A, FCGR3B | chr1 | .977 | 8 | .000 | exact | merged | merged | no |
| CF403 | RHOXF2, RHOXF2B | chrX | .998 | 139 | .730 | partial | partial + non-member | exact | no |
| CF355 | PRR23D1, PRR23D2 | chr8 | .996 | 36 | .470 | partial | partial + non-member | merged | no |
| CF236 | ZNF600, ZNF611 | chr19 | .912 | 4 | .000 | merged | merged | merged | no |
| CF54 | LYZL1, LYZL2 | chr10 | .980 | 59 | .000 | no-locus | no-locus | no-locus | no |
| CF301 | UGT2B15, UGT2B17 | chr4 | .966 | 7 | .000 | no-locus | no-locus | no-locus | no |

Pooled over the 11: members hit 17 / 23, sens .739, prec <= .944, F .829, pairwise 9 / 13. "No-locus" here means no CLUSTERED locus: LYZL1/2 and UGT2B15/17 have loci and reads but no admitted edge (the 09-25 annotation-loci arm also fails LYZL, 0 / 2). Those two are O1-definition cases, more informative than the successes.
Products: HBA1/2 comes from a HEAD run (chr16); the other ten are replays of the 09-30 e163d955 arm D, and HEAD equals it byte for byte only on chr16, chr17 and chrY.

## Class-level numbers (always print beside any subset)

| set | families | sens | prec | F | pairs |
|---|---|---|---|---|---|
| ALL scorable whole 2-4-gene families | 292 | .297 | <= .935 | **.450** | 111 / 493 |
| size 2 / 3 / 4 | 223 / 48 / 21 | .298 / .333 / .226 | .964 / .873 / .905 | .455 / .482 / .362 | 54 / 223, 34 / 144, 23 / 126 |
| the 11 that pass the input criteria | 11 | .739 | <= .944 | .829 | 9 / 13 |
| all minus the 11 | 281 | .281 | <= .934 | .432 | 102 / 480 |
| the 42 exact in all three views (selected on the OUTPUT) | 42 | 1.0 | 1.0 | 1.0 by construction | do not quote alone |

Exactness over the 292: 63 (22%) in the scorer's own view, 42 (14%) in all three views; 162 (55%) have no member in a clustered locus (63 have no locus at all, 43 have every member only in unclustered loci). 78 of the 377 families cannot be scored (55 span contigs, 23 sit on chr13, 18, 20, 21 or 22, which have no stored current-default product); 7 more are scored on a restricted piece.
Selecting on inputs helps (exact in all views 4 / 11 = 36% against 38 / 288 = 13% for the rest) but is no guarantee. Post hoc, dev, by the minimum exact-chain reads per member: all-three-exact 10 / 193 at 0 reads, 14 / 70 at 1-29, 6 / 16 at 30-99, 7 / 12 at 100-299, 5 / 7 at >= 300.
Holder-rule view (the prior-evidence sweep, A119b contigs without chr16; a different instrument from the table above): of 118 pairs with BOTH genes in a locus, 63 are complete, 40 have no admitted edge, 9 sit in one locus, 4 have an edge that is not joined; this is the edge-admission loss of the earlier §6o8 bound, not two-node MCL numerics (which explain at most 1-4).
Chain fidelity (review sweep, not rechecked by me): only 21 of the 42 all-three-exact families have every member's locus chain identical to a RefSeq transcript; of the 11, only HBA, GSTA, GTF2IRD2 and PAGE are exact in all views AND chain-exact. NOMO1/2/3 loci carry 1-2 introns.

## The two recommended cases, and what they do and do not show

- **GSTA1 / GSTA2** (CF316, chr6, tandem, same strand, 27.9 kb apart; identity .954, coverage .99; 0 bp entanglement; 111 / 115 exact-chain reads, all primary (review sweep); no read-level link between the copies, tie 0.000): exact in S1, S2 and S3, both loci chain-exact (6 introns), edge weight .929 from sequence and exon structure alone; a second human library (genome-wide testis, holder rule) gives the right family with 330 / 204 reads (both from the review sweep, not rechecked by me). chr6 is the least-exposed contig in the ledgers (quoted).
- **PAGE2 / PAGE2B / PAGE5** (CF400, chrX, same strand, copies >= 10.1 kb apart; identity .918-.973; 0 bp entanglement; tie <= .07): exact in all views, every locus chain-exact (4 introns), complete triangle (MCL density 1.0); second library right-family with 165 / 167 / 114 reads (review sweep). A119b depth is thin (minimum 12 reads).
- **They validate** a necessary condition: reads give one read-supported locus per copy, and the copy graph joins them into exactly one family with nothing else attached, on a clean, expressed, mid-divergence pair or triple.
- **They do not validate** sufficiency or recall (the class gives F .45), ties or identical copies (O2), unlisted relatives, low depth, other species, or anything blind: the outcome has been seen, so both are DEV regression anchors. Both outcomes are replay products, not HEAD runs. n = 1 each; a "recent duplication" premise is unverified.
- Not for use as "ideal": identical-copy pairs (EIF3C/EIF3CL identity 1.000, MAGEA9/9B, HSFY1/2) recover because both loci are seeded by the same tied reads (O1-trivial, useful as O2 tie controls); LYZL1/2 and UGT2B15/17 fail even with annotation as loci. HBA1/2 and GTF2IRD2/2B stay as regression anchors on the dev chromosomes.
- The literature screen (post hoc tiers, not pre-registered) proposed GPR89A/B, ROPN1/ROPN1B and WASHC2A/WASHC2C as separable twins; none passes the frozen criteria (C1 annotation quirks: partial map, a 2-exon CAT reference transcript, a non-canonical intron). The pipeline recovers GPR89A/B (CF11) and WASHC2A/C (CF62) exactly; ROPN1/1B (CF292) is partial.

## Prior evidence and held-out status

- Earlier few-copy cases (gorilla 06-20 to 08-19: RABL2A/B, AK6, CCDC196, MAGEA pairs, GSTM, DAZ, RBMY, TSPY, PCDHB, RFPL, APOBEC3; famsim SNRPB / SPAG7 simulations) all predate MCL, `--min-cov-shorter 0.70`, f1v2 and GOOD seeding; none has been re-run at the current default. No few-copy family has been pre-registered as an O1 positive control for de novo real reads at the current default (earlier pre-registrations were guided or chrY).
- **No A119b contig is unexposed** (flips were decided on products that include them: chr16 is the dev chromosome of every fusion rule, chr7 Soto-dev, chr20-22 F1v2-dev), so the 11 are DEV. The drafted held-out test uses gorilla OR6737 testis (never pooled with human), or a never-run human library (HG002: availability unverified).

## Limits and corrections the review forced

- The outcome product set: 16 of 19 contigs (chr1-12, 14, 15, 19, X) are replays of the 09-30 e163d955 arm D; replay = real run was checked only on chr17 and chrY. chr13, 18, 20, 21, 22 have no product.
- The scorer is strand blind and resolves one gene per locus (nested genes mislabel loci: the TRIM43/TRIM43B loci are labelled GPAT2); 20 families are "hidden" under S1. S3 is circular (its floor, >= 1 exonic bp, is the pipeline's own edge floor). "Exact" is unreachable where the listed size hides relatives (ZNF600/611 with ZNF808 and ZNF888 are not false merges).
- Reads: C5 uses the pipeline's own seeding pool, so for near-identical copies (identity >= .995) both loci are seeded by the same tied reads and recovery is weak evidence; the informative cases are the 8 no-tie families (cross-link <= .07). Selection on inputs is positive-heavy: an O1 score on it is a mechanism check, not a recall estimate.
- Sweep claims withdrawn or corrected: "no few-copy family was ever pre-registered" is true only for de novo real reads at the current default; the prior sweep's "best candidates (exact)" list is outcome-selected (8 of 10 fail the input criteria; RHOXF2/2B and ROPN1/1B are not exact under the scorer); the literature sweep's "93-100% of primaries also align at the other copy" does not show ties (the library was mapped with `-N 50 -p 0.1`; HBA1 to HBA2 is 398 / 407 with any secondary and 1 / 407 with a good secondary); its tier A is post hoc and none of its four families passes the criteria; the outcome headline "22%" is the inflated S1 figure (42 / 292 = 14% over all views).
- Metric traps hit and handled (reviewer): selection on outcome (F .881 on 11 chosen against .628 on all 72), denominators conditioned on outcome (80 families with all members in a clustered locus give 52 exact = 65%, against 52 / 223 = 23%), universe intersection, circular S3, unreachable "exact", threshold fragility (ZNF600/611 min reads 4 against 3; identity .912-.921 against the .90 floor), non-independence (ZNF33 and ZNF600/611 are KRAB-ZNF; PAGE and RHOXF2/2B are chrX arrays: the 11 are about 9 independent groups).
