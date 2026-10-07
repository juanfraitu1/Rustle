# Disposition of the three adversarial reviews of the T1/T2 pre-registration (2026-10-07)

Reviews: workflow wf_20046203-7d9 (three reviewers: wording 37 findings, logic 29, feasibility 26; 92 in all). The draft reviewed is the first draft of `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`. The revised file is the same path (DRAFT 2). The review results as JSON are `/mnt/linuxdisk/tmp/prereg_review_wording.json`, `prereg_review_logic.json` and `prereg_review_feasibility.json`; the numbers below are the positions in each file's `findings` list (wording W1 to W37, logic L1 to L29, feasibility F1 to F26).

**Closure review of draft 2** (workflow wf_eea75a23-1c5, 2026-10-07; four read-only reviewers on frozen copies): logic 25 findings, facts 22, wording 36, disposition audit 19 (the JSON files are `/mnt/linuxdisk/tmp/prereg_closure_{logic,facts,wording,disposition}.json`; scratch `/mnt/linuxdisk/tmp/prereg_closure/`). It found two blockers (both the Gate 7 positive control, one from the facts reviewer and one from the disposition audit) and 25 majors. The revised file is DRAFT 3; the closure dispositions are the second table below. In the first table, rows that the closure review found only partly fixed (F1, W15, W17, W22, L9) point to the closure rows that complete them.

Dispositions. ACCEPTED: fixed as suggested or in substance. ADAPTED: fixed in a different way, and the row says how and why. REJECTED: not fixed, and the row gives the evidence. NOTED: no change needed. Section numbers refer to the revised prereg. Cross-references such as 'wording 3' point to another row of this table.

## Counts by disposition

| reviewer | findings | ACCEPTED | ADAPTED | REJECTED | NOTED |
|---|---|---|---|---|---|
| wording | 37 | 27 | 10 | 0 | 0 |
| logic | 29 | 22 | 5 | 0 | 2 |
| feasibility | 26 | 19 | 6 | 0 | 1 |
| all | 92 | 68 | 21 | 0 | 3 |

## Dispositions, one row per finding

| reviewer | # | severity | location | disposition | where fixed, or why not |
|---|---|---|---|---|---|
| wording | W1 | blocker | Gate 5: G1 does not refine G0; m1c not monotone; different gene sets | ACCEPTED | Section 3 Gate 5. Gate 5a on the common universe (U1 to U2, once per famCN source), Gate 5b with the numeric bar D[S->Q] = 7 (2 if ungrouped genes are ignored); U0 dropped from the chain; the gate prints PASS or FAIL and no value. The 32 strings / 110 genes and the 10 genes were reproduced (logic/string_only_checks.out). |
| wording | W2 | blocker | T2b identity check has the wrong sign and strictness | ACCEPTED | Section 9 Identity check (h is a divergence; h(E) <= h(phi), non-strict; a violation is INVALID); age_id deleted; section 2 item 6 uses the same h; consequence C2 replaces P8. |
| wording | W3 | major | Audit not closed (deficit set, ends-only, D2, counts) | ADAPTED | Section 4 Audit: tie rule for the kept chain (lexicographically smallest maximum-weight chain, all tied chains printed), ends-only defined, D2 labels, a label never removes a gene, counts quoted in V, non-eligible defined by the S1C column. Differs from the suggestion: the audit is descriptive and is not the decisive statistic, because 'short' alone excuses 862 of 2,142 genes (finding 34) and excuses are asymmetric between data and null. Decision 3 in section 12. |
| wording | W4 | major | Gene-to-family rule, layer universes, ratios | ADAPTED | Section 4 Family layers: rules 1 to 3, V_F = 1,657, layer universes, every ratio printed with its universe. Differs: unclustered genes are outside the decisive cells (feasibility finding 6 recommends this) and enter a descriptive singleton arm, instead of being singleton families by default. |
| wording | W5 | major | Verdict algebra forks | ACCEPTED | Section 3 Verdict algebra: precedence, p only with the same cut-offs everywhere, UNINFORMATIVE defined, n_min, SPLIT rule, HOLDS removed. Class names changed (ORGANISES-ENRICHED, NOT-DISTINGUISHED). |
| wording | W6 | major | Decisive cells and expected answers | ACCEPTED | Section 3 Decisive cells table; section 10 (C1 to C5 consequences, Q1 to Q8 predictions for T1-a, T1-b, T1-c, T2a, T2b). |
| wording | W7 | major | Nulls not reproducible (rungs, singletons, strata, seeds) | ADAPTED | Section 4 Nulls and Draws: re-seed per (layer, null, universe), sorted genes and strata, draws shared only inside one universe, singletons outside the decisive cells. Strata are the homology-preserving components (logic finding 2) and not the chromosome. |
| wording | W8 | major | 'Per half' rule | ACCEPTED | Section 4 Half (gene-level column of the split file; halves A and B descriptive, never called held-out). |
| wording | W9 | major | Unit-string construction | ACCEPTED | Section 4 Exon label and unit strings (bases added, argmin key, G0 key, coordinates as stored); Gate 6 prints tie counts; the set-Jaccard sensitivity is deleted. Counts from input_checks.out. |
| wording | W10 | major | Gate 1 may have no exit (2,259 vs 2,290) | ACCEPTED | Section 3 Gate 1 requires the registered outputs only; section 2 item 9 records the erratum; section 12 item 13. Two more explanation attempts failed (c5, c5b). |
| wording | W11 | major | Gorilla p-value undefined | ADAPTED | Section 5: no control draws and no p-value (arm descriptive); family-cluster bootstrap with seed 20260930; family-weighted S-rate printed. |
| wording | W12 | major | Gorilla pair universe not closed | ACCEPTED | Section 5 Family layer and pair universe (touch rule, U(g), eligibility, 0 shared bp for every pair, distance, depth-matched, verdict set). Counts re-derived (c7_depth_matched_verdict.out). |
| wording | W13 | major | T2 staging not a partition; decisive stage underpowered | ADAPTED | Section 7 Contig roles: the hold-back unit is the expansion (LRPAP1 only), the ten contigs are a least-exposed subset without a class, accession and chr table, 'held-out' not used for gorilla, expected sizes (c8). Differs from the suggested Confirm/Run-last rule, which would leave 3 and 8 components. |
| wording | W14 | major | T2a class statistic and test undefined | ADAPTED | Section 8 Statistics, Matched controls and test, Classes. Deciding statistic is the mean R_j with absent members counted as non-co-members (logic finding 11), B, seed, two one-sided p-values, PARTIAL, size-2 stratum. |
| wording | W15 | major | T2a matching and scoring | ACCEPTED | Section 8 Members and matching and Expressed (ties, one-to-one node rule, X(n) source, overlap length, m_j, R_j). |
| wording | W16 | major | E-table definition not closed | ADAPTED | Section 7: hub cap 51, expansion = same triple (an equivalence relation, no chaining), single Liftoff pass with -copies (the second pass is dropped), strand ignored, deciding E-table is tier B only, columns. |
| wording | W17 | major | T2b decision statistic | ADAPTED | Section 9: no combined verdict; arm 1 is an estimate; arm 2 has classes with a comparator; mixed-family defined; 0.5 kept only as a majority bar. |
| wording | W18 | major | T2b arms cannot be run as written | ACCEPTED | Section 9: synteny class inlined, member classes, derived set = E minus ancestral members, arm 2 command line, rooting and UFBoot 95, minimum tips, R(D) defined, d_orth deleted. |
| wording | W19 | major | What 'go' authorises | ACCEPTED | Section 3 Run order and authorisation (table and rule); E3 deleted; one call budget. |
| wording | W20 | major | Wording changes are buried | ACCEPTED | Read first: items 1 to 7, the cell table, and the statement that the two sentences are my wording of the offer. |
| wording | W21 | major | Matched-instrument human run: no join, 'agrees' undefined | ADAPTED | Section 5: the run is optional and descriptive, the join is specified (Name and Contig with a uniqueness assertion; offsets verified by sequence), 'agrees' deleted. |
| wording | W22 | major | Gate 7 cannot fail; negative control | ADAPTED | Section 7 Gate 7: bar of 8 loci in one class, can fail, fragments are a sibling class; the negative control is a planted set of 50 fake pairs (at most 1 may be called) instead of the suggested symbol lookup, which is vacuous for a pool of copy-bearing records. |
| wording | W23 | major | Gate housekeeping | ACCEPTED | Section 3: scope of the gates, Gate 2 pass bar, Gate 5b number, Gate 6 definitions, T2 controls in Gate 4, E-table md5 recorded when section 7 closes. |
| wording | W24 | major | Symbol collisions | ACCEPTED | Renamed throughout: U0 to U2, Gate n, Q1 to Q8 (C1 to C5), Mult, pool, E-table, 'gorilla', R_J contig. |
| wording | W25 | major | Genome-wide tables exist | ACCEPTED | Section 2 item 5 corrected; section 6 R3a admits the human_testis BASE snapshot (the author's choice made, decision 5). |
| wording | W26 | minor | Numbers that do not tie out | ACCEPTED | Corrected in section 2 items 7, 8, 10, section 4 Universe V, section 5 counts (38 pairs in 33 families), section 5 tau rule, section 10 C1. |
| wording | W27 | minor | Convention gaps (KEY, results hook, G0 list, hashes) | ACCEPTED | Section 13 (names, scorers, results hook, Amendment hook, frozen inputs with sha256, atoms hash, BASE binary commit stated as not recorded); traps agent scratch path added to the header. |
| wording | W28 | minor | Jargon | ACCEPTED | Terms paragraph; 'Amendment-free' deleted; 'in-paralogs' replaced; A_viol described in section 5. |
| wording | W29 | minor | Randomisation and constants | ACCEPTED | Section 3 Gate 3 (seeds, per-cell re-seed), section 12 item 10 (constants, sweep, inherited list), call budget in section 3. |
| wording | W30 | minor | Secondary statistics | ACCEPTED | Section 4 Statistics (macro, ordered string, simple graph, Mult, m1c = m1 at U0, nesting wording). |
| wording | W31 | minor | Input handling (identical tuples, sensitivities, R3 sketch) | ACCEPTED | Section 5 identical tuple and tie rule; sensitivities defined or deleted (section 4); section 6 R3 fixed by Amendment 1. |
| wording | W32 | minor | Over-claims (P3, 'literal T1', header) | ACCEPTED | Section 10 Q2; section 1 D1 wording; header sentence on what was joined. |
| wording | W33 | minor | Markdown layout | ACCEPTED | Blank lines between paragraphs throughout. |
| wording | W34 | note | 'short' label explains 40% of V | ACCEPTED | Section 4 Audit prints the unexplained count with and without 'short'. |
| wording | W35 | note | V removes the 149 multi-family genes | ACCEPTED | Section 4 Universe V and Gate 6 (per-layer counts of clusters holding them and members deleted). |
| wording | W36 | note | Tie rule differs from the 09-30 text | ACCEPTED | Section 4 Exon label and unit strings states it. |
| wording | W37 | note | T1/T2 name collisions | ACCEPTED | Read first item 7; section 13 KEY names. |
| logic | L1 | blocker | Gate 5 false premise | ACCEPTED | Same as wording 1. Section 3 Gate 5a and 5b; D[S->Q] = 7 and 2 reproduced (logic/fq_checks.out). |
| logic | L2 | major | Tautology; both nulls destroy homology | ACCEPTED | Section 4 Nulls: governing null = components of the F_O graph. The frozen binary re-run with --dump-graph reproduces the frozen products byte for byte (c1_fo_graph.out). Chromosome and SD98-region nulls become sensitivities (C4). 30 variation-bearing strata with 416 genes (c2_fo_components.out). F_Q has no variation under its own strata and stays descriptive. |
| logic | L3 | major | Audit not tie-invariant; D2 undefined | ADAPTED | Section 4 Audit. The suggested single optimisation with excuses is not used inside the decisive statistic (asymmetric excuses; 'short' excuses 40% of V). Adopted: unique chain value, tie rule, D2 labels, ends-only reference. Section 12 item 3. |
| logic | L4 | major | HOLDS-EXACT unreachable; HOLDS-ENRICHED mislabels a refuted universal | ACCEPTED | Sections 3 and 4: T1-a is REFUTED or NOT-REFUTED and is printed first; the enrichment class is ORGANISES-ENRICHED; NOT-REFUTED only on the primary arm. |
| logic | L5 | major | Classes not mutually exclusive | ADAPTED | Section 3 Verdict algebra (precedence; 'direction differs' and '1st percentile' deleted; p only; halves descriptive; SPLIT between Phase R1 and R3a (and, from draft 3, between the two T2a libraries)). UNINFORMATIVE is defined as a null that cannot move. The independent unit is the variation-bearing stratum, so UNDERPOWERED can in principle fire on the human cells (30 against 8). Draft 3 (closure row CL10) extends SPLIT to the two T2a libraries and adds the headline rule. |
| logic | L6 | major | Gorilla arm has no unconfounded decision statistic | ACCEPTED | Section 5, option 1: descriptive report with intervals; enrichment verdict dropped; reasons stated; multi-class genes in the primary and flagged; 'agree' and 'rule frozen on human' deleted. |
| logic | L7 | major | Gorilla layer is a name join | ACCEPTED | Section 5 Family layer (one name step, declared) and section 11 Register check (lines 316, 861, 869; rows 993, 1034, 1258); matched human run optional with printed counts. |
| logic | L8 | major | T2b does not test 'younger than the root' | ACCEPTED | Section 9: arm 1 is root membership with the Dollo definition, arm 2 is a rooted clade test; d(E) < d(F) deleted; the sentence that the clock is only partly independent of the family builder is added. |
| logic | L9 | major | T2b classes conditioned on outcomes | ACCEPTED | Section 9: a non-clade counts as not younger; unit exclusion only by input properties; shares defined per arm; funnel printed; pilot unresolved rates corrected to 83% and 75%; no combined verdict replaces FAILS, PARTIAL and SPLIT. |
| logic | L10 | major | Identity check sign and strictness | ACCEPTED | Same as wording 2 (section 9). |
| logic | L11 | major | T2a estimator, classes, control | ACCEPTED | Section 8: one-to-one node rule, the deciding statistic named, the test defined, NOT-DISTINGUISHED replaces AS-CONTROLS, parent-family size reported as a covariate. The statistic is R_j with absent members as non-co-members (not E2E_j). |
| logic | L12 | major | Deciding T2 subset expected UNDERPOWERED | ADAPTED | Section 7 Expected size; section 10 Q5 and Q7. The hold-back unit is the expansion, so the verdict set is every expansion except LRPAP1 (about 12 and 21 expected) and the ten-contig subset (3 and 8 components) carries no class. |
| logic | L13 | major | Gate 1 open-ended | ACCEPTED | Same as wording 10. |
| logic | L14 | major | D1 mislabelled as the literal T1 | ACCEPTED | Section 1 (D1 up to nesting, 'literal T1' deleted, along-sequence reading stated as untested with 1,062 of 2,290); section 2 item 2 reasons corrected; 'unlabelled boundaries' deleted from Gate 6. |
| logic | L15 | major | Section 10 mixes consequences, seen quantities and guesses | ACCEPTED | Section 10 split into C1 to C5 and Q1 to Q8; the 'first evidence' sentence replaced by a pointer to register rows 732, 817, 818. |
| logic | L16 | minor | No register citations | ACCEPTED | Section 11 Register check. |
| logic | L17 | minor | Counts and wording | ACCEPTED | Section 2 items 8 and 10, section 4, section 9 (75% and 83%), header, V-based region facts (445 regions). |
| logic | L18 | minor | Free constants | ACCEPTED | Section 7 hub cap 51 (inherited), section 8 caliper sweep, section 9 UFBoot 95, seeds, section 12 item 10. |
| logic | L19 | minor | Tie rule decides many strings | ACCEPTED | Section 4 (bases are summed; sensitivities with the SET of tied labels and a 5% band). Counts from logic/neartie.out. |
| logic | L20 | minor | Single-level decision; unfrozen primary layer | ADAPTED | Section 4 level sweep at inflation 2.0 and 4.0 with the LEVEL-SPECIFIC tag (instead of a second deciding cell), size-2 stratum printed, F_O files hashed (section 13), section 6 held-out wording, layer naming. The level-C components are the strata of the governing null. |
| logic | L21 | minor | Chain tolerance broader than its rationale | ACCEPTED | Section 4: m1t (end-truncation only) printed beside m1c; wording corrected. |
| logic | L22 | minor | Gate 7 cannot fail | ACCEPTED | Section 7 Gate 7: parenthesis deleted, called a pipeline check (C3), planted negative control with a bar. |
| logic | L23 | minor | Tier B uses the best placement only | ADAPTED | Section 7 step 3: the number of qualifying placements is recorded and classes with more than one are dropped. Limit stated: Liftoff reports further copies only at identity >= 0.95, so tier C is the check. Hub cap 51. |
| logic | L24 | minor | Smaller definitional gaps | ACCEPTED | Section 4 (singleton convention, per-layer universe), section 5 (whole-family restriction), section 8 (bootstrap over pool components). |
| logic | L25 | minor | Overstatements (cnduplicon direction, 'true by construction') | ACCEPTED | Section 1 and section 2 item 4 reworded ('partial circularity', 69.7 / 15.7 / 5.2). |
| logic | L26 | minor | Same-locus duplicate exclusion defined through families | ACCEPTED | Section 4 Exclusion arms: coordinates-only definition, family-independent deletion, 10% guard from CP-5. |
| logic | L27 | note | Pigeonhole bound verified | ACCEPTED | Section 2 item 1 restated with the restricted counts and 486. |
| logic | L28 | note | Direction mapping and path defect verified | NOTED | No change needed. The suggested Gate 4 test (exception list invariant to input order) was added to Gate 4. |
| logic | L29 | note | Process note (input-only discipline) | NOTED | The same discipline was followed in this revision; every computation is listed below. |
| feasibility | F1 | blocker | Tier B cannot register same-contig paralogs | ACCEPTED | Section 7 option (a): rank-sharded runs (liftoff 1.6.3 find_overlaps checked in the source); the old design recovers at most 6 of LRPAP1's 11 loci and this is stated; Gate 7 can fail; Q5; cross-species Target code and the placement classifier named as new code. |
| feasibility | F2 | major | Gate 1 cannot currently be satisfied | ACCEPTED | Same as wording 10; two further attempts to explain 2,259 failed (c5_disc2259.out, c5b_disc2259_scan.out). |
| feasibility | F3 | major | Gate 5 premise | ACCEPTED | Same as wording 1. |
| feasibility | F4 | major | Genome-wide family tables exist | ACCEPTED | Section 2 item 5; section 6 R3a; min_cov_shorter 0 and the 09-25 date for human_testis declared. |
| feasibility | F5 | major | Deciding T2 set is empty by design | ADAPTED | Same as wording 13 (assignment rule printed with its expected size; contigs by accession and chr). |
| feasibility | F6 | major | F_O gene membership is not the 1,709 rows | ACCEPTED | Section 4 Family layers (fold map; V_F = 1,657; unclustered genes excluded with a singleton arm; half rule). Counts reproduced (c2: 76 folded V genes inherit a cluster). |
| feasibility | F7 | major | Verdict-set counts quoted for the wrong pair set | ACCEPTED | Section 5 Counts: 38 depth-matched pairs in 33 families at 0.90 (c7), source pairs.tsv named, coverage printed, 0.95 and 0.98 UNDERPOWERED, partial families out of scope. |
| feasibility | F8 | major | Controls not identity-matched; circular selection | ACCEPTED | Section 5: handled by demotion (no controls, no p, no class); identity-matched controls are infeasible because the near-distance pool has 11 and 27 pairs. |
| feasibility | F9 | major | No executable join for the matched human run | ADAPTED | Same as wording 21. |
| feasibility | F10 | major | Audit tie dependence | ADAPTED | Same as wording 3. A fixed tie rule with all tied chains printed replaces the suggested minimum and maximum over resolutions. |
| feasibility | F11 | major | Exclusion counts are over 2,334; family-dependent exclusion | ACCEPTED | Section 2 item 10 and section 4 Exclusion arms (V counts 36, 388, 105 pairs / 181 genes, 862; family-free definition). Reproduced (c10_excl_counts.out). |
| feasibility | F12 | major | E-table parameters missing | ADAPTED | Section 7: cap 51, one -copies pass, 50% outgroup overlap, unannotated copy sites listed and not scored (the recon's inheritance rule is not adopted because it cannot be verified), tier C parameters inline, the dangling 'tier C step 5' reference removed. |
| feasibility | F13 | major | Arm 1 profile of every member; synteny classifier lost; pilot rates | ADAPTED | Section 9: arm 1 uses the expansion's ancestral member only (anchor = triple); classifier written from the README with the chaining rule; unresolved rates 83% and 75%; UNDERPOWERED predicted (Q7). |
| feasibility | F14 | major | Arm 2 under-specified; identity arm input | ACCEPTED | Section 9 Arm 2 (command line, UFBoot 95, rooting) and Identity check (input = graph dump re-run, light). |
| feasibility | F15 | major | Atom builder, class rule, dedupe, counts | ACCEPTED | Section 5 Inputs and units: class builder and dedupe are new code with tests; >= 0.5 declared; counts flagged as un-deduplicated; tie rule; format argument 'gorilla'. |
| feasibility | F16 | major | T2a verdict map hole | ACCEPTED | Same as wording 14 (PARTIAL, two one-sided p-values, fewer than 20 matched sets drops the expansion). |
| feasibility | F17 | minor | Rung U2 undefined for 388 genes; 'true by construction' | ADAPTED | Section 4 U2: genes without famCN excluded with the count printed; famCN source = families_cn.json cs and co, the inputs of the frozen code, and not the famCN_sotoiv column the reviewer names; section 2 item 4 reworded. |
| feasibility | F18 | minor | Small factual slips | ACCEPTED | Section 2 items 8 and 10, section 4, '09-25/27' BASE dates, 'rows' instead of members, 'draws shared per universe'. |
| feasibility | F19 | minor | Unfrozen definitions | ACCEPTED | Section 4 (bases, G0 tie rule, ordered string); set-Jaccard and 'unlabelled boundaries' deleted; the positional reading stated as set aside (Read first 1, section 1). |
| feasibility | F20 | minor | Gate 0 inputs | ACCEPTED | Section 13 Frozen inputs (sha256, rows), interpreter pinned, BAMs and support tables named (c6_hashes.out). |
| feasibility | F21 | minor | Cost estimates incomplete | ACCEPTED | Section 7 step 2 and Tier C (one pass, RSS, serial index builds, disk); section 3 (R1 split by cell). |
| feasibility | F22 | minor | Control provenance; 'unexposed' contigs | ACCEPTED | Section 7 Gate 7 (PCDHB, TUBA, RBMY dropped) and Contig roles ('no prereg names them', ledger mentions). |
| feasibility | F23 | minor | Register pointers, KEY names | ACCEPTED | Section 11 Register check (lines and rows), section 13 KEY names, rows from 1263. |
| feasibility | F24 | minor | Overstatements and gate ordering | ACCEPTED | Section 2 item 2 (fusion sentence qualified), section 3 scope of the gates, 'Amendment-free' deleted. |
| feasibility | F25 | note | Verified facts | NOTED | Re-used. A subset was re-verified in this revision (hashes c6, V counts c0 and c2, depth-matched counts c7). |
| feasibility | F26 | note | Implementer's list | ACCEPTED | Section 13 New code to write; constants in section 12 item 10. |

## Closure review of draft 2: one row per finding (draft 3)

Counts: ACCEPTED 56, ADAPTED 7 (63 rows; severities: blocker 2, major 25, minor 32, note 4). Rows with the same issue point to one fix. The `confirmed-ok` entries of the four reviewers (39 lines) are not listed.

| id | reviewer | severity | location | disposition | where fixed, or why not |
|---|---|---|---|---|---|
| CL1 | logic | major | Gate 4 clause 'family partition equal to the strata gives UNINFORMATIVE' against algebra items 2 and 3 | ACCEPTED | Section 3 Gate 4: zero variation-bearing strata gives UNDERPOWERED; a synthetic set of at least 8 strata with two clusters and one string each gives UNINFORMATIVE. |
| CL2 | logic | major | Homology-preserving null is reached by construction when strings are nested with the cut (draft-1 L2) | ADAPTED | Class renamed CUT-ALIGNED; section 2 item 3, Read first ('What the classes mean', item 2), Q2 and section 11 state that T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity; Gate 4 gains the nested-strings case (p <= 0.01). The similarity-matched comparator (F_O at inflation 2.0 and 4.0) is not registered, because those cuts differ in granularity and change m1c and m2 mechanically; it is decision 16 and can be added by an Amendment before any run. |
| CL3 | logic | major | UNDERPOWERED counts strata while p is pooled; one stratum can decide p | ADAPTED | Section 4 'Robustness to one group': per-stratum z_s table, the most negative stratum removed, p recomputed on the same draws, the lower class issued. A sign test over strata was not adopted because it adds a second statistic with less power. |
| CL4 | logic | major | Gate 6 tie counts over the wrong universe | ACCEPTED | Section 3 Gate 6: 857 of 12,597 labelled merged exons of the genes of V (898 of 13,501 over all 2,334 genes); genes of V with at least one tie exon: 316; 'Only the counts listed in this Gate are tripwires'. |
| CL5 | logic | major | Sampling frame of the control sets not stated (draft-1 L11 d) | ACCEPTED | Section 8: control pool excludes hubs and components that contribute an expansion; a draw takes a qualifying component uniformly and then one of its matched sets uniformly; the floor counts distinct control components; the component-size class (2, 3 to 5, 6 to 11, 12 to 51) is added to the matching. |
| CL6 | logic | major | Matching variable open (median over which rows, overlap with the expansion, boundary of the caliper) | ACCEPTED | Section 8: integer identities, M = twice the median over the copy rows that join two records of the set, caliper /dM/ <= 10 inclusive (sweep 4, 10, 20), a set without a joining row cannot be a control, components that hold an E-table class are not in the control pool. |
| CL7 | logic | major | n_testable versus the drop at 20 matched sets; expected size 12 and 21 (draft-1 L12, F5) | ACCEPTED | Section 8 defines n_testable with the floor of 20 control components, evaluated before any class; section 2 item 7, section 7 and Q5 give upper estimates with the Wilson interval (8 to 17 and 14 to 30) and the Poisson chance; siblings count as separate units with a group-once sensitivity (section 7 step 4). |
| CL8 | logic | major | Arm 1 forces its own outcome; 'mixed-family' contradicts the outcome | ACCEPTED | Section 9: phi = the family holding most nodes of expressed derived members; mixed-family = phi holds a node that is not a node of a member of E; outcome, causes and funnel rewritten; the identity check uses E_nodes; the unit is named in the section 3 tables. |
| CL9 | logic | major | Arm 2 does not measure youth (D is a clade in 1//D/ of random-parent histories) | ADAPTED | Section 9 arm 2 tests E_tree (members of E in the family tree) as a clade that leaves out at least one other family member; the comparator is size-matched within the identity caliper; the text says a clade is a pipeline and gene-conversion check and that a non-clade is the informative outcome; classes are kept because T2b is not authorised and is predicted UNDERPOWERED. |
| CL10 | logic | minor | Five open cases in the algebra (T2a headline, level sweep, UNDERPOWERED sweep cell, INVALID R3a, arm 2 per library) | ACCEPTED | Section 3 verdict algebra items 4 and 5; section 4 level sweep; section 9. |
| CL11 | logic | minor | m1t, Mult, Delta not well defined; Gate 4 relabelling clause (draft-1 L21, W30) | ACCEPTED | Section 4 statistics (m1t pairwise comparable, Mult over mu >= 1, Delta over families with at least two genes); Gate 4 plurality condition. |
| CL12 | logic | minor | Gate 7 bar of 8 of 11 is a bar of 8 of 8 (draft-1 L22) | ACCEPTED | See CF1: the bar is the 6 qualifying records, a failure blocks T2a and T2b. |
| CL13 | logic | minor | n_min, floor of 20 and flags have no entry in section 12 item 10; 'the expected answer' in Amendment 1 (draft-1 L18) | ACCEPTED | Section 12 item 10 lists n_min (source: the E < 8 rule of the seed-pool prereg) and every fixed floor with its sweep; section 6 no longer lets Amendment 1 fix an expected answer (Q8 stands). |
| CL14 | logic | minor | Section 10 mixes consequences, expectations and predictions (draft-1 L15, L27) | ACCEPTED | Section 10 rewritten (C1 to C4, Q1 to Q9); C1 limited to the counted layers and rungs; old C4 became Q4 with a falsifier; old C5 became Q9 with a number; old Q4 became C4; Q3 and Q8 restated; 'number of F_S families' replaces 'cell count'. |
| CL15 | logic | note | T1-a is a count whose information is the exception list | ACCEPTED | Section 4 decisive paragraph says so and prints m1c and m2 separately. |
| CL16 | logic | note | B and the call split for the Gate 4 repeats | ACCEPTED | Gate 4: B = 10,000 per repeat in several light calls of under 180 s. |
| CF1 | facts | blocker | Gate 7 positive control unreachable (pseudogene records need -f; two long gene models cannot reach coverage 0.5) | ACCEPTED | Section 7 step 2 command now carries -f types.txt; Gate 7 bar = every annotated record of the gene-LRPAP1 component with at least 50% of its exon bases inside an LRPAP1 copy interval (6 records), the two long models are expected misses with the reason; failure blocks T2a and T2b; the fragment-in-shared-run exposure is listed among the ways the control can fail; the Amendment may change the command only on evidence from the two controls; C3 and Q5 restated. |
| CF2 | facts | major | Registered command has no -f | ACCEPTED | Section 7 step 2: full command with types.txt (681 of 2,249 pool records, 502 of 1,446 called records are pseudogenes) and the minimap2 named by the Amendment. |
| CF3 | facts | major | What counts as a placement is not stated | ACCEPTED | Section 7 step 3: primary row with coverage >= 0.5, sequence_ID >= 0.5 and neither partial flag; other rows are no placement; records whose only row is partial are counted. |
| CF4 | facts | major | The pool is blind to cis-only (tandem) expansions (second half of draft-1 F1) | ACCEPTED | Section 11: dispersed and cross-block copies only; 5 of 307 components have all records on one contig; tandem-only expansions are outside the E-table and T2 makes no statement about them. |
| CF5 | facts | major | Node join between fam.clusters.tsv and fam.copies.tsv (draft-1 W15) | ACCEPTED | Section 8: a node is a row of fam.copies.tsv; fam.clusters.tsv is not read; node counts 864 and 1,559. |
| CF6 | facts | major | Gate 6 counts with no stated source (same-locus duplicates; tie universe) (draft-1 F11, W9) | ACCEPTED | Section 4 exclusion arms name the gene span (hull of the gene's CAT v4 transcripts, families_cn.json c, s, e, st) and say at most 105 genes are deleted; tie counts as CL4. |
| CF7 | facts | minor | 'Registered OUTPUTS reproduce' overstates; 'none gives 2,259' overstates (draft-1 F2, W10) | ACCEPTED | Section 2 item 9 rewritten (cluster and unit counts reproduce; pooled delta, p and arm A shares are checked by Gate 1; scan result as measured). |
| CF8 | facts | minor | python -I ignores PYTHONHASHSEED; python3 versus python | ACCEPTED | Gates 0 and 3: python3 -B, env PYTHONHASHSEED=0 and 1 without -I. |
| CF9 | facts | minor | The graph dump is light in sections 6 and 9 and heavy in section 8; BASE binary differs | ADAPTED | One class stated in section 8 (light for human_testis; gorilla PAFs measured first under light, heavy above 2 GB); Amendment 1 requires --min-cov-shorter 0 and the BASE flags and a check that every BASE cluster lies in one component. The RSS is not measured now because the dump is not authorised. |
| CF10 | facts | minor | Five labels that do not match their source (486, 885/10,651, 157, BAM count, three contigs) | ACCEPTED | Section 2 items 1 and 7; section 5 counts (878 and 10,506 under the 100 bp rule); section 7 contig list (chr6, chr9, chrX added with their prereg names); section 13 BAM wording. |
| CF11 | facts | minor | A displaced record is not 'unplaced' in Liftoff 1.6.3 | ACCEPTED | Section 7 step 2: displaced records are recognised by coverage < 0.5 or by landing on no single in_place outgroup record; the count of records with no triple is printed. |
| CF12 | facts | note | Cost line should state 1,446 records and K = 44 | ACCEPTED | Section 7 step 2 and section 13 (split rule listed as new code). |
| CW1 | wording | major | Gate 6 tie counts (898 of 13,501; 316) | ACCEPTED | See CL4. |
| CW2 | wording | major | 'Records' undefined in the hub rule (4 and 5 against 2 and 3) | ACCEPTED | Section 7 step 1: hub = more than 51 members, members being the annotated records of any length and the copy@ sites; 4 and 5 stand. |
| CW3 | wording | major | Opening sentence contradicts step 4; control pool may share a component with an expansion (draft-1 W16) | ACCEPTED | Sentence deleted; section 7 step 5 and section 8: a component that contributes an expansion is not in the control pool. |
| CW4 | wording | major | Join from cluster rows to copies rows (draft-1 W15) | ACCEPTED | See CF5. |
| CW5 | wording | major | n_testable before or after the matched-set drop | ACCEPTED | See CL7. |
| CW6 | wording | major | E phase does not write what section 8 needs (draft-1 W16) | ACCEPTED | Section 7 step 5: member_exons.tsv, control_pool.tsv, extension.tsv, 0-based half-open ids, hashes recorded with the E-table. |
| CW7 | wording | major | Caliper boundary on a 0.001 grid (draft-1 W14) | ACCEPTED | See CL6. |
| CW8 | wording | major | Amendment 1 can be written after Phase R1 has printed (R3a hold-back) | ACCEPTED | Section 6: the gene-to-locus rule and half rule are written in the file; Amendment 1 is committed before the first Phase R1 run and may not change the rule, cells, statistics, nulls or Q8; 'expected answer' deleted. |
| CW9 | wording | major | Read first does not define the classes (draft-1 W20) | ACCEPTED | Read first: 'What the classes mean' paragraph; ORGANISES-ENRICHED renamed CUT-ALIGNED. |
| CW10 | wording | major | Read first omits links tested, V, D1/D2 mapping, nesting, verdicts that can occur (draft-1 W20) | ACCEPTED | Read first items 8 to 12. |
| CW11 | wording | minor | Algebra cases (T1-a per direction, T2a cell, lower tail only, T1-d flag, single-run labels) (draft-1 W5) | ACCEPTED | Section 3 items 4 and 5, decisive-cell table, section 4. |
| CW12 | wording | minor | Gate 2 draw rule, Gate scope, tripwires (draft-1 W23) | ACCEPTED | Gate 2 rewritten (min(20, n), one stream, at most 1 differing); scope sentence; only Gate 6 counts are tripwires. |
| CW13 | wording | minor | F_S/F_Q universes and counts (draft-1 L27, W4) | ACCEPTED | Section 2 item 1, section 4 F_S bullet, Gate 5b wording, C1. The suggested F_Q figure of 475 cells over V was not added because it is not verified. |
| CW14 | wording | minor | Section 5 sources and definitions (draft-1 W12) | ACCEPTED | Section 5: genes.tsv join, pairs.tsv join, 1-based columns, atoms versus classes, 878 and 10,506, B_viol. |
| CW15 | wording | minor | Null details: what is shuffled, string order, RNG objects (draft-1 W7) | ACCEPTED | Section 4 nulls and draws; section 8 test. |
| CW16 | wording | minor | Descriptive statistics that two implementers compute differently (draft-1 W3, W30, W31) | ACCEPTED | Section 4: m1t, audit tie key, Mult, Delta, size-2, restriction, near-tie band, 105 pairs and 181 genes, U2co versus U2cs. |
| CW17 | wording | minor | Authorisation wording (draft-1 W19) | ACCEPTED | Section 3 run order: a 'go' names its phases; light under rlock light; Amendment rows; matched human run row. |
| CW18 | wording | minor | Gate 0 hashes of 36 GB; call budgets (draft-1 W19, W29) | ACCEPTED | Gate 0: size and mtime above 1 GB, sha256 once in its own call; light calls under 180 s; heavy limits apply to heavy calls; the BASE dump is hashed when made. |
| CW19 | wording | minor | T2b-1 unit and bounds (draft-1 W17) | ACCEPTED | Section 3 table and section 9: mixed-family expansion; one-sided exact bounds; at least one ancestral member; seed named. |
| CW20 | wording | minor | Decisions list omits defaults and constants (draft-1 W19) | ACCEPTED | Section 12: items 8 and 11 corrected, 10 extended, 14 to 16 added. |
| CW21 | wording | minor | Terms and symbol collisions (draft-1 W24, W28) | ADAPTED | Terms define clean, family-less, multi-family, MCL, independent unit; R_J, D (derived set) and K renamed or defined (Dv, J); the remaining overloads (cell, E) are disambiguated by context and left. |
| CW22 | wording | minor | C4/C5 are not consequences; 'held-out' tense | ACCEPTED | Section 10; 'held-out' only after R3a runs; R3b chromosomes are 'not used to select a rule'. |
| CW23 | wording | minor | BAM count, source of 'expressed', index paths, names (draft-1 W27) | ACCEPTED | Sections 7, 8 and 13. |
| CW24 | wording | minor | R3a power count not registered | ADAPTED | Gate 6 prints the number of strata with at least two clusters among V genes of the BASE graph; Q8 applies only if it is at least 8. Not computed now because it needs the dump. |
| CW25 | wording | minor | Style sentences that hurt comprehension | ACCEPTED | Read first items 1 and 6, labels, section 5 first sentence, section 7 contig roles. |
| CW26 | wording | note | Gate 7 wording 'blind' | ACCEPTED | Removed. |
| CD1 | disposition | blocker | Row F1: Gate 7 bar unreachable, Q5 predicts a pass | ACCEPTED | See CF1. |
| CD2 | disposition | major | Row W17 and others: mixed-family contradicts arm 1 | ACCEPTED | See CL8. |
| CD3 | disposition | minor | Row W15: no join key | ACCEPTED | See CF5. |
| CD4 | disposition | minor | Registered counts over different universes (Gate 6, F_S, control pool) | ACCEPTED | See CL4, CF6, CW13, CW14. |
| CD5 | disposition | minor | SPLIT does not apply to the two libraries; no T2a headline | ACCEPTED | See CL10. The headline is the lower class of the two libraries. |
| CD6 | disposition | minor | Constants that set classes are not in section 12 item 10 | ACCEPTED | See CL13. The origin of n_min is the seed-pool prereg; the floors of 20 pairs and 8 families and the floor of 20 control components have no earlier source and are listed as fixed now with sweeps. |
| CD7 | disposition | minor | No funnel for arm 2 (draft-1 L9) | ACCEPTED | Section 9 arm 2 prints the funnel with the share. |
| CD8 | disposition | minor | Dump command classed light and heavy | ADAPTED | See CF9. |
| CD9 | disposition | minor | Row c12: orangutan outlier | ACCEPTED | Input-only table row corrected; section 7 step 2 says the estimate is from three same-genome fits and the Amendment sets the budget per outgroup. |

## Input-only computations run in this revision

Scratch directory: `/mnt/linuxdisk/tmp/prereg_revision/`. All ran under `tools/rlock.sh light` with `python3 -B -I`. None computes a deficiency, purity, containment, enrichment, S-rate or null draw. None joins a family label to a unit string beyond the marginal counts of section 2 item 1 (reproduced from the logic reviewer's scripts). The scripts `c0` and `c2` read the reviewers' strings file only for the booleans 'has an exonic duplicon' and the clean, multi-family or family-less class; they read family membership and graph structure and no unit string.

| script | output | one-line result |
|---|---|---|
| `c1_fo_graph.sh` | `c1_fo_graph.out`, `fo_graph/run.log`, `fo_graph/ours.graph.tsv` | the frozen binary (sha1 e4cc13b9) with `--dump-graph` reproduces `ours.clusters.tsv`, `ours.loci.tsv` and `ours.pairs.tsv` byte for byte; graph of 1,819 nodes and 11,031 edges (0.9 s) |
| `c0_V_list.py` | `c0_V_list.out`, `V_genes.tsv` | V = 2,142 genes (2,036 clean, 106 family-less) |
| `c2_fo_components.py` | `c2_fo_components.out`, `c2_fo_components.json` | 339 components; V_F = 1,657 genes (76 folded annotations inherit a cluster); 325 strata hold V_F genes; 30 variation-bearing strata with 416 V_F genes (29 hold at least two multi-gene clusters; 4 to 42 genes, median 11); 295 single-cluster strata with 1,241 genes; 467 V genes never connected and 18 singleton nodes |
| `c5_disc2259.py`, `c5b_disc2259_scan.py` | `c5_disc2259.out`, `c5b_disc2259_scan.out` | with either BED, 2,290 genes overlap an exon by at least 1 bp; no rule on summed, best-segment or best-exon bp from 1 to 300 gives 2,259; an exonic coverage fraction >= 0.10 gives 2,258; no per-chromosome subset explains 31 genes |
| `c6_hashes.sh` | `c6_hashes.out` | sha256 (16 hex), size and lines of the Gate 0 list; large files by size and mtime |
| `c7_depth_matched_verdict.py` | `c7_depth_matched_verdict.out` | gorilla verdict set, depth-matched eligible pairs and families: 38 / 33 (tau 0.90), 15 / 15 (0.95), 5 / 5 (0.98), identical under the 1 bp and 100 bp touch rules |
| `c8_heldout_power.py` | `c8_heldout_power.out` | supported pool components 80 (KB3781), 139 (OR6737), 157 either; of size <= 51: 76 and 134; no member on the six LRPAP1 contigs: 23 and 38; all members on the ten contigs no prereg names: 3 and 8 |
| `c9_r3a_power.py` | `c9_r3a_power.out` | 334 V genes carry a human_testis BASE cluster; 86 BASE clusters hold at least two V genes (316 genes) |
| `c10_excl_counts.py` | `c10_excl_counts.out` | 233 pairs / 339 genes over 2,334 genes become 105 pairs / 181 genes in V; the readthrough list has 76 rows, 36 in V |
| `c11_multi_duplicon_genes.py` | `c11_multi_duplicon_genes.out` | 1,062 of 2,290 duplicon-bearing genes carry two or more dominant duplicons (974 of 2,142 in V); strings only |
| `c12_shardfit.py` | `c12_shardfit.out` | Liftoff self-lift timings: gorilla, human and chimp fits give 0.05 to 0.09 s per record; the orangutan fit is dominated by one 544 s shard (a 10-fold outlier); peak RSS 10.2 to 13.8 GB. The estimate for cross-species lifts is therefore uncertain and the Amendment sets the time budget per outgroup |
| `logic/fq_checks.py`, `string_only_checks.py`, `short_counts.py`, `halves.py`, `region_strata.py`, `neartie.py`, `control_pool.py` (copies of the logic reviewer's scripts) | `logic/*.out` | D[S->Q] = 7 genes (2 if ungrouped genes are ignored); pigeonhole marginals of F_S, F_Q, F_O; 32 U1 strings / 110 genes with more than one principal duplicon; 862 genes with at most two exons; 9 straddling clusters (61 genes); 445 SD98 regions in V; 898 exact exon ties; gorilla control pool 11 / 27 / 885 / 10,651 against 13 / 5 / 15 / 21 real pairs |

**Closure round (draft 3).** The author computed nothing new in this round. The closure reviewers' scripts and outputs (synthetic scorer tests, stratum structure, tie universes, pool and direct-row counts, LRPAP1 record geometry, hash recomputations) are in `/mnt/linuxdisk/tmp/prereg_closure/{logic,facts,wording,disposition}/`; none computes a deficiency, purity, containment, enrichment, S-rate or null draw on real labels.

## Design decisions that need the user's confirmation (default chosen in draft 3; consequence of the alternative)

1. **Governing null.** Homology-preserving: labels are permuted inside the connected components of the F_O graph. Alternative: chromosome or SD98-region strata, which give p <= 0.01 by construction. Consequence of the default: T1-b and T1-c speak only for the 30 split groups (25% of the clustered genes).
2. **What T1-b and T1-c can claim.** The class is named CUT-ALIGNED and the text says it cannot separate duplicon organisation from shared sequence similarity, because strings nested with the cut give the lowest p. Alternative: register a similarity-matched comparator (the F_O cut at inflation 2.0 and 4.0 as competing sequence-only re-cuts of the same group) before any run, to test the stronger claim. It is not in draft 3 because those cuts differ in granularity.
3. **Robustness clause.** The class is the lower of the class from all strata and the class after removing the stratum with the most negative z. Alternative: a stratum-level sign test (22 of 30 for p <= 0.01), which has less power.
4. **Audit.** Descriptive; nothing is excused in a decisive statistic. Alternative: excuse 'short' and 'ends-only' inside the statistic, which makes NOT-REFUTED reachable by construction.
5. **Gorilla T1 arm.** Descriptive, no class, matched human run optional. Alternative: keep an enrichment verdict, which identity selection makes unfailable.
6. **R3a.** The human_testis BASE snapshot (pre-f1v2, RNA-derived) is a light held-out arm, powered only if the graph dump gives at least 8 strata with two clusters (Q8 is conditional). Alternative: wait for the heavy current-default run (about 17 h for A119b).
7. **T2 primary recipe.** Tier B with rank-sharded Liftoff runs and `-f types.txt` (about 45 to 60 min heavy; 1,446 records, 132 runs). Alternative: tier A or C primary (1 to 2 h each).
8. **Gate 7 bar.** The 6 LRPAP1 records with at least 50% of their exon bases in a copy interval; the two long gene models are expected misses. Alternative: use the eight 20-kb copy bodies as synthetic records, which changes the lift recipe.
9. **T2 hold-back.** By expansion: LRPAP1 only, with the ten contigs printed as a least-exposed subset. Alternative: the ten contigs decide, which is UNDERPOWERED by design.
10. **T2a controls.** Two-step draw (control component, then matched set), component-size classes added to the matching, floor of 20 control components, T2a headline = the lower class of the two libraries. Alternative: a uniform draw over all matched sets, which samples almost only large components (83% to 86% of expressed pairs lie in hubs).
11. **T2a statistic.** Mean R_j with absent members counted as non-co-members. Alternative: strict containment c_j.
12. **T2b.** Arm 1 is an estimate; arm 2 tests E_tree as a clade that leaves out a family member and has classes; no combined verdict; both predicted UNDERPOWERED. Alternative: drop T2b from this registration.
13. **T1 universe.** Genes that carry an F_O cluster (1,657). Alternative: singleton families for the 485 unclustered genes.
14. **Amendment 1 timing.** Committed before the first Phase R1 run. Alternative: after R1, which weakens the R3a hold-back.
15. **Gate 1.** It does not require the 2,259 / 2,290 explanation.

## Findings that could not be fully resolved

- The 2,259 versus 2,290 discrepancy of the 09-30 disclosure stays unexplained after five attempts; it is an erratum, not a gate.
- The 15% single-copy share behind the expected testable expansions comes from a name-linked subset (24 of 160; Wilson interval 10.3% to 21.3%) and may not transfer; the figures 12 and 21 are upper estimates and KB3781 may fall below n_min = 8.
- Rank-sharding removes conflicts inside a pool component, not between unrelated records of one Liftoff run; displaced records are detected only by coverage < 0.5 or by landing on no single in_place outgroup record. The fragments of the LRPAP1 control share runs 1 to 3 with full-length copies, which can make the control fail.
- Liftoff reports further outgroup copies only at identity >= 0.95, so the multiplicity guard of the E-table is incomplete; tier C is the check.
- The pool is blind to tandem-only expansions (5 of 307 components lie on one contig).
- T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity; a similarity-matched comparator is an option for an Amendment.
- The synteny classifier of the old pilot is lost and is rewritten from the README; T2b is predicted UNDERPOWERED.
- The R3a stratum count needs a graph dump that is not run now (Q8 is conditional), and the peak RSS of the gorilla BASE-PAF dumps is not measured.
- The time estimate for cross-species Liftoff comes from three same-genome self-lift fits; the Amendment sets the budget per outgroup.

## Can the registered cells produce each verdict class? (draft 3)

| cell | class | reachable from the inputs |
|---|---|---|
| T1-a | REFUTED | yes, and expected |
| T1-a | NOT-REFUTED | in principle (m1c = 0 and m2 = 0 on 1,657 genes), practically unreachable |
| T1-b, T1-c (Phase R1) | CUT-ALIGNED, PARTIAL, NOT-DISTINGUISHED | yes, 30 variation-bearing strata; CUT-ALIGNED is expected whenever strings and cut both follow sequence similarity |
| T1-b, T1-c (Phase R1) | UNDERPOWERED, UNINFORMATIVE | UNDERPOWERED no (30 against n_min 8); UNINFORMATIVE practically unreachable |
| T1-b, T1-c | SPLIT | only after R3a has run |
| T1-b, T1-c (R3a) | all classes | possible; UNDERPOWERED if fewer than 8 strata, known only after the BASE graph dump (86 candidate clusters); Q8 is conditional on the count |
| T1-d (gorilla) | none issued | tau 0.95 (15 pairs) and 0.98 (5 pairs) are flagged UNDERPOWERED |
| T2a, OR6737 | SPECIFIC, PARTIAL, NOT-DISTINGUISHED, LESS | yes if n_testable >= 8 (at most about 21 expected, 14 to 30 before the section 8 filters) |
| T2a, KB3781 | all content classes; UNDERPOWERED | borderline (at most about 12 expected, 8 to 17; Poisson chance of fewer than 8 about 9% at 15%) |
| T2a | SPLIT | between the two libraries, only if both are powered |
| T2a, ten-contig subset | none issued | UNDERPOWERED (3 and 8 supported components) |
| T2b arm 1 | none issued (estimate) | UNDERPOWERED expected (pilot 3 of 18) |
| T2b arm 2 | YOUNGER, PARTIAL, NOT-DISTINGUISHED | only with at least 8 informative expansions per library, which is not expected; UNDERPOWERED expected |
| any cell | INVALID | by gate or control failure; Gate 7 can fail |
