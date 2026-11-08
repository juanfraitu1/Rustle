# The advisor's standing questions — and what the codebase can answer

Compact, current version (2026-11-01). Verbatim previous text: `archive/2026-11/ADVISOR_QUESTIONS.md`; every entry ends with a pointer to it. `r####` = row of `NEGATIVE_RESULTS_REGISTER.md`. Audience: Stefan Canzar, who assumes a good number is luck or overfitting until shown otherwise; open with the discipline, then evidence, then concessions.

## Current state (2026-11-03)

- **Defaults unchanged** for families; assembly gained `--polish-tss rescue` and the `final.gtf` sub-chain drop (r1485, r1498). O2 assign judges exact AS ties by default; the ms gate and the primary-on-family guard are opt-in and measured.
- **Families stage, all FAIL or opt-in:** r1469-1476 (U1, dispersion, SKETCHFAM, semi-guided PASS not default, representative rules, chain gate, absorb), r1490/1491 (end layers break families), r1492-1497 (readthrough, truth split, fused loci PASS opt-in, shared-exon split FAIL).
- **His multimapper interest, answered 2026-11-03:** secondaries carry information down to ms ratio 0.5 of the primary (r1511; `-p 0.8` would lose 23% of the useful ones). Full family membership cannot be read off a gene graph of multimappers (r1512: pair precision <= 0.032 vs O1 0.341, silent members not found). What works is per locus: the informative-secondary pool doubles the reference-absent contigs `candidate_copies` finds (r1513, 6 vs 3). IsoCon's reference-free piece worth emulating is its significance test, already ported; unaligned reads are not where missing copies hide (they align as multimappers to the nearest paralog; o2 sample-scale memory).
- **O2 proof (his standing doubt):** on reads minimap2 could not place, O2 and the aligner share almost all verdicts (held-out 151 decidable: 124 both right, 5 O2-only, 1 primary-only, 21 both wrong; r1514). The measurable gain is correctness: without the guard 77% of held-out assignments were reads from uncatalogued parent genes assigned to retrocopies with p 1e-11 (r1510).
- **Assembly / tools:** FAMASM FAIL (r1479); IsoCon WORSE (r1478); SQANTI single-copy: StringTie, FLAIR, IsoCon worse (r1482); SQANTI3 filter loses 6.4% of exact matches (r1499).
- **Second individual (r1578/1579, 2026-11-08):** OR6737 on Jim's T2T gains about 12 expressed loci over Jim's noise floor (NBPF, IGHV-like lead), not many; details under Reference bias.
- **Strainy (his interest, 2026-11-08):** Strainy's phasing lives in `copy_phasing` (Strainy piles + label propagation/correlation clustering + refinement, IsoCon's test instead of the 0.2 allele-frequency floor, allele-slot certificate); OTHER-BINARY, not in the default pipeline. Applied 11-08: (a) OR6737 vs Jim strong loci certified 5/7 vs 0/7 (r1582, same-read caveat); (b) collapsed copies on GRCh38 (r1583 PASS as registered 10/13 but only 5/13 carry the missing gene's own reads; r1584 combined certificate + S/H held out on HG002 FAIL 2/6, null 0/7, depth-limited). Collapsed tandem copies on a non-T2T reference are visible from RNA at about a third; absent copies on a draft much better (divergence pile, r1539).
- **PSV VCF:** O2 is a pileup + caller; bcftools agrees 54.6% as run, 93.9% with MAPQ reset (r1484).
- **SPENT samples** for family rules: KB3781, OR6737, A119b, HG002.

## Status as of 2026-10-07

Status table prepended 10-07 (one row per question / section, then post-09-20 questions); answers below are folded into each entry and updated through r1484. Rows that were only about doc provenance: register then ended at r1262, `OBJECTIVES_AND_VERIFICATION.md` and `docs/REPRODUCE.md` removed 09-23 (commit 4059e4be; `REPRODUCE.md` now at repo root), cargo 1,204 tests passed (10-07). (full: archive/2026-11/ADVISOR_QUESTIONS.md:3)

Advisor questions raised after 09-20 (current answer, source):

| question (short) | current answer | where |
|---|---|---|
| 1.2 addendum: MCL size-2 families | A two-member family is clustered or split by MCL floating-point residue (w 0.99 split, 0.98 join); 57% of genome-wide families are size 2; fix candidates not applied. | memory/project_mcl_two_node_fragility.md; docs/FAMSIM_DESIGN_2026-10-02.md |
| NEW 09-20 Simple Jaccard sweep instead of MCL + conjunct | Components-of-Jaccard F 0.6044 vs shipped 0.7016 (86%); MCL and exon conjunct are a pair (r910/r911); DP chaining redundant (r912). | docs/JACCARD_SWEEP_RESULT_2026-09-20.md |
| NEW 09-21 Eichler AS-margin comparator | Agreement 97.1% where both call; his rule discards 83.7% of multimappers; we assign 145 it discards. Since 10-31 a STANDING test (r1477): decides 80-84% of family-scope reads, pick = primary on 87-93%, picks unspliced over spliced on 7-12% (wrong by r1270), undecided 16-20%. | docs/EICHLER_COMPARISON_2026-09-21.md; r1477 |
| NEW 09-22/23, 10-01 Primaries-only loci are a coin toss; use ALL alignments | GOOD seeding (AS >= 0.98 of best) default since 09-24 (r1060); at NPIP 58.7% of primaries tie, ALL adds 87 off-copy nodes; `--seed-pool` exact representatives 12/8/2 (r1261). minimap2 ranks secondaries on ms not AS (r1477). | docs/NPIP_READ_POOL_2026-10-01.md; r1477 |
| NEW 09-23 Reference bias (absent copies hide as unmapped / multimapping reads) | REFUTED arm A: extra-copy reads 100% absorbed at the template as MAPQ-60 primaries up to 5% divergence, 0 unmapped (r1089); arms B/C not run. See also `## Reference bias...` below. | docs/PREREG_o3_reference_bias_2026-09-23.md |
| NEW 09-25 Read-Jaccard + multimapper anchors as edge rule | NOT a replacement: chr16 precision .964 to .725; reads only as candidate edges gated by exon share. | docs/ADVISOR_JACCARD_ANCHORS_2026-09-25.md |
| NEW 09-25 minimap2 `-G 50k` against read-throughs | REFUTED: breaks 99.6% of real >= 50 kb introns; polyA-inside-intron rule passes held-out chr20 and gorilla. | docs/READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md |
| NEW 09-26 One locus representative loses information vs Liftoff | Four arms dev-only, PARKED (r1121-1127); at NPIP the representative IS the loss; R_J opt-in. Later rules: reads-x-junctions FAIL (r1473), chain-gate FAIL (r1475); most-reads stays default. | docs/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md; r1473, r1475 |
| NEW 09-28 You haven't replicated Soto | ARI 0.7096 literal, 0.9698 following their released code with S1C CN, 0.9277 with our copy numbers. | docs/SOTO_REPLICATION_STATUS_2026-09-28.md |
| NEW 09-30 Shared exons should exclude non-genes | Soto's code: exons of every feature, pair kept if either gene is coding or unprocessed pseudogene; 129/491 Soto families have >= 2 coding genes. | docs/SOTO_M2_MEETING_EVIDENCE_2026-09-30.md |
| NEW 09-30, 10-01 RNA famCN; some transcripts are alleles | RNA allele / haplotype-count chain pre-registered 10-01 on KB3781; became the opt-in `candidate_copies` stage (Part 3.10). | docs/PREREG_rna_allele_haplotype_count_2026-10-01.md |
| NEW 10-02 Prove it for ANY family, ANY ape | `bench/famsim` simulator, 22-rung ladder on human SNRPB and orangutan SPAG7; found the MCL size-2 fragility. | docs/FAMSIM_DESIGN_2026-10-02.md |
| NEW 10-03 Where is my hierarchy in the MCL? | The MCL cut is one level of the nested edge-test lattice; dendrogram = presentation + scoring frame; level-C reading (r1260). | docs/archive/2026-11/LEDGER_GENERAL.md:230 (ADVISOR_HIERARCHY_2026-10-03.md) |
| NEW 10-04 Analyse LRPAP1 | 11 gorilla loci (8 full + 3 5' fragments; paper: 10); O1 de novo 6+3 fibroblast, 8+3 testis; O2 0 contested; O3 maternal 7, paternal 8. | docs/LRPAP1_FAMILY_2026-10-04.md |
| NEW 10-06 Semi-guided mode will not be useful | Tab removed from the Locus Anatomy artifact (version 18). Later registered anyway: semi-guided families PASS both held-outs KB3781 + OR6737 (weaker pass: strand bug fixed between passes), not default (r1472). | memory/project_two_modes_scope.md; r1472 |
| NEW 10-07 Multi-copy families vs families within SDs / expansions | Four levels of one object (family / SD edge evidence / duplicon node class / expansion clade); T1/T2 pre-registered, NOT run. | docs/archive/2026-11/LEDGER_GENERAL.md:432 (TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md) |
| NEW 10-07 Test against the Yoo 2025 ape expansions | Gorilla only, DESCRIPTIVE: Y1 LRPAP1 10/10; Y2 JMJD7-PLA2G4B 8/9 (KB3781) / 6/9 (OR6737); MAPKBP1/SPTBN5 read-limited. | memory/project_yoo_cross_ape_expansion_test.md |
| NEW 10-29..10-31 Compare to IsoCon / single-copy tools | IsoCon WORSE both classes (multi-copy 20 of 30 copies vs ours 27; r1478); SQANTI single-copy: StringTie, FLAIR, IsoCon WORSE (r1482). | r1478, r1482 |
| NEW 11-01 Why no pileup + variant calling for PSVs? | O2 already is a pileup + caller: per-copy piles, a Poisson error test with Bonferroni, >= 6 carriers per deciding site. `bench/o3_graph/psv_vcf.py` writes the VCF (ID `family:copy:POS`, per assigned read its deciding sites). bcftools 1.25 on the same records agrees on 54.6% as run (caps base quality at MAPQ; 24% of pile records are MAPQ-0 secondaries), 93.9% with MAPQ reset to 60 (post hoc, disclosed); rest = minor alleles < 10% a diploid prior calls hom-ref. Report, no bar (r1484). | `~/o2ties_real/psv/or.psv60.*`; LEDGER_O2 `### PSV VCF RESULT`; r1484 |

| NEW 11-03 Up to what ratio are secondaries useful; is `-p 0.1` justified? | Useful (known paralog or catalog copy of the primary's family) 4.8% below ms ratio 0.5, 6.9% at 0.5-0.8, 16-33% at 0.8-0.99, 71% at >= 0.99: bound = 0.5; `-p 0.8` would drop 23% of useful secondaries. | r1511 |
| NEW 11-03 Infer full family membership (silent copies) from multimappers | Gene graph FAIL at every ratio (precision <= 0.032 vs O1 0.341; silent coding joins with a truth pair 0/27 at >= 0.99). Per-locus pool works for reference-absent copies (r1513). | r1512, r1513 |
| NEW 11-03 Emulate IsoCon's reference-free inference | Its significance test is already ported (`candidate_copies`); a silent pseudogene has no reads; reference-absent copies appear as multimappers, not unaligned reads (r1089 refuted the unaligned route). | LEDGER_O3 current state |
| NEW 11-03 Does O2 work on the reads minimap2 could not place? | Panel 284/2 vs primary 280/6; held-out after the guard 129/22 vs 125/26, 5 vs 1 where they differ; the gain is not assigning reads whose origin is outside the catalog (77% of held-out assignments before the guard). | r1509, r1510, r1514 |

### Questions with no update found
Part 3 items 3 and 12: no re-measurement since 09-20. Part 3 item 1 unchanged. (Q3 was later re-tried: FAMASM, r1479.)

Older pass notes (09-06, 09-15) and the provenance line (register 841 rows, 854 lib tests) are history; see archive lines 69-107.

## Part 0 — The object (read first)
O1's vertex set IS O2's path set: two graphs at two granularities, never "one decision rule". ~25 constants are default-reachable (not "four"); the excision abstention result is the robust-z of `de`, not the shipped α-certificate. Current O1 default: `seeded_family_definition.md` §0★ (exon-sum representatives, `mcl_families --min-shared-exon-frac 0.60`); defaults unchanged through r1484. (full: archive/2026-11/ADVISOR_QUESTIONS.md:109)

## Part 0b — What Soto 2025 is, and what it is not (READ BEFORE QUOTING ANY SOTO NUMBER)
Soto is a 98%-identity SD catalogue, not a gene-family catalogue: full S1C audit 56.5% pseudogene members, 35.9% families all-pseudogene, 35.0% bundle a fragment with a full-length member (shares later corrected to 37.3% / 31.4%). Use Soto for paired comparisons, [0.90, 1) precision and independence, never as a ceiling. Detection is class-flat (register 742). `--min-shared-exon-frac 0.30` shipped 09-14 after held-out chr5/7/21 (F 0.831 -> 0.881, precision 0.815 -> 1.000), later 0.60 in §0★. Replication: ARI 0.9698 with S1C CN, 0.9277 with our CN; 440/444 Soto families nest in one of our clusters. (full: archive/2026-11/ADVISOR_QUESTIONS.md:122)

## Part 0c — Read-throughs: what is proven, what is excluded, and where the fix lives
Positive evidence: 30 of 42 gorilla read-through junctions are used by chimp reads. Everything else excludes artefact routes (FLNC, canonical junctions, 4/46 duplicate-linked flanks, >= 3 molecules). Fix lives downstream: `-G 50k` touches 1/46 and later shown to break 99.6% of real >= 50 kb introns; polyA-inside-intron rule passes held-out; `RUSTLE_READTHROUGH_JUNCTIONS=r3` opt-in; `--bridge-regroup f1v2` default ON 09-29 (r1117-1120, r1145-1151). Do not claim the 12 non-replicating junctions are artefacts. (full: archive/2026-11/ADVISOR_QUESTIONS.md:180)

## Part 0d — What he is actually asking for: a FAMILY-AWARE assembler, and its IGV file
He wants a family-aware transcriptome assembler: isoforms, each assigned to a copy, borrowing across members, statistical resolution of tied reads. Status: isoforms `copy_assign --gtf` exists; copy tag by position, assignment by certificate (the join was the gap); borrowing live only in assignment (FAMASM assembly mode FAIL, r1479); `igv_tracks.py` moved to the attic. Hard-locus bakeoff 09-10: carried 0.835 -> 0.857, O2-assigned 0.691 -> 0.947, composite 0.846 vs flair 0.574; re-measured 09-25: all hard molecules 0.858 vs isoseq 0.732, contested stratum reverses 0.692 vs 0.813 (r1112). (full: archive/2026-11/ADVISOR_QUESTIONS.md:209)

## Part 0e — "Not a transcriptome assembler" AND "compare to flair/StringTie": the earlier reading
Earlier reading (superseded by Part 0h): the substrate is not the contribution; answer with a family-level SUBSTITUTION test (their loci into steps 6-7). Foundations near-identical (nodes 96.3% byte-identical). Substitution never run; transcript/copy-level bakeoffs instead: gorilla 10/15 cells (r1056-1058), ours >= baselines on E2 NPIP/TBC1D3 (r1189-1193), overlapping genes ours > StringTie/FLAIR not isoseq (r1256/r1259); SQANTI single-copy chr20 ours FSM+ISM 30.2% / artifact 5.0%, StringTie / FLAIR / IsoCon WORSE (r1482). Say "more copies modeled, junctions at StringTie's level", never "better junctions". (full: archive/2026-11/ADVISOR_QUESTIONS.md:258)

## Part 0f — How we guard against tuning: hold a substrate back
The deciding substrate is named before the change and was not used to build it. Killed: exon-overlap rule, three dev substrates agreed, held-out chr16+18 NPIP 26/26 -> 24/26, reverted (register 760). Shipped: `min_shared_exon_frac` 0.30 after held-out. Later instances: F1 bridge EFFECTIVE held-out (r1145); R_J FAILS held-out chr6 (r1234-1236); `candidate_copies` flip reverted after control A14 35/53; the 10-29..11-01 family rules judged on held-outs, all FAIL or opt-in (r1469-1476). Lesson recorded: MCL cluster ids are not stable across rebuilds; quarantine by coordinates. (full: archive/2026-11/ADVISOR_QUESTIONS.md:295)

## Part 0g — Yes, this is assembler-shaped work. Say so, then say what full-length changes
Concede: a StringTie clone was built and retired (41 modules dead; Rust tree now 14 files, r1249). With FLNC reads assembly is grouping + attribution: NPIP 15,960 molecules, 1,818 observed spliced structures, 703 emitted, 0 that no read shows. StringTie recovers more read-supported junctions (58.6% vs 42.1%); on >= 2-read chains Rustle is most sensitive on gorilla, 93.4% vs isoseq 82.5% (r1111). (full: archive/2026-11/ADVISOR_QUESTIONS.md:331)

## Part 1 — The cross-examination: "this is luck or overfitting"
Five moves that follow from his priors; each answered from the tree. (full: archive/2026-11/ADVISOR_QUESTIONS.md:368)

### 1.1 "You tried hundreds of things. The winners are the tail of a null."
Answer with the accounting: (at 09-15) 838 killed routes, ~330 sections, 135+ flags of which 7 ON, 3 defaults ever changed; since then the register runs to r1484 and six more defaults flipped 09-21..10-01 (junction majority r960, strict `--assemble-only` junctions r1074, GOOD seeding r1060, §0★ exon frac 0.60, cov-shorter 0.70 + f1v2, CAT annotation); none flipped by r1469-1484. Derived, not tuned: `d_max(L) = 1 − exp(−ln(L)/(T_CORE·L))` validated 20/20 on SEDEF; minimap2 k = 19 seeding model on 10,800 alignments. Concede: many hypotheses post hoc. Cleanest instance: two-sided coverage prediction committed as `5cbced4` before the run, 5/5 criteria, 9/9 edge outcomes on a human panel (rates move 2-3 events; do not say "halved"). (full: archive/2026-11/ADVISOR_QUESTIONS.md:373)

### 1.1b "You pre-registered something you already knew would pass."
`RUSTLE_ER_COVERAGE_LONGER_FLOOR` passed 5/5, 9/9 yet refused a default: D4 NULL (copy loss annotation-neutral), D5 fails (18 families deleted); named confound (length) was the wrong one, register 640. Re-priced 09-24: chr16 Compara pair precision 0.268 -> 0.529 but gorilla recall -24%, KEEP OPT-IN (r1099). (full: archive/2026-11/ADVISOR_QUESTIONS.md:434)

### 1.2 "You tuned the thresholds until the answer appeared."
Every default sits inside a measured plateau (identity 0.60-0.80, coverage 0.10-0.50, 100-500 bp unchanged; walls at 0.85, 0.60, 1 kb; §6ez); γ = 0.20 is the F1 peak of an arm trying to beat it. The 09-15 table describes the LEGACY catalog; 09-25 default: >= 300 bp, >= 70% identity, >= 30% longer / >= 70% shorter coverage, 60% shared exons, MCL 2.8; C insensitive 0.40-0.80 (r1014). Concede the one-sided coverage floor (NPIPB6 fragment 0.948 into a 38,653 bp node); fix stays opt-in (1.1b). Later relaxations of the gate (union r1435, lift r1437-1439, chains r1475) each bought recall with precision. (full: archive/2026-11/ADVISOR_QUESTIONS.md:462)

### 1.3 "Everything you have is one family."
Concede: 189/~330 sections mention NPIP. Pre-registered held-out 53 families / 8 chromosomes: F 0.7016 vs dev 0.7183 (§6s8); holds on CAT (F .7523 vs .7194). The 0.052 recall cap is a divergence cliff (protein identity < 0.60), not an edge defect. Breadth since: LRPAP1 (11 loci), few-copy class 11/377 pass (r1255), famsim, Yoo Y1/Y2; 10-29..11-01 family rules scored on four samples A119b, KB3781, OR6737, HG002 (r1469-1476). But the 10-06 audit: every NPIP/TBC1D3/YAG number at HEAD is DEV. Cross-substrate replication 87.06% of edges, shares the substrate (corroboration, not independence). (full: archive/2026-11/ADVISOR_QUESTIONS.md:512)

### 1.4 "Your validation is circular."
Present the self-inflicted retractions (ratio medians, "bases explained", prediction ⊆ truth 3 killed, conditioned denominators 7 killed, §6br zinc-finger 776/788, §6j9 projected truth 0/9 real; later: max-overlap resolver, truth-universe intersection, bipartite precision not tie-invariant; §6o1's F retracted). Non-circular O2 result: excision migrants < 50% stratum TPR 0.5066, FPR 0.0280, AUC 0.7995 vs MAPQ AUC 0.4944 (chance). Excision wall: NPIP 13 99.85% identical to copy 12, 12 of 32 reads absorbed with a valid certificate. (full: archive/2026-11/ADVISOR_QUESTIONS.md:565)

### 1.5 "You would not show me the failures."
Hand him the register (then 841 rows, now through r1484), each with its killing number, plus the deleted verdict-tag note (11/22 = 50%). (full: archive/2026-11/ADVISOR_QUESTIONS.md:616)

## Part 2 — The standing questions, with the evidence
(full: archive/2026-11/ADVISOR_QUESTIONS.md:629)

### Q1. "Do you have a method that identifies multi-copy gene families?"
Yes. Default since 09-25 (§0★): exon-sum representatives, `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`, cov-shorter 0.70, f1v2; triangle-supported leaders stay opt-in. NPIP at HEAD: 24/26 copies in one family, strict E-found 8/25 (r1250) (the 09-15 "4 clusters 25/3/2/1" is superseded). Remaining gap = annotation content (F 0.955 annotation nodes vs 0.726 de novo). All 10-29..11-01 alternatives FAIL or stay opt-in (r1469-1476, r1480); semi-guided PASS not default (r1472). Do not quote the pre-fix SD-core false-merge numbers. (full: archive/2026-11/ADVISOR_QUESTIONS.md:631)

### Q2. "Are these real families, or artifacts and overfitting?"
Narrowly defensible. External support only at >= 0.90 while 86.31% of edges sit below 0.90; 61.4% of nodes have no reciprocal gene match. Two reference annotations agree with each other only at F 0.67-0.81 (§6kl): cite before conceding any F in 0.6-0.8. Protein families F 0.977 (first pass was a superfamily artifact; at aa identity >= 0.50 sensitivity 0.867). Soto under-merges, never over-merges: pair precision 0.711 -> 0.888 (r945/946); protein edges NOT adopted (r906); nothing beats MCL held-out (r922); protein reconciler on de novo families NULL (r1480). Scoring rule: pair precision/recall only (pairscore.py). (full: archive/2026-11/ADVISOR_QUESTIONS.md:681)

### Q3. "Does the method borrow information across the family?"
In O1 assembly no: cross-copy borrowing inert/dead (§6bd). Shared multimappers are evidence; `E_c` splits large families 0.2406 -> 0.5431 (§6bi). Re-tried as FAMASM (r1479, 10-31): pool a family's reads, define chains family-wide, emit at a copy only with >= 2 certified reads. Held-out HG002 FAIL: multi-copy recall 11.3 -> 12.7 (first positive gain) but precision 14.0 -> 8.5 (~20,000 appended transcripts for 219 new matches), single-copy precision 19.2 -> 18.3; A119b recall +0.7, precision 12.7 -> 8.0. Not adopted; borrowing helps ASSIGNMENT, the 2-read certificate is not selective enough for ASSEMBLY. (full: archive/2026-11/ADVISOR_QUESTIONS.md:726)

### Q4. "What about isoforms? Could two extremely similar copies produce the same isoforms?"
Genomic read-star separates structure from origin (introns are `N`); identical isoforms tie at K = 0 (26 NPIP reads); 456 NPIP-proper reads abstain. NPIP chain comparison 10-04: only 27% of multi-junction reads are an annotated chain; sibling copies share structure, so only PSVs separate them (r1239). (full: archive/2026-11/ADVISOR_QUESTIONS.md:746)

### Q5. "Two tandem near-identical copies share a read. What happens?"
Three cases: K = 0 tie (26 NPIP; 90 of 24,462), read-through clipped locus -> `origin_rejected` (row 714), nested units special case (row 721). Tandem sim 09-24: minimap2 makes 0 cross-copy chains at identity 0.90-0.995; a cross-copy junction is a chimeric molecule or broken splice site (r1095). 21-30% of the canonical-unannotated bin replicate in both apes (§6l1); majority undetermined. (full: archive/2026-11/ADVISOR_QUESTIONS.md:756)

### Q6. "Does it port to other families, tissues, apes?"
Family/tissue yes (§4l). Apes port structurally; the 09-15 NPIP fragmentation (HSA 26/56, PTR 15/42, PPY 27/172) described the legacy catalog. Since: famsim on human SNRPB and orangutan SPAG7; F1 bridge effective on both gorilla samples (r1145); gorilla TBC1D3 held-out (r1155); family rules run on human (A119b, HG002) and gorilla (KB3781, OR6737) (r1469-1476). (full: archive/2026-11/ADVISOR_QUESTIONS.md:782)

### Q7. "Aren't the TSS/TES/UTR boundaries too convenient?"
Premise false: Soto slice 272/362 matched, median size ratio 1.00 (§6fs). Concede the coverage statistic was one-sided until `cov_longer`. 5' boundary is dispersion, not bias (d5 229 vs d3 22 bp); extension rules refuted (r973-975); `--polish-tss`/`--polish-tes` opt-in; NPIP capped TSSs 3.7-51 kb downstream of CAT first exons (r1242-1244). (full: archive/2026-11/ADVISOR_QUESTIONS.md:798)

### Q8. "1/k is a bad assignment for tied multimappers."
Agreed, never used: assign-or-abstain. Sim NPIP 95.5% right / 0.12% wrong (precision 0.9987); contested anchors 5 assigned, 5 right, 0 wrong, 28 abstain; stable-seed sim 163/163 (r1109); `--union-certificate` ties every chr16 MAPQ-0 read (r1103). At sample scale (10-31) HG002 141,265 tied reads -> 17,805 assigned, 62,158 K = 0, 67,571 origin-rejected (mostly origins outside the catalog). `--tie-ms` port FAIL as registered, opt-in (r1483); O2 levers open (r1481). Eichler 10-AS reported beside (r1477). (full: archive/2026-11/ADVISOR_QUESTIONS.md:810)

### Q9. "NPIPA and NPIPB should be distinct subfamilies."
Human CHM13 guided pipeline: gene-body unit plus intron tree recovers NPIPA|NPIPB 10/10 and every Dishuck subfamily; TBC1D3 1.000/1.000; AMY near miss (0.933 vs 0.967, annotation defect). Gorilla: 0 subfamily groups recovered (executed negative, cross-species label proxy). Unit is family-specific (concession). On the shipped single cut 21 chr16 NPIP genes are in ONE cluster (§6s8); most-junctions representatives split NPIP on dev (post hoc, r1234-1236); reads-x-junctions splits NPIP 11 -> 10 on HG002, FAIL (r1473); hierarchy answer 10-03. (full: archive/2026-11/ADVISOR_QUESTIONS.md:857)

### Q10 / Q11. Non-canonical junctions; PSV credibility.
Q10: he is right; sites recur; `RUSTLE_JUNCTION_MAJORITY` DEFAULT ON 09-21 (r960); strict junctions under `--assemble-only` (r1074); gorilla NPIP non-canonical junctions are exonic deletions/SVs vs the reference (r1152-1155). Q11: both PSV objections answered (§6aj). PSV ceiling 10-03: of 3,981 gorilla copies 639 identical, 670 beyond single-read reach. The old "do not build the VCF" is superseded: the VCF was built (r1484, see `NEW 11-01` row: bcftools agrees 54.6% as run, 93.9% with MAPQ reset, post hoc). (full: archive/2026-11/ADVISOR_QUESTIONS.md:922)

## Part 3 — What we concede before he asks
State these unprompted. (full: archive/2026-11/ADVISOR_QUESTIONS.md:930)
1. No experimental validation (no ddPCR, qPCR, FISH); DNA-side corroboration only (WGS k-mer dosage, 17/516 families differ; haplotype assemblies as O3 truth).
2. External corroboration covers the wrong end (>= 0.90 vs 86.31% of edges below); the ceiling is a protein-identity cliff at 0.60 (§6t1); chr16 Compara 80-90% loss = coverage clause (r1096).
3. 72.7% of families (88/121) without external adjudication; legacy-catalog figure, not re-measured.
4. NPIP-bound on the thesis organism: breadth added, but every NPIP/TBC1D3/YAG number at HEAD is DEV (10-06 audit).
5. Tandem "undetermined" bin: 21-30% cross-species-confirmed; r1095: the aligner never joins exons across near-identical copies.
6. NPIP fragmentation universal: superseded (guided chr16 21 NPIP in one cluster; HEAD 24/26 in one family, r1250).
7. O1's remaining losses: annotation content (F 0.955 vs 0.726; RNA-only nodes cap F 0.873-0.893); at NPIP the loss is the REPRESENTATIVE: 14 of 22 spliced-expressed copies not FOUND (r1250); R_J fails held-out (r1234-1236); later representative rules FAIL (r1473, r1475).
8. O2's wall is local identity: 99.85% copy absorbed with a valid certificate (row 711); genome-wide 639 identical + 670 beyond reach of 3,981 gorilla copies.
9. Annotation-vs-annotation ceiling F 0.67-0.81; CAT/Liftoff v2.0 human default since 10-01.
10. O3 flags a sequence, not a copy; signature failed on 5 fresh excision targets; `--flag-missing-copies` `not_tested` on both real candidates. Since: `candidate_copies` OPT-IN, held-out A13 PASSES (7,898 right = 45.7%), controls fail after Amendment 15 (r1247).
11. O3 power near zero genome-wide (10 and 6 of ~2,000 testable copies reach the test); RNA screen recall 0/17 on WGS-differing families; KB3781 53 candidates, 1 DNA-confirmed (r1090/1091).
12. No single "unit" for subfamily membership (NPIP introns, AMY exons, TBC1D3 either); not re-measured.

## Part 4 — What to put on screen, in this order
1 register and deleted verdict-tag note; 2 forking-paths accounting (count must include the six 09-21..10-01 flips; none since); 2b annotation-vs-annotation ceiling; 3 pre-registered human panel 5/5, 9/9 with commit timestamp; 4 excision O2 result with MAPQ AUC 0.4944 control; 5 cross-substrate 87.06% / 95.6% with the shared-substrate caveat; 6 NPIP coverage table; 7 `REPRODUCE.md` (now at repo root). Do not lead with catalog sizes. (full: archive/2026-11/ADVISOR_QUESTIONS.md:995)

## Part 5 — Where each question is answered

| question | current verdict | evidence |
|---|---|---|
| Q1 method exists | §0★ default (MCL, 0.60 shared exons); NPIP 24/26 in one family; 10-29..11-01 alternatives FAIL / opt-in | §0★, r1250, r1469-1476 |
| Q2 real vs overfit | narrowly defensible; annotation ceiling F 0.67-0.81 | §6kl, §6ko, r945/946, r922 |
| Q3 borrowing | no in assembly (inert; FAMASM FAIL r1479); `E_c` splits; borrowing lives in assignment | §6bd, §6bi, r1479 |
| Q4 isoforms | structure and origin separated by the read-star; chains cannot separate siblings, PSVs can | §6fd, r1239 |
| Q5 tandem reads | K = 0 / clipped / nested; aligner never chains across copies | §6fh, row 721, §6l1, r1095 |
| Q6 portability | tissue/animal yes; apes structurally; legacy fragmentation superseded | §4l, r1145, r1155 |
| Q7 boundaries | premise false; one-sided coverage conceded | §6ay, §6ba, r973-975 |
| Q8 1/k | never used; certified assign-or-abstain; tie-ms opt-in | §6fa-§6fj, r1109, r1483 |
| Q9 NPIP subfamilies | human split 10/10; gorilla executed negative | §6jm-§6js, §6kz |
| Q10 non-canonical | he is right; majority default ON | r960, r1074 |
| Q11 PSV | answered; VCF now written, O2 = pileup + caller | §6aj, r1484 |

(full: archive/2026-11/ADVISOR_QUESTIONS.md:1025)

## Part 0h — THE ADVISOR'S CLARIFICATION OF WHAT THE PROJECT IS (relayed by the user, 2026-09-09)
Supersedes Parts 0d-0g: the project IS an assembler, closer to flair / isoseq collapse than StringTie; focus = multi-copy families and multi-mapping reads; deliverable = a GTF judged on the transcripts other assemblers struggle with; O2 = all AS-tied reads, unique mappers irrelevant. (full: archive/2026-11/ADVISOR_QUESTIONS.md:1041)

### What this settles
GTF judged at hard loci (junction recovery vs StringTie 42.1% vs 58.6% is not a loss); compare against isoseq collapse, which drops 23.1% of family-region molecules (30.2% of ours that abstain); O2 scope `--as-tied-only` (09-09: human MCL0 1,118 contested, 262 / 512 / 344; "38.7% assigned" retired); `copy_assign --gtf` deliverable: phantoms at twin copies ours 2 / flair 4 / StringTie 15 / isoseq 45, coin-toss isoforms one arbitrary address ours 0 / flair 76 / StringTie 20 / isoseq 200; flair-like, no flow model. (full: archive/2026-11/ADVISOR_QUESTIONS.md:1053)

### What it changes in the standing framing
"NOT an assembler" -> "an assembler for the loci other assemblers cannot resolve" (`THESIS_OBJECTIVES.md` no longer contains "assembler"; `REPRODUCE.md` still opens "not an assembler project"); bakeoff needs a hard-locus restriction; user decides wording. (full: archive/2026-11/ADVISOR_QUESTIONS.md:1083)

## Reference bias in copy content (advisor, 2026-10-25): "the reference is an average; other individuals miss or gain copies"
Conceded at the SD level (register 1416): Jim's maternal vs paternal haplotypes differ in copy number in 44.8% of 5 kb SD windows (noise 6.6%) and 7.4% of exons (noise 0.2%), exons of 24% of our 463 families; whole expressed gene copies most stable (5-7%, register 1408). The reference is one animal's two haplotypes stitched, not an average. WGS 12x confirms 94% (SD) / 77% (exons) of differences >= 25%; confirmed core 3.0% of SD windows, 6.3% of families; 45% / 24% are upper bounds (register 1419). RNA detectors fail held-out (breakpoints precision 0.39 / sensitivity 0.09; certified tier 0.29 / 0.04; registers 1417, 1418), yet with true haplotypes 41.9% of exon gains leave >= 2 reads only the extra sequence explains (control 14%, register 1420): the limit is the detector. Multimappers find unexpressed copies IN the reference (PGM5, register 1421) but say nothing about the individual; copies NOT in the reference show only when expressed. (full: archive/2026-11/ADVISOR_QUESTIONS.md:1089)
Update 2026-11-08 (r1578/1579), on "the testis sample (OR6737, not Jim) misses many copies on Jim's T2T": against Jim's own reads at matched depth, divergence piles +5.7 points on multi-copy loci (borderline, not refuted as registered; +2.2 on loci both samples express) and the three-haplotype test +2.5 points on 7,056 shared loci (claim fails; Jim's 4.6% is the noise floor). Twelve loci show both signals; 7 fall in de novo families, led by NBPF (60 independent three-haplotype pairs) and IGHV-like: the known copy-number-variable families, about 7 of 446. Concede: per-copy claims in NBPF / IGH-type families are not exact on a non-Jim sample. Not supported: the sample is untrustworthy in general. Caveats: testis vs fibroblast (testis-only loci untested), reference copies absent in OR6737 untestable from RNA. `PREREG_crossind_missing_copies_2026-11-07.md`. **Update r1582 (2026-11-08):** copy_phasing (Strainy phasing + IsoCon test + allele-slot certificate; arm-A products of r1483, lookup only) certifies extra haplotypes at 5/7 testable strong loci in OR6737 and 0/7 in Jim (per copy 10/19 vs background 0.052): PASS as registered. Same reads chose and certify the loci (method corroboration, not independent data); IGHV-like has a V(D)J/SHM explanation; survivors NBPF, NACA-like, INTS4-like/CCT6A. `docs/PREREG_crossind_copyphasing_2026-11-08.md`.

Update 2026-11-08 (r1583/1584), on "in a non-T2T genome, can we see the over-collapsed tandem copies?": partly. copy_phasing's allele-slot certificate (a collapse = up to 4 haplotypes in 2 slots) on the A119b GRCh38 panel passes the blind bars as registered (10/13 missing vs 7/30 single-origin null) but a labelled check shows only 5/13 certificates are made of the missing gene's reads (collapsed 3/7); the rest are neighbours' haplotypes in big NPIP/GOLGA windows. It catches collapses on a heterozygous survivor (GOLGA8DP, GOLGA6L4); r1538's S/H pile catches the pure second haplotype on a homozygous survivor (GOLGA6L9), which the certificate misses (Theorem A free slot). Combined and held out on HG002 Kinnex: 2/6 missing, 0/7 null (FAIL bar 1, depth-limited). Concede: collapses of identical copies, or onto a homozygous survivor with low depth, need DNA or a better assembly. `PREREG_collapse_copyphasing_2026-11-08.md`, `PREREG_collapse_combined_2026-11-08.md`.
