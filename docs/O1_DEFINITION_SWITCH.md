# O1 — the two definitions, side by side, and the record of why a switch is warranted (2026-09-05)

**Purpose.** The advisor dislikes changes made midway without a visible trail. This file is that trail: both
definitions stay runnable from the same binary tree, both are scored on one truth with one scorer, and the
old definition remains **opt-in** (`gw_family_catalog`) rather than deleted. Nothing here is a restatement
of `THESIS_OBJECTIVES.md`; that edit is the user's decision and the proposed wording is in §5.

## 1. The two definitions
| | **E_r / γ-quasi-clique** (old, 2026-08 → 09-02) | **MCL on the annotation + SD core + read-supported units** (new, 09-03 →) |
|---|---|---|
| object | a family = connected component (≥2 loci) of the read-transcript homology graph `E_r`, refined by a max-γ-quasi-clique partition; P1 (domain-sharer exclusion) is a theorem; the exact partition is NP-hard, what ships is a verifiable certificate | a family = a set of annotated loci, clustered by genomic-sequence homology (MCL), each member carrying an SD core linked to ≥ half the family; the unit = the locus's read-supported exon chain (§6er, §6eu) |
| node | transcripts assembled from reads (de novo), collapsed to one representative per locus | annotation records → loci by exon-union overlap (§6ef); the reads shape the unit inside the annotated span |
| edge | pairwise transcript homology, identity ≥ 0.70, coverage of the shorter side (later `cov_longer` opt-in) | asm20 all-vs-all of gene spans (`minimap2 -x asm20 -c -X -N 50 -p 0.1`), identity ≥ 0.70, cov_longer ≥ 0.30, ≥ 300 bp, exonic floor 1 bp |
| partition | γ-quasi-clique certificate | MCL inflation 2.8, prune 1e-9 (size-safe; anchors invariant 2.0–4.0, §6ec) |
| constants | ~25 default-reachable (`ADVISOR_QUESTIONS.md` §1.2) | 0.70 / 0.30 / 300 bp / size 3 / "half" / 3 reads / 50 kb; 0.70, 0.30, 300 not yet justified (FRAMING_AUDIT G3) |
| binary | `gw_family_catalog --bam <bam> --fasta <fa> --out <cat> --homology-primary` (unit-tested in the suite; `soto_adj/arms.sh`) | `mcl_families --paf --gff --min-exonic-bp 1 [--sedef] [--bam --fasta]` then `copy_assign --families` |
| what it needs | reads only (no annotation) | annotation + genome (+ reads for units, + SEDEF for cores) |
| what the record says it does well | seed-free; NPIP false merges are Alu at the sensitive tier only (§6cr) | NPIP 31/31 from the first run (§6da); transports genome-wide (§6dd); duplication block ≠ family certificate (§6ef/§6eg) |
| where it failed | node construction ~58 % of the loss (§5e); NPIP fragmented into 5–6 families (§5j); 26.8 % single-exon stubs; ape catalogs not reproducible | annotation as node loses Soto members the reads found (§6ev, below); RNA admission of unannotated loci NOT implemented; three O2-facing defects §6eu |

## 2. One truth, one scorer, both definitions (§6ev; `docs/soto_two_definitions_2026-09-05.log`)
Soto 2025's 83 families / 362 members on CHM13 v2.0, HUMAN A119b, the SAME slice (`soto_adj/regions.bed`),
the same reads (`soto.bam`), `bench/soto_adjudicate.py` unchanged; MCL adapted by `bench/mcl_to_cat_copies.py`
+ `bench/mcl_edge_dump.py`. Pre-registered (`docs/PREREG_soto_two_definitions_2026-09-05.md`, md5 `507831b1…`).
Gene-overlap floor ≥ 50 % (the pre-registered unit); Wilson 95 % in brackets.
| | E_r / γ (old) | MCL units (RNA disposes) | MCL all members (DNA proposes) |
|---|---|---|---|
| Soto members detected | 128/362 = 0.354 | 238/362 = 0.658 | 254/362 = 0.702 |
| pair precision, **[0.90,1.00) band** (the only band Soto adjudicates) | **149/153 = 0.974** [0.935, 0.990] | **460/482 = 0.954** [0.932, 0.970] | 497/524 = 0.949 [0.926, 0.964] |
| pair recall, both detected | 167/191 = 0.874 [0.820, 0.914] | **560/596 = 0.940** [0.918, 0.956] | 608/646 = 0.941 |
| pair recall, ALL Soto pairs | 167/966 = 0.173 | 560/966 = 0.580 | 608/966 = 0.629 |
| family exact-set match | 21/33 | 40/56 | 42/58 |
| pairs asserted in [0.80,0.90) (Soto silent) | 21 (7 TP) | 208 (95 TP) | 219 (104 TP) |
Any-overlap floor: detection old 0.809 vs MCL 0.751 / 0.774; band precision 0.889 vs 0.904 / 0.904.
**Predictions:** P1 (band precision ≥ 0.95) held for the units arm (0.954) and missed by 0.001 for the
all-members arm (0.949); the two definitions' band precisions are within each other's CIs. P2 (recall|both
detected ≥ 0.87) held: 0.940. **P3 (detection ≥ 0.80) FAILED: 0.66–0.70** — the annotation's gene spans reach
fewer Soto members than the reads did at any-overlap (0.77 vs 0.81); the node the new definition uses is the
annotation, and the annotation does not contain everything the reads express. ⚠ Both arms are slice-conditioned
(only Soto's neighbourhoods exist, so genome-wide false merges cannot occur in either); the comparison is paired,
the levels are not absolute; Soto is SD98 and CAT-bounded (precision understated, sub-0.90 bands unadjudicable).

**Reading.** On the advisor's own benchmark the switch costs nothing measurable in precision where the benchmark
can see, gains recall through the node (3.4× at the strict floor), and asserts ten times more pairs in the band
the benchmark cannot see — which is where any objection will come from. The recall gain is a NODE effect
(annotation vs read-assembled transcripts), not a partition effect, and the annotation-node loses ~4 % of members
at any-overlap; **RNA admission of unannotated loci is the missing stage** (the pivot's own second criterion,
`project_soto_descoped_mcl_pivot`), and it is also where O3 lives.

## 3. Gorilla anchors, from the record (not re-run today)
| anchor | old | new |
|---|---|---|
| NPIP loci recovered (31 truth) | 14/31 reads-only; 30/31 with annotation seeds (§6be) | 31/31 (§6da), 29 loci one family after the locus rule |
| NPIP as one family | fragmented into 5–6 (§5j) | one family; LCR16u (SMG1P) separate; block vs family certificate 29 vs 48 loci |
| false-merge rate / corroboration | 1.33 % [0.37, 4.73] on window sets (§6bt.1) | SEDEF corroborates 39/46 clusters, 3 artefacts by the repeat library, MCL4 element-welded (§6dy–§6dz) |
| tandem array (CGB/NTF) | — | 60+58 at prune 1e-9, dissolved at 1e-5 (§6ec) |

## 4. The gaps, tracked (from `FRAMING_AUDIT_2026-09-05.md`, with owners)
| gap | what closes it | status |
|---|---|---|
| G1 two O1 definitions coexist; "family = MCL output" | restate the family as an object (§5), MCL as pre-clustering; keep E_r opt-in | **decision pending (user)** |
| G2 "at the RNA level" is true of the node only | restate O1's wording; state the five failed RNA-membership tests | with G1 |
| G3 0.70 / 0.30 / 300 bp / size 3 / 50 kb unjustified | sensitivity table on the anchors (the inflation instrument) | open |
| G4 O2 certificate conditional on the column set | chapter §4 wording; D3 pre-registration (row 691) | open |
| G5 (new, §6ev) annotation node misses expressed loci; RNA admission not implemented | admit read-supported loci outside the annotation into the unit table (different criterion ⟹ non-circular) | open |
| Q9 NPIPA/NPIPB | run the core rule on the 29 loci with the CHM13 labels | open |

## 5. Proposed restatement of O1 (for the user to accept, edit, or refuse)
> **O1.** A multi-copy gene family is a set of loci that share a duplicated core: each member's core is the
> segment of its locus linked by segmental-duplication pairs to at least half of the set, and the family's
> units are the members' read-supported exon chains within their annotated loci. Candidate sets are proposed
> by MCL over genomic homology of the annotation and certified by the core rule, SD corroboration, the curated
> repeat library and the CHM13 landing. The earlier definition — a γ-quasi-clique partition of the read-transcript
> homology graph E_r — remains available (`gw_family_catalog`) and is compared with this one on one truth in
> `O1_DEFINITION_SWITCH.md` §2.
