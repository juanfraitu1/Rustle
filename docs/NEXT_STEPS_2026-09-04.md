# Next steps after 2026-09-04 (ledger §6ds–§6ee)

Status: **O1 is two measurement sessions from a conclusion with stated limits. O2 has one honest claim on
NPIP (an abstention certificate) and needs a node fix and a new truth before more can be claimed.**
Canonical catalogs: `mcl_ann/rna_bp1_p9.clusters.tsv` (3 contigs) and `mcl_ann/gw_bp1_p9.clusters.tsv`
(genome-wide), both at `--prune 1e-9` (§6ed). Every number quoted from a catalog must say which prune.

## O1 — to conclude
| # | step | why | where it stands |
|---|---|---|---|
| O1-1 | **Non-overlapping-locus NODE rule** in `mcl_families` (overlapping annotation records on a contig = one locus; representative = longest exon-union) | nested annotations are the same DNA twice: 13 NPIP copy pairs, 607/1,221 O2 ties (§6ee); the nested MCL4 lncRNA (§6dy); NPIP giant models are real transcription units (§6ea) | **SHIPPED OFF, MEASURED (§6ef)**: exon-union-overlap loci, edges attributed to the representative; 3-contig 142→116 clusters, NPIP 44 records → 33 loci; ⚠ **NPIP reunites with MCL7 (U1) into 48 loci — the anchored boundary moves; the default flip is the user's call — §6eg: DO NOT flip with attribution semantics, it merges LCR16a with LCR16u (the duplication block); the real defect is CHIMERIC models (ABCC1+NPIP, SORL1+NPIP, EIF3C/NPIPB overlap).** O2 on it: resolved 17/43 vs 9/39, Containment 10→2 (the 2 share no exon), MAPQ-0 still 1/78 |
| O1-2 | Run the audit columns on the genome-wide catalog `gw_bp1_p9`: curated-repeat fraction, member union-exonic-coverage, size flags | the 46-cluster audit (§6eb, `docs/audit46_2026-09-04.tsv`) was on the 3-contig 1e-5 slice; the reported object is genome-wide | scripts exist: `mcl_ann/adj/audit46/audit46.py`, `adj/gwaudit/gwaudit_p9.py` |
| O1-3 | **Cut certificate**, pre-registered: "members share ≥ half their exonic ≥0.85 partners outside the cluster" (§6ed) — positive control = the tandem 60/58 split (must flag each other), negative = NPIP | the 13 cell-A cuts are uncertified (§6dy); this is the first instrument sharing no input with MCL's edge density | statistic implemented (`gwaudit_p9.py`), thresholds NOT pre-registered as a cut test yet |
| O1-4 | Library-silent element handle: a de novo-family column (`rnd-*` family share + young-mode contig spread) | MCL4's element is absent from curated Dfam and "Unknown" de novo (§6dz); the curated column misses it | scripts in `adj/rmsk/`; needs the full-genome de novo .out parsed per cluster |
| O1-5 | Two-sided exonic clause at COPY level (γ) | priced on edges only (§6dy cell B) | offline candidate; run `mcl_families` with a both-sides floor and diff clusters |
| O1-6 | Write the limitations: node boundaries are copy-inconsistent in the annotation and RNA cannot fix them (§6ea); cohesion is slice-conditioned; families >~1,600 near-identical copies exceed even prune 1e-9 at I=2.8 | | text only |
| O1-7 | Restate §6dd/§6de (inflation plateau, "zero clusters >100") under prune 1e-9 in `docs/ADVISOR_QUESTIONS.md` | they were prune artefacts (§6ec) | text only |

## O2 — to claim anything beyond NPIP's abstention certificate
| # | step | why | where it stands |
|---|---|---|---|
| O2-1 | Same as O1-1: O2 must receive non-overlapping loci | 9/39 → 16/27 resolved copies when nested copies are dropped (§6ee) | DONE via O1-1 (`o2scale/fam_NPIPloci_073242`); ⚠ O2 peaked at **20.7 GB of 25** on 43 loci — memory must be looked at before any sweep |
| O2-2 | **A truth for the correction leg's re-threaded reads** — §6ej: junction-anchored truth BUILT (`o2scale/junction_truth.py`; 52/52 on the 3-copy family, 43/43 among MAPQ<60) but **0/31 re-threaded reads are anchorable**; simulation is the remaining route | `--vg-realign-correct`'s +31 assignments are reads whose primary lies OUTSIDE every copy; unique-mapper agreement cannot score them (§6ee) | design: simulate reads from copy sequences (the existing sim harness, `-N 50`), or use copy-specific junctions as truth |
| O2-3 | Report O2 per copy WITH O1's nearest-paralogue identity as the abstention forecast (ρ −0.55; ≥0.99 ⟹ share 0) | the certificate that ties O1 to O2 | column to add to the joint catalog |
| O2-4 | Decision with the advisor: is "abstention certificate + identity forecast" the O2 result? | on NPIP in fibroblasts the multimapping pool is 592/6,545 reads and 98.6% of it is unassignable | user/advisor decision |
| O2-5 | Optional: the 86-family sweep (`o2scale/run_armA.sh`, resumable) EXCLUDING cell-B artefacts; ⚠ the 39-copy family peaked at 18.8 GB of 25 | only if O2-4 says O2 is reported catalog-wide | ready |
| O2-6 | Dispersed families: scope the claim ("assignment for local families, declared abstention for dispersed") or build a genome-wide mode | the genome-wide catalog is mostly dispersed families | user decision |

## Added 2026-09-04 late (§6ej)
| # | step | why |
|---|---|---|
| O2-8 | ~~junction/PSV conflict abstain~~ **IMPLEMENTED OFF and REFUTED (§6ek, register 681)** — O2's junction term is spliced-space, the discordance is genomic. **O2-8b DONE (§6el): read-chain copies → 5/5 agreement (was 7/11), control 57/57 (was 52/52), NPIP assigned 19→96, tied 3,347→2,766.** ⚠ 3-copy control peaked at 25.4 GB. Next: O1-10b make the read-chain node a feature (not a script); O2-5 memory. Original text: Fix the main path's identifiability bound and the copy sequences: `min_p` from the read's pairwise evidence vs its two best copies (`psv_decisive_count`, §6dg) instead of star-projected columns; copies = LOCUS sequences (core hull + read-supported exon chain), not GFF models | on NPIP cores O2 assigns 11/117 junction-anchored reads and 4/11 contradict the junction with min_p 1e-16..1e-35 while the reads differ from the candidates by 0–8 edits (register 678/680). Acceptance = the 117 anchored reads: 0 confident disagreements |
| O2-9 | State O2's NPIP result as: abstention on 3,347/3,407 molecules of the LCR16a cores (507/535 of minimap2's MAPQ-60 reads tie); junction-anchored agreement 7/11 — an abstention certificate plus a defect list | user/advisor framing |

## Added 2026-09-04 late (§6eg)
| # | step | why |
|---|---|---|
| O1-8 | **Name the O1 object**: duplication block (SEDEF-level, what attribution builds) vs gene family (LCR16a core, what the advisor means) | NPIP+MCL7 = LCR16a+LCR16u; the two are distinct families by origin (Johnson 2001) |
| O1-9 | ~~Cross-map every anchored family's copies to CHM13~~ **DONE for all 46 (§6eh, `docs/all46_chm13_landing_2026-09-04.log`)**: real families = one stem/one chromosome; artefacts scatter; 16p block families = mixed stems | make it a catalog column (needs CHM13 mapping in the pipeline or as a post-step) |
| O1-10 | **DUPLICON-FIRST node segmentation (§6eh pilot WORKS)**: per record, SEDEF partner-depth profile vs the candidate family; core = maximal segment with depth ≥ half the family; node = the core's exonic bases; re-cluster over cores. Acceptance on NPIP: 28 loci with a 19–30 kb core (median 23.4 kb = LCR16a), EIF3C (808 bp) / ABCC1-region (0) / all 14 MCL7 (0) OUT; tandem, MCL4, MCL6, MCL21 unchanged | **IMPLEMENTED (§6ei): `mcl_families --core-refine --sedef`, OFF.** TRIM arm accepted on NPIP (7 chimeras → 23–25 kb LCR16a hulls; tandem, MCL6, MCL7 as predicted; gate protects ZNF/artefacts). ⛔ DROP arm fails acceptance on MCL5/MCL25/MCL13 (register 677): `dropped` is a flag until a second pre-registration with a co-signatory. §6eh addendum over all 46: eliminates MCL0/MCL4/MCL16/MCL31 outright (core 0 for every member), trims NPIP; ⚠ **gate on SD partner depth** — old ZNF/OR families (MCL12, MCL17, MCL25) have no SEDEF pairs and must stay homology-defined. Pre-register: core as a fraction of exon-union (EIF3C 18% vs tandem 100%), the depth gate, and the NPIP/tandem/ZNF acceptance set BEFORE implementing in Rust |
| O2-7 | ~~Re-read §6ee~~ **DONE (§6eh)**: 4,810 reads on EIF3C vs 1,735 on NPIP-proper; NPIP-proper MAPQ-0 0/68, 1–59 8/524; 225 assigned over 8 copies | abstention certificate stands |

## Added 2026-09-04 late (§6el)
| # | step | why |
|---|---|---|
| O1-10b | Make the read-chain node a feature: `mcl_families` (or the bridge) emits per-locus units = SEDEF core hull + read-supported exon chain; O1 clusters over them, O2 consumes them | §6el: the GFF model was the cause of O2's confident wrong calls; the node for both objectives is the transcribed unit |
| O2-5' | **O2 memory first**: 25.4 GB on a 3-copy family with a 25 kb copy; 20.7 GB on 43 loci — the banded pairwise PSV DP grows with copy length difference | blocks every larger run |

## Housekeeping
- The whole MCL line is committed (0634a0e, dcd41d7); everything after §6ed is uncommitted.
- Repeat library on disk: `winloci_data/rmsk/` (UCSC GenArk; curated Dfam 3.8 + RepeatModeler).
- Never run O2 on the cell-B artefacts (MCL28, MCL15, NEW101/MCL0): they are L1/ERV1/L1, not families.
