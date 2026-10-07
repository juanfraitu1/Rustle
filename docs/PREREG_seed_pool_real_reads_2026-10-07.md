# Pre-registration: the seeding pool (primary / good / all) and the primary-first representative on real reads, NPIP and TBC1D3 (2026-10-07)

**Written before any arm of this study exists.** The driver switch it uses (`tools/rustle_pipeline.sh --seed-pool primary|good|all`, `--seed-as-ratio`, `--contig`, `--as-table`; commits c5cd2524, dae62624, 21 tests) is committed; no pool other than the default has been assembled with it, and no number below is read from a product of this study.
User direction (2026-10-07): "yes lets do that, however lets also make it easy to set the 'good secondaries' and 'all' in order to prove to my advisor how does this decision affect the outcome" (that = the real-read test of the primary-first representative, Rule 1 of `docs/PREREG_locus_units_2026-10-06.md`, on human NPIP and TBC1D3, then gorilla).
**Label.** Human A119b: DEV (the development block of NPIP, and the TBC1D3 audit set). Gorilla OR6737: second-species check, held out for Rule 1, spent for the pool (the default pool was decided on gorilla NC_073244.2); fixed in Amendment 1 before any gorilla arm exists. Human and gorilla are never pooled. The study is descriptive: it measures what the pool decision and Rule 1 do, with the directional predictions of section 6 stated in advance and no combined pass/fail.

## 1. The question

Which alignments should seed the assembly's loci, and what does each choice cost or buy at the multi-copy families? The three pools are nested, so they are one axis, a filtration of the read pool indexed by the width rho of the tie band:

| pool | seeds | driver |
|---|---|---|
| P | primary alignments only | `--seed-pool primary` |
| G(rho) | primaries + secondaries with AS >= rho x the molecule's genome-wide best AS | `--seed-pool good --seed-as-ratio rho` (the default is rho = 0.98) |
| A | primaries + every secondary | `--seed-pool all` |

P is the limit rho > 1 and A is rho = 0, so the pools nest: P ⊆ G(1) ⊆ G(rho) ⊆ G(rho') ⊆ A for 1 > rho > rho' > 0. The second factor is the representative of each locus (the transcript that stands for the locus in the all-vs-all, the family clusters and the copy table): MR = the transcript with the most reads (the shipped rule) or R1 = the transcript with the most reads among those whose intron chain also occurs in the primaries-only assembly (Rule 1 of the locus-unit study, implemented as `bench/entangled/locus_units.py --rule1 --primary-gtf`).

## 2. What was already seen (ledger)

- Human chr16 NPIP, the 25 CAT/Liftoff copies, **before** the bridge regroup (f1v2) and the min-cov-shorter default (`docs/NPIP_READ_POOL_2026-10-01.md`, scorer HEAD-reproduced in `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md` gate G0): own node 23 / 21 / 24, E-found within own nodes 10 / 6 / 2, locus-level 19 / 18 / 16 for P / GOOD / ALL. The current default at HEAD (G(.98) + f1v2 + cov .70): own node 24, E-found within own nodes 8, locus-level 21, U2 F .645. **P and A were not re-run at HEAD; TBC1D3 was never run by pool; Rule 1 was never run on real reads.**
- The default pool on the 09-29 copy-recovery run (`docs/PREREG_copy_recovery_tools_2026-09-29.md`): per-copy numbers for the default only.
- Ideal-read simulation windows around the same two loci (`docs/LOCUS_UNITS_LEVELS_2026-10-07.md`): P finds more reachable copies (53 of 54 against 48) at a lower cluster precision (K*_C .667 against .775); Rule 1 repairs 4 gene-runs on 3 genes with no precision cost (level C only; level M it is masked by one partition flip). Simulated reads, 39 copies x 2 draws of one annotation.
- Gorilla overlap check (`docs/GORILLA_OVERLAP_2026-10-07.md`): P costs 0-1 chains at overlapping genes and 2 of 18 at non-overlapping ones (primary denominator); the secondary-pool gain is real at N (102 of 19,953 chains).
- Not seen: any NPIP or TBC1D3 number of any pool x representative cell at HEAD other than G(.98) x MR on NPIP.

## 3. Substrates and instruments (nothing is chosen after seeing an arm)

- **Human A119b** (`winloci_data/A119b.t2t.bam`, CHM13 v2.0, `chm13v2.0.fa`): **NPIP on chr16** (25 copies, `copies.hsa.tsv` family NPIP) and **TBC1D3 on chr17** (16 copies, family TBC1D3) of the CAT/Liftoff v2.0 truth (`/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann/{copies.hsa.tsv, truth.hsa.gtf}`). Best-AS table: `human_A119b.molecules.tsv` of the 09-25 run (header `bam=` checked by the driver).
- **Copy-level instruments** (unchanged since the default re-score): `bench/default_rescore/nodes.py` (own node; NPIP with `--exons-json npip_read_pool.json` as in the re-score, TBC1D3 without) and `bench/copy_support.py --nodes` (Amendment E E-found = `tc_found`: the representative starts at a capped 5' end of the copy's reads and carries its first three introns). Read under two `PYTHONHASHSEED`s.
- **Family-level** (`bench/default_rescore/score_families.py`, `family_score` of the HEAD build): U2 (NPIP union truth, 3 families), Compara Primates, Soto 2025, on the contig, bipartite and pairwise; the all-families number is reported, and the row of the target family.
- **Own node and composition of the family's clusters** (new, `bench/seed_pool/composition.py`; its own-node flags must equal `nodes.py`'s wherever `nodes.py` is defined, gate G1b): the clusters that hold a locus overlapping a same-strand truth copy are the FAMILY CLUSTERS of the arm; their loci are the nodes; a node is ON-COPY (overlaps >= 1 exonic bp of a same-strand truth copy), IN-SPAN (inside a copy's span, off its exons), ANTISENSE (overlaps a truth copy on the other strand) or ELSEWHERE. Node precision NP = on-copy nodes / all nodes of the family clusters. Unlabelled nodes count as false (no intersection with a truth universe). The cluster file names a locus by its span only, and `nodes.py` stops when two loci share one (it did on the frozen A arm of 10-01); `composition.py` gives every locus on a shared span the cluster ids of that span (a locus is in a cluster if any row of its span is) and reports the number of shared spans per arm. This is the rule for every arm, fixed here before any arm exists.
- **Cost**: transcripts, loci, all-vs-all PAF records, wall time and peak RSS of each stage.

## 4. Arms (per contig; HEAD release binaries `rustle_target_m2/release`; every other setting the driver's default: f1v2, strict junctions, shipped polish, `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`)

Primary arms: **P**, **G98** (the default), **A**, and Rule 1 on top of the last two, **G98+R1** and **A+R1**. Descriptive arms (the filtration, section 7): **G100, G995, G95, G90** and their +R1 versions. Control: **P+R1** (Rule 1 over its own primaries-only assembly; every chain is primary-supported, so it must equal P).
Each pool arm is `tools/rustle_pipeline.sh assemble` then `families` with `--contig CONTIG --as-table W/mol/mol.tsv --no-cache --seed-pool ... [--seed-as-ratio rho]` on the full BAM (the table is genome-wide for every contig). A +R1 arm copies the pool arm's `PREFIX.gtf`, writes `PREFIX_R1.families.gtf` = `locus_units.py --base PREFIX.families.gtf --primary-gtf P.families.gtf --rule1` (a newer file than the copied `.gtf`, as the driver's guard requires) and runs the driver's `families` stage on it, so the clustering command is the driver's own. The all-vs-all of A-type arms runs through `tools/mm2_shard.sh` (byte-identical to a single run, cmp-checked 2026-09-25) in bounded calls; the others use plain minimap2, as the default re-score did. Intermediates (PAF, loci FASTA, shard caches) are deleted once an arm is scored (the work disk has 38 GB free).

## 5. Gates (an arm table is INVALID unless all hold)

**G0** the G98 arm on chr16 reproduces the registered default products of `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md`: `cmp` of the assembled GTF, `families.gtf`, `fam.clusters.tsv` and `fam.loci.gff3` against `/mnt/linuxdisk/tmp/rescore_2026-10-06/DEF/` (the driver with `--contig` equals the hand-made `--region` run). **G1** the scorers on G98 chr16 give own node 24, E-found within own nodes 8, locus-level 21, U2 F .645; **G1b** `composition.py` gives the same own-node flags as `nodes.py` on the G98 and P arms of chr16 and on a synthetic case. **G2** P+R1 equals P (`fam.clusters.tsv` identical) on both contigs. **G3** `copy_support.py` is deterministic under `PYTHONHASHSEED` 0 and 1 (same summary). **G4** sharded and unsharded all-vs-all agree: the P arm on chr17 is clustered both ways (`MM2_SHARD_MIN_BYTES=1000` forces the wrapper on a small FASTA), `fam.clusters.tsv` identical.

## 6. Outcomes, denominators and predictions (fixed now)

**Outcomes per arm and family.** N = the family's truth copies (25, 16). **E** = the copies that are spliced-expressed under Amendment E (`tc_expressed`; arm-independent, computed from reads and annotation only). M1 = copies with an own node (of N). M2 = E-found within own nodes (of E). M3 = locus-level E-found (of E): any transcript of a same-strand overlapping locus starts within 150 bp of a capped start of the copy's reads and carries the first three introns of an expressed chain (`locus_tc_found`). M4 = U2 / Compara / Soto bipartite sens, prec, F and pairwise F on the contig. M5 = NP and the node classes of the family clusters; transcripts, loci, PAF records, wall time. Beside M2, reported: `ann_found_in_npip_nodes`, `chain_found_in_npip_nodes` (Amendments A, B). A per-copy matrix (copy x arm: own node, E-found) is reported, so every difference is a named copy; counts are of copies (n = 25, 16) and a difference of one is one copy, no significance is claimed.

**UNDERPOWERED rule.** A family with E < 8 gets no verdict for the predictions that read M2 or M3; its counts are reported.

**Predictions** (each is reported as held or failed; there is no combined verdict):
- **S1** (breadth costs the representative) on NPIP: M2(P) >= M2(G98) >= M2(A), with at least one strict inequality.
- **S2** (breadth does not cost membership at the default) on NPIP: M1(G98) >= M1(P).
- **S3** (breadth costs node precision) on NPIP: NP(P) >= NP(G98) > NP(A).
- **S4** (Rule 1 repairs the representative at the default) on NPIP: M2(G98+R1) >= M2(G98) + 1, M1(G98+R1) >= M1(G98), and |F(G98+R1) - F(G98)| < .01 for U2 and for Compara (chr17, S7: Compara and Soto; U2 has no TBC1D3 family).
- **S5** (Rule 1 repairs it at A) on NPIP: M2(A+R1) >= M2(A) + 2 and M1(A+R1) >= M1(A).
- **S6** (control) P+R1 is identical to P (this is gate G2; listed so that its outcome is part of the table).
- **S7** (TBC1D3) S1-S5 are evaluated on chr17 TBC1D3 with the same wording; whichever the UNDERPOWERED rule excludes is reported with its counts.
- **S8** (the default is not a knife edge) on NPIP: M1 and M2 of G95 and G995 each differ from G98 by at most one copy.
**Summary statistic.** The non-dominated set of the arms of a family on the three counts (M1, M2, NP): X dominates Y iff X >= Y on all three and > on one. The report states the set per contig and which arm each non-dominated arm is; it does not name a winner.

## 7. The filtration (descriptive)

G100, G995, G98, G95, G90 and A, with P at the far end, give each outcome as a function of the pool width. Registered reading: the curves are reported as measured; a non-monotone curve or a step between neighbouring widths is a finding about the threshold, not a failure.

## 8. Limits declared in advance

One library per species; the human block is the development block of both families; the copy truth is a single annotation (CAT/Liftoff), its pseudogene copies are partly unexpressed (E, not N, is the found denominator); E-found rests on the 5' cap signal, which exists in the human reads only (the gorilla reading is fixed in Amendment 1); the primaries of R1 come from the same reads and the same polish as the arms they correct (a chain is primary-supported iff it occurs in P's `families.gtf`, an exact intron-chain match); NP counts real but unlabelled family members as false; the all-vs-all of A took about 20 min against 46 s for G98 on chr16 (10-01, 17 times the records), so the sharded, resumable run is the only one that fits a 10-minute call.

## 9. Amendment hook

Amendment 1 (written before any gorilla arm exists, after the human results): the gorilla contigs (NPIP NC_073241.2 / NC_073242.2 / NC_073244.2; TBC1D3 NC_073228.2 / NC_073224.2), the copy truth and family truths used, the found criterion that replaces the cap-based E-found, the table of OR6737, and the predictions carried over unchanged.
