# Pre-registration: a container for the extra pieces of fused family members, tested on simulated fusions

**Written 2026-09-28, before any simulated fusion read exists and before the container code is run on any substrate.**
User request: "can we simulate or model some sort of container for the extra pieces of the genes for the fusions?"

## 0. What was seen before this file

- Real data (dev, human A119b chr16): 11 of 26 NPIP copies sit in fused loci, 10 with co-duplicated partners (PKD1P,
  PDXDC2P, SMG1P, RRN3P1, EIF3CL). Dominant readthroughs (NPIPA1, A6, A9, B14P) are how those copies are expressed.
  Segment-level MEMBERSHIP (deciding which family a locus joins from its segments) failed: it dragged partners into NPIP
  (register rows 1121-1127, `npf_variants.md`). Families are a strict partition (§6s9).
- The gorilla ideal-expression simulation (`ggo_npip_sim`, launched 2026-09-28) supplies the genome splice index, the
  25-copy gorilla NPIP truth (T_member) and the neighbourhood gene set this simulation reuses. None of its results is
  used to set anything below.

## 1. The container (binding definition)

Run AFTER the families stage; it never changes a family. Inputs: the assembled GTF, the families outputs
(`clusters.tsv`, the locus table, `loci.paf` with CIGARs). For every clustered locus m in family F:

- **Exon blocks of m:** the union of the exons of all of m's transcripts (not only the representative), merged where they
  overlap.
- **Core:** an exon block b of m is core iff some PAF record between m and another member m' of F aligns at least one
  base of b onto an exon base of m' (projected through the CIGAR; the same exon-on-both-sides test the edge rule uses,
  `--min-exonic-bp 1`). No other constant.
- **Container (accessory pieces):** every exon block of m that is not core.
- **Relation:** for each accessory block, the same test against members of every other family F'. A hit records the
  relation (m, block, F'). A family-level relation F~F' exists iff some member carries such a block.
- Output: one row per (locus, block): coordinates, core/accessory, relation families. Unclustered loci get no rows.

## 2. The simulation

- **Substrate:** gorilla genome and the 25 T_member NPIP copies (as in `ggo_npip_sim`). Reads: 30 per transcript,
  error-free, independently jittered ends, `bench/sim.py` conventions; aligned with the same minimap2 command and splice
  index as `ggo_npip_sim`; assembled and clustered with the shipped defaults (frozen `fj_bin_frozen`, driver
  `assemble` then `families`).
- **Fusion set S:** every T_member copy that has an annotated, non-NPIP neighbour gene on the same strand within 100 kb,
  paired with its nearest such neighbour p. A seeded half of S (seed 20260928) is fused; the other half, and all copies
  without a neighbour, stay unfused as within-arm controls.
- **A fusion transcript:** the upstream gene's (in transcription order) representative transcript up to its last internal
  donor, spliced to the downstream gene's representative transcript from its second exon's acceptor. Both splice sites
  are annotated, so the junction is canonical. Assert this per fusion.
- **Expression:** the partner's own transcripts at 30 reads. The copy side at 30 reads split between standalone and
  fusion by the fusion share f: standalone 30·(1−f), fusion 30·f (rounded).
- **Arms:** f ∈ {0 (control: no fusion reads), 0.1 (minority bridge), 0.5, 0.9, 1.0 (the copy is expressed only
  through the fusion, like NPIPB14P)}. Same seed, same fused set, in every arm.

## 3. Measurements and bars (per arm)

- **M1 membership.** Fused copies placed in the NPIP family (the cluster holding the most T_member copies). Bar: at
  f ≤ 0.5 the count equals the f = 0 control. f = 0.9 and 1.0 are reported, not barred.
- **M2 partner leakage.** Partner genes (by locus overlap) in the NPIP family. Bar: at every f, no more than the f = 0
  control.
- **M3 container accuracy** (truth by construction: the partner's exon bases). For each fused member, precision and
  recall of the accessory bases against the partner bases. Bar: pooled precision ≥ 0.90 and recall ≥ 0.90 at f ≥ 0.5.
  Every partner base called core is listed with the member it aligned to (the co-duplication risk).
- **M4 relation.** For each fused member, whether its container relates to the partner's family (the family of the
  partner's own locus). Reported as a count at each f.
- **M5 false accessory.** Accessory bases on unfused members (all arms) and on fused members at f = 0. Reported; a
  non-zero value is analysed per block, not barred.

## 4. Verdict

- **Container WORKS** if M3 passes at f ≥ 0.5 and M1/M2 pass.
- **Container PARTIAL** if M1/M2 pass but M3 fails; the report says which blocks failed and why (expected cause:
  co-duplicated partner pieces aligning to other members' partner pieces).
- **Container FAILS** if M1 or M2 fail. The container cannot cause that, since it runs after the families; a failure
  therefore describes the family definition under fusion, not the container.
- Descriptive afterwards, never judged: the container applied to the real human A119b chr16 NPIP family (dev), listing
  what it puts in the container of NPIPA1, A6, A9 and B14P.

## 5. Order and machine rules

1. Implement the container as a frozen Python post-processor with unit tests on hand-made inputs (sha1 recorded in an
   Amendment before the simulation is scored).
2. Wait for `ggo_npip_sim` to finish; reuse its index and truth.
3. Build the fusions, simulate, align, assemble, run families, run the container, score (in that order, all arms).
- Heavy steps via `tools/rlock.sh heavy`, light via `tools/rlock.sh light`; foreground; never `pkill -f`; TMPDIR and
  scratch under `/mnt/linuxdisk`. Gorilla only; human numbers are never pooled with gorilla.

## Amendments

**Amendment 1 (2026-09-28 17:50, before any simulated fusion read exists).** The container is frozen:
`bench/family_container.py` sha1 e197ccb376f69f7ec8492cf1026e43626d0f0032, `bench/test_family_container.py` sha1
a2f4a13d0c0035167bb654fff5f6c257e2220b31 (copies with SHA1SUMS in `/mnt/linuxdisk/tmp/rustle_figures/container_frozen/`);
20/20 unit tests pass. Ambiguities in §1 resolved by the implementer (binding from here on):
- a block is core only when one aligned CIGAR column joins an exon base of m to an exon base of m' (insertion/deletion
  bases never count);
- the partner's exon bases are its all-transcript exon union, not the representative-only exons of the edge rule;
- a locus includes the records `loci.tsv` folds into it;
- only exons sharing >= 1 base merge into one block; touching exons stay separate;
- every PAF record counts (no identity, length or primary filter);
- relations are tested on accessory blocks only; the family relation table is directed, with a `reciprocal` column;
- output coordinates are 1-based closed, blocks numbered in genomic order.
Format smoke run on human_testis BASE families (not NPIP content): 4,824 blocks, 4,107 core, 717 accessory; 21 families
(42 loci) have no core block because their loci overlap on the genome and the shared stretch aligns to itself — the
base-level CIGAR test correctly calls those blocks accessory, where the edge rule's interval test admits them. This is
recorded as a property of the definition, not changed.

**Amendment 2 (2026-09-28 19:05, after the fusions and the read pool were built, before any read was mapped, assembled
or scored).** Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/fusion_container_sim/` (`build_pool.py` sha1 d62da791,
`ann/fusions.json` 9ba2fdf5, `ann/pool.tsv` 3e8d427d, `ann/arms.json` 0d0fd0b4, `reads/pool.fq` 06290d75). Container =
the frozen copy (SHA1SUMS verified before use). Binaries = `fj_bin_frozen/` through the driver.

*A. One read pool instead of five mappings (deviation from "map each arm", design-preserving).* minimap2 maps each read
independently, so ONE pool is simulated and mapped once, and each arm's BAM is the pool filtered by read name:
- the pool holds 30 reads for every transcript that any arm needs: per fused pair the copy representative (`C.<cid>`,
  standalone), the fusion (`U.<cid>`) and the partner's own transcripts; every other transcript (unfused copies'
  transcripts, all other neighbourhood transcripts) as in `ggo_npip_sim` NB;
- arm f keeps `C` reads round(30(1-f)) and `U` reads round(30f) (30/0, 27/3, 15/15, 3/27, 0/30). The subsets are
  NESTED: a per-transcript permutation of k00..k29 (seed stable_seed('20260928','subset',tid)), prefix of length n;
- every other read is identical in every arm. Read ends jittered 0-30 as `sim.py`, seed stable_seed('20260928',tid,k)
  (new reads, not the `ggo_npip_sim` ones); mapping = the `ggo_npip_sim` command and index, in batches.

*B. §2 ambiguities resolved (binding).*
- Neighbour = a gene/pseudogene/ncRNA_gene record of `ggo_npip_sim/ann/ggo3.gff` with >= 1 transcript, same strand,
  span gap to the copy territory <= 100 kb; overlapping spans count as gap 0 and are allowed.
- Non-NPIP = not a T_member native record, no NPIP-like description (`score.py` NPIP_LIKE: NPIP / titin / NACA /
  SRRM2), not audit class NPIP.
- "Such neighbour" must be FUSABLE, else the next nearest is taken: representative >= 2 exons (transcripts >= 120 nt,
  non-canonical introns merged as in `ggo_npip_sim`), the partner's exon union shares no base with the copy territory
  (skips EIF3C at NPIPA5, SMG1 at NPIPB11, myosin-11-like at NPIPB2, PDXDC-like at NPIPB14P, SMG1-like at NPIPB5: their
  models carry territory exons), and the donor lies 5' of the acceptor. Nearest = smallest gap, ties by start and id.
- Representative = most junctions, then longest, then id; for the copy this is its T_member primary chain (clipped to
  the territory, as `ggo_npip_sim`). Half of S = random.Random(20260928).sample(sorted cids, floor(|S|/2)).
- **S = 20 copies; fused (10):** NPIPB4~SMG1-like LOC129527585, NPIPB13~SMG1-like pseudogene LOC129527593,
  NPIPA7~LOC129527614, NPIPB8~LOC129527618, NPIPA5~PDXDC1, NPIPA9~ATXN2L, LOC124907834~CNOT3-like LOC129527626,
  NPIPB2~SNX29, NPIPB7~SMG1-like LOC129527696, NPIPB12~SMG1-like LOC101129998. **Controls (10):** NPIPB1P, NPIPB11,
  NPIPA8, LOC124907808, NPIPA2, LOC128966608, NPIPB3, NPIPB6, NPIPB5, NPIPB15. Not in S (5): NPIPB10P, NPIPA1,
  LOC124907807, NPIPA6, NPIPB14P. All 10 fusion junctions assert canonical (GT-AG).
- "Copy side at 30 reads": for a FUSED copy every simulated transcript sharing >= 1 exon base with its territory on its
  strand is removed in every arm (its isoforms, its full native records incl. chimeric ones, and 2 EIF3C + 1
  myosin-11-like transcripts that carry territory exons; 91 transcripts) and replaced by the representative alone, so
  at f = 1.0 the copy is expressed only through the fusion. Unfused copies keep all their transcripts. Pool: 545
  transcripts, 16,350 reads.

*C. §3 measurement definitions (binding, fixed before any result).*
- Holder, NPIP family and placement: `npf_audit/audit.py run()` unchanged via a config override, as `ggo_npip_sim`
  `score.py`. Unit of a locus = itself if a clusters.tsv member, else its fold target.
- **M1** = fused copies with placement NPIP-family.
- **M2** = S partner genes (all 20) whose HOLDER (the assembled locus with the most same-strand exon bp over the
  partner's exon union; ties by reads) has its unit in the NPIP family. Also reported, not barred: the count over the
  10 fused partners, and the any-overlap variant (any NPIP-family locus sharing a same-strand exon base with the
  partner).
- **Fused member** (M3/M4) = the clustered unit whose transcripts carry the fusion junction (an intron within 10 bp of
  the recorded one at both ends; the most-read carrier if several); if none, the copy's holder unit if clustered; if
  neither, the copy is listed and left out of M3 (it is an M1 failure).
- **M3** truth by construction: P_fus = the partner bases the fusion carries, P_all = the partner's all-transcript exon
  union. Per fused member, A = its accessory bases (container `blocks.tsv`); precision = |A n P_all| / |A| (every
  partner base is a correct accessory call), recall = |A n P_fus| / |P_fus| (what the construction puts in the member).
  Pooled = sums over fused members. Also reported, not barred: the literal single-set variants (P_all/P_all and
  P_fus/P_fus) and recall with the excluded copies counted as 0. Partner bases inside core blocks are listed with
  their `core_partners`.
- **M4** = fused members whose accessory blocks relate (`rel_families`) to the family of the partner's holder unit;
  reported with the reasons for "no" (partner holder unclustered / is the fused member itself / clustered but not
  related).
- **M5** = accessory bp and blocks on the holder units of the 15 unfused copies (every arm) and of the 10 fused copies
  at f = 0, each block annotated with the genes it overlaps.
- Human descriptive step: A119b chr16 regenerated with `fj_bin_frozen` (copy_assign --assemble-only --region
  chr16:0-96330374, strict + shipped polish, genome-wide best-AS seeding table; then the driver's mcl_families
  command), since the 09-26 products no longer have a PAF.

## Outcome (2026-09-28 19:50; gorilla simulation, scored with the Amendment 2 definitions, no deviation after it)

Full report: scratchpad `figs/fusion_container_sim.md`; data `/mnt/linuxdisk/tmp/rustle_figures_dev/fusion_container_sim/`
(`score/score.json` sha1 44459a90, `score.py` 8f512fb8). Mapping: all 16,350 pool reads primary; every fusion read's
primary alignment carries the exact fusion junction (300/300); every fusion is assembled with its exact intron chain in
every arm that has fusion reads (10/10, from 3 reads up). Unfused copies are stable: 14/15 in NPIP in every arm.

| arm f | copies in NPIP /25 | **M1** fused in NPIP /10 | **M2** S partners in NPIP /20 | **M3** prec | **M3** rec | M4 | M5 unfused acc bp |
|---|---|---|---|---|---|---|---|
| 0 | 23 | **9** | **0** | (0 accessory partner bp) | 0 | 0 | 10,712 |
| 0.1 | 15 | **1** | **1** | 0.317 | 0.533 (9 members) | 0 | 11,908 |
| 0.5 | 14 | **0** | 0 | **0.214** | **0.417** (9) | 0 | 11,908 |
| 0.9 | 14 | 0 | 0 | **0.214** | **0.417** (9) | 0 | 11,908 |
| 1.0 | 18 | 4 | **4** | **0.389** | **0.621** (10) | 0 | 11,386 |

- **M1 FAILS** (f = 0.1: 1 vs 9; f = 0.5: 0 vs 9). **M2 FAILS** (f = 0.1: 1 > 0; f = 1.0: 4 > 0). **M3 FAILS** at
  every f >= 0.5 (bar 0.90 / 0.90). M4 = 0 in every arm.
- **Verdict (§4): Container FAILS** — through M1/M2, i.e. the result describes the FAMILY DEFINITION under fusion, not
  the container (which runs after the families and never changes one).
- Mechanism (read from the products, not a bar): (1) from 3 fusion reads up the shipped assembler emits the copy, the
  fusion and the partner as ONE gene_id (partner "absorbed": 9/9 members at f = 0.1-0.9, 10/10 at f = 1.0), so no
  partner locus of its own exists (hence M4 = 0 by construction). (2) The families stage places that locus by its
  representative (most reads, then span): at f = 0.1-0.9 the partner's own transcript (30 reads) beats the copy (<= 27)
  and the fusion (<= 27), so the fused copy moves into the PARTNER's family (SMG1-like x4, SNX29, PDXDC, CNOT3-like,
  NSMCE-like, SAGA29-like families) or becomes unclustered (ATXN2L); at f = 1.0 the fusion ties at 30 reads and wins on
  span, so 4 loci return to NPIP carrying their partner (M2 = 4). A strict partition forces this choice (§6s9).
- The container then describes the family it was given: outside NPIP the NPIP part is the accessory piece and 6 of 9
  such members (f = 0.5) carry a relation back to the NPIP family; co-duplicated partner pieces align to paralogous
  pieces in other members and are called core (18.5-21.4 kb of partner bases per arm, listed in the report).
  Post hoc, not barred: the 4 fused members that stay in NPIP at f = 1.0 give pooled precision 0.734 / recall 0.869;
  the two with a non-co-duplicated partner (NPIPB8~SAGA29-like, NPIPA9~ATXN2L) are exact (1.0 / 1.0).
- M5: unfused accessory = the extra pieces of annotated chimeric native records (NPIPB15's 63-exon sortilin-related
  receptor-like model, 10,105 bp; LOC124907807's 77-exon model, 281 bp), a PDXDC-like exon at NPIPB14P (125 bp), and
  NPIPB10P's own exons (102-1,298 bp) because NPIPB10P sits outside NPIP (as in `ggo_npip_sim`); fused@f0: NPIPA7
  132 bp (one exon with no aligned partner, NPIPA7 already outside NPIP at f = 0).
- **Human descriptive (dev, A119b chr16, frozen binaries; clusters identical to the 09-26 run):** NPIPA1, NPIPA6 and
  NPIPA9 have EMPTY containers (21 / 33 / 24 blocks, all core): their PKD1P3 / PKD1P1 / PKD1P5 readthrough pieces align
  to one another's PKD1P pieces inside the NPIP family. NPIPB14P's container holds 19 blocks, 5,178 bp = the PDXDC2P half
  of PDXDC2P-NPIPB14P (chr16:75,820,192-75,876,271), related to no other family; its 9 core blocks (8,312 bp) are the
  NPIPB14P half.
