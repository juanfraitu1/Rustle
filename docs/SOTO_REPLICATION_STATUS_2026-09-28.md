# Soto 2025 replication: status, gaps, and the "our families are different from theirs" evidence

**Written 2026-09-28 to consolidate that day's work and let the advisor conversation resume from a fixed point;
§1, §2 and §6 rewritten 2026-09-29 after the reconciliation** (`docs/PREREG_soto_reconciliation_2026-09-29.md`,
`PREREG_soto_famcn_allwssd_2026-09-29.md`, `PREREG_soto_parcn_assembly_2026-09-29.md`; register rows 1158-1174).
Context: the advisor doubts the pipeline because we have not reproduced Soto et al. 2025 (Cell) exactly. This
file separates three things that were previously tangled: (1) how close the replication now gets on their own
inputs, (2) exactly which parts of their pipeline we still cannot rebuild and why, and (3) direct, measured
evidence that "families" in Soto's sense and in ours are related but not the same object — so an exact match
was never the right bar.

Full detail lives in the cited reports and memory files; this file is the map, not a replacement for them.

## 1. Where the replication stands today (2026-09-29)

Chain: `bench/soto/soto_replication.py` (`genesets` → `edges --exon-mapback` → `cluster --pair-mad` → `score`;
recipe and expected numbers in `REPRODUCE.md` §5a), on Soto's own genome (CHM13 v1.0), CAT v4 annotation,
2,334-gene universe and their published famCN (Table S1C, 491 multigene families). Selection was on a frozen DEV
half (225 families) and judged on the HELD-OUT half (266; `bench/soto/soto_split_2026-09-29.tsv`).

**Headline.** ~~ARI 0.7096, 235 of 491 exact families (the literal recipe, 09-28)~~ → **ARI 0.9698, 479 of 491
exact (97.6%), held-out 0.9681 / 263 of 266; pair precision 1.000, recall 0.942; bipartite MICRO 1.000 / 0.980;
0 undetected families; 36 of the 38 flagship families (all 9 zebrafish-modelled pHSD families, all six NPIP
families, TBC1D3, NOTCH2NL, FAM90A, USP17L, DUX4L) exact.** The difference from the old headline is **two named
choices of Soto's released code that their STAR Methods prose does not state**:

1. **Map SD98 exons back to the genome, not SD98 regions** (`A_SD98_regions.md` l.226-238: every CAT v4 exon fully
   inside a merged autosomal SD98 region, `minimap2 -c --end-bonus 5 --eqx -N 50 -p 0.5`, shared exon = same-strand
   ≥ 99% cover; 25,971 exons of 5,154 genes → 12,231 edges). Alone: +0.126 ARI / +103 exact from the literal recipe.
2. **Apply MAD < 1 to each shared-exon pair, then grow families through coding genes** (`B_SD98_families.ipynb`
   cells 4-11), not to connected components; a non-coding gene joins every family it pairs with, which is what
   makes S1C a cover. Alone: +0.168 / +78; both together +0.260 / +244. Their loop as released is order-dependent
   (`gene_cluster = list(set(...))` while scanning: 0.815-0.872 over 8 hash seeds); what reproduces S1C is the
   full closure it evidently intends — a repair, disclosed (register 1165).

Nothing else moves it: the MAD statistic (median vs mean) is invariant under the pair gate; SEDEF's own CIGARs, the
union with them, the B1 self-overlap fix, a threshold-free segmentation split and median∧mean add nothing once the
two choices are Soto's (best such cell held-out 0.889); minimap2 2.17 (theirs) vs 2.30 costs −0.001 / −9 exact; our
reconstructed curation rule instead of their hand-curated list −0.021 / −23 exact (register 1163).

### 1.1 The ladder: how much of the 0.97 is their published copy number

Same edges, same rule, different copy numbers (`soto_replication.py ladder`; register 1169 / 1170). Exact families
and the ARI without FAM90A are quoted beside every ARI because on this recipe the ALL-491 ARI swings ±0.035 (held-out
±0.07) on whether that one 56-gene family (ID_356, famCN 31-42) stays whole — two of its 1,533 internal pairs sit at
|ΔfamCN| = 2, the gate's edge (register 1170):

| copy numbers fed to the pair gate | ARI all 491 (DEV / HELD-OUT) | exact | ARI without FAM90A |
|---|---|---|---|
| none (sequence only) | 0.7307 (.6418 / .8693) | 345 | 0.7057 |
| our WSSD famCN, 10 SGDP samples, merged exons (`famcn_ours_all.tsv`, 08-01) | 0.9198 (.9096 / .9317) | 373 | 0.9131 |
| our WSSD famCN, all 268 SGDP samples, merged exons | 0.8855 (.9039 / .8610) | 375 | 0.9089 |
| **our WSSD famCN, 268 samples, Soto's interval (gene body ∩ SD98)** | **0.9277 (.9227 / .9343)** | **411** | 0.9251 |
| S1C famCN (Soto's published values) | 0.9698 (.9708 / .9681) | 479 | 0.9650 |

Reading: sequence alone reaches 0.73; our own copy numbers reach **0.93 / 411** once the interval is Soto's (their
gene body ∩ SD98 pieces, not all exons: +36 exact families, Pearson with S1C 0.932 → 0.977; register 1169); the
sample count is not the lever (the 10-sample table was a favourable draw that happens to keep FAM90A whole; 20
random 10-sample draws span 0.870-0.923; register 1168 / 1170). **The last ~0.04 ARI / ~68 families are agreement
with S1C's exact values** (their per-row aggregator `genotype_cn_parallel.py` is unreleased; the GFF3 gene feature
vs the transcript span; whether their median includes the outlier sample they drop in code) — unmeasured, not
searched.

**Advisor-safe sentence** (from the reconciliation's independent verification, with the 09-29 correction of its
"0.92" to the 268-sample figure): *"Using Soto et al.'s released code choices (exon map-back, a per-gene-pair
copy-number test, and their family loop completed as it evidently intends) on their genome, annotation, gene list
and published copy numbers reproduces 479 of 491 families (ARI 0.97, 0.97 held-out); sequence alone reaches 0.73
and our own copy numbers 0.93 (411 exact), so the step to 0.97 is agreement with their published copy numbers, not
an independent reconstruction."* It remains **concordance with Soto's own tables** (register 858 / 1085): the gene
universe, the 71 curated genes and, on the last rung, famCN are theirs.

### 1.2 The 09-28 "walls", re-measured — three of four were artefacts of our own readings (register 1164)

The 09-28 edition of this section said *"cannot match by construction = their cover (409 cap), 12 families that
break their own MAD rule, and 50 families with members no ≥ 98% alignment reaches"* plus an acrocentric assembly
wall. Each is now wrong or gone; the old claim is struck through:

| old claim (09-28) | what is true (09-29) |
|---|---|
| ~~"any partition scored against their full table caps at 409/491 exact because their families are a cover"~~ | A property of partition methods only. Soto's own pair-level rule MAKES the cover: EXON × PAIR puts 158 genes in ≥ 2 families, including all 149 of S1C's, and reproduces **481 of 491 families exactly with their cover genes included** (above the "cap"). |
| ~~"12 Soto families violate their own MAD < 1 rule (0/12 can be emitted by that rule)"~~ | Their code never gates a family's MAD — the notebook only reports it (cell 14); the gate is per pair. 9 of the 12 are exact under EXON × PAIR (ID_28, 63, 113 remain, in the per-pair famCN residual). |
| ~~"~50 families with a member no ≥ 98% alignment reaches (unreachable under any reading)"~~ | An artefact of the REGION reading: with region-level edges 73 families are cut at the alignment stage, with exon-level edges **2** — exactly S1C's two hand-merged families (ID_347 DUX4 / ID_482 UBTFL, `Family MAD = "Manual merge"`). |
| ~~"a genuine assembly-version wall in the acrocentric chromosomes"~~ (the v2.0 → v1.0 liftover chain) | 0 on the native v1.0 chain: no liftover, no wall. It was a property of lifting a v2.0 SEDEF, not of the replication. |

What remains non-exact (12 families): 7 the per-pair gate on S1C's one-value-per-gene famCN cannot join (ID_28, 63,
99, 113, 144, 163, 270), the 2 manual merges, 3 scoring-collapse / singleton cases (ID_62, 192, 401). The 09-28 diff
table below is kept for the record of the literal recipe; its causes no longer describe the headline.

<details><summary>09-28 literal-recipe numbers (superseded; kept for the record)</summary>

| input | ARI (median-MAD / mean-MAD) | exact families | pair P/R/F1 | bipartite MICRO P/R |
|---|---|---|---|---|
| v2.0 SEDEF + liftover (ledger §6ip) | 0.6959 / 0.6862 | 241/491 / 264/491 | .841/.595/.697 | .784/.709 |
| native v1.0 SEDEF with CIGARs (`final_v1.bed`) | 0.6985 / 0.6894 | 241/491 / 264/491 | .837/.601/.700 | .780/.715 |
| + reconstructed curation rule | 0.6983 / 0.6886 | 240/491 / 263/491 | .835/.602/.700 | .781/.715 |
| Soto's literal recipe on native v1.0 (regions mapped back, exons projected, `-f 0.99`) | 0.7096 / 0.7039 | 235/491 / 261/491 | .831/.621/.711 | .784/.718 |

Primary causes of the 256 non-exact families under that recipe (median arm): 127 over-merges the component split
could not undo (~700 real ≥ 98% cross-family edges, register 1160), 105 missing edges (62 fragments at 90-98%
identity, 49 "SEDEF row but no map-back hit" — not a `-N 50` loss, register 1161 — 17 via non-eligible members),
24 copy-number split errors; the cover was not causal. The B1 self-overlap fix was not distinguishable from its
null (register 1158); the minimum-segmentation split failed its no-harm clause on one family (register 1159).
</details>

## 2. What we still cannot rebuild, and why (five named gaps, 09-29 state)

1. **Manual curation (their step 2).** Their STAR Methods self-intersect SD98 transcripts, flag pairs with > 90%
   positional overlap, then *manually* drop redundant / readthrough-fusion transcripts; 71 genes leave between our
   reconstructed SD98 set (1,864 eligible) and their 1,793. `soto_replication.py curate` carries a one-condition
   rule in their vocabulary (drop the lower-ranked side of a same-strand overlapping SD98 gene pair) at nested
   held-out P/R ≈ 0.77 / 0.76 — but its family-score effect is not distinguishable from dropping the same number of
   genes at random (`soto_cur_critique.md`), and at the reconciled cell **their list vs our rule costs −0.021 ARI /
   −23 exact** (register 1163). We use their list (S1C `In Table S1 = Yes`); this is the one step that stays hand-made.
2. **famCN — now recomputed, and the gap is quantified.** `famcn --interval sd98 --samples all` computes our own
   WSSD famCN over Soto's interval from all 268 SGDP tracks (their outlier sample dropped, as their code does;
   `famcn_ours_allwssd.tsv`, 5,154 genes). It reaches **0.93 / 411** in the recipe (§1.1); S1C's values reach
   0.97 / 479. The remaining ~0.04 is their unreleased per-row aggregator, the GFF3 gene feature and the outlier
   sample — named, unmeasured. The 09-28 note that a WSSD recomputation "produces real false merges" (§6in) was about
   using famCN as a weak-edge lever in OUR catalog, a different use; in Soto's own recipe it is the second-best rung.
3. **parCN — QuicK-mer2 cannot run here; exact assembly counting replaces it** (`bench/soto/parcn_assembly.py`;
   register 1167, 1171-1174). QuicK-mer2's whole-genome index needs a 2^32-slot hash: `search` peaks at ~52 GB
   and `count` at ~43 GB against 25 GiB + 16 GB swap (`figs/soto_quickmer2.md`). The same k-mer rule (k = 30
   canonical, once in CHM13-noY, Hamming-1 edit depth < 100) counted exactly in complete assemblies reproduces
   S1E's population parCN on genes whose CN is fixed in humans: **HG002 within 0.5 of S1E for 321/322 Fixed genes
   (0.997; controls 297/299 = 2)**, but only 408/629 Nearly-Fixed (0.649, bar 0.70 failed: 22% of them carry < 75%
   of CHM13's "paralog-specific" k-mers — partly haplotype-specific) and 83/212 Polymorphic; Spearman 0.414. Their
   WSSD human-vs-ape calls are reproduced by the assembly famCN analogue: Duplicated in humans 105/109 (0.963),
   Expanded 105/132 (0.795), CN-called families 105/118 (0.890), non-calls 131/631 (0.208, mostly ape-individual
   differences; no bonobo); Soto's famCN > 10 exclusion hides ~234 paralogs the assemblies call human gains
   (subtelomeric DDX11L / WASH / OR4F first). Human and ape numbers are per genome, never pooled.
4. **DupMasker intersection (their final annotation step) has not been run.** It labels families with ancestral
   duplication units after clustering; it does not change membership — a completeness gap in the write-up, not a
   scoring gap.
5. ~~**A genuine assembly-version wall in the acrocentric chromosomes.**~~ Gone on the native v1.0 chain (§1.2): it
   was a liftover artefact of the v2.0 SEDEF input, "not pursued further" for the wrong reason.

**Everything else** in their stated method — the SD98 threshold, gene-to-SD98-region assignment (5,154 autosomal
genes = their count), the exon map-back, the `-f 0.99` shared-exon call, the per-pair MAD gate and the coding-gene
closure — is reproduced from their released code (`github.com/mydennislab/HSD_brain_evolution`) and checked
family for family against S1C. Two code audits are on record: the released closure loop is order-dependent
(register 1165), and the two prose-vs-code differences (regions vs exons; component MAD vs per-pair MAD) carry the
whole gap between the literal recipe and S1C (register 1166). Our 09-11 `dennislab` port read their "clusters" as
connected components — a misreading that scored 0.50-0.56 and is kept only as a diagnostic row.

## 3. Direct evidence that Soto's "families" and ours are not the same object

This is the answer to "prove their families are different," with citations to where each number was measured.

1. **Nesting, measured, not argued** (`bench/SOTO_AS_A_REFINEMENT.md`, ledger §6n7). For every Soto family with
   ≥2 genes present in our node set, we counted how many of *our* groups its members span.

   | level | Soto families (≥2 genes) | contained in ONE of ours | rate | NPIP-side (6 families) |
   |---|---|---|---|---|
   | L1a/L1b | 10 | 9 | 90.0% | 6/6 |
   | L2 | 10 | 8 | 80.0% | 6/6 |
   | L3 | 10 | 2 | 20.0% | 0/6 |

   **All six of Soto's NPIP-side families (ID_149, ID_151–ID_155) sit entirely inside our single NPIP family
   at L1/L2.** At L3 (our finer level) we split them instead — Soto's granularity sits *between* our L2 and
   L3, not orthogonal to either.
   - **Held out, never used to develop this:** chr5/7/21 MCL catalogs give 80.3%/76.4% containment (E1/E1S,
     76 and 72 Soto families) — the same ~80% band, independently.
   - **Where it breaks is a size effect, not randomness:** Soto families of ≤3 genes nest at 89.8% (53/59);
     families of >3 genes nest at only 47.1% (8/17). The disagreement concentrates in their large families.

2. **Soto's own truth is a cover, not a partition** (`SOTO_AS_A_REFINEMENT.md` header fact): 148 of their 2,333
   genes (6.4%) belong to two or more of their own published families (104 to two, 25 to three, 6 to four, 6
   to five, 7 to six). A method that emits a strict partition — ours does, by construction (§6s9) — cannot be
   graded against a non-partition truth without stating this; several early scorers got this wrong before it
   was caught.

3. **Where we disagree with Soto, the disagreement itself is measured, not assumed** (`project_cover_and_jn_
   refuted.md`, r945/r946, 118 genes both truths place). Restricted to pairs one truth calls together and the
   other does not:
   - **We do not over-merge relative to Soto:** 0 of their 343 pairs are rejected by an independent referee,
     and 0 of 24 of their families get split by our edges.
   - **We under-merge relative to Soto, and about a third of that is real:** of 247 pairs Soto keeps apart
     that we would join, 164 (66.4%) have no independent sequence-identity edge at all (ancient paralogy Soto
     is right to exclude) — but **83 carry a direct edge at median identity 0.858** (37 at ≥0.90, 30 at ≥0.95),
     against a shuffled-pair control at 0.962 (i.e. these are not noise). These are real paralog pairs Soto's
     own pipeline misses — **not via the SD98/curation gaps in §2** (superseded: attributed and independently
     re-verified below).
   - **Why Soto misses these 83, measured per-pair (attribution + independent re-check, 2026-09-28):** all 35
     genes behind the 83 pairs pass every screen in §2 (100% `in_sd98=1` on both the native-v1.0 and v2.0 SEDEF
     regions with a fully-contained exon, 100% inside the 1,793-gene eligible universe, 100% `KEPT` by the
     reconstructed curation rule) — identity floor, exon containment, and curation/universe together explain
     **0/83 (0%)**. The real causes are downstream of curation: **65/83 (78.3%)** have no shared-exon edge
     between the pair in either frozen edge file (`shared_exons_2334_finalhuman.tsv` / `..._finalv1_native.tsv`)
     — 44 chr15 (32 GOLGA6↔GOLGA8 cross-subfamily + 12 where `GOLGA8A` specifically has no edge to *any* other
     GOLGA8 paralog either, so "44/44 GOLGA6↔GOLGA8" is not quite right — 12 of the 44 are a GOLGA8A-specific
     isolation, not a clean subfamily boundary) plus 21 chr16 NPIP A-clade/B-clade pairs; and the remaining
     **18/83 (21.7%)**, all chr16 NPIP, do have an edge but land in different S1C Family IDs, split by Soto's
     own sequential famCN-MAD grouping at a modest famCN gap (e.g. `NPIPA1`↔`NPIPA7`: famCN 9.05 vs 47.58).
     Category "same-family artifact" is 0/83 by construction. **Verified independently**: a seeded random
     sample of 20/83 pairs, re-derived from raw inputs with freshly written code (not reusing the attribution
     script), matched the original per-pair calls 20/20 (0% mismatch); a full independent re-classification of
     all 83 pairs reproduced the exact 65/18/0/0/0/0 breakdown with 0 mismatches, and the first-applicable-
     cause rule is mutually exclusive by construction and sums to 83 in both scripts. Detail:
     `scratchpad/soto_attr/attributed_pairs.tsv` (per-pair), `scratchpad/figs/soto_um_attribution.md` (report),
     `scratchpad/figs/soto_um_verify.md` (independent re-check).

**The composite statement for the advisor** (09-28 wording; its "ARI≈0.70 / 49%" and gap decomposition are
superseded by §1 — use §1.1's sentence for the replication half, this paragraph for the "different object" half):
*"Our pipeline reproduces Soto's own numbers to ARI≈0.70 / 49% exact families using only their inputs. The
remaining gap is not evidence the pipeline fails to find real structure — it decomposes into (a) an unreproducible
manual curation step we've now characterized to ~0.77 precision/recall, (b) one non-algorithmic assembly-version
wall, and (c) a genuine granularity difference: Soto's families are measurably a coarser, non-partition cover that
nests inside ours 80–90% of the time (100% on the NPIP case the advisor cares about) except at their largest, most
heterogeneous families, where we can independently show 83 real paralog pairs their own truth misses."*

## 4. What was corrected today (say this too, don't let it surface as a gotcha)

- The earlier memory claim "71/71 flagged by their own command" was wrong: their literal self-intersect
  command flags 117 of 1,864 candidate genes (52 real matches, 65 false positives against the 71) — precision
  0.44, not a clean match. Corrected in `soto_cur_critique.md`.
- A "fully native chain, no S1C gene list" framing was found still reading famCN from S1C for one arm
  (`curate_F_c03.log`); corrected — nothing in the native-v1.0 chain is famCN-independent yet (see §2 item 2).
- The curation rule's originally reported held-out P/R (0.86/0.72, later 0.78/0.80) is optimistic: the rank
  order, tie-breaks and 0.3 threshold were chosen after seeing all 71 answers before the held-out split was
  run. A fully nested refit (rank *and* threshold both re-chosen per half) gives the more honest 0.766/0.760 —
  use this number, not the higher one, if quoted again.

## 5. Where the work products live (to resume from)

- Session reports (scratchpad, **not durable** — copy anything still needed into `bench/soto/` or here before
  the session ends): `soto_v1_cigar.md`, `soto_v1_literal.md`, `soto_cur_label.md`, `soto_cur_rule.md`,
  `soto_cur_critique.md`.
- Durable: `bench/soto/soto_replication.py` (now has `--native-v1` and the `curate` subcommand),
  `bench/SOTO_AS_A_REFINEMENT.md`, memory `project_soto_full_replication.md`, `project_cover_and_jn_refuted.md`,
  `project_soto_family_pseudogene_fragment_audit.md`, `reference_soto_2025_hsd_brain.md`.
- Data: `/mnt/c/Users/jfris/Desktop/final_v1.bed` (native v1.0 SEDEF with CIGARs — user-supplied today; not yet
  copied anywhere durable — consider moving a copy under `winloci_data/soto_replication/` before it's lost).
- Labelled gene table: `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_cur_label/genes_labelled.tsv` (scratch disk,
  regenerable from `final_v1.bed` + `cat_v4.bed` via `soto_replication.py genesets`/`curate`).

## 6. Where things stand and what is left (2026-09-29)

**Durable now (reproducible from `bench/`; `REPRODUCE.md` §5a):** `soto_replication.py edges --exon-mapback`,
`cluster --pair-mad`, `score --split/--half/--drop-family`, `famcn --interval exons|sd98 --samples all` and
`ladder`; the frozen edge table `bench/soto/shared_exons_5154_exon_mapback.tsv` and split
`bench/soto/soto_split_2026-09-29.tsv`; `bench/soto/parcn_assembly.py` (+ `test_parcn_assembly.py`) for the
assembly parCN. Each was re-run once on 09-29 and reproduced its frozen product byte for byte (edges, exon BED /
FASTA, the reconcile partition, `famcn_ours_allwssd.tsv`, the 269-sample matrix, all 113 parCN summary values).
The literal 09-28 chain is unchanged and byte-identical without the new flags.

1. **Nothing is left to search in the family recipe.** Both DEV / HELD-OUT halves are spent for edge sources,
   famCN splits and the curation rule (register 1158-1163). Do not re-open the MAD statistic, B1, SEDEF-CIGAR,
   DP-split or curation arms; their numbers are in the reconciliation prereg's `cells.tsv`.
2. **famCN independence, if the advisor wants it closed further:** the three named, unmeasured causes of the
   0.93 → 0.97 step (§2 item 2). The first would need Soto's `genotype_cn_parallel.py`, which is not public.
3. **parCN:** QuicK-mer2 itself needs a ≥ 64 GB node or a hash-partitioned re-implementation validated
   bit-for-bit on one chromosome (register 1167); the assembly route is what we have. Its open finding is the
   Nearly-Fixed class (0.649): "paralog-specific" k-mers of one haploid reference are partly haplotype-specific,
   so a single individual's exact count disagrees with a population median — and the same haplotype effect makes a
   k-mer presence gate for O2 worse than random (register 1180-1182). Judge copy presence by synteny, not k-mer
   survival.
4. **DupMasker** (item 4) only if the downstream human-specific-expansion labels are asked for by name.
5. **The nesting result (§3.1) is still the single most advisor-facing artifact** — `SOTO_AS_A_REFINEMENT.md`'s table
   as one slide, now beside the ladder of §1.1 (sequence 0.73 → our CN 0.93 → their CN 0.97).
