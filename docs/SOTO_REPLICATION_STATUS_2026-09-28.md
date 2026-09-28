# Soto 2025 replication: status, gaps, and the "our families are different from theirs" evidence

**Written 2026-09-28 to consolidate today's work and let the advisor conversation resume from a fixed point.**
Context: the advisor doubts the pipeline because we have not reproduced Soto et al. 2025 (Cell) exactly. This
file separates three things that were previously tangled: (1) how close the replication now gets on their own
inputs, (2) exactly which parts of their pipeline we still cannot rebuild and why, and (3) direct, measured
evidence that "families" in Soto's sense and in ours are related but not the same object — so an exact match
was never the right bar.

Full detail lives in the cited reports and memory files; this file is the map, not a replacement for them.

## 1. Where the replication stands today

Chain: `bench/soto/soto_replication.py` (`genesets` → `edges` → `cluster` → `score`), on Soto's own 2,334-gene
universe and their published famCN truth (Table S1C, 491 multigene families).

| input | ARI (median-MAD / mean-MAD) | exact families | pair P/R/F1 | bipartite MICRO P/R |
|---|---|---|---|---|
| v2.0 SEDEF + liftover (previous best, ledger §6ip) | 0.6959 / 0.6862 | 241/491 (49.1%) / 264/491 (53.8%) | .841/.595/.697 | .784/.709 |
| **native v1.0 SEDEF with CIGARs** (`final_v1.bed`, today) | **0.6985 / 0.6894** | 241/491 / 264/491 (unchanged) | .837/**.601**/**.700** | .780/**.715** |
| + reconstructed curation rule (below) | 0.6983 / 0.6886 | 240/491 / 263/491 | .835/.602/.700 | .781/.715 |
| **Soto's literal recipe** on native v1.0 (getfasta → minimap2 `-c --end-bonus 5 --eqx -N50 -p0.5` → exons projected, `-f 0.99`) | **0.7096 / 0.7039** | 235/491 (47.9%) / 261/491 (53.2%) | .831/**.621**/**.711** | .784/.718 |
| literal ∪ CIGAR-walk edges (diagnostic) | 0.7098 / 0.7053 | **243/491 / 269/491** | — | undetected 100 / 89 |

Full numbers: `soto_v1_cigar.md`, `soto_v1_literal.md`, `soto_cur_rule.md` (session scratchpad, paths in §5).

**Reading:** switching to the real, native CHM13 v1.0 SEDEF calls with their actual CIGAR strings (rather than
a lifted-over v2.0 SEDEF or a v1.0 SEDEF file without CIGARs) closes the residual gap from missing/liftover-
dropped rows and gives a small, consistent recall gain (+.006 pair recall, +.006/.007 bipartite recall) at a
small precision cost. The exact-family set does not change. **This is concordance with Soto's own numbers on
their own inputs, not an independent replication** — the famCN truth we score against is theirs (register
rows 858, 1085 flag this explicitly).

### 1.1 The literal recipe and the family-by-family diff (2026-09-28, `soto_v1_literal.md`, `soto_v1_diff.md`)

- **Their recipe, run literally on their genome, is our best replication: ARI 0.7096 (median-MAD).** It finds
  11,935 edges (truth precision .879) against the CIGAR walk's 4,384; 81% of the extra edges are ≥98%-identity links
  SEDEF never reported as pairs.
- **Which SD catalog seeds it does not matter:** the published Vollger v1.0 SD set gives a byte-identical edge set and
  byte-identical families. The SD catalog only chooses what gets mapped back; the map-back makes the edges.
- **Their methods sentence read literally (region-level intersect) gives ARI 0.19.** Whatever they ran kept
  exon-to-exon correspondence, so the exon-projection reading is ours and is disclosed as such.
- **More edges recover more genes but fewer whole families** (235 exact vs 241): about 700 real ≥98% edges join
  different Soto families that the copy-number split then cannot separate again.

**Why the remaining families are not exact** (median arm, 256 non-exact, primary cause):
| cause | families | what it is |
|---|---|---|
| over-merge we cannot split | 127 | the joined true family is a minority of its predicted group in 59; a median-MAD cannot see a sub-half population |
| edges we lack | 105 | 62 fragments whose closest SD row is 90–98% identity with no map-back hit (unreachable under any reading); 49 where SEDEF has a ≥98% row but minimap2 finds nothing; 17 linked only through non-eligible members |
| copy-number split | 24 | 12 Soto families violate their own MAD < 1 rule (0/12 can be emitted by that rule); 20 avoidable cuts by our greedy split |
| their cover / universe | 0 | not causal: removing the 149 multi-family genes changes nothing |

**Ceilings:** a perfect split of our own components would reach ARI 0.946 / 420 exact; any partition scored against
their full table caps at 409/491 exact because their families are a cover, and two readings of their own table agree
only at ARI 0.94.

**Advisor sentence:** *replicated* = the alignment-and-graph half of their recipe on their genome, annotation and
universe (ARI 0.71, 48–53% exact); *concordance only* = everything that consumes their tables (famCN split, universe,
the 71 curated genes); *cannot match by construction* = their cover (409 cap), 12 families that break their own MAD
rule, and 50 families with members no ≥98% alignment reaches.

## 2. What we still cannot rebuild, and why (five specific gaps)

Each is a *named* gap, not a general "it doesn't work" — this is the list to hand the advisor directly.

1. **Manual curation (their step 2).** Their STAR Methods say they self-intersect SD98 transcripts, flag
   pairs with >90% positional overlap, and then *manually curate* to drop redundant/readthrough-fusion
   transcripts. 71 genes are removed this way between our reconstructed SD98 set and their published 1,793.
   - Today's work (`soto_cur_label.md`, `soto_cur_rule.md`, `soto_cur_critique.md`) found a one-condition rule
     in their own vocabulary — drop the lower-ranked side of a same-strand overlapping SD98 gene pair (rank =
     biotype > name-quality > transcript count > length) — that matches their 71 removals at **precision ≈0.77,
     recall ≈0.77** when the rank order and threshold are validated on a held-out half of the genes.
   - **Caveat, found by the same pass's hostile critique:** the family-score effect of this rule is *not*
     distinguishable from dropping the same number of genes at random from the same candidate pool (random
     subsets reach the same ARI about half the time). So the rule explains their curation step reasonably
     well; it does **not** explain, or improve on, the family-recovery numbers. Say both halves of this.
   - ~9 of the 71 are coin-flips between identical duplicate copies (e.g. DUX4L28 vs DUX4L29) and a handful
     are protein-level readthrough judgments (e.g. ISY1-RAB43) that no coordinate rule can see.
2. **famCN is not independently recomputed** — we use their published Table S1C values. famCN comes from WSSD
   (whole-genome shotgun read-depth of multi-mapping short reads across SGDP n=269 + 4 archaic/great-ape
   genomes). We built a real WSSD-based recomputation (`bench/soto/famcn_from_wssd.py`) and found it produces
   real false merges (§6in: 4/7 accepted weak-edge components were false merges across unrelated true
   families). **Not adopted; this stays concordance, not re-derivation**, and is the main reason full
   independence from their pipeline isn't yet claimed.
3. **parCN (QuicK-mer2 paralog-specific k-mer copy number) is not built.** It feeds their human-vs-ape
   expansion calls, not family clustering itself, so it hasn't blocked anything so far — but reproducing those
   downstream calls would need a fresh k-mer pipeline.
4. **DupMasker intersection (their final annotation step) has not been run.** It labels families with
   ancestral duplication units after clustering; it does not change membership, so it's a completeness gap in
   the replication write-up, not a scoring gap.
5. **A genuine assembly-version wall in the acrocentric chromosomes (chr13/14/15/21/22).** One small window in
   CHM13 v1.0 corresponds to many separate, tandem-periodic segments in v2.0 — real new assembly content, not
   a missing liftover anchor (checked directly: realigning wide windows around every gap edge re-confirms
   existing anchors with no headroom). Several v2.0 "extra" copies have no v1.0 counterpart to compare against
   at all, independent of any liftover fix. Ledger §6iq: "not pursued further — a real ceiling."

**Everything else** in their stated method — the SD98 threshold, gene-to-SD98-region assignment, the
minimap2 map-back, the `-f 0.99` shared-exon graph, the MAD-based family split (both the mean-vs-median-MAD
and the discard-then-bridge-merge discrepancies between their paper prose and their released code) — has been
reproduced from their own released code (`github.com/mydennislab/HSD_brain_evolution`,
`section_I&II/B_SD98_families.ipynb`) and checked byte-for-byte or family-for-family against their published
tables. See `project_soto_full_replication.md` for the full step-by-step trace.

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

**The composite statement for the advisor:** *"Our pipeline reproduces Soto's own numbers to ARI≈0.70 / 49%
exact families using only their inputs. The remaining gap is not evidence the pipeline fails to find real
structure — it decomposes into (a) an unreproducible manual curation step we've now characterized to ~0.77
precision/recall, (b) one non-algorithmic assembly-version wall, and (c) a genuine granularity difference:
Soto's families are measurably a coarser, non-partition cover that nests inside ours 80–90% of the time
(100% on the NPIP case the advisor cares about) except at their largest, most heterogeneous families, where
we can independently show 83 real paralog pairs their own truth misses."*

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

## 6. Next steps if resumed

1. Decide whether to pursue DupMasker (item 4) and parCN (item 3) — neither blocks the family-recovery
   argument above; both are needed only if the advisor asks about the downstream human-specific-expansion
   calls specifically.
2. Consider a WSSD recomputation attempt with a stronger corroboration requirement (item 2) if the advisor
   wants full independence from Soto's famCN — flagged as unsolved, not attempted further today (§6in's
   generalizable lesson: a weak signal correlating with a real property still needs a structural admission
   guarantee, not just corroboration from a second imperfect signal).
3. The nesting result (§3.1) is the single most advisor-facing artifact here — consider turning
   `SOTO_AS_A_REFINEMENT.md`'s table into a one-slide figure.
