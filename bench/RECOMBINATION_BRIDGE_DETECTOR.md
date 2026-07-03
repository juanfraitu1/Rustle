# Recombinant / exon-shuffled bridge genes: do they cause real over-merges, and does VG path-level colinearity split them?

**Scripts / outputs**
- Detector: `bench/recombination_bridge_detector.py` (build stage, deterministic `PYTHONHASHSEED=0`, RNA-only).
- Per-bridge table: `bench/recombination_bridge_detector.tsv` (90 data rows; 17 cols incl. new `bridge_kind`, `confidence`).
- Machine summary: `bench/recombination_bridge_detector.json` (byte-identical across runs; JSON md5 `ebde024a…`, TSV md5 `9481eb55…`).

---

## The question

The pairwise family-edge definition (E_r transcript-homology component → connected-components / γ-quasi-clique)
can be fooled by a **recombinant / exon-shuffled bridge gene**. Scenario from the theory:

- gene1 (exon1+exon2) at locus A, gene2 (exon1+exon2) at locus B — A and B share *nothing* directly.
- gene3 is a **mosaic**: its exon1 matches gene2's exon1 and its exon2 matches gene1's exon2.
- gene3 makes edges G1–G3 (shared exon2) and G2–G3 (shared exon1), so connected-components merges
  {G1,G2,G3} **transitively**, and a density γ-quasi-clique (density 2/3 on the 3-node barbell > γ=0.20)
  does **not** split it.

This is the **recombination obstruction** of the copy-assignment theory — the **K-frontier**: for K≥3 copies a
cross-copy recombinant path defeats any purely local/pairwise separator (see `bench/copy_assignment_theory.md`,
K-frontier / Prop). The **VG framing** should catch it: gene3 is a **recombinant path** whose exon-chain
traverses nodes from two otherwise-disjoint clusters, so a **path-level colinearity** criterion sees the chain
is a mosaic and does *not* colinearly match any single family.

**The distinction that matters** (reshuffling ≠ domain-sharing):
- **RESHUFFLING / recombinant:** *different* exons of the bridge match *different* families (bridge exon1 →
  family B, bridge exon2 → family A) — a **colinear mosaic**.
- **DOMAIN-SHARE:** the *same* exon/domain of the bridge matches *multiple* families (one shared domain
  bridges) — the GSTM2 case, **not** recombination.

---

## Method

For each over-merged / E_p-impure block (74 E_p-impure blocks from `family_level_pr_current.json` + the named
`fam17` hub, GSTM2, FOXO1, MAGE, MPHOSPH8, NLGN1 oversize):

1. Reconstruct per-locus **exon coordinates** from `denovo_skeletons.tsv` (exons = complement of introns) and
   extract each exon's genome sequence from `GGO.fasta` (same spliced-exon provenance as `aln_frac`).
2. Align every locus's exons to every other member's exons (**edlib** HW/infix, fwd+revcomp, **id≥0.70**, length
   prefilter ratio≥0.5, MIN_EXON=30 bp).
3. Build the intra-family **locus graph** (edge iff ≥1 shared exon at id≥0.70); find **bridge loci** = articulation
   points whose removal splits the family into ≥2 sub-families.
4. Map each bridge exon → the sub-family it best matches, and **classify**:
   - **RECOMBINANT** — different exons → different sub-families (`single_exon_subfams ≥ 2`, checked *before*),
   - **DOMAIN-SHARE** — the same exon → multiple sub-families (`any exon hits ≥2 sub-families`),
   - else NEITHER; no-bridge blocks → CARDINALITY (near-identical copies) or DISCONNECTED (sub-exonic bridge).
5. **VG path-level test:** best single-neighbour **colinear cover** (longest order-preserving matched-exon chain /
   n_exons via LIS). A colinearity criterion **splits** the block iff the bridge is RECOMBINANT **and** no single
   neighbour colinearly covers ≥ `COLINEAR_COV` (0.60). Recall cost estimated on pure single-gene multi-copy
   control families.
6. **New (this revision):** for each recombinant bridge, a **genomic co-location** check separates a true
   **DISTRIBUTED-MOSAIC** (sub-families Mb apart / different chromosomes) from a **READTHROUGH-CHIMERA**
   (sub-families on genomically *overlapping* adjacent loci = an assembly artifact, not reshuffling), plus a
   **confidence tier** (HIGH = distributed, both sub-families annotated, minor block ≥2 exons).

---

## Taxonomy of the 83 over-merged / E_p-impure blocks

| mechanism | blocks | share | what it is |
|---|---:|---:|---|
| **DISCONNECTED** — sub-exonic repeat/low-complexity/Alu bridge (merged members share **no full exon** at id≥0.70; cross-component best-exon-id **median 0.586**) | **34** | 41% | not recombination — a partial-sequence E_r leak |
| **DOMAIN-SHARE** (16 no-bridge + 4 bridge) — the **same** exon/domain bridges ≥2 families | **20** | 24% | one shared domain (GSTM2-style) |
| **CARDINALITY** — near-identical over-counted copies of ONE family (med exon-coverage ≥0.7) | **19** | 23% | recent tandem dup / MAGE-style count |
| **RECOMBINANT-BRIDGE** — **different** exons → **different** families (colinear mosaic) | **10** | 12% | reshuffling / K-frontier |

**Named placements** (all reproduced):
- `fam17` 20-locus hub (fid 17) and FOXO1 hub (fid 332) → **DISCONNECTED** (repeat bridge, cross-comp id 0.636 /
  0.573 — *not* recombination).
- GSTM2 hub (fid 9, med_cov 0.875 / med_id 0.995) → **CARDINALITY** (recent tandem dup; E_p splits the protein
  families but RNA correctly groups — *not* domain-share).
- MAGE (fid 549, 547) and MPHOSPH8 (fid 374) → **CARDINALITY**.
- The reshuffling / domain-share / cardinality distinction holds: **no** domain-share block is mis-called
  recombinant (a pure domain-share has `single_exon_subfams<2` because its non-domain exons have empty hits).

---

## Is reshuffling a REAL over-merge mechanism? YES — but it is the smallest bucket, and rare

**10 / 83** blocks are RECOMBINANT-BRIDGE. Of these, VG path-colinearity **splits 5** (blocks **94, 110, 187, 210,
321** = 7 bridges); the other 5 (37, 75, 288, 339, 443) are spared because one family colinearly covers ≥0.60.

The 5 split blocks are **not all** distributed exon-shuffles — the new co-location check breaks them down:

| VG-split block | bridge_kind | confidence | what it actually is |
|---|---|---|---|
| **fid 210** `LOC115935264` | DISTRIBUTED-MOSAIC | **HIGH** | **flagship** — the exact G1–G2–G3 barbell |
| **fid 187** `LOC115935376` | DISTRIBUTED-MOSAIC | **HIGH** | 3-exon internal block → a different family inside a 23-exon gene |
| fid 94 `LOC134757313` | DISTRIBUTED-MOSAIC | LOW | 2-exon fragment, minor exon id 0.799 (borderline) |
| fid 110 `LOC134758618`/`NA` | DISTRIBUTED-MOSAIC | LOW | both sub-families **unannotated (NA)** → distinctness unconfirmable |
| fid 321 `RDH14` | **READTHROUGH-CHIMERA** | LOW | NT5C1B / RDH14 loci **overlap by ~6 kb** → adjacent-gene readthrough, an assembly artifact, **not** reshuffling |

So the honest accounting of the VG-split set is **4 distributed exon-shuffles + 1 readthrough chimera**, and of the
distributed set only **2 are high-confidence** (fid 210, fid 187).

**Flagship — fid 210 (the textbook barbell).** On `NC_073230.2`: bridge `LOC115935264` (7 exons) at 81.09 Mb;
its exons 2–5 → `LOC101126070` (76.43 Mb, id 0.978–0.993), its exons 6–7 → `GALNT17` (80.58 Mb, id 0.983–0.994),
exon1 unmatched. The two flanks are **~4.1 Mb apart** and share **no exon directly** (best cross-flank exon-id
0.598 < 0.70 both directions) — the bridge is the **sole articulation**. The protein oracle independently confirms
two distinct families (**PRFAM2156** vs **PRFAM86**). Connected-components + γ-quasi-clique (density-2/3 barbell)
merge {LOC101126070, GALNT17, LOC115935264}; VG splits it (colinear cover 0.571 < 0.60). This is the
recombination obstruction realised on real gorilla loci.

---

## Would VG path-level colinearity split them where the pairwise-edge graph merges? YES — at zero recall cost

**Splits** the 5 recombinant blocks the pairwise-edge/γ-quasi-clique graph transitively merges — a fix **no density
or edge criterion catches**, because the barbell is dense (2/3) and connected. This is exactly the K-frontier the
theory predicts a *path-level* (not vertex/edge-local) criterion is needed to break.

**Recall cost = 0.** `0/24` pure single-gene multi-copy control families wrongly split at colinear≥0.60, and **0
across the entire 0.50–0.80 sweep**. There is a **clean valley**: the max colinear among **split** bridges =
**0.571**, the min among **spared** bridges = **0.75**; the 0.60 threshold sits in the gap.

**Threshold sweep (robustness of 0.60):**

| colinear thresh | blocks split | bridges split | control false-splits |
|---:|---:|---:|---:|
| 0.50 | 3 | 3 | 0 |
| 0.55 | 4 | 6 | 0 |
| **0.60–0.75** | **5** | **7** | **0** (robust plateau) |
| 0.80 | 10 | 12 | 0 |

Thresholds 0.60–0.75 give the **identical 5-block / 7-bridge split at zero recall cost**; the count only jumps to
10 at 0.80, where it starts eating legitimate alt-first-exon / terminal-exon variants.

**Does it spare real families with legit alternative first exons?** Yes — demonstrated directly by the **spared**
recombinant bridges: `fid 288` (DDX11 **alt-first-exon**, colinear 0.75), `fid 339`/`fid 443` (terminal-exon
variants, colinear 0.75–0.94) are all correctly spared at 0.60 and would only be lost above their own cover.

**Honest caveats (disclosed):**
- The pure control contains **0 alt-first-exon families** (`n_alt_first_exon_families=0`), so "0/24" is measured on
  a set with no at-risk members. The real recall-protection evidence is therefore the **spared alt-first / terminal
  bridges above + the 0.571/0.75 valley**, which is empirically adequate for *this* dataset but is dataset-specific,
  not a guarantee — carry this caveat with the recall number.
- **Honest miss — fid 75** (`LOC129526398` = KMT2C×6 exons + TMEM128×2 exons; distributed across 3 chromosomes;
  PRFAM2125 vs PRFAM457) is a **genuine high-confidence distributed mosaic** but colinear = 0.75 (KMT2C covers
  6/8 exons), so VG **spares** it at 0.60. Catching it needs threshold >0.75, which starts splitting the legit
  alt-first / terminal-exon families — so 0.60 is the correct **precision-preserving** choice, at the cost of
  missing this one mosaic.

---

## Connection to the theory (K-frontier / recombination obstruction)

The copy-assignment theory's K-frontier states that at **K≥3** a cross-copy **recombinant** defeats any purely
pairwise / vertex-local separator: a mosaic path stitches two otherwise-disjoint copy-clusters into one dense,
connected component. That is precisely what a recombinant bridge does to the family graph — connected-components and
the density γ-quasi-clique are **local** oracles and cannot see it. A **VG path-level colinearity** criterion is the
*non-local* separator the theory calls for: it requires a member's **exon-chain to colinearly match one family**,
and a recombinant chain fails that by construction. fid 210 is the concrete K=3 witness; the VG split is the
theory's predicted resolution.

---

## Verdict

**Reshuffling is a REAL over-merge mechanism in the gorilla catalog, but the smallest bucket and rare.** Only
**10/83** over-merged blocks are recombinant; VG cleanly splits **5**, of which **4 are distributed exon-shuffles +
1 readthrough chimera (fid 321)**, and only **2 are high-confidence distributed mosaics** (fid 210, fid 187), plus
one high-confidence distributed mosaic (fid 75) that VG conservatively spares. The residual over-merge is dominated
by **sub-exonic repeat bridges (34 DISCONNECTED, cross-comp id median 0.586) + cardinality (19) + domain-share
(20)** — **none** of which recombination-detection addresses.

**Does VG recombination-detection earn wiring in?** **On its niche, yes — but as a precision top-up, not a headline
lever.** It is the *only* criterion that catches the K-frontier recombination obstruction (density/edge criteria
provably cannot), and it splits true mosaics while sparing real alt-first-exon families **at zero measured recall
cost**, sitting in a clean 0.571/0.75 valley robust over 0.60–0.75. But the net precision gain is
**~5/606 blocks ≈ 0.8 pt** (realistically ~2–4 blocks once the readthrough artifact and NA-vs-NA / 2-exon-fragment
LOW-confidence calls are set aside). So it is worth wiring in as a targeted, theory-grounded filter for the K≥3
mosaic case, with the honest framing that **the residual over-merge is overwhelmingly domain-share + cardinality +
sub-exonic-repeat, and reshuffling is real but rare**, with **fid 210 (LOC101126070 + GALNT17 → LOC115935264)** the
one textbook barbell.

### Numbers to quote
- **Recombinant over-merge count:** 10/83 blocks classified RECOMBINANT-BRIDGE; VG splits **5** (= **4 distributed
  exon-shuffles [94, 110, 187, 210] + 1 readthrough chimera [321]**; **2 high-confidence distributed** [187, 210]).
- **Path-colinearity result:** splits all 5 (7 bridges) that connected-components + γ-quasi-clique merge; clean
  valley max-split 0.571 < min-spared 0.75; robust plateau over thresholds 0.60–0.75.
- **Recall cost:** 0/24 pure control families wrongly split (0 across the 0.50–0.80 sweep); alt-first-exon families
  (DDX11 fid 288) and terminal-exon variants (fid 339/443) correctly spared. Caveat: control has 0 alt-first
  families; recall evidence rests on the spared recombinant bridges + the valley.
- **Verdict:** reshuffling is real but rare (smallest bucket); VG recombination-detection earns its keep on the
  K-frontier niche at zero recall cost, but is a ~0.8 pt precision top-up — the residual over-merge is dominated by
  domain-share + cardinality + sub-exonic-repeat bridges.
