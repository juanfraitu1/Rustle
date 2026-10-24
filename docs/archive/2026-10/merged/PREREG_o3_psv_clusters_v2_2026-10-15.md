# PREREG — O3 PSV clustering v2 (linkage test) in Rust (`copy_assign --psv-clusters`), held-out on base and hard worlds (2026-10-15)

**Written after development and before any held-out world below is built.** Follows `docs/PREREG_o3_unexplained_2026-10-15.md`
(v1, PASSED). Design: `docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`.

## What changed since v1, and why (development only: seeds 32, 33, 51, 52)

- **Harder worlds** (`world.py --design hard`): the heterozygous allele background (H = A + 0.1 / 0.2 / 0.4%), a copy plus a background
  allele on one host, 8-copy families (intact and with one copy removed), gene conversion (an absent copy with two exons converted back
  from its host; present copies with converted exons), a mosaic absent copy, 5x sequencing error (a pessimistic stress case; Kinnex
  HiFi is ~0.1%), expression imbalance (host 200 reads, absent copy 6), exon-skipping isoforms, an absent copy that lost an exon,
  A-to-I RNA editing (8 sites, each in 40% of a present copy's reads), an absent copy at 4 reads, and an absent copy whose reads split
  between two hosts 0.4% apart. A 10x-error scenario was built and dropped (1% error is not a Kinnex regime).
- **v1 broke on editing** (false unexplained clusters of a present copy's reads: 3-9 reads, A->G at every site) **and on 1% error**
  (one 5-read cluster). **v2 adds the linkage test**: a missing sequence is a haplotype, so its alts travel together; edits and errors
  hit sites independently. Over every pair of a cluster's unexplained sites, P(alt at both | alt at either) among the host's reads
  covering both; median >= 0.8 or the cluster is `unlinked` (editing signature reported as the A->G / T->C fraction). Dev: every true
  missing sequence 0.95-1.00, editing 0.27-0.40, the error cluster 0.20.
- **Dev result of v2** (arm P, 4 worlds x E0/E2): absent copies >= 0.6% 207/208; control families 0/80; spurious unexplained clusters
  0; multi-group separation 16/16; alleles labelled `allele_or_copy` 70/70; copies >= 2% labelled `copy` 120/124 (the 4 are CONV_ABSENT,
  whose effective divergence is lowered by the converted exons). v1's base worlds unchanged (35/35, 43/43).
- **Rust port** (`rustle::missing_copy::psv_clusters`, `copy_assign --psv-clusters`; a post-pass over the run's own
  `<out>.assignments.tsv`, read sequences from one BAM sweep in their original orientation). **P0 parity: 1,241 / 1,241 clusters
  identical** to the prototype (membership by source, status, site counts, explaining copies, est_div / linkage / editing within
  rounding) on all 8 dev runs; the assignments it reads are byte-identical to the prototype's input; every other `copy_assign`
  product is byte-identical with and without the flag. 4 unit tests (cs parsing, Counter tie rule, a missing haplotype, unlinked edits
  and an explaining copy).

Frozen (md5): `psvclust.py` 17888625c3d31245d03e2221e5086dc0 (the reference), `score.py` 51b602dab31809c3743fb4d7aaacbcb8,
`world.py` b9b9cc581cdf845eb65001ed748a94d9, `parity.py` 858a8c4282b22fb2278a71831ba31b86; `copy_assign` sha1 f36d4e7a.

## Held-out run (the RUST product is scored: `score.py rows ... --rust`)

Base design seeds **61, 62, 63** and hard design seeds **71, 72, 73, 74, 75** (`--strands`), reads E0 and E2, arm P
(`copy_assign --families --flag-missing-copies --no-as-tied-only --psv-clusters`). Parity with the prototype is re-checked on every
held-out run (`parity.py`) and reported.

## Bars (pooled over all held-out runs)

| id | measure | bar | dev v2 |
|---|---|---|---|
| **H1** | absent copies with d >= 0.6% recovered (pure unexplained cluster) | >= 0.90 | 207/208 |
| **H2** | control families with an unexplained cluster (base: NULL_TWIN, NULL_MULTI, NULL_CLOSE, ADDBACK, DECOY; hard: MANY_NULL, NULL_CONV, HIGHERR_NULL, EDIT_NULL) — 152 units | <= 3 | 0/80 |
| **H2b** | spurious unexplained clusters (majority of a PRESENT copy's reads), all families | <= 3 | 0 |
| **H3** | families with two or three missing sequences on one host (COMBO_EH, DOUBLE, HET_COPY): all recovered separately | >= 0.80 | 16/16 |
| **H4** | recovered alleles labelled `allele_or_copy`; recovered copies with d >= 2% labelled `copy`, excluding CONV_ABSENT | >= 0.90 each | 70/70; 116/116 |
| **P0'** | held-out Rust-vs-prototype parity | every cluster identical | 1,241/1,241 |

## Reported beside (no bar)

Every scenario; the editing families (EDIT_NULL / EDIT_COPY: their unlinked clusters, linkage and A->G fraction); the 5x-error
scenarios; copies at 0.3% and 1%; CONV_ABSENT's labels; the existing detector (`--flag-missing-copies`) as the baseline.

## Reading, fixed now

- **H1-H4 and P0' pass** → `copy_assign --psv-clusters` (with `--no-as-tied-only`) is O3's method at simulation level, robust to the
  allele background, editing, conversion, many-copy families, imbalance and isoforms. Next level: real families on the main machine
  (remove copies from the GENOME, re-map, add back) with a single-copy-gene control for the real allele background, pre-registered.
- **P0' fails** → the Rust port has a defect: fixed in the port, never in the rule; the bars are read only after P0' passes.
- **H1, H2 or H2b fails** → v2 is not adopted; the failing scenarios are named; no constant is re-tuned on held-out worlds.
- **H3 or H4 fails alone** → detection is claimed, separation / labelling is not.

# OUTCOME (2026-10-15) — all bars PASS; the Rust product scored

`pixi run bash bench/o3_unexplained/heldout_v2.sh ~/rustle_data/o3_unexplained/heldout_v2` (base 61-63, hard 71-75, E0 / E2, arm P;
`copy_assign` sha1 f36d4e7a; frozen md5s as above).

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| P0' | Rust vs prototype, every held-out run | **16 / 16 runs, 0 cluster differences**; assignments byte-identical 16/16 | identical | ✅ PASS |
| H1 | absent copies >= 0.6% recovered | **425 / 426 = 0.998** (existing detector 414/426) | >= 0.90 | ✅ PASS |
| H2 | control families with an unexplained cluster | **0 / 152** | <= 3 | ✅ PASS |
| H2b | spurious unexplained clusters of a present copy's reads | **0** | <= 3 | ✅ PASS |
| H3 | two or three missing sequences on one host, all separated | **41 / 44 = 0.932** | >= 0.80 | ✅ PASS |
| H4 | alleles labelled `allele_or_copy`; copies >= 2% labelled `copy` (CONV_ABSENT excluded) | **135 / 135; 262 / 262** | >= 0.90 each | ✅ PASS |

**Reported beside.**
- RNA editing: EDIT_NULL 0/20 false flags, EDIT_COPY's absent copy 20/20; the editing produced 68 `unlinked` clusters (11 in
  EDIT_NULL, 57 in EDIT_COPY), none unexplained — the linkage test is what keeps them out.
- The one missed absent copy: HIGHERR@0.01 (5x error, a stress case), hard73 E2 f26 (its cluster was `unlinked`). HIGHERR@0.006
  20/20, labelled `allele_or_copy`.
- The 3 H3 misses are HET_COPY's background allele at 0.2% (17/20 recovered); its absent copy 20/20.
- The allele background: HET_BG at 0.1% recovered 14/30 (one or two variant sites: below what a cluster can call), 0.2% 28/30, 0.4%
  30/30, all labelled `allele_or_copy`; never a false `copy`.
- CONV_ABSENT labelled `copy` 20/20 here (4/8 `allele_or_copy` on dev: the label follows the conversion-reduced divergence).
- Labels at delta: copies at 1% split (TWIN 8/4, UNIQ 8/4, MULTI 8/4, MANY_REMOVED 12/8, DOUBLE 6/6, SPLIT 8/12 copy / allele_or_copy).
- Copies at 0.3% recovered 24/24 (existing detector 3/24).

**Reading (as fixed).** `copy_assign --families --no-as-tied-only --psv-clusters` is O3's method at simulation level: every missing
copy from 0.3% upward is recovered (425/426 at >= 0.6%), no intact family is flagged across the allele background, editing,
conversion, 8-copy families, imbalance, isoforms and 5x error (0/152), two or three missing sequences on one host are separated, and
the label states which ones RNA cannot tell from an allele. Next level: real families on the main machine (remove copies from the
GENOME, re-map the real reads, add back) with a single-copy-gene control for the real allele background, pre-registered before it runs.
