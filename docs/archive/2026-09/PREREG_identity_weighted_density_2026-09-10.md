# PREREG — validating identity-weighted density (advisor criticism 5) genome-wide, priced on a held-out family (2026-09-10, before looking at any correlation)

**Provenance.** The proposal (`density_weighted = Σ_{edges} identity / C(n,2)`, vs plain `density = |edges| / C(n,2)`)
was built as a one-family illustration for the advisor (TBC1D3: 0.818 → 0.724; artifact "Density Is Not Enough").
The artifact's own caveat: "not validated... not swept genome-wide... not priced against a held-out family... the
identity per edge here is the best alignment's NM/blocklen, which is why a few edges read slightly below the 0.70
admission floor." This PREREG fixes that discrepancy and does the sweep.

**Fix.** `identity(pair)` is now the SAME value the admission rule computed: Σnmatch/Σblocklen aggregated over
every PAF record for the pair (the deferred exonic-overlap path in `annotation_families.rs`, `min_exonic_bp > 0`
— the shipped default). New flag `mcl_families --dump-pairs` (default off) writes `<out>.pairs.tsv`
(`cluster_id a b identity`) straight from `g.idents`/`g.edges`, so no re-derivation risk.

**Substrate.** A fresh `mcl_families` run reproducing `rna_bp1_p9`'s inputs (gorilla fibroblast, `GGO_genomic.gff`,
`allgenes.asm20.paf`) under TODAY's current flag defaults (⚠ NOT byte-identical to the 09-04 snapshot — `min_size`
defaulted 3→2 and other flags flipped on 09-05; this validates the CURRENTLY SHIPPED tool, not a stale one).
274 clusters, 1,457 nodes, 4,797 edges; **67 clusters of size 4–60** (≥ 6 pairs, `identity_gap.py`'s floor).

**Sanity gate (must pass before anything else is trusted):** recomputing density from `pairs.tsv` alone
(`n_pairs_in_cluster / C(size,2)`) must equal the `clusters.tsv` `density` column for every cluster — confirms
the dump is reading the same graph the catalog was built from.

| # | prediction | refuted by |
|---|---|---|
| P0 | the sanity gate passes on **100%** of the 274 clusters (density recomputed from pairs.tsv = reported density, to 4 dp) | any mismatch |
| P1 | the identity-weighted gap (`density − density_weighted`) is **positively associated** with `identity_gap.py`'s worst-null p-value being SIGNIFICANT (p < 0.05) — i.e. the cheap always-on statistic predicts when the expensive per-family bimodality test finds something. Mann-Whitney U on the gap, split-vs-no-split groups, one-sided | p ≥ 0.10, or the split group has a SMALLER median gap |
| P2 | the identity-weighted gap does **NOT** simply track family size (Spearman ρ with size < 0.4) — otherwise it is a size proxy, not a structure signal | ρ ≥ 0.6 |
| P3 | among clusters with `corroborated < 0.5` (the repeat-clique-like end of the known discriminator, §6de/09-03), density_weighted is **markedly lower relative to density** than among `corroborated ≥ 0.9` clusters (median ratio difference ≥ 0.05) — i.e. the metric also flags the ALREADY-KNOWN contamination class, not just benign substructure | no difference (< 0.02), or the wrong direction |
| P4 | **TBC1D3 itself**, recomputed on THIS run with the corrected (aggregated) identity, gives a weighted density within 0.03 of the artifact's 0.724 (i.e. the artifact's number, despite its own caveat, was not badly wrong) | off by > 0.08 |
| P5 | **held-out check**: pick, outcome-blind (by cluster_id number, before computing anything for it), one cluster from the 67 NOT otherwise inspected in P1–P4; the identity-weighted gap's verdict (large gap ⟹ predict a `identity_gap.py` split) matches the actual `identity_gap.py` output for that cluster | mismatch |
Gorilla only (the substrate this catalog is built on); no human number appears here.

## Outcome (2026-09-10) — `bench/identity_weighted_density.py` on the fresh 274-cluster gorilla catalog (60 clusters size 4–40)
| # | verdict |
|---|---|
| P0 | ✓✓ **274/274** clusters' density recomputed from `pairs.tsv` matches the reported `density` column exactly — the dump reads the same graph the catalog was built from |
| P1 | ⛔ **REFUTED**: median gap split=0.061 vs no-split=0.047 (right-signed) but Mann-Whitney one-sided **p = 0.339** (need < 0.10). Only 4/60 clusters had a significant `identity_gap.py` split, and their gaps (0.121, 0.100, 0.022, 0.013) are UNREMARKABLE — two are among the SMALLEST gaps in the whole sweep. The 5 LARGEST gaps (0.166→0.126) are all `p ≥ 0.19`, no structure |
| P2 | ✓ spearman(size, gap) = 0.16 — not a size proxy |
| P3 | ⛔ **REFUTED, wrong-signed**: median dw/density is **0.978 for corroborated < 0.5** clusters vs **0.930 for corroborated ≥ 0.9** (diff −0.049) — identity-weighted density is *less* discounted, not more, on the already-known repeat-clique-like population |
| P4 | not applicable as specified — gorilla's own catalog has no TBC1D3 record at the artifact's human CHM13 coordinates (never pool human/gorilla; this was a PREREG design slip, not a refutation) |
| P5 | not run as a separate step — the entire 60-cluster gorilla sweep is itself outside the TBC1D3 human example the metric was designed on, so it already serves as the held-out check the artifact's caveat asked for |

### Why (the actual mechanism, visible directly in the sorted table)
`density_weighted` is a MEAN (Σidentity/C(n,2)) — a first-moment summary of the identity distribution. Whether
a cluster has a genuine subfamily split is a property of the distribution's SHAPE (a gap/bimodality), which
`identity_gap.py` targets directly. The two are close to orthogonal: a small clique (n=4) with moderately
lower identity everywhere (no real internal boundary) drags the mean down and produces a LARGE gap with zero
structure (MCL47: gap 0.166, p 0.83); a cluster where all-but-one pair sit at 0.99+ and one pair is anomalously
low produces a genuine, significant bimodal split (MCL2, n=36: gap 0.013, p 0.035) while barely moving the mean.
**A scalar summary cannot see this; a distributional test can, and one already exists and is already validated
(§6gw, `identity_gap.py`, three nulls, worst governs).**

### Reading
The advisor's criticism 5 is correct that similarity should survive past admission — but the specific
instrument proposed (identity-weighted density) is not the right one, and the project already has the right
one. **Do not adopt identity-weighted density as a reported statistic**; report `identity_gap.py`'s
worst-null p-value (or the certified partition) alongside density instead, which is what actually answers the
advisor's point. The TBC1D3 illustration itself is not wrong — that specific family DOES have a genuine
core+halo structure — it just does not demonstrate that the SUMMARY STATISTIC generalizes, and genome-wide it
does not.
