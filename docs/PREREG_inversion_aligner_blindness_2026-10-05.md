# Inverted copies and the strand-aware spliced aligner — mechanism + a parked, pre-registered measurement (NOT RUN)

2026-10-05. Parked by the user's cost decision ("I don't have much usage — just record it"). No data was looked at beyond
what is cited; every number below is a *claim to be tested*, not a result.

## The mechanism (why this is expected, not speculative)

`minimap2 -x splice:hq` — the aligner every BAM in this project was made with — is **strand-aware**: it forces the
transcript strand from the canonical splice sites and does not emit a spliced alignment whose strand disagrees with them.
A wholly-inverted gene copy's transcripts are reverse-strand relative to the forward copies of its family, so:

- the reads of a whole-copy inversion do not get *weakly* aligned back to the family — the splice-strand constraint
  **refuses** the spliced chain outright; only short intra-exon local matches can survive. The reads land in the
  **unmapped** pool (or as poorly-placed fragments). Seeding is not the first wall; the strand constraint is.
- an **exon-internal inversion** partially survives: the flanking canonical splice sites still anchor the read's strand,
  and the inverted exon body aligns with mismatches/indels.
- an **intron-span inversion** is intermediate: breakpoints lose spliced support, interiors may still align locally.

Consequence for the pipeline as shipped: inverted copies are **invisible to O1 and O2 by construction** — their reads are
not in the mapped pool at all, so there is nothing to assemble into a locus and nothing to assign. This is a property of
every strand-aware spliced aligner, so it holds equally for the *real* gorilla/human BAMs: inverted copies, if any exist
in these gene families, are already sitting in the real unmapped pools today.

## The designed mitigation already in the tree (to be tested, not assumed)

`o3_candidates` pass B (`src/rustle/vg_family/o3_candidates.rs`, Amendments 13/13b) consumes **unmapped and poorly-placed
reads ≥ 300 bp** and attributes them by **`map-ont`** — which is *not* strand-constrained. A wholly-inverted copy is
exactly the class of object that mechanism exists for; whether it actually fires for inversions is an open empirical
question (the inversion changes the sequence locally — snp-level divergence — and `map-ont` needs ≥ 50% read coverage at
`de ≤ 0.20` against the family's net or copies).

## The parked experiment (famsim; ~1 h of compute on the 5-core box when run)

Substrate: `bench/famsim` synthetic background (self-contained, no real data needed; `examples/synthetic_smoke.json` is
the template). Runs: one 100%-identity control (`three_identical`: copies A, B, C, no ops) and one arm per inversion
type, each rung verify-PASS before scoring:
- `inv_whole` (whole-copy strand flip — the wall case),
- `inv_exon` (one exon inverted inside an otherwise-identical copy),
- `inv_intron` (one intron inverted).

Claims, written before any run (adjudicate on `famsim`'s own products, not the pipeline's):
- **C1 (aligner loss):** the fraction of each inverted copy's reads that are unmapped in `align/`s BAM is ~1.0 for
  `inv_whole`, intermediate for `inv_intron`, low for `inv_exon`; the control's unmapped fraction is the ordinary
  error/simulation floor. Measured from the BAM + `families.tsv`/reads truth, independent of rustle.
- **C2 (O2 abstention):** among the inverted copy's reads that *do* map, O2's assignment rate does not exceed the
  control's (they map to the forward copies; PSV/junction evidence should abstain or misassign — reported, not barred).
- **C3 (O3 net):** `o3_candidates --candidates` flags the inverted copy as a candidate for at least one of
  `inv_exon`/`inv_intron`, and for `inv_whole` **iff** pass B's map-ont attribution recovers ≥ `min-support` (6) reads
  into its family's net. Reported beside: the pass-B attributed-read count per inversion type, and whether the flagged
  candidate's nearest locus is the family's own interval.
- **C4 (advisor framing):** the deliverable number is the pair (C1 loss, C3 recovery) per inversion type — "inverted
  copies are an O3/unmapped-read phenomenon by construction of strand-aware spliced alignment; the aligner loses X, the
  net recovers Y" — never a claim that the primary mapping sees inversions.

Not in scope: real-data inversion screening (no annotation of inverted gene copies in these genomes is assumed);
`sim.py`-style standalone inversion support (famsim covers it; nothing else in bench/ makes inversions — see the
2026-10-05 inventory in the session record).
