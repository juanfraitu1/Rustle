# Annotation-free two-pass: read-coherence preserves per-read PSVs → recovers copies flow loses

Prototype of the architecture: **don't use StringTie flow** (which collapses reads and discards the
per-read PSV linkage); use **read-coherence / direct transfrags** as Pass 1 so per-read PSVs survive to
Pass 2, where families are identified and reads assigned to copies. Demonstrated on the *collapsed*
5-copy regime (all copies' reads piled onto a single-copy reference), where the contrast is starkest.
No annotation, no StringTie. `bench/twopass_denovo.py`.

- **PASS 1 (read-coherence):** map all 5 copies' reads to the single-copy reference, group reads by
  exact intron chain → transcript skeletons, keeping each read's sequence. → **ONE skeleton** (the
  copies share the chain), with all reads + their PSVs.
- **PASS 2 (no annotation):** call PSVs *de novo* from the skeleton's collapsed pileup → split reads by
  PSV allele-vector into copies (copy_split) → declare the multi-copy family → assign each read to a copy.

## Result (PSV ladder K = identifiability axis)
| K | Pass-1 skeletons | de-novo PSVs | copies recovered (of 5) | read→copy assignment acc |
|---|---|---|---|---|
| 0 (identical) | 1 | 0 | 0 | 0.0 — UNASSIGNABLE (no info) |
| 1 | 1 | 1 | 4 | 0.80 (1 column < 5 copies) |
| 2 | 1 | 2 | **5** | **1.00** |
| 4 / 8 | 1 | 2 | 5 | 0.99 |

(de-novo PSVs saturate at 2 = ⌈log₄5⌉, the true # of distinguishing columns for 5 copies.)

## What this shows
- **The decisive contrast:** Pass 1 yields ONE skeleton at every K — exactly where **StringTie flow
  would stop and emit ~1 transcript, losing the 5 copies entirely** (their differences collapsed away).
  Read-coherence keeps the skeleton *with its reads*, so Pass 2 can call PSVs and **recover the 5 copies**.
- **The pipeline is annotation-free end-to-end** — skeletons and PSVs come from the reads, not a GTF.
- **It hits the same identifiability boundary** as the annotation-anchored version: 5/5 copies at K≥2,
  partial at K=1, impossible at K=0 (identical) — confirming the substrate change doesn't move the ceiling,
  it just *preserves the signal* needed to reach it.

## Honest scope
- This is the COLLAPSED regime (the copy_split target) — the case where read-coherence's advantage over
  flow is clearest. For dispersed copies (separate loci) the coordinate already resolves them and flow's
  collapse is less harmful.
- Identifiability ceiling unchanged: a read spanning no distinguishing PSV (tied) stays unassignable.
- Prototype on one synthetic family with truth; a production read-coherence Pass 1 would also handle
  noise (chimeras, minor splice variants) that this clean simulation doesn't stress.

## Verdict
Confirms the architectural point: **read-coherence (not flow) is the correct Pass 1 for a copy-aware,
annotation-free transcript caller** — it preserves the per-read PSVs that let Pass 2 recover and assign
copies that flow would have collapsed. The two-pass `(read-coherence skeleton) → (family + PSV split)`
is `copy_split`'s `(intron chain ⊕ PSV)` realized as a full pipeline.

## Genome-wide: not restricted to a handful (bench/twopass_genomewide.py)
The synthetic prototype above is one locus only because it needs a constructed collapsed reference +
ground truth. The pipeline itself is genome-wide — every stage already ran at genome scale this session
(read-coherence: 25 chroms; family graph definition: 1,337 families; de-novo PSV / co-segregation:
10,178 loci in the hidden-collapse scan). Composed over ALL graph-defined families:

- **1,337 families processed genome-wide** (603,267 reads at family loci) in ~5 s.
- **848 (63%) dispersed** — copies at separate loci → coordinate-resolved (Pass-1 separate skeletons).
- **489 (37%) co-located** — copies share/near one frame → de-novo PSV split applies.
- **hard multimappers (MAPQ-0): 3,803 (0.6%)**; **46 co-located families have ≥5** — those are where the
  PSV-split is genuinely decisive.

So: genome-wide-capable and fast; most families resolve by coordinate, the PSV-split is decisive at the
~46 co-located+hard families (the sparse hard regime — consistent with every prior finding; abundant
only in deep co-located data like testis HiFi).

## Reproduce
- `MINIFORGE python bench/twopass_denovo.py` ; `python3 bench/twopass_fig.py`  (synthetic, with truth)
- `MINIFORGE python bench/twopass_genomewide.py`  (genome-wide tally over all families)
