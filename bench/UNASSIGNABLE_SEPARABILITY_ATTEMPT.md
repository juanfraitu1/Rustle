# Can ANYTHING separate the unassignable reads? — an exhaustive BAM + scikit attempt (negative)

To show the K-frontier ("tied", `n_decisive = 0`) reads are unassignable by *evidence*, not by lack of
trying, we threw every BAM-derivable feature and scikit-learn at them. **Nothing separates them by copy
except the PSV alleles themselves.** Documented so the impossibility is empirical, not just asserted.

## What we tried

**17 BAM-record features per read** (`dsfam42_separability_features.py`): exonic length, # exons,
genomic span, soft-clip fraction + 5′/3′ soft-clip, strand, MAPQ, **AS / NM / de** (alignment score,
edit distance, divergence), GC%, base-quality mean/min, sequence entropy, and 5′/3′ offset relative to
the copy locus. Then (`dsfam42_separability_sklearn.py`): RandomForest (supervised) + KMeans
(unsupervised), with chance baselines and permutation/CV.

## Results

### Test 1 — supervised, DSFAM42 (real): CONFOUNDED (and the confound is itself the floor)
DSFAM42 confidently resolves **only one** copy (copy 8 = 245/250 assigned; the other 8 are the
collapsed near-identical ones), and the de-novo family has a 33 kb **container** copy that dominates
read-overlap — so the labels are ~single-class and cannot benchmark separability. *That only the
divergent copy is assignable is exactly the floor.* (Honest: not a clean test, reported as such.)

### Test 1b — supervised control, sim5x (balanced, TRUE labels): AT CHANCE
sim5x K8: 5 tandem copies, **40 reads each (balanced)**, true copy in the read name. RandomForest
predicting the true copy from sequence/quality features (GC, QV mean/min, entropy, length, # exons):
**accuracy 0.245 vs chance 0.200** — at chance. **No sequence/quality feature separates the copies;
only the PSV alleles do.** The clean, balanced, ground-truth version of the negative.

### Test 2 — unsupervised, the actual tied reads (DSFAM42): structure exists, but NOT copy structure
KMeans on the 1,051 tied reads (non-AS features): silhouette **0.47** (so they *do* cluster) — but the
**adjusted Rand index vs the minimap2 copy placement = 0.00** (random). The clusters are driven by
**read length / technical** structure (the container vs small-copy bimodality), *not* copy identity. No
feature groups the unassignable reads by copy.

## Verdict

Across 17 BAM features, RandomForest, and KMeans: **nothing recovers copy identity for the
`n_decisive = 0` reads.** On balanced ground-truth (sim5x) a classifier is at chance (0.245 vs 0.200);
on real tied reads the only structure is read-length, orthogonal to copy (ARI 0.00). This is the
information-theoretic floor, now **empirically checked with ML** rather than asserted: a read that spans
no copy-distinguishing site carries no copy signal that any feature, classifier, or EM can extract —
because the molecule is invariant under relabeling the copies. The PSV alleles are the *only* separator;
where a read has none, the honest output is abstention.

(Caveat: DSFAM42 is size-heterogeneous, which degenerated the supervised real-data test; the balanced
sim5x control + the unsupervised real-data ARI together carry the negative. A fully clean *real* balanced
multi-copy supervised benchmark would strengthen it further, but every axis tested here came back null.)

Artifacts: `dsfam42_separability_features.py` · `dsfam42_separability_sklearn.py` ·
`/tmp/dsfam42_features.tsv` · `/tmp/sim5x_feat.tsv`.
