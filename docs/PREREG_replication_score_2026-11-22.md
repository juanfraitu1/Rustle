# Row 1898: a one-read ranking fitted to replication instead of to the annotation (pre-registered 2026-11-22)

User goal (2026-11-22): "... better rules for the ones only covered by 1 long read or the junctions shared by 1 read."

State. The one-read tier is ranked by `oneread_score` (row 1806): 23 read features, among them the weakest junction's read
support (`min_isup`, `min_isup_rel`), the partial reads on the chain, its share of the locus, read-end concentration. So the
support of single-read junctions is already a feature. Row 1806's verdict was "not separable": on real libraries no fit
reaches the precision of two-read chains, and the label was the ANNOTATION, which cannot judge a real transcript seen once.
Rows 1890 / 1895: two hand-made orders are no better.

New evidence, looked at before this registration (`r1876_whats_left/repl.py`): WTC11 has three replicates of one sample.
Of the 52,756 one-read chains of replicate 2, 15,257 (29%) are seen again in replicate 1 or 3; those are 31.3% annotated, the
rest 5.3%. The frozen score hardly predicts it (best decile 41% seen again, worst 20%). Replication is a label that does not
depend on the annotation, and it is the same quantity row 1808's partner promotion uses with a second library in hand.

**Rule tested:** the same 23 features and the same design as `oneread_r1806.json`, refitted to "seen again in another
replicate" on WTC11 r2. The fitted ranking needs ONE library when applied; the replicates are training data only.

## Fit and data
- Features: `chain_rules_feat.py` over the floor-1 assemblies that exist (rows 1880, 1887), design of `oneread_tier.py`,
  standardised inside each library; dumped by `r1898_replication_score/feats.py`. No new assembly.
- Label on a library L of WTC11: a one-read chain is REPLICATED when its intron chain is in `x.gtf`, `final.gtf` or the
  one-read tier of either other replicate.
- **Model R:** logistic fit (`chain_rules.fit`) on r2's one-read chains, label = replicated.
- **Model A (control, same data, other label):** the same fit with label = annotated (RefSeq on CHM13). It tells whether a gain
  comes from the label or from refitting on this library.

## Measures
1. Inside WTC11, on r1 and r3 (not fitted): share replicated and share annotated in the best 10% and 30% of the tier, frozen
   score against R and A.
2. **Transfer, the bar:** on the datasets of row 1880 where StringTie is larger than `final.gtf` (17), correct chains of the
   ranked list at StringTie's size with tier 2 ordered by R, against the frozen order (`rank.py`'s count). Same for A.
3. The `list_sensitive` cut (best 10% of the tier): annotated chains under each order.

## Reading, fixed now (bars of rows 1890 / 1895)
- **ADOPTABLE** = R gives at least as many correct chains as the frozen order at StringTie's size on at least 15 of the 17
  datasets, on at least 3 of the 4 non-primate ones, and loses no AHEAD of row 1880. WTC11 r2 is fit data and is counted but
  marked; the verdict must also hold without it.
- **LABEL EFFECT** is claimed only if R is ADOPTABLE and A is not.
- If R is not adoptable the frozen ranking stays, and the reading is that single-library features cannot predict which
  one-read chains recur, whatever the label.
- An ADOPTABLE R is a candidate replacement of `tools/data/oneread_r1806.json` (the user's call); it changes only the order of
  the one-read tier, hence `PREFIX.ranked.gtf` and `PREFIX.oneread.ranked.gtf`, never `final.gtf`.

## Limits, stated in advance
- Replicated is not "real": an artefact of the aligner or of the library chemistry recurs too, and a real transcript at one
  molecule per library often does not (replicates 1 and 3 are larger than 2, which helps).
- One sample (human iPS), three replicates. Transfer to the other nineteen datasets is the only protection against fitting it.
- The judging truth on the transfer bar is the annotation again; a ranking that finds real unannotated chains is not rewarded.
