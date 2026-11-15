# Row 1680 — a second sample as witness for a refused family pair (2026-11-13)

Goal: "find empirical rules to improve edge recall in the multi-copy gene families definition". User's choice after
row 1679 stopped ("Cross-sample witness"). Families layer (O1), opt-in only if found, not for `advisor_share`. No read
is assigned to a copy; families and `fam.copies.*` are not changed by this row.

**Mechanism tested.** A refused pair that scores high in one sample is often an older relative or noise. If the false
picks differ from sample to sample while young duplicates score high in every sample that expresses them, requiring
two samples of the same species to agree should raise precision. Caveat registered now: older relatives are a property
of the genome, not of the sample, so they may agree across samples just as well; that is what the test measures.

**Rule (frozen, `bench/denovo_eval/edgewitness.py`).** Row 1675's model unchanged (L2 logistic, the 20 features, truth
positives, fitted on A119b + HG002). A pair is keyed by its two gene names. Witness score = the lower of the pair's
scores in the two samples; pairs refused (not linked, stage not `edge`) in only one sample are out of the population.
Out-of-fold folds are by gene-name pair, so a pair's two rows always share a fold.

**Step 0 (ceiling, before any fit).** Truth pairs refused in both human samples. Fewer than 20 -> STOP.
**Go/no-go.** Threshold = lowest witness score with out-of-fold precision >= 0.75 and >= 20 true pairs (pairs counted
once). None -> STOP, gorilla not scored.
**Test (frozen threshold; rebuilt OR6737 + KB3781, same GGO reference):** >= 20 new true pairs at precision >= 0.60.
**If it passes:** row 1674's step 2 (family level through `--extra-edges`, same bars). It would be a rule only for
species with two samples, and is reported as such.
**If it stops or fails:** the row closes. No further variant of the row-1675 score is proposed by me.

## Result (2026-11-13, workspace `r1680_witness/witness.out`) — STOP

Step 0 GO: 52,937 pairs are refused in both human samples, 60 of them truth pairs. No witness threshold reaches 0.75
with 20 true: out of fold, 20 true come with 43 false (0.317) and 40 with 165 false (0.195). Row 1675's single-sample
score had 0.606 at 20 true, so asking a second sample to agree lowers precision. The registered caveat holds: the false
high scorers are properties of the genome (older relatives) and repeat in every sample, while only 60 of the ~280
refused truth pairs are refused in both samples to begin with. Gorilla not scored. Row closed.
