# Row 1678 — among the relative tier's picks, a shortlist of young duplicates (2026-11-13)

User choice after row 1676: "Young-duplicate links only". Families layer (O1); opt-in only if found; not for
`advisor_share`. This reopens a model on the row 1674 features, which row 1675 had closed, at the user's request.
Prior stated before the run: row 1675 asked nearly the same question in one stage and reached 0.61 at 20 true pairs
out of fold; a second stage changes the negatives (relatives and unrelated picks, 5 to 1, not all refused pairs,
1,300 to 1), so a modest gain is the most I expect.

**Stage 1 (frozen, row 1676):** the relative-tier model and threshold 0.0511; its picks are the population here.
**Stage 2 (frozen, `bench/denovo_eval/edgeyoung.py`):** the same 20 features and the same L2 logistic; positives =
simiiformes truth pairs, negatives = every other pick. Fit on the A119b + HG002 picks pooled. Threshold = the lowest
stage-2 score whose 5-fold out-of-fold precision is >= 0.75 with >= 20 truth pairs.
**Go/no-go:** no such threshold -> STOP, gorilla not scored.
**Test (frozen, comparable tables only):** OR6737 rebuilt (`rel.OR2.tsv`) + KB3781: pooled >= 20 truth pairs at
precision >= 0.60 against the young-duplicate truth.
**If it passes:** row 1674's step 2 (family level through `--extra-edges`, same bars) decides whether a rule for the
family definition is found. **If it stops or fails:** the row closes and I do not propose another model on these
features.

(Numbered 1678: `r1677_chainrules` in the workspace belongs to another session; my check and the write ran in one command.)

## Result: FAIL held out (2026-11-13, `r1674_edgerule/young.out`)

Stage-1 picks on the training tables: 1,328, truth pairs 207. Out-of-fold threshold: stage-2 score >= 0.261 with
23 true / 7 false (0.767): go/no-go GO. Further down the out-of-fold ranking: 40 true at 0.494, 80 at 0.374.

Gorilla, frozen: rebuilt OR6737 3 true / 9 false; KB3781 1 / 7; pooled 4 true / 16 false, precision 0.200.
Bars: >= 20 true at >= 0.60. **FAIL.** The row closes; no further model on these features.

Reading: the few young duplicates the stage-2 score ranks first on human (long alignments, few exons, high partner
rank) are specific to the human samples; on gorilla the same score picks older relatives. This agrees with rows
1674 and 1675: nothing in these features dates a duplication.
