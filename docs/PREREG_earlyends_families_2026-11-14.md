# PREREG row 1740 (2026-11-14): families non-regression of the early-end chains

User question: can the chains of the `earlyends` stage (rows 1721 / 1722) go into the main assembly file, which the families
stage reads? Rows 1490 / 1491 failed this kind of check for other assembly layers, so it is measured before anything is merged.

**Arms** (same binary, frozen copy in the workspace `r1728_earlyfam/bin`; same driver, default settings, SEDEF damping on):
- base: `families` on the sample's `PREFIX.gtf` + `PREFIX.families.gtf` as they are.
- arm: `families` on `PREFIX.earlyends.gtf` and on `PREFIX.families.gtf` with the added transcripts inserted after the last record
  of their gene (A119b 2,893 of 2,930 chains; 37 have a gene absent from the families input and are left out. OR6737 842 of 842).

**Samples:** A119b (assembly of `r1662_edgemiss/asg_def`, the quoted default) and gorilla OR6737 (`r1628/orrun`). Both are seen samples;
this is a non-regression check, not a held-out test.

**Score:** `pairscore.py` at its defaults on `x.fam.copies.tsv` (with the extra representatives), `--gtf` = the BASE families input
for both arms so the expressed truth set is the same; gorilla with `--species gorilla` (annotation symbols, no protein-hit map).

**Bar (set before the run):** PASS = on both samples precision and recall of the arm are each no more than 0.005 below base.
Reported beside it, no bar: families, copies, predicted pairs, hits, false.
A PASS makes the merge eligible (the user decides); a FAIL keeps the chains in the side file `PREFIX.earlyends.gtf`.

## RESULT (job 4471971, 2026-11-14): PASS

| sample | arm | families | copies | predicted pairs | hit | false | precision | recall |
|---|---|---|---|---|---|---|---|---|
| A119b | base | 1,380 | 6,636 | 1,960 | 1,768 | 192 | 0.902 | 0.808 |
| A119b | arm | 1,380 | 6,637 | 1,960 | 1,768 | 192 | 0.902 | 0.808 |
| OR6737 | base | 382 | 1,358 | 35 | 18 | 17 | 0.514 | 0.200 |
| OR6737 | arm | 382 | 1,358 | 35 | 18 | 17 | 0.514 | 0.200 |

Both samples inside the bar (no change at three decimals). Clusters tables byte-identical on both samples; copy rows differ in a few representatives (see below). The gorilla score rests on
35 coding pairs (annotation symbols only), so it is a weak check on its own. Nothing merged: the user decides.
`r1728_earlyfam/` in the workspace.

Row number: registered as 1728, renumbered 1740 on 2026-11-15 because another session registered the shifted-site
reads under 1728 at the same time. Content unchanged; the workspace directory keeps its name (`job.sbatch` holds its paths).

### Copy tables, arm against base (compared by family and gene id, not by line)

| sample | rows base -> arm | representative replaced by an early-end chain | rows added / removed | other differences |
|---|---|---|---|---|
| A119b | 6,636 -> 6,637 | 1 (MCL152: `DN_chr17_5022171_14` -> `END_chr17_5033542_14_10`, 6 reads both) | +2 (MCL1160: pieces `DN_chr2_130455581_5.p2` and `DN_chr2_131733204_6.p2`, each represented by an early-end chain) / -1 (MCL152: piece `DN_chr17_28369204_2.p2`) | 18 rows with a new `copy_idx` only (MCL56 13, MCL152 3 of which one also changes `max_family_identity`, MCL639 2); 4 gene ids with a changed suffix |
| OR6737 | 1,358 -> 1,358 | 1 (MCL152 copy 3: `DN_NC_073237.2_84404798_18`, 2 reads -> `END_NC_073237.2_84403254_17_6`, 9 reads) | 0 | 0 |

So 3 early-end chains are representatives on A119b (1 replaced, 2 in new rows) and 1 on OR6737. The added and removed
rows are extra representatives (`.pN` pieces, row 1662 famreps); the clusters table is byte-identical. Suffix changes, same
transcript in each: MCL639 `DN_chr2_130455581_5` and `DN_chr2_131733204_6` gain `.p1` (their loci gain the two MCL1160
`.p2` pieces), MCL56 `DN_chr16_29048209_2` gains `.p1` and MCL1160 `DN_chr2_131112717_5` loses it (cause of these two
not examined).

Not tested: `assign` on the arm (it reads the copy table, so the changed representatives could move read assignments);
`all` end to end with the merged file. The merge itself is not done.
