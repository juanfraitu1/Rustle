# PREREG r1630: tie-group loci by correlation clustering plus a merge of duplicate models (2026-11-11)

User (2026-11-11, AskUserQuestion): resume the tie-group layer as the way family members share information. Its open
point is records -> loci: no rule passed the duplicate-model bar on A119b (r1620: span 0.058, cc 0.068, bar 0.05).
Written before the rule was run on either sample.

## Rule (one change, `tie_groups.py build ... ccm`)
`cc` as registered in r1620, then: loci of the same chromosome and strand whose MODELS attract (aligned blocks overlap
by >= 50% of the smaller model, the r1612 criterion) are merged, the merged locus is modelled again, and this repeats
until no two models attract (at most 10 rounds). Everything else frozen.

## Honest scope
Both samples have been seen under cc, and the merge is aimed at the bar that failed (D1). So D1 passing proves little
by itself. What the run tests is whether the merge costs anything: D2 (swallowed loci) and the r1607/r1611 bars must
still pass.

## Bars (A119b `r1603/recs.slim.tsv`, gorilla OR6737 `r1611/recs.slim.tsv`), as r1620
- D1 <= 0.05, D2 <= 0.01; B1g grouped reach (scorer's bar), B2s >= 0.90, B3 TAF11L >= 20 (A119b), B5r, B4 <= 1 h.
- PASS on both: ccm is the recipe's locus rule, figures replace r1607/r1611, and the next step is the genome-wide
  run and the product form (group table + GTF as an opt-in pipeline stage). FAIL: reported, locus rule stays open.

## Outcome (2026-11-11, job 4466248, `r1630/`): PASS as registered, with the scope stated above

| sample | loci | D1 duplicate models | D2 swallowed | grouped reach | purity B2s | TAF11L <= 5 kb | largest group | build |
|---|---|---|---|---|---|---|---|---|
| A119b | 11,356 | 0 | 33 = 0.003 | 69 of 72 | 0.999 | 30 of 31 | 155 loci, 1 chromosome | 3,324 s |
| gorilla OR6737 | 1,554 | 0 | 3 = 0.002 | 127 of 128 | 1.000 | - | 33 loci, 1 chromosome | 227 s |

- Every bar passes on both samples. D1 = 0 is by construction (the merge stops when no two models attract) and
  proves nothing; the test was the cost, and there is none on these bars: reach, purity and swallowed loci are where
  cc left them (cc: 69 of 72, 0.999, 0.003 on A119b; 127 of 128, 1.000, 0.003 on gorilla).
- The merge removed 900 of cc's 12,256 loci on A119b and 25 of 1,579 on gorilla.
- B4 passes by 276 s on A119b; the build is close to its hour.
- Per the registration, ccm is the recipe's locus rule. Changed since r1607: A119b largest group 155 loci (86 before),
  still on one chromosome. Not done: the genome-wide run and the product form.
