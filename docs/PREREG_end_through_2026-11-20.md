# Row 1841: does an end peak judged against the reads running through separate true contained chains? (pre-registered 2026-11-20)

User (2026-11-20, after row 1840): do the sub-chain question. Rows 1839 / 1840: neither the read ratio nor a peak among a
contained chain's OWN read ends tells a true contained chain from a truncated copy. Left open: the peak judged against the
reads that run through the position (IsoQuant's form, `isoquant_lib/model_construction`, terminal position cutoff 0.05), at
the 3' side as well. **Measurement only, no assembler code.** A build is a separate row and only follows a PASS here: rows
1064, 1486 and 1715 already found no read-end lever for nested chains on the apes, so the larger build is not started on
hope. Workspace `r1841_end_through` (`ends.py`, `ends.sbatch`). Committed before the script was written or run.

## Candidates
The chains switching the collapse off adds before the final filter: in `r1839_subchain_ratio_001/NAME/r0/x.gtf`, not in
`NAME/base/x.gtf` (the set row 1840's separation was read on). A candidate Y is used when it is a contiguous part of a longer
chain X of either file on its strand; a side is *internal* when some X has introns beyond it. Candidates without a container
are counted and left out. CORRECT = contig, strand and chain equal a truth transcript's (truths of row 1839).

## Evidence (primary records of the BAM, no strand filter: an intron fixes the strand)
For an internal 5' side, with J = Y's first intron and p = the commonest 5' end of the reads whose chain equals Y's:
- **S** = reads whose outermost intron on that side is J and whose end lies within 10 bp of p;
- **C** = all reads containing J (the reads running through are C minus those ending in J's flanking exon);
- **F** = reads whose outermost intron on that side is J and whose end lies within 100 bp of p.
- *share* = S / C. *contrast* = S / F.
The same at an internal 3' side with Y's last intron and the 3' ends. A candidate's share and contrast are the minimum over
its internal sides; its S likewise. No polyA signal test (row 1064 measured it: it separates weakly).

## Signatures, fixed now
| | holds when |
|---|---|
| T05 | share >= 0.05 and S >= 2 (IsoQuant's cutoff) |
| T20 | share >= 0.2 and S >= 2 |
| T50 | share >= 0.5 and S >= 2 |
| K50 | contrast >= 0.5 and S >= 3 |
| T20K | T20 and K50 |

## Reading, fixed now
- A signature SEPARATES on a dataset when the candidates holding it are at least twice as precise as those not holding it
  (row 1840's wording) and it holds at least 20% of the correct candidates.
- **PASS** = one signature separates on at least 6 of 7 non-primate and 4 of 5 primate datasets. A PASS makes that signature
  the rule of a build row (opt-in switch, judged by "pays" as in rows 1839 / 1840); the user decides.
- **NON-PRIMATE ONLY** = at least 6 of 7 non-primate and fewer than 4 primate: reported as such; a rule that acts on mouse and
  not on apes is the user's call.
- Anything else = FAIL: read ends at the collapse give no lever, and the question is closed with rows 1064, 1486, 1715, 1840.
- Reported beside it, no bar: precision of the candidates by share in the bins < 0.05, 0.05-0.2, 0.2-0.5, >= 0.5 and by
  contrast in < 0.25, 0.25-0.5, 0.5-0.75, >= 0.75; the same split by which side is internal (5', 3', both). This shows
  whether ANY cutoff could separate, without choosing one after the fact.

## Expectation, written before the script
FAIL. The share is close to the read ratio times the peakedness, and each failed alone (rows 1839, 1840); row 1715's `own`
feature on the apes was the ratio rule itself. If anything separates I expect it at the 3' side (row 1715: NEST3 0.03-0.04
annotated against NEST5 0.006) or as K50 on the mouse libraries. On the simulated sets the starts are exact, so true
contained chains there should have contrast near 1; the admitted chains of row 1840 were LESS precise there, which I cannot
explain yet and expect this table to show.

## Limits
- Candidates come from the polish output with the collapse off, not from its input; a contained chain another polish step
  removes is not seen.
- The reads are not deduplicated as the polish does, and the junction snap is not applied to them: a candidate whose
  junction was re-placed has no reads of its own chain; such candidates are counted and scored as holding no signature.
- Real-library truths are incomplete annotations; simulated reads carry no degradation.
- Twelve reused datasets (rows 1837-1840), said.
