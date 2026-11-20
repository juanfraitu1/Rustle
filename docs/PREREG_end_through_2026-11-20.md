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

## Result (2026-11-20): FAIL; no end signature separates on the non-primate side

Jobs 4489918-4489929, all exit 0 (24-167 s); `r1841_end_through/ends.NAME.txt`, per candidate `cand.NAME.tsv`. Every
candidate has a container; 2 of 255,804 have no reads of their own chain.

Precision of the candidates holding the signature / of the rest (share of the correct candidates held); **bold** = SEPARATES:
| | fly sim | mouse sim | zebrafish sim | mouse ES | heart | adrenal | cortex | OR6737 | A119b | chimp | bonobo | KB3781 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all candidates | 0.182 | 0.533 | 0.155 | 0.052 | 0.086 | 0.087 | 0.081 | 0.005 | 0.006 | 0.006 | 0.009 | 0.006 |
| T05 | 0.151 / 0.238 | 0.335 / 0.632 | 0.098 / 0.244 | 0.185 / 0.046 (0.15) | **0.150 / 0.074** | 0.129 / 0.080 | **0.181 / 0.066** | **0.027 / 0.004** | **0.031 / 0.005** | **0.021 / 0.004** | **0.026 / 0.007** | **0.023 / 0.005** |
| K50 | 0.129 / 0.389 | 0.474 / 0.684 | 0.111 / 0.329 | **0.098 / 0.040** | 0.116 / 0.072 | 0.103 / 0.079 | 0.111 / 0.069 | **0.010 / 0.004** | **0.012 / 0.004** | **0.010 / 0.004** | 0.013 / 0.007 | **0.010 / 0.004** |
T20, T50 and T20K hold 0-3% of the correct candidates everywhere (a share of 0.2 is reached by 1,097 of 255,804 candidates)
and separate nowhere.

- **T05** separates on 2 of 7 non-primate (heart, cortex; on mouse ES it is four times as precise and holds 15% of the
  correct, under the 20% bar) and 5 of 5 primate. **K50**: 1 of 7 and 4 of 5. Bars were 6 and 4: **FAIL**, and not
  NON-PRIMATE ONLY either. The pattern is the reverse of the one hoped for: the signatures separate on the apes.
- Where they separate the admitted chains are still far too poor to keep: 0.02-0.03 on the apes, 0.13-0.19 on the mouse
  libraries, against `x.gtf` precision of 0.33 (mouse ES) to 0.55 (heart, adrenal).
- On the three simulated sets both signatures pick the WRONG candidates (as row 1840's did). No cutoff repairs it: precision
  falls with the contrast on all three (mouse simulated 0.75 at contrast 0.25-0.5, 0.47 at >= 0.75).
- Expectation (FAIL) held. The part "if anything, the 3' side or K50 on mouse" did not hold as a signature.

**Observation, not registered (post hoc): the internal SIDE carries more than any end signature.** Precision of the
candidates by which side is internal (5' = cut at the start, 3' = ends early, both):
| | fly sim | mouse sim | zebrafish sim | mouse ES | heart | adrenal | cortex | OR6737 | A119b | chimp | bonobo | KB3781 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5' | 0.132 | 0.341 | 0.078 | 0.037 | 0.071 | 0.065 | 0.060 | 0.004 | 0.005 | 0.005 | 0.007 | 0.005 |
| 3' | 0.955 | 0.984 | 0.401 | 0.236 | 0.335 | 0.222 | 0.338 | 0.015 | 0.023 | 0.011 | 0.017 | 0.010 |
| both | 0.091 | 0.563 | 0.350 | 0.055 | 0.068 | 0.092 | 0.085 | 0.002 | 0.002 | 0.001 | 0.002 | 0.004 |
| 3' candidates, correct / all | 21 / 22 | 1,054 / 1,071 | 71 / 177 | 449 / 1,901 | 124 / 370 | 314 / 1,414 | 230 / 680 | 57 / 3,884 | 301 / 13,086 | 34 / 3,228 | 56 / 3,234 | 70 / 7,240 |
A contained chain that only ends early (no container reaches upstream of it) is 3-6 times as precise as a 5'-cut one on
every non-primate dataset, and nearly always true on fly and mouse simulated. On the apes it is 0.010-0.023: rows 1064 and
1486, which closed this class, were read on A119b alone. This also explains the simulated reversal: there the true candidates
are the early-ending ones, and the 5' signatures are computed on the other side. Not a result: the split was looked at
after the tables. A rule "the collapse spares a chain that is internal at the 3' side only" would need its own registered row
(judged by "pays"; mouse ES would gain 2.5% correct chains at about 3 wrong per correct).

**Decision:** read ends at the collapse give no lever (with rows 1064, 1486, 1715, 1840); no build. Open, the user's call: the
3'-only exemption as a registered row.
