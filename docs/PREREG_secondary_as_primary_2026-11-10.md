# PREREG r1613: does StringTie recover the multimapper members when secondary alignments are relabelled primary? (2026-11-10)

Advisor's objection (2026-11-10): StringTie and FLAIR score 0 on family members with no primary read only because they skip
secondary alignments; clear the secondary flag and they would find them, so the result is not about the assembler. Written before
any run. Quick test bed first (gorilla OR6737, `GGO_mm.bam`); A119b only after reading this one. FLAIR is not in this row (it
realigns the reads to its own candidates; a separate row if this one leaves the question open).

## Arms (StringTie 3.0.3 `-L`, no annotation, installed here; all three run here so the version is the same)
- `base`: the BAM as it is.
- `all`: flag 0x100 cleared on every secondary record.
- `tied`: flag 0x100 cleared only on secondary records with AS >= 0.98 x the molecule's best AS (`or.molecules.tsv`).
- Rustle: `~/runs/gorilla_OR6737/or.gtf` (2026-10-29), unchanged.

## Measures (`bench/famcmp/flip_score.py`; panel members and classes from `r1611/or.genes.tsv`, exons `r1611/ggo.genes.tsv`, RefSeq GTF)
- M1 loose: members with a model on their strand overlapping their exons, by read class (the plotted measure).
- M2 exact: members with such a model whose intron chain equals a RefSeq chain (multi-exon), by read class.
- P family: of the distinct multi-exon chains overlapping any panel member's exons, the share equal to a RefSeq chain.
- P genome: `chains.py` precision genome-wide; transcripts in total.

## Reading (registered)
The objection HOLDS (the flag is all that separates the tools on these members) if some flipped arm has, versus Rustle:
M1 no-primary-read >= Rustle's, M2 no-primary-read >= Rustle's, and P family >= Rustle's - 0.02, with P genome not below the
`base` arm by more than 0.02. Otherwise it does NOT hold as stated, and the failing measure says why (recall not reached, or
reached at a precision cost). Either way the plot's caption is corrected to say what the control showed.

## Outcome (2026-11-10, `r1613/`, job 4456503, 1,642 s)
| arm | transcripts | no-primary M1 / M2 (of 110) | ambiguous M1 / M2 (of 49) | ordinary M1 / M2 (of 552) | P family | P genome |
|---|---|---|---|---|---|---|
| Rustle | 73,912 | 72 / 35 | 37 / 30 | 392 / 265 | 0.325 | 0.352 |
| StringTie base | 68,851 | 0 / 0 | 26 / 21 | 439 / 267 | 0.295 | 0.350 |
| StringTie all secondaries flipped | 117,863 | 109 / 58 | 46 / 37 | 524 / 268 | 0.158 | 0.222 |
| StringTie tied secondaries flipped | 70,763 | 75 / 34 | 40 / 29 | 455 / 271 | 0.289 | 0.344 |

**As registered: does NOT hold** (`tied`: M2 34 < 35 and P family 0.289 < 0.305; `all`: P family 0.158, P genome -0.128).
**In substance the advisor is right about the members' recovery:** with the tie filter (AS >= 0.98 x best, our own rule)
StringTie models the no-primary members as often as Rustle (75 vs 72 loosely, 34 vs 35 exactly) at a small precision cost in
family loci (0.289 vs 0.325) and none genome-wide. Flipping every secondary recovers more (109, 58 exact) but halves precision
(49k extra transcripts). So the plotted gap is the input rule (tied secondaries are read), not the assembly. What the trick does
not give: each tied molecule is counted at every copy, and nothing says which copy is expressed; that is the assign step and
the tie groups, and it is where the claim has to rest. Not run: A119b; FLAIR.
