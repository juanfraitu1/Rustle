# PREREG: admit to the seeding pool the alignments that tie on `ms` and lose on AS only by their introns (row 1953, 2026-11-25)

Finding (gorilla OR6737, `docs/PREREG_support_tier_2026-11-25.md`, last sections): 207 near-tied reads sit on intronless
copies of ribosomal-protein genes and tie EXACTLY on `ms` with a spliced alignment at the parent gene; on AS the copy is
ahead by 24 per intron, and all 207 parent alignments fall outside the assembler's pool (AS >= 0.98 x best). Elsewhere
the AS pool and the `ms` band agree (91%). A pool on `ms` alone was tested and not adopted (rows 1942, 1942b).

## The rule (emulated, no change to the assembler)
`bench/o2_families/pool_bam.py --rule asintron:0.98:10:30`: a secondary alignment is kept when AS >= 0.98 x the read's best
AS (today's pool), OR when its `ms` is within 10 of the read's best `ms` and AS >= 0.98 x best AS - 30 x (its introns).
The filtered BAM is assembled with `--seed-pool all`, as in row 1942 (whose control showed the emulation reproduces
today's pool exactly).

## Data
- Real: gorilla OR6737 (the sample the class was described on; so this is its development sample) and, unseen for this
  rule, gorilla KB3781.
- Planted human named families, the six 60-read sets of row 1952 (seeds 1967-1969, E and S): families of multi-exon
  copies, where the rule should change almost nothing.

## Claims
- **H1, targeted:** on OR6737 the rule admits at most 5% more secondary alignments than today's pool, and the parent
  alignment of at least 150 of the 207 reads of the class is among them.
- **H2, no harm on real samples:** on OR6737 and on KB3781, at the named families the mode's annotated copies are not
  lower and its junction precision is not lower by more than 0.01; and at least 0.99 of the base assembly's transcripts
  (by exon structure) are shared with today's.
- **H3, no harm on planted reads:** in each of the six sets the mode's exact isoforms are not lower by more than 1 and its
  transcripts on planted copies not higher by more than 2.
- Reported: the transcripts only the new arm has (number, reads, the annotated genes they sit on), and for the genes the
  class ties with (RPL21, RPL12, RPL11, RPL9, MICU2) the reads of their transcripts before and after.

## Predictions (mine)
H1, H2, H3 hold. The new arm adds fewer than 100 transcripts on OR6737; the parents' read counts rise slightly. There is
no truth for whether those reads belong to the parent or to an expressed copy: this test can show that the rule is
narrow and harmless, not that it is right. A default is the user's call and would need the rule inside the assembler.

## Result (job 4504876)
- **H1 HOLDS.** OR6737: 111,552 secondary alignments admitted against today's 110,012 (+1,540, +1.4%; bar 5%), and the
  parent alignment of all 207 reads of the class is among them (today: none).
- **H2 HOLDS on both gorillas.** Named families, the mode: OR6737 37 = 37 annotated copies at junction precision 0.67 =
  0.67; KB3781 28 = 28 at 0.63 = 0.63. Base assembly shared with today's: 0.9998 (OR6737: 11 transcripts lost, 16 new)
  and 0.9999 (KB3781: 5 lost, 8 new).
- **H3 FAILS, narrowly.** In all six planted sets the exact isoforms are unchanged (197, 103, 197, 105, 197, 104), but the
  transcripts on planted copies rise by 3, 2, 3, 3, 6, 3 (bar: at most 2): the extra alignments seed two to six chains per
  set that are no planted isoform, and nothing correct. My prediction (almost no change) was wrong in direction.
- **What the rule changes on OR6737:** the read counts of 42 transcripts rise by 770 reads in total (RPL21 +73, RPL12 +39,
  MICU2 +19, RPL11 +14, and a group of `LOC1295297xx` genes with +147 to +169); 16 spliced transcripts appear (reads 157,
  29, 17, then under 10) and 11 disappear.
- **Reading.** The rule is narrow (+1.4% alignments, 99.98% of the assembly untouched) and reaches exactly the class it
  was aimed at, but what it buys is a few hundred reads on transcripts that already exist, and on planted families it
  adds only wrong chains. Not adopted, and not worth a change inside the assembler. The class is better handled where it
  costs nothing: in the layer summary, by counting a near-tied read as carried when ANY of its tied alignments sits in
  a transcript (55% of the "unexplained" reads on OR6737).
