# PREREG: where minimap2's primary is wrong on reads of known origin, and can the alignments' own context overrule it (row 1922, 2026-11-23)

User: "focus on understanding how minimap2 operates and how can we overwrite wrong calls, for example minimap might be
blind to substitutions or indels but we have the contexts".

## How minimap2 decides (what the measurement is built on)
- Seeds -> chains -> base-level alignment per chain. With `splice:hq` a match is +1, a mismatch -4 (5 points lost), a
  gap -6 to open and -1 per base (an indel of one base about 7 lost), a non-canonical splice site is penalised.
- The primary is the chain with the highest `ms` (the score of the best-scoring segment, introns free); `AS` (introns
  and splice sites charged) is reported but does not rank. Exact `ms` ties: the primary is picked by a hash of the read
  name, a coin toss.
- So minimap2 is NOT blind to substitutions and indels inside the aligned part: the score gap between two copies is
  the count of positions where the read fits one copy and not the other, weighted 5 / about 7. Where it can be blind:
  (a) the ends: bases it soft-clips are not scored; (b) introns and splice sites: free on `ms`; (c) exact ties;
  (d) a candidate it never reports (chain filter `-p`, `-N`); (e) it weighs a sequencing error like a true difference.
  Closed before (do not repeat): any re-placement on real reads with truth lost to the primary (rows 1642-1651, 1777,
  1781, 1783, 1784, 1844); the near-tie census (row 1829) put reference PSVs at 82% of the decidable.

## Measurement (planted human families, ragged ends; every read's copy of origin is in its name)
Development: seeds 1920-1924, scenarios E and S (`r1918_rerun_ragged/lab`). Check: seeds 1927-1929, S
(`r1921_weighted_em/c`), opened only if a rule qualifies. A placement is RIGHT when it lies on a copy of the origin's
unit (copies with identical planted isoforms cannot be told apart by any method and are not errors).
1. Primary right / wrong. Wrong ones by cause: origin tied on `ms` (coin toss); origin reported at a lower `ms` (error
   in the read outweighs the copies' differences); origin not reported.
2. Among the `ms`-tied candidates of a read, does anything the aligner already wrote separate them: `AS`, soft-clipped
   bases, `NM`.

## Rules (fixed now; each replaces the primary only among candidates tied with it on `ms`)
R1 highest `AS`; R2 fewest soft-clipped bases; R3 lowest `NM`; R4 fewest clipped bases, then highest `AS`.
A rule QUALIFIES on development when it fixes at least twice what it breaks AND removes at least 10% of the primary's
wrong calls. A qualifying rule is then scored once on the check seeds (same two conditions).

## Limits stated beforehand
Planted reads carry uniform errors (0.1% substitutions, 0.2% insertions, 0.2% deletions), not the HiFi profile
(homopolymer indels), and the sample equals the reference. So this row can show blindness of kinds (a)-(d); kind (e)
with real error context needs real reads with truth (HG002 on its own haplotypes) and is not tested here.

## Predictions (mine)
The primary is right on more than 95% of reads. Most wrong calls are coin tosses between copies that differ outside
the read or at its clipped ends. R2 / R4 (clips) fix some; I expect none of the four to remove 10% of the wrong calls
with a 2:1 ratio, because for most tied candidates the alignments are identical in every tag.

## Result (2026-11-23, workspace `r1922_minimap_wrong_calls`, `wrong.py`; ten development sets, 60,900 reads)

| | reads | share |
|---|---|---|
| primary on the origin's unit | 53,220 | 87.4% |
| primary wrong | 7,680 | 12.6% |
| - origin tied with the primary on `ms` (coin toss) | 7,409 | 96.5% of the wrong |
| - origin reported at a lower `ms` (gap 1-5: 121; 6-10: 45; 11+: 105) | 271 | 3.5% |
| - origin not reported | 0 | |

Reads with two or more `ms`-tied candidates on different units: 16,282 (26.7% of reads); the primary is right on 8,859
of them (54%: a coin toss). What separates the tied candidates in the aligner's own output: `AS` in 573 reads (3.5%),
soft-clipped bases in 6, `NM` in 19. In 96% of the ties the candidates are identical in score, edit distance and
clipping: the read fits both copies equally.

| rule (among `ms`-tied candidates) | moved | fixed | broke | share of the wrong calls removed |
|---|---|---|---|---|
| R1 highest `AS` | 88 | 75 | 11 | 0.8% |
| R2 fewest clipped bases | 2 | 0 | 1 | 0 |
| R3 lowest `NM` | 4 | 0 | 2 | 0 |
| R4 clips, then `AS` | 90 | 75 | 12 | 0.8% |

- **No rule qualifies** (R1 / R4 meet the 2:1 ratio, 7:1, but remove under 1% of the wrong calls against the 10% bar).
  The check seeds were not opened.
- **My prediction was half wrong:** I expected the primary right on more than 95% of reads; it is 87.4%. The other
  half held: the tied candidates carry nothing that overrules the primary.
- **Reading.** minimap2 is not blind inside the read: no origin went unreported, and only 0.4% of reads have the origin
  at a lower score (an error in the read outweighing the copies' differences; 121 of them by one substitution, the case
  the weighted EM of row 1921 addresses). Its wrong calls are almost all exact ties in which the read itself holds no
  information: it fits two copies identically, so no re-scoring, clip rescue or context of THAT read can decide. Units
  are sets of isoforms, so two copies of different units can share one isoform exactly; reads of that isoform are among
  these ties.
- What can decide such a read is evidence outside it: the other reads of the locus (abundance, which is the EM of the
  `o2` labels) or nothing. That is why Rustle labels such copies `shared` and why per-read assignment was closed.
- R1 (highest `AS` among `ms` ties, i.e. the spliced / canonical candidate) is real but small: 75 fixed, 11 broken in
  60,900 reads. Not built.
- Not tested, as stated: sequencing-error context of real HiFi reads (homopolymer indels weighed like true
  differences). It can only bear on the 0.4% class, and needs real reads with truth.
