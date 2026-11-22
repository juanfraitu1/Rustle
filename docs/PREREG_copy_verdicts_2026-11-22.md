# PREREG: use the near-tied multimappers with evidence, not all of them and not none (row 1904, 2026-11-22)

User (2026-11-22): the advisor's point is that near-tie multimappers are a gap no other tool uses well: they are ignored
or all taken; "what we want is to prove that in comparison with other tools ours does things more intelligently and uses
real info to decide, like real PSVs or other statistical tests that actually test provenance".

## Seen before writing (row 1903's planted reads; `r1903_o2_family_benchmark/room.py`, legacy `assign` on scenario E)
- Reads by tie class, E / D / S: clear 2,925 / 3,815 / 1,486, primary on the true copy 99.6-99.7%; near tie (rival 1-10 ms
  below) 1,746 / 2,135 / 1,020, primary right 95.9-99.0%; exact tie 3,309 / 5,410 / 1,694, primary right 31-41% (a toss
  among 3-5 copies).
- In an exact tie nothing in the read decides: fewest edits single out the true copy in 19 of 3,309 (E) and a wrong one
  in 46. The legacy `assign` stage (PSV certificate per read) assigns 1 of 3,125 tied reads in E.
- So the choice per read is empty for exact ties; what differs between tools is which COPIES end up with a transcript:
  S, silent copies called, primaries only 24, all near ties 41-42, StringTie 25 (of 63).

## Idea: a provenance test per copy, from the reads that do distinguish two copies
A read with records on copies A and B whose `ms` differ is DECISIVE for the better one: sequencing errors cost the same
on both records, a score difference is a base where the copies differ (a PSV the read carries), no catalogue needed.
Per locus of the assembly (`tools/copy_verdicts.py`, loci and read-to-locus rule of `eq_classes.py`):
- SUPPORTED: >= 2 reads whose best record is here alone.
- REFUTED: not supported, and some SUPPORTED copy A has >= 3 reads decisive for A over it and none decisive the other way
  (the reads that could tell the two apart all say A).
- UNDECIDED: reads reach it, nothing tells it from its partners (identical for these reads): kept, flagged with partners.
The extension: assemble with near-tie secondaries (the mode), then drop the transcripts of REFUTED loci.

## Test (planted reads of row 1903, same BAMs, same scorer; nothing realigned)
Arms: Rustle primaries / secondaries / mode, each with and without the verdict filter, against the six other rows of
row 1903's table (StringTie, FLAIR, IsoQuant, bambu unchanged).
- V1 (fewer false copies): in S the mode + verdicts calls at most 25 silent copies (StringTie's number; the mode alone 42).
- V2 (no expressed copy lost): copy recall of mode + verdicts within one copy of the mode's in E, D and S.
- V3 (the claim against other tools): in S, mode + verdicts has copy precision > every other tool's with copy recall >=
  StringTie's minus one copy.
All three hold = the extension is VALIDATED on planted reads and proposed as O2's default; otherwise said which failed.
Thresholds (2 own reads, 3 decisive) are fixed here and not tuned on these BAMs.

## Predictions (mine)
V1 holds (most silent copies of NPIP, TBC1D3, GOLGA8 differ from an expressed copy and are refuted; identical TSPY copies
stay). V2 holds in E and S; D is the risk (a 20-read copy next to a 600-read one). V3: precision highest, holds.
Limits known in advance: sample = reference; identical copies cannot be judged by any read; a real sample's private
variants would look like decisive reads for no reference copy.

## Result of row 1904 (2026-11-22, `r1904_copy_verdicts/{E,D,S}/score.tsv`): V1 FAILS, V2 holds, V3 FAILS
- V1: S, silent copies called, mode 42 -> mode + verdicts 40 (bar 25). Only 4 of 131 loci are REFUTED; 82 are UNDECIDED.
- V2: copy recall mode / mode + verdicts: E 129 / 128, D 129 / 128, S 68 / 68 (one expressed copy lost in E and in D).
- V3: S copy precision 0.630, below every other tool (StringTie 0.734, FLAIR 0.772).
Why (looked at after): 73 of the 82 undecided loci have decisive reads both ways. A silent copy that matches expressed
copy A1 at one place and expressed copy A2 at another beats A1 on A2's reads and A2 on A1's reads, so no single partner
refutes it; pairwise counts cannot tell it from an expressed copy. The question is whether the reads can be explained
WITHOUT the copy, which is a joint (likelihood) question, not a pairwise one.
Predictions: V1 wrong, V2 right, V3 wrong.

# PREREG: the joint provenance test, a copy must keep reads under the EM (row 1905, 2026-11-22)

Seen: row 1904 above; rows 1874 / 1877 (on the annotated loci of an earlier planted set, half silent: the EM with one
no-locus member per class leaves >= 2 reads on 16 of 63 silent copies and under 2 on 1 of 70 expressed ones).
Rule (`bench/o2_families/em_filter.py`): on the assembly's loci, reads tied between loci are shared by the EM of
`eq_classes.py` (ms, exact band; one `outside` member per class as row 1877); a locus is KEPT when its EM count
(placed + share) is >= 2, otherwise its transcripts are dropped. Nothing else changes. Threshold fixed here.
Same planted BAMs, arms and claims as row 1904: V1 (S: mode + EM calls <= 25 silent copies), V2 (copy recall within one
copy of the mode's in E, D, S), V3 (S: copy precision above every other tool's with recall >= StringTie's minus one).
All three = VALIDATED on planted reads. Prediction (mine): V1 holds (about 16), V3 holds, V2 holds in S and D and is the
risk in E (identical expressed copies with no placed read keep an equal share, so they should stay).

## Result of row 1905 (2026-11-22, `r1905_em_filter/{E,D,S}/score.tsv`): V1 HOLDS, V2 FAILS, V3 FAILS narrowly. NOT VALIDATED
| S (63 silent, 70 expressed) | silent copies called | copy precision | expressed copies found | isoform recall / precision |
|---|---|---|---|---|
| Rustle, primaries | 24 | 0.739 | 68 | 0.719 / 0.575 |
| Rustle, primaries + EM test | 19 | 0.776 | 66 | 0.711 / 0.611 |
| Rustle mode (all near ties) | 42 | 0.618 | 68 | 0.750 / 0.409 |
| Rustle mode + EM test | 20 | 0.767 | 66 | 0.742 / 0.522 |
| StringTie | 25 | 0.734 | 69 | 0.648 / 0.606 |
| FLAIR | 18 | 0.772 | 61 | 0.648 / 0.656 |
| IsoQuant | 22 | 0.750 | 66 | 0.594 / 0.628 |
| bambu | 24 | 0.739 | 68 | 0.719 / 0.503 |
- V1 holds: 42 -> 20 silent copies called (StringTie 25).
- V2 fails: expressed copies found, mode -> mode + EM: E 129 -> 128, S 68 -> 66, D 129 -> 119. In D the EM gives the tied
  reads to the 600-read copy and leaves 10 expressed 20-read copies under 2 reads.
- V3 fails narrowly: S copy precision 0.767 against FLAIR 0.772 (above StringTie 0.734, IsoQuant 0.750, bambu 0.739), and
  copy recall 66 against the bar of 68.
Predictions: V1 right (I said about 16; 20), V2 wrong (D loses 10; I named E as the risk), V3 wrong by 0.005 and on recall.
What stands: taking every near tie doubles the false copies (24 -> 42); ignoring them or tossing gives 18-25; the joint
test brings the mode back to 20 while keeping its isoform recall (0.742 against 0.59-0.72), at the price of 2-10
expressed low-count copies. Not a default; the rule that keeps a low-expressed copy beside a dominant identical-in-part
one is the open problem, and any next rule needs a fresh seed (these three BAMs are now development data).
