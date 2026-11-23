# PREREG: are the unconfirmed sequences explained by a haplotype hit that scores lower (row 1930, 2026-11-23)

Row 1929: on KB3781, 47 of 108 recurrent non-reference sequences among tied copies have variants that are not in the
animal's maternal or paternal assembly at the read's BEST hit. Two explanations were left: the assemblies lack the
sequence there, or the read's true origin is another copy that minimap2 scores lower. No genomic reads of the animal
are available (user), so the second is tested: it only needs the assemblies.

## Measurement (`r1930_lower_hits/explain.py`)
The reads of every kept sequence are aligned again to both haplotypes with the candidate list opened wide
(`-N 200 -p 0.2`, otherwise the usual recipe). For each read of a non-reference sequence its variant positions (the
substitutions of its sequence's pattern, as read positions) are looked up in EVERY hit on either haplotype.
- A read is EXPLAINED BY ONE PLACE when some single hit matches >= 0.8 of its variant positions.
- A sequence is explained when >= 0.8 of its reads are.
Reported beside it: the cost (edit distance + clipped bases) of the explaining hit minus the best hit's, and the
weaker per-site figure (a variant position matching at ANY hit), which does not make an origin.

## Claim
- E1: at least half of the sequences that row 1929 found NOT in the DNA at the best hit are explained by one place at
  a lower-scoring hit. Yes = the second explanation: the sequences are real copies of the animal that the aligner
  ranks second, and the tool's splits are mostly right after all. No = they are in neither assembly anywhere the read
  aligns, which leaves the assemblies (or something in the library not thought of) and cannot be settled here.

## Prediction (mine)
E1 fails. A hit that matches all of a read's variants and still scores lower must lose more elsewhere, and reads tied
between copies are full-length over nearly identical sequence, so I expect few such hits. My three predictions on
row 1929 were all wrong, so this one deserves little weight.

## Result (job 4497069, `explain.py`; KB3781, 107 non-reference sequences with a substitution and reads on a haplotype)

| sequences | n |
|---|---|
| variants in the DNA at the read's best hit | 47 |
| not at the best hit, explained by ONE lower-scoring place | 5 |
| not at the best hit, partly explained (0.2-0.8 of the reads) | 8 |
| not at the best hit, explained NOWHERE on either haplotype | 47 |

Read-sites: 8,495; matching at the best hit 6,072 (0.71); matching at some hit 6,692 (0.79). Hits per read: median 4,
up to 329, so the candidate list was open.
- **E1 FAILS:** 5 of the 60 sequences not confirmed at the best hit are explained by a lower-scoring place (the claim
  needed 30). The second explanation accounts for a twelfth of them.
- Prediction right this time.
- **What stands for KB3781:** of 107 recurrent non-reference sequences among tied copies, 52 are in the animal's
  assembled DNA (47 where the read fits best, 5 at another copy the aligner ranks lower), 8 partly, and 47 are in
  neither haplotype assembly at any place the read aligns. Those 47 are independent molecules (row 1929) carrying
  substitutions with the spectrum of related copies. With mosaics, duplicates, sequencing-error context and a
  lower-ranked origin all tested and removed, what is left is a sequence the two assemblies do not contain, or a
  property of the library not thought of. It cannot be settled without the animal's genomic reads.
