# PREREG r1772: are the siblings of an expressed family member reached by its multimapping reads? (2026-11-17)

Descriptive measurement, written before any count was looked at. Workspace `r1772_sibling_reach`,
script `bench/famcmp/sibling_reach.py`.

## Question (advisor, as relayed by the user)
"If one member of a family is found, the others should be covered by multimapping reads, so their loci can be found
whether or not they are expressed." Expression of the sibling is deliberately not part of the question: the unit is the
sibling LOCUS and whether reads of the other members leave alignment records on it.

## Definitions
- Truth families: Compara paralog pairs at the pairscore default level (simiiformes) among protein-coding genes,
  connected components (`famtruth`). Pseudogene siblings are not in this truth and are not counted.
- Per gene, records overlapping its exons on its strand (the r1521 rule); supplementary records skipped.
- Source molecules of a sibling = molecules with a PRIMARY record on any OTHER member of its family.
- A sibling is judged when its family has >= 10 source molecules (at fewer, absence of records is depth).
- Reached(K, band) = >= K distinct source molecules with a SECONDARY record on the sibling, the record's AS at least
  band x the same molecule's primary AS. Bands: any (what minimap2 `-N 50 -p 0.1` printed), 0.5, 0.98 (the seeding pool).
- Headline cell: K = 2, band any.
- Identity of a sibling = best identity of its spliced exon sequence to a family member with >= 2 primary molecules
  (minimap2 `-c -x asm20 -N 500 -p 0.05`, shorter sequence covered >= 0.5; seqscore's rule). Bins: >= 0.98, 0.95-0.98,
  0.90-0.95, 0.85-0.90, < 0.85, no alignment.
- Split by the sibling's own primaries: < 2 molecules ("no own expression by the P rule") and >= 2.
- Also recorded, not judged: all secondaries at the sibling from any molecule; whether our default `PREFIX.gtf` has a
  model on the sibling's exons and strand.

## Samples
Amended 2026-11-17 before any full result was read (only a 15-family gorilla smoke test had been seen): the user dropped
A119b as too complex a sample; its job was cancelled unread. Samples are now apes only:
- OR6737 (gorilla testis) and KB3781 (gorilla fibroblast, the T2T individual), `DENOVO_GGO_ENS_MAP` truth of r1611;
  models `artifact_refresh_2026-11-16/g.def/or.gtf` and `apetools/GGO_KB3781_fib_mm/ours/ours.gtf`.
- PTR (chimpanzee testis), RefSeq annotation with symbol-inherited human ids only (no protein-hit map exists for
  chimpanzee, so its truth is the smaller one); models `apetools/PTR_mm/ours/ours.gtf`.
All are spent samples; nothing is tuned here.

## Reading rule (per identity bin, siblings with < 2 own primaries, n >= 20 to read a bin)
- reach(K = 2, any) >= 0.80: the advisor's statement HOLDS in that bin;
- < 0.50: it FAILS in that bin;
- between: partial.

## Predictions
Holds at >= 0.98; partial at 0.95-0.98; fails below 0.90 and at "no alignment". Most judged siblings fall in the bins
where it fails.

## Limits stated beforehand
- The BAM holds at most 50 secondaries per read at `-p 0.1`; reach is of what this alignment printed.
- AS is used for the bands; minimap2 ranks secondaries on `ms`.
- A record on a sibling shows sequence similarity, not expression and not origin (r1671).
