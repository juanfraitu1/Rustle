# gffcompare on chr20: ours vs StringTie vs FLAIR, with today's assembly settings

Same substrate and reference as `bench/CHR20_ASSEMBLER_COMPARISON.md` (2026-09-15): T2T-CHM13 human testis
IsoSeq restricted to chr20, reference `chr20_ref.gtf`, all tools at defaults. Re-run because every
assembly change of 2026-09-18/19 postdates that report.

## Result

| tool | query mRNAs | **matching intron chains** | intron-chain Sn / Pr | transcript Sn / Pr |
|---|---|---|---|---|
| **ours — today's settings** (`--assemble-only --read-isoform-k 3`, `RUSTLE_JUNCTION_MAJORITY=1`) | 1,275 | **358** | **8.4** / 33.4 | **7.9** / 28.2 |
| ours — baseline flags | 976 | 345 | 8.0 / 44.6 | 7.6 / 35.6 |
| **StringTie** | 712 | 331 | 7.7 / **47.4** | 7.3 / **47.1** |
| FLAIR | 820 | 264 | 6.2 / 35.1 | 5.8 / 32.3 |

⭐**On "does it find real transcripts" — the count of reference intron chains correctly reconstructed —
ours leads: 358 vs StringTie's 331 (+8.2%) and FLAIR's 264 (+35.6%).** Today's settings add 13 chains over
our own baseline (345 → 358).

**StringTie wins precision decisively** (intron-chain 47.4 vs our 33.4, transcript 47.1 vs 28.2), which is
the expected shape: its network-flow model emits far fewer transcripts (712 vs our 1,275) and is much more
conservative. We are not claiming to beat it there.

**The trade today's flags make is explicit:** +13 true chains for +299 emitted transcripts, i.e. precision
44.6 → 33.4. That is a recall-oriented setting. **For a precision-oriented run, use the baseline flags** —
they still beat StringTie on sensitivity (345 vs 331 chains) at much closer precision (44.6 vs 47.4).

## Validation of `--assemble-only`

The baseline arm run through `--assemble-only` reproduces the 2026-09-15 published numbers **exactly** —
976 mRNAs in 456 loci, 345 matching intron chains, base 11.6/69.3, intron-chain 8.0/44.6, transcript
7.6/35.6. The mode changes no assembly output.

## Known precision cost, from the 09-15 SQANTI3 breakdown

Our extra transcripts are dominated by **incomplete-splice_match (ISM): 217 vs StringTie's 79** — i.e.
5′-truncated fragments of real transcripts, which §6p4 localised precisely (TBC1D3 5′UTR only 75.3%
covered while its 3′UTR is 100%). That is the same 5′-truncation signature, and it is what the precision
gap is made of. It is a fragment problem, not a false-junction problem: intron-level precision stays high
at 85.3-86.1%.

⚠ chr20 is an ORDINARY chromosome chosen specifically because it was not tuned for multi-copy families;
this measures the assembler, not the project's multi-copy contribution.
