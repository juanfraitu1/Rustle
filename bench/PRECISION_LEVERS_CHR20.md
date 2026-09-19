# Getting more transcripts, or better precision — measured on chr20; and what `-R`/`-Q` actually do

## The two levers, on one table (gffcompare vs `chr20_ref.gtf`, NO `-R`/`-Q`, shared denominators)

| arm | query mRNAs | intron-chain Sn / **Pr** | transcript Sn / Pr | **matching intron chains** |
|---|---|---|---|---|
| **MAX RECALL** — `--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1` | 1,275 | **8.4** / 33.4 | 7.9 / 28.2 | **358** |
| ours, baseline flags | 976 | 8.0 / 44.6 | 7.6 / 35.6 | 345 |
| balanced — recall flags **+ ISM collapse** | 856 | 7.7 / 47.3 | 7.3 / 38.7 | 329 |
| **MAX PRECISION** — baseline **+ ISM collapse** | 784 | 7.6 / **51.9** | 7.1 / 41.6 | 324 |
| StringTie | 712 | 7.7 / 47.4 | 7.3 / **47.1** | 331 |
| FLAIR | 820 | 6.2 / 35.1 | 5.8 / 32.3 | 264 |

### More transcripts
`--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1` gives **358 matching intron chains, the most of any arm
or tool** (StringTie 331, FLAIR 264). The cost is precision: 44.6 → 33.4, because it emits 1,275 transcripts.

### Better precision — collapse the ISM fragments
`bench/ism_collapse.py` drops any transcript whose intron chain is a **contiguous sub-chain** of another
transcript's chain on the same contig/strand — i.e. a 5′-truncated fragment of something longer we already
emit. (Single-exon transcripts are dropped only when contained in a multi-exon span.) It is a pure GTF
post-filter; nothing in the assembler changes.

⭐**On the baseline arm it takes intron-chain precision 44.6 → 51.9 — the best of any arm OR tool,
StringTie included (47.4) — at essentially unchanged sensitivity (8.0 → 7.6).** It drops 192 of 976
transcripts and costs 21 matching chains (345 → 324).

This is the fragment problem §6p4 localised: our excess is incomplete-splice_match (217 vs StringTie's 79),
5′-truncated pieces of real transcripts, not wrong junctions — intron-level precision was already 85-86%.

⚠**Transcript-level precision still trails StringTie (41.6 vs 47.1)** even after the collapse, because
transcript-level scoring requires the ENDS to match too, and our 5′/3′ boundaries are less exact. That is a
separate problem from fragment collapse and is not fixed here.

## Do `-R` / `-Q` help? Yes for reading one tool, NO for ranking tools

(The lowercase `-r` is just the reference; the flags meant are uppercase.)
- **`-R`** counts only reference transcripts that overlap the query — removes unexpressed genes from the
  sensitivity denominator.
- **`-Q`** counts only query transcripts that overlap the reference — removes novel/intergenic predictions
  from the precision denominator.

With `-R -Q` every tool's sensitivity roughly doubles (ours 8.4 → 19.2, StringTie 7.7 → 16.6) because chr20's
unexpressed genes leave the denominator. **That is a much more honest read of "how well did this tool do on
what is actually expressed."**

⚠**But both flags make the denominator TOOL-DEPENDENT** — each tool gets a different reference subset (`-R`)
and a different query subset (`-Q`) — so **cross-tool Sn/Pr under `-R -Q` is not like-for-like.** Use them
to characterise one tool; use the shared-denominator run (or the absolute **matching intron chain count**)
to compare tools.

**Recommendation:** report the `-R -Q` numbers for "how good is our assembly on expressed genes", and the
matching-intron-chain count for "how do we compare to StringTie/FLAIR". Both are in this file.
