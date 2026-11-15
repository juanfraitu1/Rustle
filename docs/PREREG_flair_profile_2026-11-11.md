# PREREG r1638: a FLAIR-like chain rule from read support with merged subset reads (2026-11-11)

Why: r1634 found no arm moving toward FLAIR. Seen before this registration (so the gorilla figures below are post hoc):
every FLAIR chain on gorilla has at least one full-chain primary read (`r1633/or.reads.tsv`), half of them exactly one, and
9,616 of r1637's 15,464 merged one-read chains are FLAIR's. FLAIR counts partial reads as support; our arms did not.

Rule family (prototype `bench/famcmp/subset_merge.py`, no Rust change): emit every distinct primary read chain with
full-chain reads + merged subset reads >= S, S in 2..6, merged = unique container (U) or any maximal container (A); each
alone and in union with the default assembly. 20 arms.

Metric = r1634's: Jaccard with FLAIR, share of FLAIR's chains reached, share of the arm inside FLAIR; price = RefSeq-correct
chains and precision against the default.
- "Almost exactly" stays 0.80 / 0.90 / 0.90, on gorilla and on A119b.
- Stop rule: no arm reaches Jaccard 0.60 on gorilla -> A119b is not run; the closest arm is reported with its price.
- Otherwise the arm with the highest gorilla Jaccard is frozen and run once on A119b (unseen for this rule).

## Result (2026-11-11, `r1638/or.profile.txt`): STOP at the registered rule (best Jaccard 0.456 < 0.60); A119b not run
Gorilla OR6737; FLAIR has 138,125 chains (24,293 correct, 0.176), all of them some primary read's chain; 549,019 distinct
read chains. Default: 73,348 chains, 25,794 correct (0.352), Jaccard 0.309, reached 0.361, inside 0.680.
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| U2 | 180,247 | 28,212 | 0.157 | 0.574 | 0.440 | 0.331 |
| U3 + default | 134,976 | 27,131 | 0.201 | 0.487 | 0.498 | 0.327 |
| A2 | 327,326 | 32,321 | 0.099 | 0.912 | 0.385 | 0.371 |
| A3 | 264,632 | 29,652 | 0.112 | 0.863 | 0.450 | 0.420 |
| A5 + default | 225,132 | 29,856 | 0.133 | 0.814 | 0.499 | 0.448 |
| A6 + default (closest) | 210,726 | 29,512 | 0.140 | 0.791 | 0.518 | 0.456 |
Counting merged subset reads as support is the first lever that moves toward FLAIR (0.309 -> 0.456, 79% of its chains
reached, +3,718 correct chains) and its price is precision 0.352 -> 0.140, below FLAIR's own 0.176. Half the arm stays
outside FLAIR: read support alone does not say which of the 549,019 read chains FLAIR keeps (its realignment of reads to
the first-pass isoforms does). Full 20-arm table in the result file. Not ported.

## r1640 (registered 2026-11-12, user chose "Prototype FLAIR's realignment"): chain-level emulation of FLAIR collapse
FLAIR 3.0.1 as the benchmark ran it (`collapse --trust_ends`, no annotation, support 3), read from its source:
first-pass isoform = every distinct read chain; `filter_collapsed_isoforms` default drops a chain that is a contiguous
subset of a longer one with a matching end exon (always when < 4 exons, else unless its reads > 3 and > 1.2 x the
container's); reads are realigned to the first-pass sequences and each goes to the best alignment score, ties to the
shortest transcript; isoforms with >= 3 reads are printed.
Emulator (`bench/famcmp/flair_like.py`, no minimap2): same first pass and subset filter on chains (ends = median read
ends); a read goes to its own chain if kept, else to the shortest kept chain that contains it (ends within 10 bp), else
nowhere; print chains with >= 3 reads. E0 = the same without the subset filter (ablation). Multi-exon chains only.
- PASS: E1 meets Jaccard >= 0.80, reached >= 0.90, inside >= 0.90 with FLAIR on gorilla, then on A119b (run once, frozen).
- Stop rule: E1 Jaccard < 0.60 on gorilla -> A119b not run; real minimap2 realignment is the only step left and equals
  running FLAIR itself, so the FLAIR side closes with the price table.
- Between 0.60 and 0.80: the residual is split by cause (ends, realignment ties, reads with no container) and reported;
  one amendment may be registered from it before A119b.
Price in every case: RefSeq-correct chains and precision against the default (25,794 at 0.352 on gorilla).

### r1640 result, E1 as registered (2026-11-12, `r1640/or.txt`): STOP (Jaccard 0.516 < 0.60)
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| default | 73,348 | 25,795 | 0.352 | 0.361 | 0.680 | 0.309 |
| E0 (no subset filter) | 105,980 | 23,269 | 0.220 | 0.380 | 0.495 | 0.273 |
| E1 | 91,038 | 23,302 | 0.256 | 0.564 | 0.856 | 0.516 |
| E1 + default | 115,094 | 27,296 | 0.237 | 0.603 | 0.724 | 0.490 |
What E1 prints is FLAIR's (0.856 inside); it misses 60,166 FLAIR chains, 59,516 of them kept in the first pass with 1 or 2
own reads and too few assigned reads.

### Amendment E2 (registered 2026-11-12 after seeing E1's residual; gorilla is post hoc for it, bars unchanged)
E1 left out a step FLAIR runs: single-exon reads are realigned too, and a single-exon first-pass isoform inside an exon
of a multi-exon one is dropped, so those reads land on multi-exon isoforms (ties to the shortest transcript). E2 = E1 + every
unspliced primary read assigned to the shortest kept chain with an exon that contains it (ends within 10 bp).
Same bars (0.80 / 0.90 / 0.90); E2 >= 0.60 on gorilla -> frozen and run once on A119b (unseen); below -> FLAIR side closed.

### E2 result (2026-11-12): STOP (Jaccard 0.547 < 0.60)
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| E2 | 101,072 | 23,545 | 0.233 | 0.612 | 0.837 | 0.547 |
| E2 + default | 124,520 | 27,445 | 0.220 | 0.648 | 0.719 | 0.517 |
Still missed: 53,544 FLAIR chains, 52,894 of them first-pass chains with 1-2 own reads. Chain containment does not find
the reads FLAIR's realignment gives them.

### E3 (registered 2026-11-12 before running): the realignment itself, the arm the user chose
E1's first pass (kept chains, median ends) written as transcript sequences; every primary read realigned with FLAIR's own
command (`minimap2 -a -N 4`, default preset) and assigned by FLAIR's rule as read from `count_sam_transcripts.py`: drop
an alignment with an indel > 25 bp, then best alignment score, least clipping, shortest transcript; print chains with
>= 3 reads. Same bars (0.80 / 0.90 / 0.90 on gorilla, then A119b once). Below 0.60 on gorilla -> FLAIR side closed.
This arm costs a second minimap2 pass over all reads; its run time is reported with the result.

### E3 result, gorilla (2026-11-12, `r1640/or.e3.txt`): above the 0.60 stop, below "almost exactly"
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| first pass (no read cut) | 386,385 | 33,218 | 0.086 | 0.995 | 0.356 | 0.355 |
| E3 | 154,769 | 25,097 | 0.162 | 0.896 | 0.800 | 0.732 |
| E3 + default | 176,704 | 28,872 | 0.163 | 0.919 | 0.718 | 0.675 |
FLAIR itself: 138,125 chains, 24,293 correct, 0.176. Run time 1,927 s on 64 threads (minimap2 1,857 s, 16 GB).
Bars: Jaccard 0.732 (< 0.80), reached 0.896 (< 0.90), inside 0.800 (< 0.90): not "almost exactly", but past the stop rule,
so E3 is frozen and run once on A119b.

### E3 result, A119b (2026-11-12, `r1640/a.e3.txt`, run once, frozen): FAIL "almost exactly" (Jaccard 0.636, inside 0.674)
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| first pass | 1,645,470 | 49,218 | 0.030 | 0.993 | 0.378 | 0.377 |
| E3 | 855,000 | 41,822 | 0.049 | 0.919 | 0.674 | 0.636 |
| E3 + default | 942,243 | 47,649 | 0.051 | 0.931 | 0.620 | 0.592 |
Default: 219,517 chains, 38,715 correct, 0.176. FLAIR itself: 626,983 chains, 39,230 correct, 0.063. Run time 9,866 s on 64
threads, 33 GB. E3 prints 228,000 more chains than FLAIR here; the first pass differs from FLAIR's (ends, and whatever
FLAIR's own collapse drops), not checked further. Verdict: FLAIR's selector is first pass + subset filter + realignment;
copying it reaches Jaccard 0.732 / 0.636 and costs precision 0.352 -> 0.162 (gorilla) and 0.176 -> 0.049 (A119b), for
-698 / +3,107 correct chains. Not ported.

### E3 residual, gorilla (`r1640/or.e3.counts.tsv`) and amendment E4 (registered 2026-11-12 before running)
E3's 31,013 chains outside FLAIR: 14,522 with exactly 3 assigned reads, 5,865 with 4, 2,976 with 5, 7,650 with 6+.
FLAIR chains E3 misses: 8,783 with 2 assigned reads, 4,495 with 1, 441 with 0, 650 dropped as subsets. The error sits at
the 3-read cut and in the subset filter, where E3's first pass differs from FLAIR's (`collapse_isoforms_precise.py`):
- FLAIR writes up to 2 starts x 2 ends per chain (best site by distance-weighted read count in a 100 bp window; a second
  site only if >= 4 reads sit exactly on it and it takes at least the average share), each its own isoform that needs 3
  reads on its own; E3 wrote one isoform per chain at the median ends.
- the support the subset filter compares is the read count exactly at the chosen end (max of start and end), not the
  chain's read count.
E4 = E3 with both copied. Post hoc for gorilla; same bars (0.80 / 0.90 / 0.90). Passing on gorilla -> A119b once.
Not passing -> the FLAIR side closes with E3/E4 as the price table; no further amendment.

### E4 result, gorilla (2026-11-12, `r1640/or4.e4.txt`; the file labels the arm "E3"): closest, still short of the bars
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| E4 | 153,772 | 24,705 | 0.161 | 0.934 | 0.839 | 0.792 |
| E4 + default | 176,232 | 28,855 | 0.164 | 0.940 | 0.737 | 0.704 |
Jaccard 0.792 (< 0.80), reached 0.934 (>= 0.90), inside 0.839 (< 0.90). As registered: A119b not run for E4, no further
amendment, FLAIR side closed. Price of the closest FLAIR match on gorilla: 25,795 correct at 0.352 -> 24,705 at 0.161.

### Amendment E5 (2026-11-12; the user chose "One more FLAIR amendment" over the no-further-amendment rule; registered before running)
E4's 24,759 chains outside FLAIR: 13,363 with exactly 3 assigned reads, 4,922 with 4, 2,259 with 5, 4,215 with 6+.
Remaining known difference from FLAIR's source: its first pass also holds single-exon isoforms (unspliced reads grouped
by overlap, one best start/end per group; dropped when inside an exon, +-100 bp, of any multi-exon first-pass isoform, or
inside a larger single-exon one unless reads > 3 and > 1.2 x its reads). They compete for reads in the realignment; E4
had none, so unspliced and partial reads could only add to spliced isoforms.
E5 = E4 + these single-exon candidates (not scored themselves). Same bars (0.80 / 0.90 / 0.90). Gorilla post hoc;
passing -> A119b once. Not passing -> closed; any further step is the user's call.

### E5 result, gorilla (2026-11-12, `r1640/or5.e5.txt`; the file labels the arm "E3"): PASS all three bars (post hoc sample)
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| E5 | 137,784 | 24,257 | 0.176 | 0.931 | 0.933 | 0.873 |
| E5 + default | 161,297 | 28,545 | 0.177 | 0.937 | 0.803 | 0.762 |
93,264 single-exon candidates. FLAIR itself: 138,125 chains, 24,293 correct, 0.176. Frozen; A119b run once next.

### E5 result, A119b (2026-11-12, `r1640/a5.e5.txt`; run once with the rule frozen on gorilla; A119b had been scored under E3 before E4/E5 were written): PASS all three bars
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| E5 | 630,228 | 39,213 | 0.062 | 0.931 | 0.926 | 0.866 |
| E5 + default | 729,771 | 46,619 | 0.064 | 0.935 | 0.803 | 0.760 |
FLAIR itself: 626,983 chains, 39,230 correct, 0.063. Default: 219,517 / 38,715 / 0.176. 10,242 s, second minimap2 pass.
Against the default, scored separately: A119b E5 has 7,904 correct chains the default lacks and lacks 7,406 the default has (shared 31,309); gorilla 2,750 gained, 4,288 lost (shared 21,507). Union of both: A119b 46,619 correct at 0.064, gorilla 28,545 at 0.177. Chains 2.9x the default on A119b (precision 0.176 -> 0.062). Scope: spliced intron chains only; single-exon isoforms and ends not scored. Bench script only, not ported.

### Port (2026-11-12)
PORTED opt-in 2026-11-12 (user chose 'Port FLAIR-like mode opt-in'): driver stage `flairlike` -> PREFIX.flairlike.gtf (bench/famcmp/flair_like.py modes firstpass + gtf; never part of `all`). Gorilla parity run through the driver (`r1640/port`): 137,897 spliced chains, 128,633 shared, Jaccard 0.873 (bench 137,784 / 128,586 / 0.873; same candidate set as the bench run, different candidate order moved ~0.1% of chains; order then fixed (imap), two first passes byte-identical); 154,317 isoforms with single-exon ones; 1,953 s.
