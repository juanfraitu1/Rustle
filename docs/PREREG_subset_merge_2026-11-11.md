# PREREG r1637: merging reads that are subsets of one longer read (strict containment, unique container; 2026-11-11)

User's question (2026-11-11): "if we have reads that are clearly subsets of one another, can we merge them into one long read
that preserves the information?" Existing modes (`--compat-support`, `--direct-transfrags`, r1618) failed their RefSeq bars
but are looser than the question. Prototype `bench/famcmp/subset_merge.py`, written before any data run.

Rule (primary records only, spliced, per contig; a read has no strand, the default assembly's chains do):
- Read chain = its introns. A read S is **contained** in an observed read chain C (longer than S's) when S's introns are a
  contiguous run of C's, S starts inside or after the start of the C exon before that run, and S ends inside the C exon
  after it (C's outer ends = the widest read carrying C; inner exon bounds = C's introns; 10 bp tolerance).
- **Maximal container** of S = a container not itself contained in another container of S. S is merged only when it has
  exactly one maximal container; then S's read counts as one read of that container. No stitching: a container is always a
  chain some read carries in full.
- **New chain** = a container with exactly 1 full-chain read and >= 1 merged read (support >= 2), absent from the default GTF.

Bars (gorilla OR6737 development, A119b held out, same bars):
- M1 share RefSeq-correct among new chains >= the sample's default 2-read chains (0.170 gorilla, 0.066 A119b).
- M2 new correct chains >= 200 (gorilla) / >= 500 (A119b).
PASS = both on both samples. Reported: new chains in StringTie's and FLAIR's output; the share of merged reads that end
within 10 bp of their container's 3'-side end (either end, strand unknown); the same numbers without the uniqueness rule.

## Result (2026-11-11, job 4466601, `r1637/`): FAIL on gorilla (M1); A119b not run
| rule | new chains | correct | share | in StringTie | in FLAIR |
|---|---|---|---|---|---|
| unique maximal container | 15,464 | 762 | 0.049 (bar 0.170) | 1,933 | 9,616 |
| any maximal container | 162,456 | 4,861 | 0.030 | 7,768 | 56,298 |
343,764 reads merged under the unique rule; 67.5% end within 10 bp of an outer end of their container.
Reading: the reads can be merged without ambiguity, but a chain that then rests on one full read plus its fragments is
correct 5% of the time, well under a chain two full reads carry (17%). The merged support is not equivalent to a second read.
