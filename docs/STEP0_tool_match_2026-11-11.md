# Step 0, r1633: how far is the initial assembly from StringTie's and FLAIR's, and what would matching cost? (2026-11-11)

Descriptive, no bar, no new assembly run. Gorilla OR6737 (StringTie 3.0.1 `-L`, FLAIR collapse, the benchmark GTFs), multi-exon
intron chains on contig + strand; correct = exact RefSeq chain. Scripts `bench/famcmp/tool_match.py`, `tool_match_reads.py`,
`nascent_class.py`; outputs `r1633/`.

## 1. Agreement of chain sets (gorilla)
| assembly | chains | correct | precision | share of StringTie's chains | Jaccard | share of FLAIR's chains | Jaccard |
|---|---|---|---|---|---|---|---|
| default | 73,348 | 25,794 | 0.352 | 0.639 | 0.442 | 0.361 | 0.309 |
| StringTie | 67,491 | 23,248 | 0.344 | 1 | 1 | 0.277 | 0.228 |
| FLAIR | 138,125 | 24,293 | 0.176 | 0.566 | 0.228 | 1 | 1 |
| 1-read floor | 127,014 | 29,006 | 0.228 | 0.729 | 0.339 | 0.412 | 0.273 |
| polish off | 155,147 | 27,392 | 0.177 | 0.680 | 0.260 | 0.461 | 0.277 |
| extend-transfrags | 76,139 | 26,424 | 0.347 | 0.650 | 0.440 | 0.364 | 0.306 |
| compat-support | 74,907 | 26,111 | 0.349 | 0.648 | 0.443 | 0.369 | 0.315 |
| transfrags union | 79,296 | 26,726 | 0.337 | 0.669 | 0.444 | 0.400 | 0.341 |
| junction snap 10 | 68,240 | 24,196 | 0.355 | 0.651 | 0.478 | 0.327 | 0.280 |
The two tools agree with each other less (0.228) than either does with ours. No existing switch moves agreement much.
A119b: default shares 0.408 of StringTie's chains (Jaccard 0.251), 0.191 of FLAIR's (0.165); extend and compat within 0.01.

## 2. Why: tool-only and ours-only chains by full-chain primary records (gorilla)
| class | 0 reads | 1 read | 2 | 3-4 | 5+ | correct share |
|---|---|---|---|---|---|---|
| StringTie-only (24,341) | 3,944 | 15,463 | 2,424 | 1,110 | 1,400 | 0.07-0.13 |
| ours-only vs StringTie (30,198) | 718 | 314 | 10,130 | 7,864 | 11,172 | 0.12-0.25 |
| FLAIR-only (88,236) | 0 | 68,459 | 9,637 | 6,259 | 3,881 | 0.03-0.05 |
| ours-only vs FLAIR (23,459) | 720 | 331 | 11,897 | 4,596 | 5,915 | 0.13-0.28 |
StringTie-only: 64% have one full-chain primary record, 16% none (causes not measured); StringTie omits 29,166 chains
that have two or more full-chain reads. FLAIR's side: 78% one-read chains.

## 3. The price of matching (gorilla; sets built from the tools' own output, so each is an upper bound for any rule)
| set | chains | correct | precision |
|---|---|---|---|
| ours | 73,348 | 25,794 | 0.352 |
| exactly StringTie | 67,491 | 23,248 | 0.344 |
| ours + StringTie-only with >= 2 reads | 78,282 | 26,208 | 0.335 |
| ours + StringTie-only with 1 read | 88,811 | 27,778 | 0.313 |
| ours + StringTie-only with 0 reads | 77,292 | 26,192 | 0.339 |
| ours or StringTie | 97,689 | 28,590 | 0.293 |
| exactly FLAIR | 138,125 | 24,293 | 0.176 |
| ours + FLAIR-only with >= 2 reads | 93,125 | 26,654 | 0.286 |
| ours or FLAIR | 161,584 | 28,559 | 0.177 |
| ours and StringTie | 43,150 | 20,452 | 0.474 |
| in at least two of the three | 64,974 | 25,182 | 0.388 |
On gorilla, matching either tool exactly is a loss on both axes (A119b not yet split; FLAIR has more correct chains there, r1619). Taking a tool's extra chains buys at most +2,796 (StringTie) or
+2,765 (FLAIR) correct chains, at precision 0.293 / 0.177.

## 4. Nascent (`~/stringtie_modded`, `-N`: `mark_if_nascent`, rlink.cpp)
A prediction is a nascent of a higher-coverage one when it shares its first introns from the 5' end and then runs into the next
intron. The guide-based synthetic nascents need `-G`; this link rule does not. In our default it marks 450 of 73,348 chains
(0.6%), 68 correct (0.151); 259 of them are in StringTie's `-L` output, where the rule marks 1,474 (0.165 correct).

## 5. Merging reads that are subsets of a longer read
Already in the binary: `--direct-transfrags` (reads of contained transfrags are transferred to the maximal one) and
`--compat-support` (truncated reads support the chain that contains them). On RefSeq (r1618) their new chains are correct
0.150 / 0.094 (gorilla), 0.046 (A119b compat) against 0.170 / 0.066 for 2-read chains; above they add 575 / 2,015 of
StringTie's chains. A merged read cannot make a StringTie-only chain with 0 full reads by containment alone.

## 6. r1636: what StringTie does with the chains we have and it lacks (decision log, gorilla NC_073242.2)
`~/stringtie_modded` rebuilt from a workspace copy (`r1636/src`; the shipped binary needs a newer glibc). Its `-L` output on
NC_073242.2 equals the 3.0.1 baseline (3,122 of 3,122 transcripts, strand + exons), with and without the log. Read-level
fate (`bench/famcmp/st_fate2.py`, `-p 1` log with per-read events): the reads whose alignment gives exactly the chain, the
transfrag StringTie built from each (rebuilt from its graph nodes), then the predictions with that chain. A first
chain-string version (`st_fate.py`) put 52% in a catch-all class and is superseded. Check on the 1,836 chains both tools
share: 78% equal, 22% "junction moved", so that label carries about that much reconstruction noise.
| fate of the 1,346 chains ours-only on the contig | chains | share | RefSeq-correct |
|---|---|---|---|
| the reads' transfrag lacks an end intron of the read | 366 | 0.272 | 0.197 |
| the reads' transfrag has a junction moved (<= 25 bp; noisy label) | 331 | 0.246 | 0.199 |
| transfrag equals the chain, never predicted (172 are sub-chains of a StringTie chain) | 197 | 0.146 | 0.076 |
| transfrag differs otherwise | 150 | 0.111 | 0.180 |
| predicted, killed as retained intron | 81 | 0.060 | 0.185 |
| no read with that chain in StringTie's input log | 73 | 0.054 | 0.068 |
| predicted, killed by isofrac | 43 | 0.032 | 0.209 |
| predicted, killed as included / deleted / flagged off | 67 | 0.050 | |
| all introns lost / no transfrag | 38 | 0.028 | |
Terminal exons: of 87,882 spliced reads on the contig, 6,808 lose their first and 6,091 their last exon in the transfrag;
75% of those exons are inside a node of the read's own graph, 15% only in another strand's graph, 10% in none. Their length
(median 225 bp) is close to kept terminal exons (315 bp). The cause is not identified.
StringTie-only chains on the contig (1,253): 332 have no read whose transfrag is the chain (0.081 correct), 576 one such
read (0.113), 345 two or more (0.107).
Reading: no single rule separates the two outputs. Junction correction and sub-chain absorption (profile S3f, r1634) cover
about 40% of the drop side; an end-intron loss of unidentified cause covers another 27%, on chains that are 0.197 correct.

## 7. A119b price table (`r1633/price.a.tsv`)
| set | chains | correct | precision |
|---|---|---|---|
| ours | 219,517 | 38,715 | 0.176 |
| exactly StringTie | 211,865 | 31,781 | 0.150 |
| ours or StringTie | 344,926 | 43,072 | 0.125 |
| exactly FLAIR | 626,983 | 39,230 | 0.063 |
| ours + FLAIR-only with >= 2 reads | 381,502 | 42,383 | 0.111 |
| ours or FLAIR | 726,693 | 46,651 | 0.064 |
| in at least two of the three | 179,854 | 37,330 | 0.208 |

## 8. Oracle: is read trimming the lever toward StringTie? (registered 2026-11-11 before running; r1639)
On NC_073242.2, replace each read's own chain by the chain of its longest StringTie transfrag (from the `-p 1` log, so it
carries StringTie's trimming and junction correction) and keep every chain with >= k such reads (k = 1, 2, 3). Jaccard /
reached / inside against StringTie's output on the contig, beside the default's on the same contig.
- Bar: some k reaches Jaccard >= 0.60 -> trimming is the lever, find its rule in `rlink.cpp`.
- Below 0.60 at every k -> StringTie's output is decided after the read level (flow); the StringTie side closes with this
  ceiling. Caveat: log-derived chains carry junction noise (22% of reads on shared chains), so the figure is a lower bound.

Result (`r1636/fate2.oracle.txt`; StringTie has 3,089 chains on the contig): **below the bar at every k.**
| set | chains | shared | reached | inside | Jaccard |
|---|---|---|---|---|---|
| default | 3,182 | 1,836 | 0.594 | 0.577 | 0.414 |
| transfrag chains, >= 1 read | 14,434 | 2,360 | 0.764 | 0.164 | 0.156 |
| >= 2 reads | 6,520 | 1,777 | 0.575 | 0.273 | 0.227 |
| >= 3 reads | 4,254 | 1,294 | 0.419 | 0.304 | 0.214 |
| default + (>= 2 reads) | 7,890 | 2,181 | 0.706 | 0.276 | 0.248 |
Even with StringTie's own trimmed read chains, a read-count rule is further from StringTie than the default, and 24% of
StringTie's chains are no read's transfrag chain at all: a perfect selector over read-level chains tops out at 0.764
reached. What StringTie prints is decided by its flow over the graph, after the read level. StringTie side closed here.

## 9. What IsoQuant does differently, and is any of it portable? (2026-11-12, descriptive, gorilla OR6737; r1646 de novo `--polya_trimmed all`)
IsoQuant 4.1.0 `default_pacbio` (isoquant.py): intron graph from corrected reads; an intron within 10 bp of one with twice
its count is merged into it; a novel isolated intron needs 2 reads and 2% of the neighbourhood; novel model needs 2 reads
and 2% of the gene's; only canonical splice sites; ends from clustered read ends (polyA when present).
| set | chains | correct | precision |
|---|---|---|---|
| ours | 73,348 | 25,794 | 0.352 |
| IsoQuant | 37,089 | 16,677 | 0.450 |
| both | 31,285 | 16,239 | 0.519 |
| ours only | 42,063 | 9,555 | 0.227 |
| IsoQuant only | 5,804 | 438 | 0.075 |
Ours-only: 36,865 carry an intron IsoQuant reports nowhere (0.232 correct, 8,570 correct chains), 5,198 are other
combinations of its introns (0.189). None of ours has a non-canonical intron. A relative floor on our output (reads >=
x of the locus maximum): 0.05 keeps 23,507 at 0.411 (plain reads >= 3: 22,270 at 0.423; about +0.01 over the read-count
curve at that count); 0.10 keeps 20,682 at 0.451 (reads >= 4: 20,053 at 0.474; below the curve). Reading: IsoQuant is a
stricter subset of ours; its intron filter would cost a third of our correct chains at a precision the read-count cut
already gives (r1646 matched-count table); intron merging = the junction snap that failed r1527 / r1634. Nothing ported.

## 10. IsoQuant's ends and read assignment (2026-11-12, descriptive, gorilla OR6737; scratch script, RefSeq ends as truth)
Ends on the 16,239 correct chains in ours and IsoQuant (share within 10 bp of the nearest RefSeq transcript with the chain):
| ends | 5' within 10 / 50 | 3' within 10 / 50 |
|---|---|---|
| ours emitted (`or.gtf`) = most extended full-chain read | 0.351 / 0.603 | 0.532 / 0.609 |
| IsoQuant emitted | 0.321 / 0.579 | 0.661 / 0.737 |
| rule on our reads: mode | 0.117 / 0.373 | 0.668 / 0.741 |
| rule: median | 0.054 / 0.323 | 0.651 / 0.727 |
| rule: FLAIR window peak | 0.073 / 0.343 | 0.680 / 0.750 |
On the 13,139 chains all four tools share: 3' within 10 bp ours 0.527, StringTie 0.623, IsoQuant 0.673, FLAIR 0.695;
5' ours 0.361, IsoQuant 0.327, StringTie 0.248, FLAIR 0.061; `or.final.gtf` 5' 0.056 (median 80 bp short), 3' 0.535.
Reading: RefSeq 5' ends sit at the longest read, RefSeq 3' ends at the commonest read end. Our most-extended rule is the
best 5' rule and the worst 3' rule; the 3' mode / window peak on our own reads reproduces IsoQuant's 3' accuracy.
The layered `final.gtf` moves the 5' end to the median, which by RefSeq is 80 bp short (r1487 judged it on another measure).
Read assignment: IsoQuant skips secondary and supplementary records (alignment_processor.py; `--use_secondary` off), so
every read stays at its primary locus; it only chooses the isoform inside the gene. Of 4.40 M reads: 1.94 M unique,
0.08 M unique with a minor difference, 0.47 M ambiguous between isoforms, 1.94 M inconsistent with every model, 0.52 M
non-informative. Nothing there bears on copy placement (O2). Not tested: a 3' rule held out on A119b; strand of read ends.

## 11. Why the end rules differ between gorilla OR6737 and A119b (2026-11-12, descriptive, `bench/famcmp/end_profile.py`, `r1660/prof.{or,a}.tsv`)
Per RefSeq-correct chain of the default: signed error of the most extended and of the commonest full-chain read end
(+ = beyond RefSeq), by depth. Truth classes differ: gorilla RefSeq is 96% XM (predicted), A119b 67% NM (curated).
| | gorilla 5' | A119b 5' | gorilla 3' | A119b 3' |
|---|---|---|---|---|
| median error of the mode, all | -59 | 0 | -3 | -8 |
| ... at 100+ reads | -127 | -2 | -3 | -78 |
| median error of the most extended, 100+ reads | 0 | +141 | +122 | +4 |
| most extended beyond RefSeq by > 10 bp, 2-3 / 100+ reads | 0.16 / 0.28 | 0.43 / 0.94 | 0.13 / 0.63 | 0.13 / 0.36 |
| share of reads within 5 bp of the mode (median) | 0.36 | 0.38 | 0.73 | 0.69 |
- 5': on A119b the commonest read start is RefSeq's TSS to the base (NM 0, XM -5) and the most extended read overshoots
  more the deeper the chain. On gorilla the commonest start is 60 bp inside RefSeq's (XM) and only the longest read
  reaches it. Human XM also sits at -5, so the record class alone does not explain gorilla; either gorilla reads are
  5'-shortened or gorilla XM models are extended from other evidence (not separated here).
- 3': commonest end > 100 bp upstream of RefSeq's in 0.205 of gorilla chains and 0.428 of A119b chains. A-rich genome
  just downstream (internal-priming signature; 0.02 at on-target ends): 0.451 of those gorilla ends, 0.234 of the A119b
  ones. So most of A119b's short 3' ends are not internal priming: the commonest 3' end is an upstream site of its own
  (hypothesis: proximal polyA use in this sample; tissue not checked). Soft clips are 0 at both ends in both samples.
- The most extended end degrades with depth at whichever side has outliers (gorilla 3', A119b 5'): it is a maximum.
Reading: no single statistic is right for both because the two truths and the two libraries differ at opposite ends.
A depth-robust rule would have to be chosen per end from the read-end distribution itself; not registered.

### 11b. PolyA signal at the commonest 3' read end (2026-11-12; AATAAA or ATTAAA 5-40 bp upstream; background 0.034 in the same window 300 bp further upstream, both samples)
| commonest 3' end vs RefSeq | gorilla: share / signal / A-rich downstream, no signal | A119b: share / signal / A-rich downstream, no signal |
|---|---|---|
| on RefSeq (within 10 bp) | 0.644 / 0.881 / 0.003 | 0.459 / 0.879 / 0.003 |
| short 11-100 bp | 0.076 / 0.534 / 0.356 | 0.059 / 0.748 / 0.096 |
| short > 100 bp | 0.205 / 0.423 / 0.407 | 0.428 / 0.581 / 0.205 |
| beyond > 10 bp | 0.074 / 0.631 / 0.194 | 0.054 / 0.525 / 0.270 |
Ends > 100 bp short with a signal and no A-rich stretch: 0.078 of all gorilla chains, 0.236 of all A119b chains; A-rich
without a signal (internal priming): 0.083 and 0.088. Reading: internal priming is equally common in both samples; the
A119b excess is upstream sites with a real polyA signal, i.e. the sample mostly uses a polyA site upstream of the one
RefSeq annotates in about a quarter of its correct chains. The r1660 measure counts those as misses for the mode rule.

### 11c. Reported, not removed (user 2026-11-12): opt-in driver stage `altends`
`tools/rustle_pipeline.sh altends --bam B --fasta G --out PREFIX` -> `PREFIX.altends.gtf` = `PREFIX.gtf` plus, on multi-exon
transcript lines whose commonest full-chain read 3' end is > 10 bp from the emitted end and has >= 2 reads: `alt_tes`,
`alt_tes_reads "n/total"`, `alt_tes_polya_signal`, `alt_tes_a_rich` (`bench/famcmp/alt_ends.py`, selftest). No coordinate
change, `PREFIX.gtf` untouched. Gorilla NC_073242.2: 684 of 3,182 models get one (434 with a polyA signal, 134 A-rich).
