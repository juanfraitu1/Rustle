# Difference model: why a chain is in our output and not StringTie's or FLAIR's, 2026-10-11

**Pre-registered after the development numbers (chr20, both species) and before any held-out chromosome is run.** User
(2026-10-11): establish a rule that distinguishes our assembler from StringTie (and FLAIR), as a base for an operating point
between the two tools.

## Instrument (frozen)

`bench/tool_difference/diff_model.py` (sha1 `a78ee425fd2874a04f67f1601c21168ef0c0c6f4`). Inputs per chromosome: the primary
alignments (`samtools view -F 2308`, columns 1-6), the figure 1 reference restricted to the chromosome (gorilla
`GGO.GCF_029281585.2_RefSeq.gtf.gz`, human `A119b.chm13v2.0_RefSeq.gtf.gz`), our shipped default assembly
(`copy_assign --assemble-only --region <chrom> --assembly-junctions strict` + the shipped polish, primary pool), and the lab's
annotation-free StringTie 3.0.3 / FLAIR 3.0.1 GTFs restricted to the chromosome. Read chains use our default de-duplication
key. Every multi-exon intron chain exclusive to one side gets the first class that applies:

- **other-only**: `ours_gate` (>= 2 reads carry it exactly; we lose it in gate / polish / motif), `single_shared` (1 read,
  every junction also carried by another read), `single_novel` (1 read, a junction only it carries), `corrected` (no read
  carries it; a read chain matches within 25 bp at every junction), `stitched` (no read carries it or a near copy).
- **ours-only**: `near_stronger` (a junction within 25 bp, donor and acceptor, of one more reads carry), `two_reads`
  (exactly 2 reads), `subset` (a contiguous sub-chain of another of our transcripts), `other`.

`reference-exact` = share of a class whose intron chain a reference transcript has exactly.

## Development numbers (chr20; seen)

| | gorilla NC_073244.2 | human chr20 |
|---|---|---|
| StringTie-only chains (reference-exact) | 1,274 (11.5%) | 2,954 (3.3%) |
| single_shared share / ref-exact | 63.1% / 14.4% | 37.6% / 5.0% |
| stitched share / ref-exact | 13.7% / 5.2% | 23.2% / 2.2% |
| ours-only vs StringTie (ref-exact) | 1,412 (24.3%) | 2,853 (10.4%) |
| near_stronger share / ref-exact | 40.2% / 34.9% | 21.8% / 17.8% |
| FLAIR-only chains (ref-exact) | 3,217 (3.5%) | 10,474 (1.7%) |
| FLAIR-only single_novel / single_shared ref-exact | 2.2% / 8.5% | 0.5% / 2.8% |
| FLAIR-only corrected + stitched | 0 | 0 |
| ours-only vs FLAIR (ref-exact); two_reads ref-exact | 1,180 (24.7%); 19.8% | 2,051 (10.4%); 5.1% |

## Held-out (decides)

Chromosomes fixed before running: gorilla OR6737 **NC_073242.2** and **NC_073245.2**; human A119b **chr18** and **chr22**.
(These samples' other chromosomes served as held-out for register 1264-1265; these four chromosomes are reused only as the
genome-wide remainder there, and this test adopts no assembler rule.) Counts are summed over the two chromosomes of a species.
Claims, each in both species:

- **D1** `single_shared` is the largest StringTie-only class.
- **D2** StringTie-only `single_shared` is reference-exact more often than StringTie-only `stitched`.
- **D3** ours-only (vs StringTie) `near_stronger` is reference-exact more often than StringTie-only chains overall.
- **D4** FLAIR-only chains are all read-carried (`corrected` + `stitched` = 0), and `single_novel` is reference-exact less often
  than `single_shared`.
- **D5** ours-only (vs FLAIR) `two_reads` is reference-exact more often than FLAIR-only chains overall.
- **D6** ours-only chains are reference-exact more often than other-only chains, against both tools.

**Verdict.** D1-D6 in both species => the difference model HOLDS and is the stated basis for the in-between operating point;
any failure is named and the model is restated without that claim. Class shares are reported per chromosome with their dev
values (no bar).

**Reading the model (fixed now).** StringTie's lead is mostly single-read chains whose junctions other reads share (D1, D2);
its 25 bp junction correction removes chains that are right more often than its own additions (D3). FLAIR's lead is
single-read chains, mostly with a junction no other read has (D4), and its support floor of 3 drops two-read chains that are
right more often than its additions (D5). An operating point between the two keeps our chains and adds StringTie's
`single_shared` class (our `--read-isoform-k 2`), not FLAIR's `single_novel` class.
