# PREREG: counting expressed copies from variants that recur across reads, on a planted sample that differs from the reference (row 1924, 2026-11-23)

Rows 1922 / 1923: nothing inside ONE read overrules minimap2; its wrong calls are exact ties. User: "IsoCon can
cluster errors". IsoCon's principle is recurrence across reads (an error is in one read, a real difference in many).
Proposed and accepted ("yes, can we try that"): use that principle at the locus level, on a planted set where the
sample differs from the reference, IsoCon as comparator. Every planted set so far had sample = reference, where copies
with identical reference sequence produce identical reads and nothing could be counted.

## Built before any run
- `plant.py --private RATE`: every copy gets its own variants (by genome position, shared by its isoforms; 80%
  substitutions, 10% one-base deletions, 10% one-base insertions); reads come from that sequence; `sample.fa`,
  `private.tsv` written. RATE = 0.001 (one per kb, the order of a human genome against a reference).
- `bench/o2_families/sample_psv.py`: loci joined by exactly tied reads form a group; the tied reads are read at one
  locus of the group; a difference from the reference seen in >= 3 reads and >= 10% of the covering reads is a site of
  the sample; reads are grouped greedily into haplotypes that agree at every shared site; haplotypes of >= 3 reads are
  counted. It says HOW MANY different sequences are expressed among copies the reference cannot separate, never which
  copy. Thresholds fixed now.
- Scorer `r1924_sample_psv/score_psv.py`: units = sets of >= 2 copies with identical reference isoforms and >= 1
  expressed copy; truth = number of distinct expressed sample sequences in the unit (copies that drew no variant in
  their transcripts stay identical and count once). A haplotype is credited to the unit most of its reads come from
  (the read names are used for that credit only). IsoCon (pooled per family as in `isocon_arm.sh`): the number of those
  sequences one of its candidates reproduces exactly.

## Data
Human named families, fresh seeds 1930, 1931, 1932, scenarios E (all expressed) and S (half silent), `--private 0.001`.
Control: seed 1930 with no variant (`--private 1e-9`): sample = reference.

## Claims
- V1: with variants, the haplotype count equals the truth in >= 70% of units (three seeds, E and S pooled) and beats
  both constants (one sequence per unit; every copy of the unit).
- V2: control: the count is 1 in >= 95% of units (errors do not invent copies).
- V3: exact rate not more than 0.10 under IsoCon's.
All three = worth wiring as a column of the `o2` report (expressed sequences among `shared` copies), the user's call.

## Predictions (mine)
V2 holds. V1 is uncertain, 60-80%: partial reads that cover no site, and copies whose only variants lie in an exon of
one isoform, will be merged (under-count). IsoCon is at home here and I expect it ahead of the tool, within 0.10 or not
I cannot say. Limits known beforehand: uniform simulated errors, not HiFi's; the unit structure is the truth's; a real
sample also carries allelic (heterozygous) variants, which this count would report as extra copies.

## Result (job 4496800, workspace `r1924_sample_psv`; units with the count equal to the truth)

| | with private variants (3 seeds, E + S, 67 units) | control, sample = reference (22 units) |
|---|---|---|
| haplotypes from recurrent variants (`sample_psv.py`) | 29 (0.43); over 4, under 34 | 13 (0.59); over 0, under 9 |
| constant: one sequence per unit | 18 (0.27) | 22 (1.00) |
| constant: every copy of the unit | 35 (0.52) | 0 |
| every EXPRESSED copy distinct (needs the truth on expression) | 60 (0.90) | 5 (0.23) |
| IsoCon (sample sequences reproduced exactly) | 66 (0.99) | 22 (1.00) |

- **V1 FAILS** (0.43 against 0.70; below "every copy of the unit"). **V2 FAILS** (0.59 against 0.95), though not in the
  direction feared: no copy is invented (over 0), units are left without a haplotype. **V3 FAILS** (0.43 against 0.99).
- **IsoCon does this almost perfectly:** 66 of 67 units with every distinct expressed sequence reproduced exactly, and
  every control unit as one sequence. Its cost: 6-13% of its candidates match no sequence of the sample (14-54 per run).
  It was given the reads pooled per family by the truth, as in every earlier IsoCon arm.
- Why the tool fails (read off the counts, not investigated further): it under-counts. 30% of its haplotypes with
  variants mix reads of different sequences (171 of 243 pure): greedy grouping lets partial reads that cover no site
  bridge two sequences. And groups chained by partial reads span several units, so a unit can be left with nothing.
- Predictions: V2 wrong, V1 below my range (I said 60-80%), IsoCon ahead as expected, by far more than 0.10.
- **Not built into `o2`, and not tuned on these seeds.** The principle (recurrence across reads) is confirmed by IsoCon,
  not by my tool. The supported next step is a hybrid, the user's call: Rustle supplies the groups of tied loci and
  their reads, IsoCon is run on each group's reads and its candidates give the count and the sequences. It needs its
  own registered test with reads pooled by Rustle's groups, not by the truth's families.
- `plant.py --private` stays (the planted set the benchmark lacked); `sample_psv.py` stays as a failed arm, called by
  no stage.
