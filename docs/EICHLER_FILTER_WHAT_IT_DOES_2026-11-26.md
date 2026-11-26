# The Eichler lab's 10-point filter: what it is, what it needs, and what it does not do (2026-11-26)

Written for one recurring statement: "the correct copy is the one whose primary is 10 AS above the others, so the rule
tells you quickly which copy a read came from." This page puts the paper's text, the three ways we ran it, and the
numbers side by side. Short version:

- The filter is real and it works **in the setting the paper used it in**: reads aligned to the haplotype assemblies of
  the genome, compared between paralog clusters.
- It is a **filter that discards**, not an assignment method: reads that fail it are "marked as ambiguous and ignored".
- It **cannot be carried to a single reference**. The same reads that it decides on an individual's own haplotypes are
  ties on CHM13: for TBC1D3, 96% of them.
- At an exact tie nothing in the read decides, for this rule or any other; the aligner's primary there is a draw.

## 1. What the paper says (Guitart et al. 2024, Methods, "Iso-Seq and transcript analyses")
> "Human FLNC reads from fetal brain and iPSCs were mapped to all validated human haplotypes. Next, we compared these
> primary alignments to one another and considered the cluster paralog from which they were derived. Any Iso-Seq read
> with primary minimap2 alignment scores of 10 or greater for a given paralog cluster relative to all other cluster
> mappings were retained, while other mappings were marked as ambiguous and ignored."

Mapping in the same section: `minimap2 -ax splice --sam-hit-only --secondary=yes -p 0.5 --eqx -K 2G -G 8k -N 20`, reads of
at least 1,000 bp and 99.9% quality. Source in the vault: `~/ThesisVault/Fulltext/Guitart_2024_TBC1D3_fulltext.md`.

Read closely, the rule has four parts that the short form drops:
1. **Many targets, one per haplotype.** The read is mapped to each validated haplotype assembly separately.
2. **Primaries are compared with primaries**, one per haplotype. It is not "primary against its secondaries".
3. **The unit is a paralog CLUSTER** (a group from the paper's phylogeny), not an individual copy.
4. **Failing reads are dropped.** The output is a set of confident reads per cluster, used to say which clusters have
   expression support. It does not place every read.

## 2. Three ways we ran it
| | what was run | where | rows |
|---|---|---|---|
| A | our ADAPTATION to one reference: the primary leads every other alignment of the same mapping by 10 AS | HG002, KB3781, OR6737, A119b | 1775-1786, 1954 |
| B | the PAPER'S PROCEDURE: reads mapped separately to the individual's maternal and paternal assemblies with the paper's parameters, primaries compared per paralog, lead of 10 = retained | HG002 (hg002v1.1) | 1972 |
| C | what the aligner does when nothing leads | planted reads, renamed and re-seeded | 1954 addendum, 1969 |

## 3. Results
### B. The paper's procedure on the individual's own haplotypes (row 1972)
| | all families | TBC1D3 | NPIP |
|---|---|---|---|
| family reads | 116,099 | 97 | 1,463 |
| retained for a paralog | 98.9% | 53.6% | 92.0% |
| ambiguous, ignored | 1.1% | 46.4% | 8.0% |
With the genome's own copies as targets nearly every read has one best home. Even so, in the family the rule was
written for, **almost half of the reads are ambiguous and are thrown away.**

### B against A: the same reads on one reference
| reads the paper's procedure RETAINS | all families | TBC1D3 | NPIP |
|---|---|---|---|
| near or exact tie on CHM13 (`ms`) | 5.3% | 96.2% | 44.9% |
| our single-reference adaptation abstains | 6.2% | 94.2% | 40.4% |
| CHM13 primary is on the paralog the haplotypes name | 91.5% | 26.9% | 50.9% |
For TBC1D3 the procedure names a paralog for 52 reads; on CHM13 the primary is on that paralog for 14 of them. Among
about ten near-identical reference copies that is what a draw gives.

### A. The adaptation on one reference (row 1954)
| | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| reads that are near ties on `ms` | 16,245 | 53,712 | 316,155 | 198,767 |
| of them, the rule fires (primary >= 10 AS above every other) | 11.7% | 9.6% | 15.4% | 8.7% |
Where it fires and the read is clear on `ms`, the primary is right in 99.89% (HG002) and 99.995% (KB3781): the rule is
correct where it speaks. It speaks for about one near-tied read in ten. And taking "the best AS above 10" instead of the
primary is WRONG in 8-42% of reads at leads of 11-50 (row 1786), because AS charges about 24 per intron and `ms`,
which the aligner ranks by, does not: a processed copy beats its parent gene on AS for that reason alone.

### C. When nothing leads
At an exact tie minimap2 orders the alignments by a hash of the read's name, its length, a seed and the first anchor.
Renaming the reads moves 29% of the tied primaries (all of them mapping quality 0); twins of identical copies receive a
copy's reads 343 to 344. After renaming, 6-25% of the transcripts of every primary-only assembler change (row 1969).

## 4. Why it is not a quick way to know a read's copy
1. **It needs the haplotype assemblies of the genome the reads come from** (or a panel of them). That is the expensive
   part, and it is what does the work: with them 99% of reads are decided, without them the TBC1D3 reads are ties.
2. **It discards what it cannot decide.** 46% of TBC1D3 reads and 8% of NPIP reads even on the individual's own
   haplotypes. A method that must give every read a copy cannot use it as its answer.
3. **Its unit is a cluster of paralogs**, coarser than a copy. Run per copy, as an assembler needs, more reads fall
   between units.
4. **The number 10 is not a threshold with a meaning on one reference** (row 1786: no break at 10 in the error curve).
5. **On AS it favours intronless copies** over the spliced gene they come from.
6. **At exact ties it has nothing to compare.** Of the TBC1D3 reads the paper's procedure retains, 42% are exact ties on
   CHM13 and another 54% near ties.

## 5. What is right about it (so that it is not undersold)
- Where a read leads clearly, the primary is where the read came from: 99.9%. That is most reads of most genes.
- In the paper's setting it is a sound, conservative filter for "which paralog clusters have expression support".
- It is the right picture of the problem: the evidence that decides a read is the true sequence of the copies.

## 5b. The part that adapts to one reference: the cluster as the unit (row 1974)
The paper compares paralog CLUSTERS. Asked that way on CHM13 (clusters of copies at 99% sequence identity or the
pipeline's families; truth = the paralog the haplotype procedure names), the same 10-point lead does this to the 7,051
HG002 reads that are near ties between copies:

| lead asked between | reads decided | decided right | the answer names |
|---|---|---|---|
| individual copies | 14.7% | 79.3% | 1 copy |
| clusters at 99% identity | 78.4% | 98.0% | 2 copies in the median |
| the pipeline's families | 81.0% | 99.2% | 2 copies in the median |
TBC1D3 is one cluster of 9 copies under every definition: every read is decided for the array and none is placed in it.
NPIP at 99% identity: 77.5% decided, 98.5% right, clusters of 6. So the adaptation is "which cluster, and how many
copies it holds", with the rest discarded, exactly the paper's output. Per copy the answer does not exist on one
reference. Record: `docs/PREREG_cluster_filter_2026-11-26.md`.

It generalizes (row 1975, `docs/PREREG_cluster_filter_denovo_2026-11-26.md`): with clusters made from the run's own copy
sequences and no annotation, HG002 80.0% decided and 96.6% right, KB3781 (gorilla) 74.3% and 99.3%, planted reads with
private variants 82-85% and 100%; per locus the rule decides 8-16%. The retained cluster holds 2 loci in the median on
real samples and 6-7 in the large planted families, so a verdict always states how many copies it covers.
In the pipeline: stage `o2` with `RUSTLE_O2_CLUSTERS=1` (opt-in) writes `PREFIX.o2.clusters.tsv` and
`PREFIX.o2.cluster_loci.tsv`; the built step decides 83.7% of HG002's tied reads at 97.3% right.

## 6. Words to use
- "The Eichler lab's 10-point filter (Guitart 2024)", and say which was run: "as the paper describes it, on HG002's
  own haplotypes" (row 1972) or "adapted to a single reference" (rows 1775-1786, 1954). Never "Eichler's rule".
- "We ran the filter as the paper describes it. It retains 99% of family reads on the individual's own assemblies and
  about half of the TBC1D3 reads. On one reference the same TBC1D3 reads are ties in 96% of cases. The filter works
  because the copies' true sequences are the targets; an assembler works where they are not available."
- Not: "the rule fails". Not: "the primary is random" (it is a reproducible, uninformed draw at exact ties only).

## 7. Limits of row 1972
One individual's two haplotypes (the paper: all validated human haplotypes, reads from other individuals). A paralog
cluster here is the nearest reference copy inside the read's span, not the paper's phylogenetic clusters, so MORE reads
fall between clusters here than would there. No read-quality filter (these reads carry none). 97 TBC1D3 reads in this
sample. The paralog label was changed once after the first scoring; both versions are in the pre-registration and the
verdicts are the same.

## 8. Where everything is
Page (private): https://claude.ai/artifact/CCNiwfa8YVJfKmsWGG4Gif (source `pages/eichler_filter/index.html` in the workspace).
Pre-registrations with results: `docs/PREREG_eichler_haplotypes_2026-11-26.md` (row 1972),
`docs/PREREG_as10_baseline_2026-11-25.md` (row 1954, the hash on one read, ties by family),
`docs/PREREG_rename_and_abundance_2026-11-26.md` (row 1969). Tools: `bench/o2_truth/eichler_haplotypes.py`,
`as10_census.py`, `as10_truth.py`, `stability.py`. Workspace: `r1972_eichler_haplotypes` (`score_window.txt`),
`r1776_hg002_room` (haplotype assemblies and the family reads), `coin_toss_2026-11-25`.
