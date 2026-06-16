# Read-coherence + PSV copy-split on REAL co-located/tandem paralog families (GGO)

Driver: `bench/copy_split_colocated.py` (adapted from `bench/copy_split_realdata.py`, the
RABL2 driver). Same algorithm port of `split_readchain_by_psv`
(`src/rustle/vg_family/copy_split.rs`), same CIGAR/allele extraction, same
`Bio.Align.PairwiseAligner` PSV discovery. Generalised from 2 copies on separate contigs to
N tandem copies on one contig; reads pulled once over the whole cluster span (incl
secondaries); each alignment assigned to the copy it overlaps most.

Inputs: BAM `/home/juanfra/winloci_scratch/GGO.bam`, FASTA
`/home/juanfra/winloci_scratch/GGO.fasta`, annotation
`bench/copy_recovery_eval/results/ref.gtf`. Contig NC_073235.2 (autosomal, ~100x library
baseline). DAZ (collapsed single locus, 3.4x) and RBMY (unannotated here) are unusable as
noted in the task; these two autosomal clusters are the best available substitutes.

> **Annotation caveat (important).** The task framed these as "6 tandem copies of
> LOC101144552" and "3 tandem copies of LOC101132628". In this RefSeq annotation the
> *gene* `gene-LOC101144552` actually sits at 144779454-144823106 (OUTSIDE the given region
> 144691596-144756221), and `gene-LOC101132628` is a single transcript at the START of its
> given region. The given regions instead contain a CLUSTER of distinct paralog gene IDs.
> I took the multi-exon protein-coding paralogs inside each cluster as the "copies":
> LOC129532018 / LOC129532017 / LOC101144552 and LOC101132628 / LOC115933736 / LOC115933737.
> There is no annotated 6-copy tandem array here.

---

## Family 1 — LOC101144552 cluster (NC_073235.2:144691596-144756221)

| copy | span | len | pairwise %id to ref |
|---|---|---|---|
| LOC129532018 (ref) | 144622575-144683955 | 61381 bp | — |
| LOC129532017 | 144691596-144721773 | 30178 bp | 68.9 % |
| LOC101144552 | 144779454-144823106 | 43653 bp | 59.1 % |

Direct pairwise genomic identities: 18-vs-17 **68.9 %**, 18-vs-552 **59.1 %**, 17-vs-552
**73.2 %**.

- **#copies:** 3 (putative; see caveat)
- **#PSVs:** 23,828 columns, but this number is **meaningless as a "family"** — it reflects
  ancient/divergent paralogs over full genomic spans (incl. introns), not a tight tandem
  array. At ~60-73 % identity these are not a recent-duplication family the PSV machinery is
  designed for.
- **Total alignments over the cluster span:** 212 (29 primary, 182 secondary, 1
  supplementary).
- **Identifiability budget:** of the 212 alignments, **182 observe 0 PSV columns** and only
  **30 span >= 2 (and >= 3) PSV columns.** The 30 informative ones are the genuine
  transcripts of the three genes; the 182 zero-PSV alignments are a pile of MAPQ-0
  secondaries (many stacking at identical start coords 144655231 / 144816108 with
  11-19 introns) — cross-mapped multimappers from a well-expressed paralog elsewhere, not
  evidence for these copies.
- **Shared intron chain across copies (the co-located signal):** **0.** 83 distinct spliced
  chains, none assigned across >= 2 copies.
- **Multimapper molecules placed on >= 2 copies:** 3 — all MAPQ-0 / mostly 0-intron stubs
  (not informative).
- **Identifiable splits found:** **0** at K=2 and K=3 (min_reads=2). Every chain is a
  singleton or all-one-copy; no chain carries two competing >= 2-read allele vectors.

**Verdict (Family 1): sparsity- AND structure-limited.** These are divergent paralogs
(~60-73 % id), not a tight tandem family; each copy's own transcript is supported by only a
handful of reads, and the bulk of the "coverage" is uninformative cross-mapped MAPQ-0
secondaries. No phasing is possible or meaningful here.

---

## Family 2 — LOC101132628 cluster (NC_073235.2:135066522-135115177)

This is the genuine high-identity tandem trio.

| copy | span | len | pairwise %id to ref |
|---|---|---|---|
| LOC101132628 (ref) | 135066522-135069510 | 2989 bp | — |
| LOC115933736 | 135094425-135098162 | 3738 bp | 97.2 % |
| LOC115933737 | 135111301-135115177 | 3877 bp | 97.2 % |

Direct 736-vs-737 genomic identity = **99.06 %** (the two real co-located copies are nearly
identical, incl. introns).

- **#copies:** 3 (high-identity tandem)
- **#PSVs:** 99 (ref-anchored). Direct 736-vs-737 = **35 PSVs**, of which **18 are exonic**
  (coverable by a spliced read) and 17 fall in introns (skipped by spliced reads).
- **Total alignments over the cluster span:** **10** (7 primary, 3 secondary).
- **Identifiability budget:** **6 / 10** alignments span >= 2 (and >= 3) PSV columns;
  4 observe 0. The informative ones observe 34-75 PSV columns each.
- **Shared intron chain across copies:** **0** (5 distinct chains; none crosses copies in the
  by-coordinate sense — see below).
- **Multimapper molecules placed on >= 2 copies:** **3** — and these are the real signal:
  molecules 234128, 596395, 668711 each align as **PRIMARY (MAPQ 60) on LOC115933737** and
  **SECONDARY (MAPQ 0) on LOC115933736**, both with an identical 10-intron chain.
- **Identifiable splits found:** **0** at K=2 and K=3.

### Does PSV evidence actually resolve the multimappers? Yes — but not via this split.
For all 3 molecules, the PRIMARY-on-737 alignment observes 34 ref-anchored PSV columns and
matches **737 perfectly (34/34) vs 736 imperfectly (26/34)** — i.e. **8 exonic columns
decisively assign them to 737.** minimap2's MAPQ-60 primary is correct; the PSVs confirm it.

### Why the copy-split emits 0 anyway (the co-located structural point).
`split_readchain_by_psv` groups reads **by intron chain** and splits a chain iff >= 2 allele
vectors each have >= min_reads support and pairwise differ at >= K observed columns. Here the
copies sit at **different genomic coordinates**, so reads on 737 and reads on 736 live in
**disjoint column frames** — a read on one copy is `None` at the other copy's PSV columns.
The 3 real 737 reads share one vector (no second competing vector with >= 2 reads on the same
coordinate frame), and the 736 secondaries observe **0** of the 35 direct 736-vs-737 PSVs
because those reads' spliced exons + the column anchoring don't coincide. So there is never a
single chain carrying two >= 2-read competing copies → no split fires.

This is exactly the **RABL2 contrast**: RABL2A/B also live on separate contigs/coordinates,
and the algorithm was designed for reads **piling onto ONE coordinate frame** and being
phased apart, not for tandem copies whose reads land at distinct coordinates. Co-located
tandem copies do **not** collide into one chain here.

**Verdict (Family 2): sparsity-limited (decisively).** This is a clean, near-identical
(99 %) tandem pair with informative exonic PSVs, but only **10 alignments total / 7 primary**
over the whole cluster and **6 PSV-spanning alignments**. The one real multimapper event
(3 molecules, 737-primary / 736-secondary) is correctly assignable by PSVs (8 decisive
columns), but 3 reads on a single copy is far below the depth needed to *phase* a second copy
out of a shared frame, and the split algorithm's by-chain/one-frame design does not apply to
tandem copies at distinct coordinates regardless.

---

## Overall honest verdict

| | LOC101144552 | LOC101132628 |
|---|---|---|
| copies (putative) | 3 | 3 |
| pairwise identity | 59-73 % (divergent) | 97-99 % (tight tandem) |
| PSV columns | 23,828 (not a real family) | 99 (35 direct 736/737; 18 exonic) |
| total alignments | 212 (182 are MAPQ-0 cross-map stubs) | 10 |
| primary molecules | 29 | 7 |
| span >= 2 PSV / >= 3 PSV | 30 / 30 | 6 / 6 |
| shared intron chain across copies | 0 | 0 |
| multimapper molecules (>=2 copies) | 3 (uninformative stubs) | 3 (real: 737-pri / 736-sec) |
| identifiable splits (K=2, K=3) | 0 | 0 |

**Enough coverage to phase? No.** Both families are far too sparse. LOC101132628 has only
7 primary molecules / 6 PSV-spanning alignments over the whole cluster; LOC101144552's
apparent depth (212) is 86 % uninformative MAPQ-0 cross-mapped secondaries, with just 29
real primaries spread over three divergent genes.

**Splits recovered? None (0/0).** No chain in either family carries two competing >= 2-read
allele vectors, so the ported `split_readchain_by_psv` correctly emits no splits.

**Sparsity-limited — and additionally a structural mismatch.** This is a **valid negative
result**: (1) these families are lowly expressed exactly as the task anticipated, so the
identifiability budget is tiny (6 informative alignments in the clean 99 %-id pair); and
(2) the read-coherence+PSV split, like the RABL2 case, operates on reads piling onto a single
coordinate frame — co-located **tandem** copies put their reads at distinct coordinates, so
there is no within-chain collision to split. The genuine multimapper event in LOC101132628
(3 molecules, primary-737/secondary-736) *is* PSV-resolvable (8 decisive exonic columns,
34/34 vs 26/34), confirming minimap2's primary; but that is per-read **assignment**, not a
copy **split**, and 3 reads is well below phasing depth.

Net: with the available GGO long-read data, co-located/tandem paralog copy-splitting is
**not demonstrable** — the data is identifiability-starved (sparsity) and the test families
are either too divergent to be a real family (LOC101144552) or too shallow to phase
(LOC101132628). Consistent with the broader memo that the multimapping well is exhausted and
real recall levers are non-multimapping.
