# Copy Assignment And Gate (consolidated)

> Merged from 17 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [copy_split_colocated](#copy-split-colocated) · [copy_split_realdata](#copy-split-realdata) · [copy_split_sim](#copy-split-sim) · [core_gate_atscale](#core-gate-atscale) · [core_gate_gw](#core-gate-gw) · [core_gate_gw_downstream](#core-gate-gw-downstream) · [core_gate_pipeline](#core-gate-pipeline) · [psv_llr_vs_votes_killtest](#psv-llr-vs-votes-killtest) · [psv_tau_sweep](#psv-tau-sweep) · [copy_specific_junctions](#copy-specific-junctions) · [copy_assign_phasing](#copy-assign-phasing) · [copy_resolution_boundary](#copy-resolution-boundary) · [identifiability_bound](#identifiability-bound) · [eval_significance_gate](#eval-significance-gate) · [resolution_improvement_bound](#resolution-improvement-bound) · [stringtie_vs_copy_headtohead](#stringtie-vs-copy-headtohead) · [quant_benchmark](#quant-benchmark)


---

## copy_split_colocated

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


---

## copy_split_realdata

# Read-coherence + PSV copy-split on the REAL RABL2 family — honest result

Run: `python bench/copy_split_realdata.py` (deterministic).
Mirrors the rust port `split_readchain_by_psv` in `src/rustle/vg_family/copy_split.rs`.

**Verdict up front: NEGATIVE / already-resolved-by-minimap2.** On the real RABL2A/RABL2B
family the PSV copy-split fires **0 identifiable splits** at every tested K and min_reads.
The two copies are separated *before* the PSV step — by the structural (intron-chain) axis
and by minimap2's primary placement — so the PSV axis has nothing left to resolve. This is a
useful confirmation that the mechanism is correctly conservative, not a recovery win.

## Inputs
- RABL2A `NC_073235.2:15131653-15147533` (+), 15,881 bp
- RABL2B `NC_086018.1:48818440-48832011` (-, reverse-complemented for the A-vs-B alignment), 13,572 bp
- BAM (HiFi long reads, with secondaries): `/home/juanfra/winloci_scratch/GGO.bam`
- genome FASTA: `/home/juanfra/winloci_scratch/GGO.fasta`

## (a) PSV budget — copies are richly distinguishable in PRINCIPLE
- A-vs-B global alignment (Bio.Align.PairwiseAligner): **189 PSV columns** (fixed-base diffs).
- Spacing on copy A: min 1 bp, median 47 bp, max 1472 bp, mean 79 bp → ~1 PSV per 84 bp.
- So the reference sequences are NOT the limit — there is ample fixed divergence.

**But the per-read budget is thin.** Of 397 alignments pulled over both loci, only
**77 span ≥ 2 PSV columns** (and the same 77 span ≥ 3 — the K=2 vs K=3 distinction is moot
here). The other 318 alignments span **0** PSV columns. Those 318 are almost all secondary
alignments that minimap2 placed at the *other* copy: they cover the locus interval but their
CIGAR consumes little reference inside the small RABL2 exons, so they touch no PSV. The
identifiability budget that actually reaches the split logic is 77 alignments / 77 spanning
reads, ~29 on A and ~48 on B.

## (b) Copies the method resolves
| K | min_reads | distinct intron chains | identifiable (split) copies |
|---|-----------|------------------------|-----------------------------|
| 2 | 2 | 245 | **0** |
| 3 | 2 | 245 | **0** |

Every emitted unit is a `merged` (non-identifiable) chain-group. No chain ever splits into
two PSV-haplotype copies.

## (c) The contrast — do reads that share a junction chain get split by PSVs? **No — they never share a cross-copy chain.**
This is the crux. The split can only fire when ≥ 2 reads share one intron chain AND carry
allele vectors that differ at ≥ K PSV columns. On real RABL2 that situation does not occur:

- **245 distinct exact intron chains** for 397 alignments. Chain-group sizes:
  191 chains have 1 read; the largest spliced chains have 7–12 reads, **all single-locus**.
- **Exactly 1 intron-chain mixes both loci: the monoexon (0-intron) group** (13 A + 11 B).
  But within it, only 4 reads span ≥ 2 PSVs and **all 4 are RABL2B**; the RABL2A monoexon
  reads span 0 PSVs. So there is no A-vs-B allele contrast to act on, and 0 candidate vectors
  clear `min_reads = 2`. The split correctly stays merged.
- For every spliced chain, RABL2A and RABL2B reads live on **different contigs with different
  intron coordinates**, so they are already in different chain-groups before PSVs are even
  consulted. The structural axis fully separates the copies on its own.
- Even after a generous counterfactual fuzzy-merge of intron chains (round splice sites to
  10 bp to absorb HiFi jitter), the 11 multi-read fuzzy chains are **still each single-locus**,
  and within each the reads agree with their own copy's reference at the PSVs
  (`disagree_with_own_ref` ≈ 0–3 of ~12–16 observed columns; the few mismatches are sporadic
  sequencing errors, not a second haplotype). There is no hidden second copy inside a chain.

## (d) Comparison to minimap2's own placement — the split agrees, adds nothing
minimap2 already assigns each read molecule to exactly one copy via the primary flag:
**29 primaries at RABL2A, 48 at RABL2B** (77 distinct molecules — the same 77 spanning reads).
For every multi-read chain, the reads' PSV alleles match the copy minimap2 chose
(`disagree_with_own_ref` ≈ 0). So where the PSV signal is observable it **confirms** the
primary placement; it never contradicts it and never reassigns a read. The 170 + 148 secondary
alignments span no PSVs, so the PSV axis cannot even weigh in on them — it has strictly less
information than the primary flag, not more, for this family.

## (e) Honest verdict
On the real RABL2 family the read-coherence + PSV copy-split **does not recover copy structure,
because there is no residual copy structure to recover**:

1. **Already resolved by the structural axis.** RABL2A and RABL2B sit on different contigs with
   different intron coordinates, so reads never share a cross-copy intron chain. The copies fall
   out of the read-coherence grouping for free; PSVs are never the deciding factor.
2. **Already resolved by minimap2.** The one chain that does mix both loci (monoexon) has no
   spanning-read allele contrast, and elsewhere the primary flag already places every molecule
   correctly — the PSV alleles only echo it.
3. **Identifiability is thin where it would matter.** 189 reference PSVs exist, but only 77/397
   alignments span ≥ 2 of them and the 318 cross-mapped secondaries span 0, so the PSV channel
   has no leverage on exactly the reads (secondaries) where copy ambiguity would arise.

**Implication.** The PSV copy-split is correctly conservative (it does not invent spurious
copies on a genuinely two-locus family) but it is **inert as a recovery lever on RABL2**: the
work is already done upstream by intron-chain coordinates and by minimap2's primary flag. The
mechanism would only pay off in the harder regime where two copies are **co-located / share an
intron chain** (true tandem duplicates, or a collapsed assembly where both copies map to one
locus) AND reads span multiple PSVs — that is where a shared chain could split on alleles.
RABL2A/RABL2B, living on separate chromosomes, is not that regime. This matches the prior
finding that rustle-VG's secondary advantage is confined to *tied, co-located* copies, not to
cleanly separated paralogs like RABL2.

## Key numbers
- 189 reference PSVs; 1 per ~84 bp (median spacing 47 bp).
- 397 alignments over both loci (200 at A: 30 primary / 170 secondary; 197 at B: 49 primary / 148 secondary).
- 77 alignments span ≥ 2 PSV columns (= ≥ 3 too); 318 span 0.
- 245 distinct exact intron chains; exactly 1 mixes both loci (monoexon), and it carries no
  A-vs-B spanning-allele contrast.
- minimap2 primary placement: 29 RABL2A / 48 RABL2B (77 molecules).
- Identifiable PSV splits: **0** at K ∈ {2,3}, min_reads = 2.

## Verdict

**NEGATIVE — minimap2 (and the structural intron-chain axis) already separated the copies; the joint read-coherence + PSV split recovered NO additional copy structure on RABL2, and this is not an identifiability ceiling.**

Reproduced end-to-end (`python bench/copy_split_realdata.py`, exit 0) and the rust bridge is green (`cargo test --lib copy_split` → 19 passed / 0 failed: 11 pre-existing + 8 new `bridge_` tests).

Which of the three possible outcomes happened:
- *Did the method recover copy structure?* **No.** 0 identifiable splits at K=2 and K=3, min_reads=2; every emitted unit is `merged`.
- *Was it identifiability-limited (too few PSVs / reference indistinguishable)?* **No.** The reference is richly distinguishable: **189 fixed PSVs**, ~1 per 84 bp (median spacing 47 bp). The limit is not reference divergence.
- *Did minimap2 already separate the copies?* **Yes — together with the structural axis.** RABL2A and RABL2B sit on different contigs with different intron coordinates, so reads **never share a cross-copy intron chain** (245 distinct chains; the spliced multi-read chains are all single-locus; exactly 1 chain mixes both loci — the monoexon group — and it carries no A-vs-B spanning-allele contrast). minimap2's primary flag already places all 77 molecules (29 A / 48 B), and where PSVs are observable they only **confirm** that placement, never reassign.

Key numbers: **189 PSVs** (1 per ~84 bp) · **397 alignments** (A 200 = 30 pri / 170 sec; B 197 = 49 pri / 148 sec) · **77 alignments span ≥ 2 PSV cols** (= ≥ 3; the K-distinction is moot here) and **318 span 0** · **245 distinct intron chains**, exactly **1** mixes both loci · minimap2 primary **29 A / 48 B** (77 molecules) · **identifiable splits = 0** at K ∈ {2,3}, min_reads = 2.

Implication for the thesis:
1. **The method is correct and conservative**, not broken: it does not invent spurious copies on a genuine two-locus family, and the PSV alleles it does observe agree with minimap2's placement (`disagree_with_own_ref` ≈ 0–3 of ~12–16 observed cols = sequencing error).
2. **RABL2 is the wrong family to demonstrate the win.** It defines the boundary: the PSV-split has leverage *only* when two copies **share one intron chain** (co-located / collapsed tandem) AND reads span ≥ K PSVs. On separated-chromosome paralogs the structural axis + primary flag do the work for free, and the 318 cross-mapped secondaries — exactly the reads where copy ambiguity could arise — span **0** PSVs, so the PSV channel has *strictly less* information than the primary flag here.
3. **Next family to try** = the co-located / collapsed-tandem regime where copies map to one locus and reads carry multiple PSVs: DAZ / RBMY (tied, co-located paralogs) — consistent with the prior finding that rustle-VG's secondary advantage is confined to tied co-located copies, not cleanly separated ones like RABL2. That is where read-coherence + PSV-split would actually be the deciding axis.


---

## copy_split_sim

# Copy-split in-silico validation (collapsed-tandem regime + identifiability boundary)

`bench/copy_split_sim.py` -> `bench/copy_split_sim.png`

Definitive controlled simulation that the joint **read-coherence + PSV copy-split**
separates co-located / collapsed-tandem paralog copies that read-coherence alone
cannot, and an empirical map of the **identifiability boundary** (the thesis
identifiability theorem). The split logic is a faithful line-by-line port of
`split_readchain_by_psv` (+ `select_identifiable_copies`, `read_consistent_with`,
`distinct_columns`, `consensus_haplotype`) from
`src/rustle/vg_family/copy_split.rs`.

## Model

- **N copies share ONE intron chain** (the collapsed-tandem regime). Read-coherence
  (the structural axis) keys transcripts by exact ordered intron chain, so the whole
  family collapses to a single transcript — this is the adversarial case where
  structure alone is blind.
- **P PSV columns** carry the only separating signal. Each copy gets a fixed allele
  vector, rejection-sampled so any two copies differ at >= K columns (the
  identifiability premise).
- **Reads** are simulated per copy (`C` reads/copy). Each read covers a random
  **contiguous sub-interval** of the transcript (models finite read length L /
  partial spanning): it observes only the PSVs inside its span (others = `None`).
  Each observed allele is corrupted at per-base error `e` (HiFi ~ 0.003).
- A copy is **recovered** iff enough error-free reads spanning >= K of its columns
  land on its exact allele vector to reach `min_reads`, and it stays >= K-distinct
  from every already-admitted copy (greedy `select_identifiable_copies`).

## Determinism

Fully deterministic. `random.seed` / `numpy.random.seed` are set from **fixed
integer constants only** (`MASTER_SEED=20260616`, `DEMO_SEED=11`,
`SWEEP_SEED_BASE=700`); every sweep cell re-seeds from a constant + cell offset.
Never seeded from the wall clock. Verified byte-identical across reruns:
- stdout: identical
- `copy_split_sim.png` md5 stable = `1e1e2a32b9032f6f0ea46e757ef5603a`

## (1) Single-config demo  (N=3, P=8, C=10, e=0.003, K=2, min_reads=2)

```
read-coherence ALONE -> 1 transcript  (the whole family collapses)
read-coherence + PSV -> 3 identifiable copies
truth = 3 copies; fraction recovered = 1.000
  copy reads: 10, 10, 9   (all identifiable=True)
```

The PSV axis recovers all 3 copies where structure alone gives 1.

## (2) Sweep — identifiability boundary

N=3, P=8, K=2, e=0.003, min_reads=2, 40 seeded replicates/cell.
Sweep over read span (12.5% .. 100% of PSVs) x coverage (1 .. 32 reads/copy),
plus an error-`e` sweep (0 .. 0.20) at fixed generous geometry.

- **Panel A** heatmap: fraction of copies recovered. Sharp navy(merged) ->
  green(recovered) transition; white contour = 0.99 boundary.
- **Panel B**: recovery vs the theorem's order parameter = per-copy clean
  (>= K-spanning, error-free) support for the candidate vector. A clean sigmoid:
  ~0 below `min_reads`, snapping to 1.0 once support clears the candidate
  threshold (green dashed). This is the empirical identifiability boundary.
- **Panel C**: error floor. Flat at 1.0 through HiFi e=0.003; degrades as `e`
  rises (errors fragment support off the exact true vector below `min_reads`).

### Where the boundary sits (key numbers)

- Full recovery at **full span** needs **coverage >= 3 reads/copy**.
- Full recovery at **max coverage** needs **span >= 62% of the PSVs**
  (a read must reliably span >= K columns above the error floor).
- **Error floor**: recovery = **1.000 at e=0.003 (HiFi)**; holds to ~e=0.01,
  then degrades, falling to **0.550 at e=0.20**.

## Interpretation (the theorem, empirically)

Recovery -> 1.0 exactly when the expected count of error-free reads per copy that
span >= K distinguishing columns crosses `min_reads`; below that the chain stays
**merged** (non-identifiable). Both the geometric axis (coverage x span must put
>= min_reads clean reads on each copy's >= K columns) and the error axis (per-base
`e` erodes that clean support) reproduce the predicted boundary. This is the
in-silico validation of the read-coherence + PSV copy-split and its
identifiability bound.

---

## Verdict

**The method works in-silico in its target regime, the recovery boundary matches
the identifiability theorem, and the real-data win is blocked only by coverage —
not by the algorithm.**

**(1) Real GGO co-located families — coverage-limited, valid NEGATIVE.**
Both real autosomal co-located/tandem clusters were tested (`bench/copy_split_colocated.md`)
and recovered **0 splits (K=2 and K=3, min_reads=2)** — but the cause is sparsity,
not a method failure. The two clusters bracket the failure modes:
- **LOC101144552 cluster**: putative copies are **divergent paralogs (59-73 % pairwise
  id)**, not a tight family, with only **29 real primaries** (212 alignments, but
  **182 = 86 % are uninformative MAPQ-0 cross-mapped secondaries**), and **0** shared
  intron chains. Structure- AND sparsity-limited.
- **LOC101132628 cluster (736/737)**: the genuine clean tandem pair — **99.06 %
  genomic id, 35 direct PSVs (18 exonic)** — but only **10 alignments total / 7 primary /
  6 PSV-spanning** over the whole cluster. The one real multimapper event (3 molecules,
  PRIMARY-MAPQ60 on 737 / SECONDARY-MAPQ0 on 736) **is** PSV-resolvable (**8 decisive
  exonic columns, 34/34 vs 26/34**) — but 3 reads is far below phasing depth, so this is
  per-read **assignment**, not a copy **split**.
This matches the task's prior: the deep-tandem families that need this method (DAZ
collapsed at 3.4x; RBMY unannotated) are exactly the ones GGO does not cover at depth.
GGO marks the boundaries of the multimapping well, it does not contain the target case.

**(2) Simulation — the split DOES recover what read-coherence merges, and the
boundary matches the theorem.** In the collapsed-tandem regime (N copies sharing ONE
intron chain), read-coherence alone collapses to **1 transcript**; read-coherence + PSV
recovers all copies (demo N=3 -> **3 copies, fraction 1.000**). The faithful port of
`split_readchain_by_psv` recovers a copy **iff >= min_reads error-free reads span >= K of
its distinguishing PSV columns** — exactly the identifiability condition. The sweep
confirms the boundary empirically:
- full span needs **coverage >= 3 reads/copy**;
- max coverage needs **span >= 62 % of PSVs** (a read must reliably span >= K columns);
- **error floor: recovery = 1.000 at HiFi e=0.003**, holding to ~e=0.01, degrading to
  **0.550 at e=0.20**.
Panel B is the clean order-parameter sigmoid: recovery snaps to 1.0 once clean
>= K-spanning support per copy clears the `min_reads=2` candidate threshold, and is ~0
(merged) below it. **Recoverable iff enough reads span >= K PSVs above the error floor** —
the theorem, reproduced.

**(3) Thesis implication.** The read-coherence + PSV copy-split is **validated in its
target regime in-silico**, with a recovery boundary that **matches the identifiability
theorem on both the geometric (span x coverage) and error axes**. The real-data
demonstration is **gated on coverage, not on the method**: it needs a **deep co-located
long-read dataset** in the tight-tandem regime (e.g. **testis HiFi for DAZ/RBMY**, where
the high-identity arrays are expressed at depth). GGO usefully **marks the boundaries** of
where the method can and cannot fire — RABL2 copies are too separated (separate coordinate
frames, no within-chain collision), and its co-located tandems (LOC101144552/LOC101132628)
are either too divergent to be a real family or too shallow (6-7 informative reads) to
phase. The honest position: **method correct and boundary-matched in silico; real win
awaits a depth-adequate co-located dataset.**


---

## core_gate_atscale

# Contiguous-core family gate AT SCALE

> **Of 406 annotated within-family gene pairs across 62 universe families, the contiguous-core gate (T=0.13) KEEPS 361 (88.9%) as true cores and DROPS 45 (11.1%) as domain-sharing false-family memberships.**

> On the Compara-labeled subset the gate is fully consistent: it KEEPS 5/5 confirmed pairs and DROPS 7/7 domain-sharers.

## What this measures (and what it does NOT)

This is the **universe-annotated-family scale**: it takes the universe TSV's family assignment as GIVEN and asks, pair by pair, which of those declared within-family memberships the contiguous-core gate would uphold (KEEP) versus reclassify as domain-sharing (DROP). It is a **proxy for the full de-novo pipeline** -- the production gate actually runs on de-novo assembled loci, not on the RefSeq-annotated family universe -- so the genome-wide rate here is an annotation-scale estimate, not the exact de-novo rate. The Compara labels exist only for a small subset (12 named pairs); the rest are UNLABELED, so for them KEEP/DROP is the gate's verdict, not a verified ground truth.

## The gate

```
core(a,b) = largest_ungapped_equal_block(a,b) / min(len(a), len(b))
KEEP  iff core >= T   (T = 0.13)
DROP  iff core <  T
```

computed from a GLOBAL Needleman-Wunsch alignment (match=+2, mismatch=-1, gap-open=-5, gap-extend=-1; `Bio.Align.PairwiseAligner`, the robust aligner reused verbatim from `poa_family_definition.py`). The 'largest ungapped equal block' is the longest single contiguous gap-free run of paired columns in the global alignment (one block between two gaps), divided by the shorter gene -- the contiguous homologous core, robust to the gappy chance-match filler a global aligner pads non-homologous pairs with. This is exactly `poa_family_definition.py`'s `biggest` metric.

## (1) Coverage / evaluability

| quantity | value |
|---|---|
| universe families | 62 |
| universe distinct genes | 195 |
| genes present in gene_rep.fa | 195 (100.0%) |
| genes MISSING from gene_rep.fa | 0 |
| of present genes, LOC* (provisional) loci | 154 (79%) |
| families evaluable (>=2 genes present) | 62 of 62 |
| within-family pairs evaluated | 406 |

**Match rate is 100%** -- every universe gene (including all 154 LOC* provisional loci) is present in `gene_rep.fa`, so all 62 families and all 406 within-family pairs are evaluable; nothing is dropped for missing sequence.

## (2) The at-scale effect

| verdict | pairs | % of pairs |
|---|---|---|
| KEPT (true core, core>=T) | 361 | 88.9% |
| DROPPED (domain-sharer, core<T) | 45 | 11.1% |
| total | 406 | 100% |

**Genome-wide (annotation-scale) false-family reclassification rate = 11.1%** of declared within-family pairs would be split off by the gate as domain-sharers.

### Per-family fraction-of-pairs-kept distribution

| keep-fraction | families |
|---|---|
| 1.00 (all pairs kept) | 44 |
| mixed (0 < frac < 1) | 8 |
| 0.00 (all pairs dropped) | 10 |

Median per-family keep-fraction = 1.000; mean = 0.801. 44 of 62 families are fully upheld (every within-family pair kept); 10 are fully reclassified (every pair dropped); 8 are mixed (the gate splits the family into a true-core subset + domain-sharing outliers).

Quartiles of the per-family keep-fraction: q25=0.797, q50=1.000, q75=1.000.

## (3) Compara-labeled subset (the gate is correct where we can check)

The 12 Compara-checkable named pairs are the only within-family pairs with external ground truth. The gate's verdict on each:

| pair | family | Compara label | contiguous-core | gate verdict | correct? |
|---|---|---|---|---|---|
| RFPL1 <-> RFPL3 | LOC134758217 | confirmed | 0.608 | **KEEP** | yes |
| RFPL1 <-> RFPL2 | LOC134758217 | confirmed | 0.606 | **KEEP** | yes |
| RFPL2 <-> RFPL3 | LOC134758217 | confirmed | 0.520 | **KEEP** | yes |
| APOBEC3D <-> APOBEC3F | APOBEC3D | confirmed | 0.302 | **KEEP** | yes |
| RABL2A <-> RABL2B | RABL2A | confirmed | 0.174 | **KEEP** | yes |
| CASP8 <-> FLACC1 | CASP8 | domain-sharer | 0.055 | **DROP** | yes |
| ASDURF <-> ASNSD1 | ASDURF | domain-sharer | 0.031 | **DROP** | yes |
| GPR39 <-> LYPD1 | GPR39 | domain-sharer | 0.020 | **DROP** | yes |
| CCDC188 <-> ZDHHC8 | CCDC188 | domain-sharer | 0.014 | **DROP** | yes |
| CDPF1 <-> PPARA | CDPF1 | domain-sharer | 0.012 | **DROP** | yes |
| CREB1 <-> METTL21A | CREB1 | domain-sharer | 0.010 | **DROP** | yes |
| GCA <-> KCNH7 | GCA | domain-sharer | 0.008 | **DROP** | yes |

**The gate KEEPS all 5/5 confirmed pairs and DROPS all 7/7 domain-sharers** -- it is perfectly consistent with Compara on the labeled subset. Confirmed contiguous-core range [0.174, 0.608]; domain-sharer range [0.008, 0.055]. The shipped T=0.13 sits in the gap between them.

## Per-family table (all evaluable families)

Sorted by keep-fraction (most-reclassified first), then size.

| family | genes | pairs | kept | dropped | keep-frac | has Compara label |
|---|---|---|---|---|---|---|
| ASDURF | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| CASP8 | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| CDPF1 | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| COPS8 | 2 | 1 | 0 | 1 | 0.000 | - |
| CREB1 | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| GCA | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| GPR39 | 2 | 1 | 0 | 1 | 0.000 | domain-sharer |
| LOC129529456 | 2 | 1 | 0 | 1 | 0.000 | - |
| LOC134756662 | 2 | 1 | 0 | 1 | 0.000 | - |
| LOC134756677 | 2 | 1 | 0 | 1 | 0.000 | - |
| CCDC188 | 4 | 6 | 3 | 3 | 0.500 | domain-sharer |
| LOC134758217 | 4 | 6 | 3 | 3 | 0.500 | confirmed |
| LOC101150852 | 5 | 10 | 6 | 4 | 0.600 | - |
| LOC129529434 | 4 | 6 | 4 | 2 | 0.667 | - |
| LOC101144552 | 6 | 15 | 11 | 4 | 0.733 | - |
| LOC101129569 | 10 | 45 | 34 | 11 | 0.756 | - |
| LOC101123878 | 14 | 91 | 84 | 7 | 0.923 | - |
| LOC101126655 | 11 | 55 | 54 | 1 | 0.982 | - |
| GGTLC2 | 11 | 55 | 55 | 0 | 1.000 | - |
| LOC129532044 | 8 | 28 | 28 | 0 | 1.000 | - |
| LOC129529666 | 6 | 15 | 15 | 0 | 1.000 | - |
| LOC129529667 | 5 | 10 | 10 | 0 | 1.000 | - |
| DDT | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC101132221 | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC101132628 | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC101142890 | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC115931911 | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC129529548 | 3 | 3 | 3 | 0 | 1.000 | - |
| LOC129529560 | 3 | 3 | 3 | 0 | 1.000 | - |
| APOBEC3D | 2 | 1 | 1 | 0 | 1.000 | confirmed |
| AQP12A | 2 | 1 | 1 | 0 | 1.000 | - |
| DGCR6 | 2 | 1 | 1 | 0 | 1.000 | - |
| FAM246A | 2 | 1 | 1 | 0 | 1.000 | - |
| GGT1 | 2 | 1 | 1 | 0 | 1.000 | - |
| GP1BB | 2 | 1 | 1 | 0 | 1.000 | - |
| IGLL1 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101127159 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101127811 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101131416 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101134642 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101136027 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101138607 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101144123 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101145825 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101145885 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101146886 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101146937 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101147293 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101147656 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC101151758 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC109024534 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC109025943 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC115931965 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC115931973 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC115933728 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129529430 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129529513 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129529592 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129529611 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129530539 | 2 | 1 | 1 | 0 | 1.000 | - |
| LOC129532202 | 2 | 1 | 1 | 0 | 1.000 | - |
| RABL2A | 2 | 1 | 1 | 0 | 1.000 | confirmed |

## Honest caveats

- **Annotation-scale proxy, not the de-novo pipeline.** This runs on the RefSeq-annotated universe families (`universe.tsv`), taking their assignment as given. The production gate runs on DE-NOVO assembled loci; the genome-wide rate here is an annotation-scale estimate of the gate's behavior, not the exact de-novo reclassification rate.
- **Most pairs are UNLABELED.** Only the 12 Compara-checkable named pairs (5 confirmed + 7 domain-sharer) have external ground truth. For the other 394 pairs, KEEP/DROP is the gate's own verdict, validated only indirectly by its perfect agreement with Compara on the labeled subset.
- **LOC families included.** All 154 present LOC* (provisional/computationally-predicted) loci are included if in `gene_rep.fa`. Many of the largest, most-fragmented families are LOC tandem arrays (e.g. GGTLC2, LOC101126655), whose internal pairs the gate splits aggressively; these dominate the DROP count and may include both genuine sub-family structure and annotation noise.
- **Gene-representative sequences.** One representative sequence per gene (`gene_rep.fa`); a different representative isoform could shift a borderline pair's coverage.
- **Block vs exact-run metric.** The gate metric is the largest ungapped ALIGNED block / shorter gene -- identical to `poa_family_definition.py`'s `biggest`, the metric that separates the labeled classes at T=0.13. A stricter variant (largest ungapped base-IDENTICAL run) was also computed but is NOT used as the gate: it does not separate the labeled set at T=0.13, because genuine recent-duplicate cores carry scattered SNPs that break the exact run (confirmed exact-run range [0.046, 0.085], all below T). The aligned-block metric tolerates those SNPs and is the gate-faithful one.
- **Single threshold, no per-family tuning.** A single shipped T=0.13 is applied uniformly; the separation was established on the 12-pair labeled set in `poa_family_definition.py` and is re-used here unchanged.
- **Determinism.** Alignment is deterministic; only the figure jitter uses a fixed seed (1729). Every reported number is reproducible.

---

## Verdict

The contiguous-core family gate is **real, robust, and lands in the expected direction at every level tested** — Rust unit, annotation-scale Python, and a real `rustle --vg` pipeline run. (Note: `core_gate_atscale.py` regenerates everything above this line; this Verdict is appended by the synthesis review and is re-added by hand if the script is re-run.)

**(1) RUST — the poasta robustness fix is real and verified.** `cargo test --lib` is fully green (451 passed, 0 failed, 3 ignored). The robustness regression guard `contiguous_core_coverage_divergent_flanks_still_high` PASSES: two true copies with a 400 bp shared core and divergent 5'+3' flanks that scored ~0.005 under poasta's weak default gap-open (`GapAffine::new(1,1,2)`) now score high under the dedicated strong gap-open (`GapAffine::new(1,1,32)`). The fix is clean: `poa_msa` is untouched (still 1,1,2 everywhere else); a new `poa_msa_with_costs` takes explicit affine costs and `contiguous_core_coverage` is the only caller using the strong gap-open. Default-off byte-identical (`family_min_core_coverage_default_off_merges_like_jaccard`) and the domain-sharer split (`family_min_core_coverage_splits_domain_sharer_keeps_copies`, with its precondition `..._off_jaccard_groups_all_three`) all PASS. The merge-site gate (family_graph.rs:438-447) is strictly additive: it only runs when `min_core_cov > 0.0` (or the trace env is set), can only PREVENT merges, and is inert at the default 0.0.

**(2) PYTHON at-scale — computed correctly, deterministic, honestly scoped.** Re-running `core_gate_atscale.py` reproduces the report byte-for-byte. Coverage is honest and complete: 195/195 universe genes present in `gene_rep.fa`, all 62 families and all 406 within-family pairs evaluable. At-scale effect: KEEP 361 (88.9%) / DROP 45 (11.1%). The Compara subset check is perfectly consistent: KEEPS 5/5 confirmed (core 0.174–0.608), DROPS 7/7 domain-sharers (core 0.008–0.055), T=0.13 cleanly in the gap. Verified independently that the labeled pairs (e.g. APOBEC3D/APOBEC3F, CASP8/FLACC1) really are co-assigned to one `family_id` in `universe.tsv`, so the gate's verdicts on them are genuine. The report is candid that this is an **annotation-scale proxy** for the de-novo pipeline, that 394/406 pairs are UNLABELED (gate's own verdict), and that LOC tandem arrays dominate the DROP count — appropriate, not overclaimed.

**(3) PIPELINE — a REAL `rustle --vg --vg-layer2` on/off run was achieved on NC_073235.2; every quantitative claim corroborated by the run artifacts (`/tmp/coregate/`), and the gate changed family formation in the expected direction.** Confirmed against the actual GTFs and trace logs: 2623 = 2623 transcripts (no gain/loss); 67 jaccard-passing cross-copy pairs, sharply bimodal with a genuine empty gap (5 pairs at 0.004–0.011, then 62 pairs at ≥0.314); gate ON @0.13 fires `would_gate=true` on exactly those 5 sub-0.013 domain-sharer pairs (all in one family, cid0/cid2, mismatched lengths incl. 1140 vs 2846 bp) and leaves all 62 true-copy merges intact. Downstream the diff is exactly 3 genes: RSTL.447 and RSTL.508 change ONLY the `family_id` token (verified: every other field byte-identical → cosmetic enumeration renumber), and RSTL.647 is the real effect — same 4 transcripts/coords but the spurious `rescue_class "strand_pure_minority"` / `copy_status "novel"` tags drop (2→0 each), cleaning the contaminated novel-copy attribution. This is NOT overclaimed: it is honestly scoped to one contig where domain-sharers happen to be rare (5 of 67 pairs, all one family), the architecture finding is correct and verified in source (the gate lives in `merge_singletons_by_sequence`, the dispersed-paralog Stage-1b path; co-located `cluster_by_position` families bypass it), and a genome-wide measurement is correctly deferred as future work. The earlier-stage `/tmp/off.gtf` / `/tmp/on.gtf` files are stale artifacts from an unrelated experiment (1905/1989 tx) and are NOT the basis of these numbers — the load-bearing artifacts are in `/tmp/coregate/`.

**Bottom line:** the robustness fix is genuine and well-isolated; the at-scale measurement is correct, deterministic, and honestly caveated; and the pipeline result is a real (not partial, not asserted) on/off comparison whose every headline number checks out against the run logs. Nothing in the three workstreams is overclaimed.


---

## core_gate_gw

# Genome-wide contiguous-core family-gate measurement
**Site measured:** rustle's de-novo cross-copy singleton-exon merges (`merge_singletons_by_sequence`) — the gate's actual call site.
**Threshold:** T = 0.13 on contiguous-core coverage (`core_cov`).
**Gate status:** DEFAULT-OFF. This is a *measurement* of what the gate would do, not a behavior change.

## Headline
- **25 chromosomes** parsed (22 carry cross-copy merge pairs; 3 are empty).
- **1876** cross-copy merge pairs total.
- **155 (8.3%) DROP** at core_cov < 0.13.
  - 145 have core_cov < 0.05 (deep lobe).
  - 10 sit in the valley [0.05, 0.13).
- **1721 (91.7%) KEEP** at core_cov ≥ 0.13.
- **616** pairs have core_cov ≥ 0.95 (near-identical true copies).
- `NC_073235.2` = **5 / 67** drops (matches the prior pipeline test).

## Bimodality — why 0.13 is a principled threshold
The core_cov distribution is **bimodal**: a heavy mass near 1.0 (near-identical true copies; 616 pairs ≥ 0.95) and a second lobe down near 0.0 (145 pairs < 0.05). Between them is a sparse **valley** spanning roughly [0.05, 0.13) holding only 10 pairs. T = 0.13 is placed inside that valley, so the decision boundary does not slice through a dense region — small threshold perturbations move very few pairs. This is exactly the kind of separated, non-arbitrary cut a gate wants: the data, not the knob, picks the split.

## Why jaccard alone fails and core fixes it (the danger zone)
At the same numeric cut, **jaccard alone keeps all 1876 pairs** (every pair has jaccard ≥ 0.13) — i.e. a jaccard-only gate at 0.13 would never drop anything here. Yet **106 pairs have high jaccard (≥ 0.50) but low contiguous-core (< 0.13)**, and 64 of those have jaccard ≥ 0.80. These are **domain-sharers**: two copies that share a large fraction of *k*-mers/segments (high jaccard) but whose shared material is *not contiguous core* — it is a repeated domain, a shared exon cassette, or scattered homology rather than a single colinear copy body. Jaccard cannot tell a tandem-repeated domain from a real copy; contiguous-core coverage can, because it demands the overlap form one run.

Worst offenders (jaccard, core_cov) — jaccard says "merge", core says "don't":

| jaccard | core_cov |
|---|---|
| 0.997 | 0.004 |
| 0.995 | 0.004 |
| 0.995 | 0.004 |
| 0.984 | 0.004 |
| 0.971 | 0.008 |
| 0.971 | 0.008 |

## Per-chromosome
| chrom | pairs | drop (core<0.13) | drop %% |
|---|---|---|---|
| NC_073228.2 | 135 | 54 | 40.0% |
| NC_073224.2 | 334 | 24 | 7.2% |
| NC_073246.2 | 28 | 17 | 60.7% |
| NC_073244.2 | 69 | 15 | 21.7% |
| NC_073238.2 | 135 | 14 | 10.4% |
| NC_073240.2 | 219 | 6 | 2.7% |
| NC_073235.2 | 67 | 5 | 7.5% |
| NC_073236.2 | 86 | 3 | 3.5% |
| NC_073229.2 | 9 | 3 | 33.3% |
| NC_073248.2 | 395 | 2 | 0.5% |
| NC_073230.2 | 85 | 2 | 2.4% |
| NC_073247.2 | 42 | 2 | 4.8% |
| NC_073243.2 | 9 | 2 | 22.2% |
| NC_073242.2 | 7 | 2 | 28.6% |
| NC_086017.1 | 6 | 2 | 33.3% |
| NC_073237.2 | 48 | 1 | 2.1% |
| NC_086018.1 | 11 | 1 | 9.1% |
| NC_073241.2 | 114 | 0 | 0.0% |
| NC_073233.2 | 67 | 0 | 0.0% |
| NC_073227.2 | 6 | 0 | 0.0% |
| NC_073231.2 | 2 | 0 | 0.0% |
| NC_073234.2 | 2 | 0 | 0.0% |
| NC_073232.2 | 0 | 0 | 0.0% |
| NC_073239.2 | 0 | 0 | 0.0% |
| NC_073245.2 | 0 | 0 | 0.0% |

## Honest scope & caveats
- **What this measures.** rustle's *de-novo cross-copy singleton-exon merges* (`merge_singletons_by_sequence`). That is the gate's literal call site, so these counts are the population the gate would actually act on.
- **Per-chrom / within-chrom families only.** The OOM-safe genome-wide protocol runs one chromosome at a time, so only *within-chromosome* families are grouped. Cross-chromosome paralogs (e.g. RABL2A/RABL2B) are not co-grouped here — but those are the coordinate-separated case this gate does **not** target anyway; the gate is about distinguishing true copies from domain-sharers inside a candidate family, not about linking copies across chromosomes.
- **`would_gate` trace field is vacuous — ignored.** Every trace line reads `would_gate=false` because it was computed against the active threshold 0.0 (gate-OFF during the trace run). We do **not** trust that field; every drop/keep number above is recomputed directly from `core_cov` against T = 0.13.
- **Gate is default-off.** Nothing here changes rustle's emitted output; this is a measurement of the gate's *would-be* effect.

## Reproduce
```
/home/juanfra/miniforge3/bin/python bench/core_gate_gw.py
```
Inputs: `/tmp/gw/cg_NC_*.trace` (deterministic parse).

## Genome-wide verdict

**What is solidly established (independently re-verified).** Reparsing all 25
`/tmp/gw/cg_NC_*.trace` files (1876 `[CORE_TRACE]` lines) reproduces every
headline exactly: **1876 pairs, 155 (8.3%) drop at core<0.13** (145 deep-lobe
<0.05 + 10 in the [0.05,0.13) valley), 616 near-identical (≥0.95),
**NC_073235.2 = 5/67**. The `core_cov` distribution is genuinely bimodal and
T=0.13 sits in a sparse valley (only 10 of 1876 pairs in [0.05,0.13)), so the
cut is data-placed, not knob-tuned. Core does separate domain-sharers from true
copies that jaccard cannot: **106 pairs (64 at jaccard≥0.8) have high jaccard
but core<0.13** — these would be merged by jaccard, dropped by core. The
`would_gate` trace field is confirmed vacuous (all 1876 = `false`, computed
against gate-OFF 0.0) and is correctly ignored; all decisions are recomputed
from `core_cov`.

**Downstream effect is benign-cleaning, NOT recall loss — confirmed by my own
coordinate-keyed diff.** On the worst-case contig NC_073228.2 (gate fires on
exactly 54/135 pairs), an independent re-keying of every transcript by
(strand, exon-chain) finds **0 chains lost, 0 chains gained, 4989=4989 in both**
(genuine identity, not cancellation). The only output change is **2 attribute
edits at one locus (RSTL.262)**: two guide-backed transcripts (guide:STRG.440.1
and .5) shed a spurious `copy_status "novel"` + `rescue_class
"strand_pure_minority"` tag while keeping identical exon coordinates and
identical cov/FPKM/TPM (verified byte-equal). The transcripts stay emitted.
gffcompare is unchanged both ways: **FSM 1733=1733, intron-chains 1716=1716,
loci 945=945, Sn/Pr 23.3/34.7=23.3/34.7**. So the answer to the downstream
question is the same as NC_073235.2: **benign-cleaning, no real isoform lost** —
the gate removes a false novel-copy attribution and nothing else.

**Honest deflation — do not overclaim.**
- **The downstream payoff is essentially negligible in magnitude.** Even on the
  worst-case contig with 54 gated merges (11x NC_073235.2's 5), the entire
  output effect is **2 attribute tags on 1 of 1497 genes**. Drop count does not
  scale to output churn — the gated cross-copy fusions simply did not survive to
  the final flow-extracted transcripts. The gate's value is precision *hygiene*
  on copy attribution, **not a headline recall/precision metric move** (Sn/Pr
  literally do not change).
- **The "jaccard alone drops nothing" line is true but tautological.** Every one
  of the 1876 pairs already passed an upstream jaccard prefilter (observed
  **min jaccard = 0.300** across all pairs; 0 pairs below 0.13). So this is not a
  fair head-to-head of two free-standing gates at 0.13 — it is "among
  jaccard-selected candidates, core rejects 155 that jaccard's own scale would
  never reach." The substantive finding (core catches 106 high-jaccard
  domain-sharers) stands; the rhetorical framing slightly oversells it.
- **Scope is narrow and must stay stated.** This measures only de-novo
  *within-chromosome* cross-copy singleton-exon merges (`merge_singletons_by_
  sequence`), per the OOM-safe one-chrom-at-a-time protocol. Cross-chromosome
  paralogs are not co-grouped (and are not this gate's target). The genome-wide
  drop measurement is from traces; the downstream-effect proof is **one contig
  (NC_073228.2) + the prior NC_073235.2** — not all 22 contigs run gate-on/off.
- **Net.** The gate is **strictly safe** (provably loses nothing real on the two
  contigs tested) and **principled** (valley-placed, separates a class jaccard
  cannot), but its demonstrated benefit at scale is **small precision hygiene**,
  not a metric win. Keep it default-off as a measured, optional precision lever;
  do not pitch it as a recall or F1 improvement.


---

## core_gate_gw_downstream

# Contiguous-core gate: downstream effect AT SCALE on the biggest-drop contig

Date: 2026-06-16. Confirms the (robustness-fixed) contiguous-core family-merge
gate's REAL downstream output effect on **NC_073228.2** — the chromosome with
the most would-be drops (54 of 135 cross-copy merge pairs have core_cov < 0.13,
40% drop rate; the biggest in the genome-wide measurement `core_gate_gw.md`).
Gate OFF vs ON inside rustle's actual `--vg --vg-layer2` pipeline.

This is the decisive companion to `core_gate_pipeline.md` (which ran the same
test on NC_073235.2 = only 5/67 drops). Question: does 11x more dropped merges
translate into 11x more downstream churn, or worse, lost real transcripts?

## Setup (OOM-safe, reproducible)

- Binary: `/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle`, built
  2026-06-16 16:30 from HEAD `9452964` (has the gate + the robust-POA fix).
- Contig: **NC_073228.2** (gorilla autosome, 195.3 Mb, **371,526 reads** — the
  largest contig; more reads than NC_073235.2's 188k).
- Slices from `/home/juanfra/winloci_scratch/GGO.{bam,fasta}`; guide
  `/tmp/gw/st_NC_073228.2.gtf` (StringTie -L, 4777 tx). `RAYON_NUM_THREADS=4`,
  serial, single process (~0.7 GB peak RSS observed; well under the 17 GB free).
- Gate env: `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (default 0.0 = OFF). ON value: **0.13**
  (the bimodal-valley threshold from `core_gate_gw.md`).
- `RUSTLE_VG_CORE_GATE_TRACE=1` on both runs (additive, never changes decisions).

Commands:
```
samtools faidx GGO.fasta NC_073228.2 > NC_073228.2.fasta; samtools faidx NC_073228.2.fasta
samtools view -b GGO.bam NC_073228.2 > NC_073228.2.bam; samtools index NC_073228.2.bam
# OFF
RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 rustle --vg --vg-layer2 \
  --genome-fasta NC_073228.2.fasta -G st_NC_073228.2.gtf -L NC_073228.2.bam -o cgon_off.gtf
# ON @0.13
RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
  rustle --vg --vg-layer2 --genome-fasta NC_073228.2.fasta -G st_NC_073228.2.gtf \
  -L NC_073228.2.bam -o cgon_on.gtf
```

## TL;DR

The gate **fires on exactly 54 merge pairs** (the predicted count), yet the
downstream output change is **smaller than on NC_073235.2** despite 11x more
dropped merges. **Zero transcripts gained or lost** (4989 = 4989). The ONLY real
change is at **one locus (RSTL.262)**: two transcripts lose their spurious
`copy_status "novel"` + `rescue_class "strand_pure_minority"` tags while keeping
identical coordinates and identical abundance. **Recall and precision vs the
reference are byte-identical** (FSM 1733 = 1733; transcript Sn/Pr 23.3/34.7 both
ways). This is benign cleaning — the gate's intended precision behaviour. **No
real isoform is lost.**

## Result 1 — the gate is active: 54 merges dropped

- Fresh OFF trace: 135 Jaccard-passing cross-copy singleton pairs;
  **54 have core_cov < 0.13** (52 deep-lobe < 0.05, 2 in the [0.05, 0.13)
  valley) — reproduces `core_gate_gw.md` exactly.
- ON @0.13: `would_gate=true` count = **54** — fires on precisely those pairs.
  The 81 true-copy pairs (core_cov >= 0.13; 48 in [0.31, 0.95), 29 >= 0.95) are
  untouched. Direction exactly as designed: drop short-block domain-sharers,
  keep real copies.

## Result 2 — downstream output delta (full contig)

| metric | OFF | ON | delta |
|---|---|---|---|
| transcripts | 4989 | 4989 | **0** |
| genes | 1497 | 1497 | 0 |
| `family_id`-tagged tx | 49 | 49 | 0 |
| `copy_id`-tagged tx | 49 | 49 | 0 |
| `copy_status "novel"` | 111 | 109 | **-2** |
| `rescue_class "strand_pure_minority"` | 62 | 60 | **-2** |
| FSM (=) vs reference | 1733 | 1733 | **0** |
| transcript Sn / Pr | 23.3 / 34.7 | 23.3 / 34.7 | 0 / 0 |

Method note: the raw line diff is large (~57k lines) but that is GTF
line-ORDER nondeterminism (rayon), not content. The numbers above come from a
**coordinate-keyed diff** (`/tmp/gw/diff_keyed.py`): each transcript keyed by
(strand, exon-chain). Keyed result:
- chains only in OFF (lost): **0**
- chains only in ON (gained): **0**
- chains in both: **4989**
- chains with a REAL attr change (copy_status / rescue_class / copy_id /
  source / reference_id / family_size): **2**
- chains with ONLY a cosmetic family_id renumber: **0**

So the 4989 = 4989 is NOT cancellation — it is literally the same transcript set
by coordinate, with 2 attribute edits.

## Result 3 — the one real change: RSTL.262 (48,724,206-48,742,369, + strand)

Same 5 transcripts in both runs, **identical exon coordinates and identical
cov / FPKM / TPM / abundance_min** values. The only difference: the two
transcripts sourced from `guide:STRG.440.1` and `guide:STRG.440.5` carry
`rescue_class "strand_pure_minority"` + `copy_status "novel"` with the gate OFF;
with the gate ON those two tags are **dropped**. (Transcript_id .1-.5 ordering
reshuffles, but the set keyed by source/coords/cov is identical.)

Interpretation: removing the spurious cross-copy domain-sharer merges that were
contaminating this family's variation graph cleaned up its false novel-copy
attribution. The transcripts are real and guide-backed (STRG.440.x) — they
**stay emitted**; only the unwarranted "novel copy" / rescue tags go away. This
is the exact analogue of the NC_073235.2 RSTL.647 result in `core_gate_pipeline.md`.

## Verdict

- **Benign-cleaning (expected), not harmful.** 0 transcripts lost, 0 real
  isoforms dropped (FSM 1733 = 1733; Sn/Pr unchanged). The only output change is
  removing 2 spurious `copy_status "novel"` / `rescue_class` tags at one locus —
  a precision improvement on copy attribution, not a recall cost.
- **Drop count does NOT scale to output churn.** 54 dropped merges (11x the 5 on
  NC_073235.2) produce FEWER downstream edits (2 tag cleanings, 0 family
  renumbers) than NC_073235.2 (1 locus cleaned + 2 cosmetic renumbers). The 54
  gated merges sit on family graphs whose final flow-extracted transcripts did
  not depend on those specific cross-copy exon fusions — confirming the gate
  prunes contaminating edges without disturbing the emitted assembly.
- **Honest caveat — the at-scale downstream effect is essentially negligible in
  magnitude** (2 tags on 1 of 1497 genes), even on the worst-case contig. The
  value of the gate is precision hygiene (it removes a real false novel-copy
  call), not a headline metric move. It is also strictly safe: nothing real is
  lost. Per-chrom / within-chromosome families only (OOM protocol); cross-chrom
  paralogs are not the gate's target.

## Artifacts

- GTFs: `/tmp/gw/cgon_off.gtf`, `/tmp/gw/cgon_on.gtf`
- Logs (with CORE_TRACE): `/tmp/gw/cgon_off.log`, `/tmp/gw/cgon_on.log`
- Keyed diff script: `/tmp/gw/diff_keyed.py`
- gffcompare: `/tmp/gw/gc_off.*`, `/tmp/gw/gc_on.*`
- BAM/FASTA slices removed after the run (regenerate via the commands above).


---

## core_gate_pipeline

# Contiguous-core gate: real-pipeline measurement (Phase C)

Date: 2026-06-16. Measures the (robustness-fixed) contiguous-core family-merge
gate inside rustle's ACTUAL `--vg --vg-layer2` family formation on a real
chromosome, gate OFF vs ON.

- Gate env: `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (default 0.0 = OFF). Operating
  value tested: **0.13**.
- Binary rebuilt from `src/rustle/vg_family/family_graph.rs` (the prior-phase
  Rust fix is in `contiguous_core_coverage` / `merge_singletons_by_sequence`).
- Contig: **NC_073235.2** (gorilla autosome, 150.2 Mb, 188,134 reads), full
  contig. Inputs: `/home/juanfra/winloci_scratch/GGO.{bam,fasta}`, guide
  `/tmp/gw/st_NC_073235.2.gtf`. RAYON_NUM_THREADS=4, serial, OOM-safe.

## TL;DR

YES — the gate changes real-data family formation, in the **expected direction**.
At threshold 0.13 it drops exactly **5 domain-sharing exon merges** (contiguous
core 0.004–0.011 of the shorter exon) while leaving all **62 true-copy merges**
(core 0.31–1.00) intact. Downstream the change is small and benign: one real
locus (RSTL.647) gets its copy/novelty attribution cleaned; two others only get
cosmetic family-id renumbering. No transcripts gained or lost (2623 = 2623).

## Where the gate actually lives (architecture)

The gate is NOT at the cross-locus family-GROUPING level (`discover_family_groups*`,
k-mer Jaccard 0.20 over bundles — that decides which loci join a `FamilyGroup`,
and the gate never touches it). It is at the **cross-copy singleton-exon merge**
inside `merge_singletons_by_sequence` (family_graph.rs:~416-445), invoked from
`assemble_family_graph_from_copy_exons` (family_graph.rs:538-540) only when
`n_singletons >= 2 && total_exons <= 2000`. It governs which homologous exons
from DIFFERENT copies fuse into one shared family-graph node (ExonClass / PSV
column).

Consequence: families whose copies overlap in genomic POSITION merge via Stage-1
`cluster_by_position` and never see the gate. Only the dispersed-paralog
sequence-merge path (Stage-1b) is gated. The `[DIAG] family ... member_spans`
dump (gated by `RUSTLE_DIAG_TARGET`) reports family→loci membership but is
therefore the WRONG instrument for this gate; the right instrument is the
per-pair core-coverage at the merge loop.

## Instrumentation added (additive, default-off)

Added `RUSTLE_VG_CORE_GATE_TRACE` in `merge_singletons_by_sequence`: for every
Jaccard-passing cross-copy pair it prints
`[CORE_TRACE] cid_i= cid_j= jacc= core_cov= len_i= len_j= would_gate=`.
Computes `core_cov` for the trace even when the gate is off (so the full
distribution is visible); never changes merge decisions. Family-graph unit tests:
17 pass / 0 fail.

## Result 1 — core-coverage distribution of Jaccard-passing pairs (gate OFF)

67 Jaccard-passing cross-copy singleton pairs on the contig. Bimodal:

  core_cov < 0.05 : 5     <- domain-sharers (0.004, 0.008, 0.011, 0.011, 0.011)
  core_cov < 0.13 : 5     (no pairs between 0.011 and 0.314 — clean empty gap)
  core_cov >=0.31 : 62    true copies (mean of all 67 = 0.739)

The threshold 0.13 sits squarely in the empty bimodal gap.

The 5 low-core pairs (all in ONE family, copies cid=0 vs cid=2):

  jacc=0.386 core_cov=0.004 len_i=1140 len_j=2846
  jacc=0.636 core_cov=0.008 len_i=808  len_j=662
  jacc=0.707 core_cov=0.011 len_i=180  len_j=233   (x2)
  jacc=0.605 core_cov=0.011 len_i=180  len_j=233

These are textbook domain-sharers: enough shared 15-mer minimizers to clear the
Jaccard bar, but the longest contiguous homologous run covers <1.1% of the
shorter exon, and lengths are mismatched (e.g. 1140 vs 2846 bp). Not paralog
copies.

## Result 2 — gate ON @0.13 fires on exactly those 5 pairs

With `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13`, `would_gate=true` count = 5 —
precisely the 5 sub-0.013 pairs. The 62 true-copy merges are untouched. Direction
is exactly as designed: drop short-block domain-sharers, keep real copies.

## Result 3 — downstream output delta (OFF vs ON), full contig

Transcript count identical: 2623 vs 2623. GTF differs (md5 differs); the diff is
3 genes:

- **RSTL.447** (FAM_13 -> FAM_11) and **RSTL.508** (FAM_0 -> FAM_4): ONLY the
  `family_id` integer label changed. All coordinates, copy_id, abundances,
  family_verdict, copy counts byte-identical. Cosmetic renumbering caused by
  removing 5 merges upstream (family enumeration order shifts). Not a real
  formation change at these loci.

- **RSTL.647** (NC_073235.2:112,808,333-112,845,893, the domain-sharer family,
  cid0/cid2): REAL effect. Same 4 transcripts, same exon coordinates, but
  transcripts are re-ranked and the `rescue_class "strand_pure_minority"` +
  `copy_status "novel"` tags present with the gate OFF are DROPPED with it ON.
  Removing the spurious cross-copy merge that was contaminating this family's
  graph cleaned up its novel-copy attribution. This is the gate's intended
  precision behaviour on a real locus.

## Verdict

- Did the gate change family formation on real data? **Yes** — 5 merges removed
  on this contig.
- Expected direction (drop domain-sharers, keep copies)? **Yes, exactly** — all
  5 removed pairs are <0.013 contiguous core (domain-sharers); all 62 true-copy
  pairs (>=0.31) are preserved. Bimodal gap at 0.13 is clean.
- Output cost? Minimal/benign: 0 transcripts gained or lost; 1 locus
  (RSTL.647) cleaned; 2 cosmetic family-id renumbers.
- The gate is narrow on this contig because (a) only the dispersed-paralog
  Stage-1b merge path is gated (co-located families bypass it), and (b)
  genuine domain-sharers passing Jaccard are rare (5 of 67 pairs here, all in
  one family).

## Reproduce

  # rebuild
  cargo build --release
  # contig slice
  samtools faidx GGO.fasta NC_073235.2 > NC_073235.2.fasta; samtools faidx NC_073235.2.fasta
  samtools view -b GGO.bam NC_073235.2 > NC_073235.2.bam; samtools index NC_073235.2.bam
  # OFF (trace shows full distribution)
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 rustle --vg --vg-layer2 \
    --genome-fasta NC_073235.2.fasta -G st_NC_073235.2.gtf -L NC_073235.2.bam -o off.gtf
  # ON @0.13
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
    rustle --vg --vg-layer2 --genome-fasta NC_073235.2.fasta -G st_NC_073235.2.gtf \
    -L NC_073235.2.bam -o on.gtf
  # compare: grep -c would_gate=true; diff off.gtf on.gtf

## What a genome-wide run would need

Per-contig serial (JOBS=1, the OOM protocol), 24 contigs, ~25-35 min total,
mirroring `bench/copy_recovery_eval/results_genomewide/gw_run.sh`. Run each
contig twice (gate OFF/ON) with `RUSTLE_VG_CORE_GATE_TRACE=1`; aggregate
`would_gate=true` counts and the core_cov histogram across contigs, and diff the
two GTF sets per contig (separating cosmetic family-id renumbers from real
structural/attribution changes, as done here). Expected scale: a handful of
domain-sharer merges per contig, concentrated in a few families.


---

## psv_llr_vs_votes_killtest

# Kill-or-validate: what actually beats the production vote-rule on copy assignment

The workflow's first recommendation was "wire the LLR (log-likelihood) engine into production." A
labeled-data decomposition (`psv_llr_vs_votes_killtest.py`, ground-truth sim copies) shows that the
**likelihood engine *per se* buys nothing** — the real lever is the **admission gate**. Three effects,
isolated on identical reads:

| effect | comparison | result |
|---|---|---|
| **SCORING** (votes → likelihood) | `VOTE_gate1` vs `LLR_flat` | **identical in 16/16 configs.** With flat per-base error the LLR margin = vote_margin·ln(3(1−e)/e), so likelihood ≡ vote-count. Switching scorers alone = 0 gain. |
| **GATE** (`min_psv=3` → `n_decisive≥1`) | `VOTE_prod` vs `VOTE_gate1` | **the dominant lever.** A family with only 1–2 PSVs can never reach 3 votes, so the production rule **discards 100% of the near-identical tail** (recall 0.000). The `n_decisive≥1` gate resolves it: +1151–2901 reads recovered at **94.7–100% precision**. |
| **QV-INFO** (flat vs calibrated per-base QV) | `LLR_flat` vs `LLR_cal` | the **only** thing votes structurally cannot do. Real but **modest and high-error-only**: at normal HiFi error it removes ~0 misassignments; under stress (err×3) it ~halves them (e.g. C5K1 97→38 wrong, prec 0.947→0.978) at a small recall cost. A no-op at realistic QV. |

## The reframed lever

The headroom over Eichler's AS≥10 gate (measured: ~28% of discarded GGO reads recoverable) AND over
Rustle's own production rule comes from the **same place** — the **1–2 PSV near-identical tail** that
both gates throw away:

- Eichler discards it because the *aggregate* alignment-score margin is <10 (signal diluted over the
  full molecule).
- Rustle-production discards it because `min_psv=3` (a 3-concordant-PSV hard floor) can't be met.

The fix is one principled gate: **`n_decisive≥1` (identifiability) + a single calibrated decisive
margin τ = ln((1−p)/p)** that sets the precision/recall operating point by a *stated* misassignment
bound P(misassign)≤p — the Canzar-clean restatement of Eichler's conservative AS≥10, with the magic
integers (`min_psv=3`, `margin=1`) collapsed into one interpretable knob.

The LLR engine (`copy_assign.rs::assign_read`) is the right **vehicle** for τ (it makes τ a probability
bound and can fold in QV + junction columns), but the **mechanism that recovers reads is the gate, not
the likelihood**. Sell it that way.

## What this does NOT cross

The `n_decisive=0` floor is untouched: reads spanning no column where candidate copies differ
(exonically identical silent copies; single-shared-exon reads) stay discarded under every arm. Every
lever here re-adjudicates the **resolvable** tail more correctly; none manufactures a decisive feature.

## Operating-point evidence (the one knob)

Resolving the 1–2 PSV tail is a precision/recall trade, governed by τ:
- K=2 families: recovered at ~99–100% precision (cheap, near-free win).
- K=1 families: recovered at ~95–98% precision (one PSV → one error can flip it; this is where τ and QV
  earn their keep). Production's `min_psv=3` is just the τ→∞ corner (discard the whole tail).

Artifacts: `psv_llr_vs_votes_killtest.py` · `psv_llr_vs_votes_killtest.png` · `.json`.


---

## psv_tau_sweep

# Operating point: precision/recall vs the decisive margin τ

τ is now a real parameter of the shipped engine (`copy_assign --margin <τ>`, default 2.0; the doc-string
gives τ = ln((1−p)/p), τ=6.9 ≈ p 1e-3 = the Eichler AS≥10 analog). The identifiability gate
(`n_decisive≥1`) is independent of τ and always applied. This sweeps τ on two substrates using the
production engine (`copy_assign.rs::assign_read`, mirrored by `bench/copy_assign.py`).

## Sim5x K-ladder — TRUE labels (read name = source copy)

| regime | behaviour |
|---|---|
| **K ≥ 2** | trivially perfect (recall 1.0, precision ~1.0) for every τ up to ~7. Multi-PSV reads are decisive. |
| **K = 1** | the only τ-sensitive case. argmax (τ→0): recall 1.0 but **precision 0.80** (genuine ties tie-broken wrong). τ ≥ 0.5: recall 0.60 / **precision 1.0** (ties correctly discarded). |

A single clean HiFi PSV gives margin = ln(3·0.997/0.003) = **6.90**, so **τ=6.9 (p=1e-3) sits exactly at
the one-confident-PSV knife-edge** — an independent calibration check that τ = ln((1−p)/p) is the right
parameterisation.

## Real GGO co-located families — 47,732 reads / 70 families (silver-standard agreement)

| τ | recall (assigned/all) | silver-standard agreement |
|---|---|---|
| 0.0 (argmax) | 0.976 | 0.9340 |
| 0.5 | 0.964 | 0.9411 |
| **2.0 (default)** | **0.964** | **0.9412** |
| **6.9 (p=1e-3, Eichler-faithful)** | **0.964** | **0.9412** |
| 12.0 | 0.961 | 0.9416 |
| 25.0 | 0.958 | 0.9419 |

**The real-data curve is nearly flat.** Two readings:

1. **τ=2 and τ=6.9 give the *identical* operating point** (recall 0.964, agreement 0.9412). The
   principled conservative choice — τ=6.9, a *stated* per-read misassignment bound p≤1e-3, the PSV-space
   restatement of Eichler's AS≥10 — costs **nothing** in recall over the loose default. Adopt it for free.
2. **The only material decision is τ=0 vs τ>0.** Argmax (assign ties) buys +1.2% recall but *loses*
   ~0.7% agreement; any τ>0 discards the ties and is strictly better. Beyond that, the margin distribution
   on real data is **bimodal** — reads are either decisive (margin ≫ 25) or genuine ties (margin ≈ 0),
   with little thin-margin middle — so τ is a robust, insensitive knob. The lever is the **gate**
   (`n_decisive≥1`, already in this engine), not fine-tuning τ.

## Recommendation for the meeting

Default **τ = 6.9 (p_target = 1e-3)** + **discard ties (τ>0)**: a principled, Canzar-clean operating
point with a stated misassignment guarantee, free on real data, and the exact PSV-space analog of the
AS≥10 criterion your advisor cites. The curve is the evidence; the point is robust.

## Honest caveats

- "Silver-standard agreement" treats each read's best-overlap copy as truth — exactly what is unreliable
  for hard multimappers — so ~94% is a *proxy*; some of the ~6% disagreement is PSV *correcting* the
  overlap call, not error. The sim panel (true labels) is the rigorous anchor.
- 70 families (≥2 valid copies after the spliced-length filter), not genome-wide; capped at 150
  candidates. Representative of the co-located tandem regime, not a full census.
- Flatness is consistent with the GGO headroom measurement (Eichler-discard only 3.4% of multimappers —
  most gorilla paralogs are divergent enough that the margin is decisive): the thin-margin tail where τ
  would matter is genuinely small in this family set.

## Rust cross-check

`copy_assign --margin {2.0,6.9,12.0}` on `NC_073242.2:3771193-3799186` (a 2-copy de-novo family): runs
end-to-end on real GGO, 18 reads all ASSIGNED at every τ (margins 20–330, the decisive regime) — confirms
the engine builds, honors τ, and applies the identifiability gate.

Artifacts: `psv_tau_sweep.py` · `psv_tau_sweep_fig.py` · `psv_tau_sweep.png` · `.json` ·
`src/bin/copy_assign.rs` (--margin/--error-rate).


---

## copy_specific_junctions

# Copy-specific junctions — do paralog COPIES splice differently?

The unification: (1) discovered cross-chromosome copies → (2) reads assigned per copy (by locus here; PSV-based copy_split applies in the collapsed regime, data-limited in GGO) → (3) per-copy junction usage compared between copies (the ASJ machinery, allele→copy).

- cross-chrom recent-dup pairs (core_ident≥0.9, both ≥3 exons, BOTH copies expressed ≥15 reads): **248**
- homologous junctions compared between copies: **696**
- copy-specific junctions raw (q<0.05, |dPSI|≥0.3): **191**, of which:
  - **DIFFERENTIAL splicing (both copies splice): 146** across **66** copy pairs
  - one copy UNSPLICED/retrocopy (trivial dPSI): 45
- copy-private exon junctions (one copy includes an exon the other lacks): **1408**

## Top DIFFERENTIAL copy-specific junctions (both copies splice, differ at this junction)
| copy A | jA | PSI_A | copy B | jB | PSI_B | dPSI | q |
|---|---|---|---|---|---|---|---|
| LOC101147259 | 93534519-93538824 | 0.0 | LOC115933254 | 32884967-32889244 | 1.0 | 1.00 | 1.1e-19 |
| HERC2 | 29471916-29472174 | 1.0 | LOC101149126 | 68795812-68796072 | 0.0 | 1.00 | 1.8e-03 |
| HERC2 | 29473320-29476000 | 1.0 | LOC101149126 | 68791948-68794664 | 0.0 | 1.00 | 4.3e-03 |
| LOC115932779 | 28925295-28926477 | 0.0 | SORL1 | 129494627-129495811 | 0.995 | 0.99 | 1.1e-29 |
| LOC115933039 | 21007327-21008509 | 0.0 | SORL1 | 129494627-129495811 | 0.995 | 0.99 | 6.7e-28 |
| LOC101124243 | 6830621-6832115 | 0.0 | NOM1 | 165843104-165844598 | 0.974 | 0.97 | 2.0e-25 |
| LOC101147259 | 93526275-93527224 | 0.0 | LOC115933254 | 32876845-32877795 | 0.964 | 0.96 | 2.0e-19 |
| LOC101125415 | 17453004-17453086 | 0.961 | LOC101139163 | 141285207-141285366 | 0.0 | 0.96 | 2.3e-17 |
| LOC101125415 | 17453004-17453086 | 0.961 | LOC101152938 | 75725353-75725435 | 0.0 | 0.96 | 4.8e-12 |
| FRG1 | 210904802-210905583 | 0.955 | LOC115933219 | 42792484-42793264 | 0.0 | 0.95 | 2.1e-32 |
| LOC129529226 | 15687227-15687359 | 0.951 | PLA2G4B | 42658550-42658680 | 0.0 | 0.95 | 1.1e-17 |
| LOC129531205 | 4579716-4582007 | 0.0 | SAFB2 | 11378953-11381252 | 0.949 | 0.95 | 2.9e-17 |
| LOC101137900 | 32597303-32597305 | 0.0 | LOC129531862 | 31668721-31671325 | 0.941 | 0.94 | 3.5e-07 |
| LOC101142904 | 2232115-2232899 | 0.935 | LOC115933219 | 42792484-42793264 | 0.0 | 0.94 | 3.4e-21 |
| LOC129531982 | 29073750-29073882 | 0.931 | PLA2G4B | 42658550-42658680 | 0.0 | 0.93 | 6.3e-15 |
| LOC101125415 | 17453352-17454181 | 0.962 | LOC101139163 | 141284114-141284941 | 0.038 | 0.92 | 6.4e-16 |
| DPY19L2 | 37744041-37746734 | 0.922 | LOC115935345 | 37139557-37143806 | 0.0 | 0.92 | 5.3e-09 |
| DAZ1 | 42826894-42827637 | 0.074 | DAZL | 26430807-26431538 | 0.992 | 0.92 | 2.6e-151 |
| LOC101142904 | 2223093-2224079 | 1.0 | LOC129526395 | 24713570-24714553 | 0.087 | 0.91 | 2.3e-17 |
| LOC129531205 | 4578243-4579576 | 0.071 | SAFB2 | 11381383-11382713 | 0.974 | 0.90 | 8.7e-17 |

## Honest notes
- BOTH copies expressed is rare (~3% of recent-dup pairs) — most paralogs have one dominant copy (consistent with the headroom findings). So this is a small, real set, not a large catalog.
- read→copy assignment here is by genomic locus (cross-chrom copies map to distinct coords). The PSV-based copy_split assignment (interest #2) is exercised only in the COLLAPSED regime, which is essentially empty in GGO (task c) — it needs a deep co-located dataset (testis HiFi for DAZ/RBMY).
- homologous-junction mapping is via exon-length alignment; copy-private exons flag structure divergence (a copy gained/lost an exon).
- coverage-limited; deeper data finds more both-expressed pairs.

## Reproduce
- `MINIFORGE python bench/copy_specific_junctions.py`


---

## copy_assign_phasing

# Facultative long-read phasing (`copy_assign --phase`) — dependency-free

The PSV harness already phases reads into copy-haplotypes internally (`copy_split::split_readchain_by_psv`
+ `psv_linkage::assign_read_to_copy`). `--phase` surfaces that as a first-class, **opt-in,
dependency-free** phasing output — no external phaser (WhatsHap/HapCUT2), no neural variant caller
(DeepVariant/Clair3), no CNN. It is the N-copy generalisation of read-backed phasing: a long molecule
spans multiple PSVs, the linked pattern assigns it to a haplotype, and the formal object is the
minimum path-cover of the PSV graph (MCC = χ(H), the identifiability theorem).

## Outputs (written only under `--phase`; all other outputs unchanged → additive/default-off)

- **`<out>.phase.gfa`** — the phasing AS A VARIATION GRAPH (self-contained GFA1, loadable in
  Bandage / `vg`): each PSV column is a **bubble** (one `S` segment per allele, carrying the base +
  `PO:i:<genomic_pos>`); each copy is a **path** (`P` line) threading the bubbles. Copies that share an
  allele share the node (a bubble anchor) and fork where they differ — so the graph topology *is* the
  copy structure. PSV-identical copies collapse to **one path** = the K-frontier rendered as graph.
- **`<out>.phase_blocks.tsv`** — one **phase set (PS)** per family:
  `block_id  chrom  n_haplotypes  n_psv_sites  n_reads_phased  n_unphased`
- **`<out>.phased_haplotypes.tsv`** — each haplotype's phased alleles:
  `block_id  haplotype  copy_tid  n_support_reads  variants` where `variants` = `pos:allele;pos:allele;…`
- **`<out>.phased_reads.tsv`** — the **haplotag (HP)** per read:
  `read_name  block_id  haplotype  n_psv_spanned  margin  status`
  (`haplotype = -1` when the read is Ambiguous/Tied = unphaseable = the K-frontier).

`block ≙ PS`, `haplotype ≙ HP`, `bubble ≙ PSV`, `path ≙ copy` — standard phasing/VG terms. The `.gfa` +
`.phased_reads.tsv` are jointly self-contained: the graph holds copies-as-paths-through-PSV-bubbles, the
haplotag maps each read to its copy/path.

### GFA demo (NC_073242.2 family: 3 copies, 43 PSV columns)

```
phase.gfa: 86 bubble-nodes, 84 links, 3 copy-paths
S  CAFAM0_c0_C  C  PO:i:3790785     S  CAFAM0_c0_G  G  PO:i:3790785   <- the bubble (2 alleles)
P  CAFAM0_copy0  c0_C+,c1_G+,...  ┐ identical path (PSV-identical copies collapse = K-frontier)
P  CAFAM0_copy1  c0_C+,c1_G+,...  ┘
P  CAFAM0_copy2  c0_G+,c1_C+,...    the divergent copy forks at every bubble
```

## Demo (NC_073242.2:3771193-3799186, a clean tandem family)

```
phase_blocks:  CAFAM0  NC_073242.2  3 haplotypes  43 PSV sites  18 phased  0 unphased
haplotypes:    hap2  2 reads  3790785:G;3790787:C;…   (PSV-distinguishable; margin 330)
               hap0  8 reads  3790785:C;3790787:G;…   ┐ identical PSV alleles, resolved by
               hap1  8 reads  3790785:C;3790787:G;…   ┘ LOCUS/read-coherence axis (margin 20)
phased_reads:  m64076…/161284984/ccs  CAFAM0  hap0  45 spanned  margin 20.0  assigned
```

Two axes at once: PSV linkage (hap2) **and** read-coherence/locus (hap0 vs hap1, PSV-identical copies).
On a pathological size-heterogeneous family (PCDHB container-locus, 116 kb backbone + ~4 kb paralogs) it
phases **0/15,386** reads — correctly **abstaining** rather than guessing (DAZ3/Canzar discipline).

## Properties

- **Dependency-free**: pure Rust over data already computed (`copy_psv_alleles`, `psv_col_pos`,
  `assignments`); no new crates, no external tools, no ML (grep-confirmed: zero torch/onnx/sklearn/CNN in
  the harness).
- **Opt-in / additive**: gated behind `--phase`; without it, the families/assignments/quant outputs are
  byte-identical.
- **Confidence-gated**: a read is phased iff it clears the decisive-margin gate τ = ln((1−p)/p)
  (default 6.9 ≈ p 1e-3, the Eichler AS≥10 analog); unphaseable reads are emitted as `-1`, never guessed.
- **Phases copies, not alleles**: it phases into paralog copy-haplotypes (the thesis problem), not the
  two parental allelic haplotypes; allele-level het is the separate ASJ machinery.

Implementation: `src/bin/copy_assign.rs` (`--phase` flag + emission). Build: `cargo build --release --bin copy_assign`.


---

## copy_resolution_boundary

# The copy-resolution boundary: which multi-copy gene families can RNA resolve, and which it provably cannot

*The empirical capstone to the copy-assignment identifiability theory (`bench/copy_assignment_theory.md`).
Measured on GGO HiFi IsoSeq (FLNC, full-length mature mRNA) against the T2T reference. A clean two-tier
boundary: a resolvable majority (exonic PSVs), and a K=0 residual that is provably unresolvable per-read —
splice-site divergence is real in the reference but aligner-masked, so the entire K=0 residual co-quantifies.*

## Bottom line

Per-read copy assignment **works for the large majority** of multi-copy families and **provably cannot work**
for a small, stereotyped residual — and we can tell them apart in advance. The hard cases are not a broad
failure of RNA; they are recent inverted segmental duplications (largely X-linked cancer-testis-antigen genes)
whose entire mature transcript is identical between copies, so no spliced read carries copy information.

## The census (definitive prevalence)

Across all 1130 de-novo multi-copy families we formed all co-located copy pairs and kept the
**assignment-relevant** ones (≥3 cross-mapping MAPQ-0 reads, so copy assignment is actually at issue):
**236** pairs across **97** families (**42,313** classified reads). Each pair is K=0 (per-read unresolvable)
iff its cross-mapping reads have `NM_A == NM_B` (the spliced read aligns equally to both copies' exons).
Script + per-family results: `bench/copy_resolution_census.py`, `bench/copy_resolution_census.tsv`.

| class | pairs | fraction |
|---|---|---|
| **resolvable** (exonic PSVs, `NM_A ≠ NM_B`) | 206 | **87%** |
| K=0 (≥95% reads NM-identical) | 30 | 13% |
| K=0 strict (100% identical) | 7 | 3% |

**84%** of all 42,313 cross-mapping reads carry resolving evidence. K=0 is monotonically confined to the
highest-identity tail: K=0 pairs have **median 0.000% exonic divergence** vs **0.42%** for resolvable pairs.

**What the K=0 residual *is* (30 pairs / 11 families):** overwhelmingly the X-chromosome (NC_073247.2)
cancer-testis-antigen region — `MAGEA9` and the `LOC129529xxx` / `LOC1011xxxxx` tandem clusters (DNFAM1,
DNFAM14, DNFAM160, …), plus a few autosomal recent dups (`ZNF793`, scattered LOC pairs). These are recent
**inverted** segmental duplications. The big co-located families one might fear — protocadherins, keratins,
tubulins, APOBEC, APOL, BTN2A — are all in the *resolvable* majority.

## The two tiers

### Tier 1 — Resolvable majority (~87% of pairs / ~84% of reads)

The copies carry exonic PSVs; a read covering one fits its true copy better (`NM_A ≠ NM_B`). This is the regime
of **Theorem 2/3** (`copy_assignment_theory.md`): under the identifiability condition the true copies are the
unique minimum cover and are recovered in polynomial time. The shipped per-read PSV + junction assigner handles
this tier (sim K-ladder K≥2 → 100%, GGO silver-standard 100%). **This is the headline: the method works.**

### Tier 2 — K=0 residual: splice divergence is real in the reference but per-read-masked (~13% of pairs)

The exons are identical, so per-base PSVs give 0. The *introns* do diverge (4–60 SNPs/indels across the
duplicated block), and where that lands at a splice site it shifts a junction (e.g. MAGEA pair0's 3 bp
intronic indel makes copy A's intron 3702 nt and copy B's 3705 nt). One might hope a read using that junction
betrays its copy — but it does **not, per read**: minimap2's spliced aligner independently snaps every junction
to the nearest canonical `GT-AG` at *each* copy, so the read is fully canonical at **both** copies and, because
the indel lives inside the long intron (never in the spliced mature read, and absorbed without an exonic edit),
NM stays tied. Measured on the panel (`bench/splice_divergence_resolver.py`): per-read splice **resolution = 0
for every pair, including pair0** — the divergence is *detectable* (pair0's reads carry a different intron-length
chain at A vs B) but **non-directional** (it cannot say which copy a read came from). The earlier "pair0 ~33%"
was a reference-level count of junctions whose homologous site differs, not a per-read assignment.

So there is **no Tier-2 per-read rescue**: the K=0 residual collapses into Tier-3. (The reference-level splice
distinguishability remains relevant to copy-*model* / family-graph work — interest #3 — but not to assigning an
individual mature read.)

### Tier 3 — The K=0 core: co-quantify (no per-read assignment)

Exons identical **and** splice sites identical: the byte-identical block extends through the splice signals and
into the proximal intron, so every junction is `GT-AG` at both copies (`MAGEA pair2/pair3`-class). The 4–60
distinguishing SNPs all sit in introns/flanks a spliced FLNC read never covers; verified **0 of the 372
junction-reads at the splice-identical pairs (pair2 300 + pair3 72) resolvable**
(`bench/splice_divergence_resolver.py`). Because FLNC reads are already full-length mature transcripts,
**no read-length or chemistry improvement helps** — the entire mature transcript is identical, so the information is simply not in the RNA.
This is the identifiability theorem's K=0 floor instantiated: the distinguishing-column set is empty.

For this core the correct output is **not** a forced per-read assignment but a **co-quantified ambiguity set** —
report the copies as a shared-evidence bundle under the family graph and estimate per-copy expression *ratios*
(with a DNA/copy-number prior), rather than pretending to assign individual reads. With Tier-2 shown to be
per-read-masked, the *entire* K=0 residual (not only the strict core) lands here: report a co-quantified
ambiguity set, not forced per-read assignments.

## Why this is sound (not an artifact)

The K=0 copies are genuinely distinct loci, verified against the T2T reference: distinct annotated genes
(`MAGEA9` vs `LOC129530018`), independent expression (own primary reads at each), **inverted** orientation
(B = revcomp of A — an isoform cannot be its own reverse-complement), and — decisively — each is a **bounded
44–75 kb homology block embedded in divergent unique genomic context**, the signature of a real segmental
duplication rather than an assembly false-duplication or same-locus comparison.

## The boundary as the result

```
multi-copy family with cross-mapping reads
        │
        ├─ exonic PSVs (NM_A≠NM_B)         → Tier 1  ~87%  → per-read PSV assignment (Thm 2/3)        [interest #2]
        │
        └─ K=0 (exons identical)            ~13%
              → splice divergence is reference-real but per-read-masked (minimap2 snaps junctions)
              → Tier 2 has NO per-read rescue → co-quantified ambiguity set (Tier 3)                  [interest #3]
```

This unifies the three advisor interests under one measured boundary: detection (#1) forms the family, PSV
assignment (#2) resolves the majority, the K=0 residual is shown per-read-unresolvable and co-quantifies
under the family graph (#3 — allele-specific junctions remain relevant at the copy-*model* level, not per-read
assignment), and the same identifiability quantity (the distinguishing-column count) governs all three.

## Reproducibility

- Census: `python3 bench/copy_resolution_census.py` (classifier) → `bench/copy_resolution_census.tsv` (per-family).
- Theory: `bench/copy_assignment_theory.md` (+ `bench/copy_assignment_theory_checks.py`).
- K=0 frontier attack + splice-site mechanism: workflows `wf_e36e46e4-aa1`, `wf_0c7be571-7ed`;
  `bench/resolution_improvement_bound.md`.

## Open / next

- Extend the census beyond the ≥3-read assignment-relevant subset for a fully exhaustive per-family K=0 label.
- Splice-divergence is per-read-masked by aligner junction-snapping (`bench/splice_divergence_resolver.py` — negative-result probe); the K=0 residual co-quantifies (Tier-3).
- Formalize the Tier-3 co-quantification (per-copy expression ratio under the family graph + copy-number prior).


---

## identifiability_bound

# The identifiability bound on copy detection (3 independent confirmations)

The method recovers every copy the read data can DISTINGUISH; the rest are unidentifiable by definition.
Confirmed on three independent fronts this session:

1. DIVERGENT PARALOGS (Phase 1, divergent_phase1_features.py): on a gold Compara set (781 paralog /
   1294 non-paralog), no RNA sequence feature separates true divergent paralogs from domain-sharers
   (best AUC 0.629, protein; DNA <=0.56). Distinguishable only with phylogenetic-level input.

2. COLLAPSED COPIES (--recover-copies, recover_collapsed_copies): feeding AS-tied secondary reads +
   running copy_split PAST the family gate recovers 0 NEW copies in no-family collapsed arrays
   (DAZ1/RABL2A/GOLGA6L7/PRAMEF). Distinguishable copies form families (in-family split already gets
   them); collapsed copies are not distinguishable (that is why they collapsed). Exhaustive + complementary.

3. PRIOR VG work: "decisive secondary signal self-contradictory"; multimapper advantage confined to
   tied co-located copies which are by definition non-separable.

THESIS STATEMENT: copy detection/assignment is COMPLETE up to identifiability; identifiability is
information-theoretic (does a copy carry a copy-specific PSV/junction?), not algorithmic. Engineering
(secondary reads, longer k-mers, masking, protein-ORF) cannot cross it. The in-family PSV split realizes
the bound. Transcripts themselves are validated real (bench/transcript_validation.md: 86% intron precision,
98.9% gene-associated), so the bound is about copy SEPARABILITY, not transcript quality.


---

## eval_significance_gate

# Evaluation — IsoCon significance gate vs the τ-margin gate

The copy-assignment decision was changed from a fixed log-LR margin τ to IsoCon's real-vs-error significance
test (`copy_assign::assign_read`, default gate). This note records the validation. Reproduce with the
`copy_assign` binary flags `--alpha <α>` (significance gate, default 1e-3) and `--margin-gate --margin <τ>`
(legacy τ-gate).

## 0. The guarantee (what α controls)

`best` is the argmax over `n` copies, so a per-competitor p-value only bounds the *pairwise* error; selecting
the winner inflates the family-wide rate to the union bound `(n−1)·α`. The gate applies a **Bonferroni
correction** — each competitor must clear `α/(n−1)` — so the **family-wide misassignment rate over ASSIGNED
reads is bounded by α** (a theorem, validated in `sig_gate_bonferroni_controls_familywide_error`: in a
10-competitor 1-PSV regime the corrected gate keeps realized ≤ α by abstaining, while the un-corrected
per-pair level inflates to ~K·(e/3) > α). For `n=2`, `thr=α` (unchanged).

## 1. Unit calibration (Rust test `sig_gate_is_calibrated_realized_error_tracks_alpha`)

Simulated reads from a known copy over 6 Q30 PSVs (n=2), 20k reads per α ∈ {1e-2, 1e-4}: the gate assigns
~100% with **0 realized misassignments** at both α (the bound holds with slack). Confirms the certificate is
conservative and high-recall when PSVs are present.

## 2. sim5x labeled-truth ladder (ground truth, NOT the circular silver standard)

`bench/build_sim5x.py` builds a 5-copy tandem gene with a private-exonic-PSV ladder K=0..8 and reads named
with their true copy (`K{K}_c{copy}_r{n}`). Per-read assignment scored against the label (detected→true copy
mapped by majority vote), significance gate (α=1e-3) vs legacy τ-gate (τ=6.9):

| K | assigned% (both gates) | tied% | accuracy of assigned | misassigned |
|---|---|---|---|---|
| 0 (exonically identical) | **0%** | **100%** | — | 0 |
| 2 | 100% | 0% | 100% | 0 |
| 4 | 100% | 0% | 100% | 0 |
| 8 | 100% | 0% | 99% | 2 / 200 |

**Findings:**
- **The significance gate is identical to the validated τ-margin gate** at every K — same recall, same
  accuracy, same 2/200 errors at K8. So it is a drop-in, no-regression replacement: it reproduces the
  validated behaviour while adding the calibrated per-read certificate and a cleaner identifiability rule.
- **K=0 is correctly 100% Tied.** Exonically-identical copies have no distinguishing PSV, so every
  competitor has an empty distinguishing set → `min_p_value = 1.0 ≥ α` → `Tied`. The identifiability
  certificate flags the unresolvable regime *by construction*, replacing the heuristic `n_decisive ≥ 1`
  with the power-calibrated `min_p_value < α`.
- The 2/200 misassignments at K8 are produced by **both** gates (including τ=6.9 ≡ p=1e-3), so they are a
  property of the data (reads whose few spanned PSVs match a sibling after error), not a miscalibration of
  the significance gate.

## 3. Default α

**α = 1e-3 confirmed** as the default: it reproduces the production τ=6.9 (p=1e-3 / Eichler-AS≥10) operating
point on the labeled ladder, with the precision behaviour the assign-or-abstain philosophy wants. Tunable
via `--alpha`.

## Status / follow-ups
- The labeled sim5x ladder (ground truth) is the load-bearing validation here; a genome-wide GGO re-run
  would only add circular silver-standard numbers and is deferred (runnable via `copy_assign --regions … `
  with/without `--margin-gate`).
- **L17 closed:** `<out>.assignments.tsv` now emits `p_value` and `min_p_value` (the per-read certificate +
  identifiability bound) as the final two columns.
- **Honest limits:** the εⱼ = e/3 model is uniform-error (anti-conservative for biased substitutions, esp.
  unfiltered A→G **RNA editing** at A/G PSVs — the deferred Clair3-RNA editing filter addresses this).


---

## resolution_improvement_bound

# Would resolving multimapper ambiguity improve families and their isoforms?

*Bridging advisor interests #1 (family definition) and #2 (copy assignment): the achievable improvement from
resolution is bounded by RNA-level identifiability, and on the GGO panel the bound bites to ~0.*

Measured on `GGO.bam` over a 7-family panel + the chrY mega-array (workflow `wf_ad5ebaa7-0fb`, 4 measurements
+ synthesis). The answer separates into two questions that resolve differently.

## Families: detection unchanged, completeness +3 copies

- **Detection (membership) is unimproved — by construction.** The conflict graph is built from the ambiguity
  itself; resolution decides *where* an ambiguous read goes, not *which* reads are ambiguous. Family membership
  is defined prior to resolution.
- **Completeness improves via starved-copy recovery, but modestly.** A copy that loses the arbitrary MAPQ-0
  primary-flag battle to a tied sibling has ~no primary reads and is invisible to primary-only assembly.
  Reclaiming its secondaries (= resolution) recovers it. Found **3 starved copies** (gate: primary < 3, secondary
  ≥ 3): `LOC115934278` (1 prim / 23 sec; sibling AK6 hoards 848 prim) + 2 chrY-array loci
  (`22450601` 2/57, `22474453` 2/91). The 2 chrY copies are the strongest case — their reclaimed secondaries are
  100% MAPQ-0, ~90% spliced, with the arbitrary primary parked at a tied in-array sibling, so resolution recovers
  **both the copy and its multi-intron isoform**. (Lower bound: gap-clustering understates the ≥12-copy array as
  8 loci, and minimap2's N=5 cap drops secondaries past the 6th placement.)

## Isoforms: two clean regimes

| regime | families | ambiguous (primary-MAPQ-0) mass | resolution corrects |
|---|---|---|---|
| **separate-chrom / well-separated** | RABL2, AK6, CCDC196, MAGEA_p2 | 0–12% | **~0** — minimap2's `de` already de-ties them (RABL2: 0/195 MAPQ-0, primary on the better-fit copy 195/195) |
| **co-located tandem** | MAGEA_p1/p3/p4 | **95–99%** | the coin-flip label on **474/495 spliced reads** |

In the co-located MAGEA pairs, 95–98% of family reads are spliced **and** primary-MAPQ-0 **and** exactly
de-tied (de-gap = 0.0000) — nearly every isoform is pinned to a coin-flip copy. That is the headroom resolution
would, in principle, correct.

## The bound — and why it is the whole story

**0 of the panel's 494 ambiguous reads lie in a K≥2 (PSV-bearing) family.** The two levers are **disjoint**:

- The ambiguous mass that **needs** resolving (488/494 = 98.8% MAGEA co-located) is **K=0**: per-read edit
  distance is *identical* against both copies (NM_A == NM_B, 1–3 mismatches = HiFi error floor; 0/311 reads
  differ between copies on MAGEA_p3). The copies are **sequence-identical over the transcribed exonic region** —
  the 67–72% genomic divergence is entirely intronic/intergenic. No PSV column can adjudicate; attribution is
  intrinsically arbitrary.
- The copies that **can** be resolved (chrY AS-gap median 259; RABL2; AK6; CCDC196) are already MAPQ>0 → already
  correctly placed by minimap2. No ambiguity remains to resolve.

> **Where resolution is possible it is unnecessary; where it is needed it is impossible.**
> Achievable-improvement ceiling on this panel's ambiguous mass = **0%** (resolvable 0/494; stuck 494/494).

## What this means (not a defeat — a characterization)

This is the identifiability theorem made concrete, and a precise statement of **when copy-resolution helps**:
**co-located AND exonically-divergent (K≥2)** families. The GGO panel contains no such family — co-located ⇒
exonically identical (K=0), divergent ⇒ already separated (MAPQ>0) — so the two conditions never co-occur here.

- The improvement is real and was demonstrated in the regime it targets: the **sim5x K-ladder** (synthetic
  co-located copies, PSV ladder K=0..8) shows K≥2 → 100% correct assignment. The win exists; it is gated by K.
- The value for the thesis is the **bound itself**: the achievable gain from resolution is a measurable function
  of RNA-level identifiability (K = exonic PSV count), predictable a priori per family. The same K governs
  detection, resolution, and isoform separation — the through-line to the identifiability theorem.
- The **joint-inference loop** (detect → resolve → recover copies/isoforms → refine loci → re-detect) is virtuous
  only in the co-located-AND-divergent regime; on GGO it stalls at the resolve step because the reclaimable mass
  is K=0. The one place it closes (the 2 chrY starved copies) is granularity- and cap-limited — a trickle, not a
  cascade.

**Bottom line.** Yes, in principle: family completeness gains 3 starved copies, and 474/495 co-located spliced
reads carry a coin-flip label resolution would fix. But the achievable gain on this real panel is ~0%, because
the ambiguity that needs resolving is RNA-level-unresolvable (K=0) and the resolvable copies are already
resolved. A real win requires a co-located-AND-divergent family — present on the sim ladder, absent from this
GGO panel.


---

## stringtie_vs_copy_headtohead

# StringTie vs. copy-aware de-novo at multi-copy families (GGO) — what the copy layer adds, and what it doesn't

**Date:** 2026-06-21 · **Data:** T2T gorilla (`GGO`) IsoSeq · **Branch:** `vg/flow-capacity-apportionment`
**Reproduce:** `bench/stringtie_vs_copy_headtohead.py` (+ `copy_assign.py real`)
**Status:** revised after a 3-lens adversarial review that correctly deflated an
earlier, overstated draft. Read the "What this is NOT" section.

**Question (the advisor's):** does the copy-aware pipeline do *better* than StringTie
on GGO *because it captures copy information*? Standard isoform recall (gffcompare
FSM) can't answer it — the annotation collapses paralogs, so a correct copy split
looks redundant.

**Answer (honest, narrow):** the genuine, unique contribution is **per-copy read
attribution** — assigning each read to *which paralog copy* it came from — a
capability StringTie structurally lacks. It is **not** transcript recovery:
**in 60/63 collapse cases (95 %) StringTie already models both copies' transcripts**
(distinct intron chains), it just labels them one `gene_id` and offers no per-copy assignment.
The set of genuine tandem families where this matters is **small and must exclude
domain-sharer "families"** (the ZNF mega-family DSFAM0 accounts for 28/63 collapses
and is a false family by our own detection validation). Net: a real *capability* gain,
**bounded** by identifiability and by family-definition validity — not a headline
copy-count.

---

## Setup

- **StringTie:** `genome_st.gtf` = `stringtie -L -p 4 GGO.bam` (v3.0.1, de-novo) → 68,166
  transcripts (19,416 gene_ids). StringTie *does* read secondary alignments (verified:
  removing them changes output), but it has **no per-copy concept** — paralog copies
  whose transcripts overlap get one `gene_id` and no read-to-copy assignment.
- **Copy-aware:** the de-novo family + copy-assignment pipeline (`denovo_families*`,
  `copy_assign.py`; `bench/denovo_family_pipeline.md`).
- **Unit:** co-located multi-copy family = ≥3 *distinct* same-chrom loci within 5 Mb
  (78 such families; collapses are counted per shared StringTie gene_id).

## The over-split correction (and a real pipeline bug it exposed)

Spot-checking the top hit — **"PRNP, 5 copies"** — exposed that all five spans are
14,600,00x–14,615,51x: **one locus** split into 5 pseudo-copies (the locus union-find
merges isoforms only on an *identical* junction, so near-identical-junction variants stay
separate and then look like a multi-copy family). This was not rare: **495 of the 1,190
de-tie "families" (≈42 %) were entirely one over-split locus** (verified: ASTN, SDHAF,
PRNP, novels — all isoform-fragments of a single gene), and ~21 % of co-located "copies"
were fragments.

**Fixed at source** in `denovo_family_split.py` — an *output-level* dedup (merge family
members whose spans reciprocally overlap ≥50 %), applied **after** family detection so the
POA detection/decomposition is untouched: **1,190 → 695 family-class families** (495
over-split single-locus families dropped), MAGEA 15→**11**, **APOBEC3 recovered**,
RFPL/RABL2/KRT preserved, copies 3,636→2,362. (A first attempt merging loci *before* POA
was reverted — it perturbed detection and dropped genuine copies, MAGEA 15→6, APOBEC3
broken.) **All numbers below are on the fixed list**, where the head-to-head's own
over-split guard is a 0 % no-op — confirming the fix subsumes it.

## What the copy layer actually adds

**1. Per-copy read attribution (the real, unique capability).**
`copy_assign.py real` over 25 co-located families / 30,709 reads:
- **Unique-mapper agreement = 28,704/28,726 = 99.9 %** — this is the *silver-standard*,
  but it is on the **confident unique-mappers** (minimap2's MAPQ>0 reads, ~97 % here). It
  validates that our per-copy labels agree with minimap2 where minimap2 is already confident.
- **Hard reads:** 95 % of all reads get a confident copy assignment; these are *not*
  independently validated (no orthogonal ground truth) — reported, not claimed as proven.
- **Copy-specific junctions** resolve **+167 reads** PSVs alone cannot — e.g. **DSFAM42**
  (5 copies, 95 % MAPQ-0): **10 % → 99 % resolvable**.
StringTie produces none of this — it has no read-to-copy assignment at all.

**2. Copies StringTie merges under one gene_id (attribution, not recovery).**

| | Count | Honest reading |
|---|---|---|
| Collapse instances (≥2 distinct loci share one StringTie gene_id) | **63** | |
| └ **StringTie emits distinct intron chains at both loci** | **60/63 (95 %)** | → it *models* the copies; only the *label* is merged. **Attribution, not recovery.** |
| └ in **DSFAM0** (178-copy/18-chrom ZNF mega-family = false domain-sharer family) | **28/63** | excluded as a non-genuine family |
| Collapses in **genuine** (non-DSFAM0) families | **35** (mostly disjoint-tandem) | many still other ZNF |
| └ in genuine families **copy-assignment validated** | **12** | the defensible core |
| └ in genuine **non-ZNF** families (PDPR, GSTM, TRIM, MYH, RFPL, novels) | **~16** (≈6 validated) | cleanest, free of domain-sharer doubt |

Named, defensible wins (StringTie merges the copies' transcripts under one label; we
attribute reads per copy): **RFPL** (incl. RFPL2, a family-rescued single-read copy),
**GSTM, TRIM, MYH, PDPR, PCDHB, MAGEA/MAGEB, KRT**.

**3. Family-aware rescue (soft, overlapping, mostly thin).**
94 co-located distinct loci are **family-rescue-flagged** (POA-confirmed against a family
member — a mechanism StringTie lacks), and 55 distinct loci get **no StringTie model** at
all. **These two sets overlap heavily (53 of 55) and must NOT be summed** with the
collapses or each other. Most are single/two-read assemblies (the rescue stage is 85
single-read + 25 two-read; a handful have high support). Honest framing: a thin,
real-but-soft tail, not a headline.

## What this is NOT (corrections to the first draft)

- **Not transcript recovery.** StringTie models 60/63 collapsed copies; the gain is
  per-copy labeling/attribution.
- **Not 40–60 "recovered" copies.** That conflated overlapping win-rows and a
  domain-sharer family. The defensible collapse core is **~12 validated** (≈6 non-ZNF).
- **The 99.9 % is unique-mapper consistency** (easy reads), not proof on the hard reads
  the method exists to resolve.
- **StringTie is not primary-only.** It reads secondaries; it simply has no copy model.

## Why standard recall can't see even this

gffcompare scores against a paralog-collapsed annotation, so per-copy attribution reads
as redundant. The (separate) read-coherence recall win (+1,735 FSM) is an *assembler*
lever, orthogonal to copy-awareness.

## Q: more multimappers in GGO.bam?

GGO.bam was aligned `minimap2 -ax splice:hq -uf` (default cap `-N 5 -p 0.8`). Re-aligning
array reads uncapped (`-N 50 -p 0.1`, `arrayfix/`) turned **737 records → 10,513 (9,757
secondaries) across ~413 loci** — the multimapper graph the copy-aware pipeline uses to
*enumerate and assign* array copies (the 5→11 array-core recovery). StringTie reads
secondaries but has no per-copy model, so more of them give it no copy resolution; the
copy-aware side is the beneficiary. Caveat: 737 reads cost 18.5 GB RAM (full-genome
index) — do it **targeted at family loci**.

## Verdict for the advisor

*"On GGO, StringTie usually models the paralog copies' transcripts but merges them under
one gene_id with no read-to-copy assignment. The copy-aware pipeline's genuine, unique
contribution is **per-copy attribution** — 99.9 % agreement with confident unique-mappers,
plus copy-specific-junction resolution (DSFAM42 10→99 %). The set of genuine tandem
families where this is demonstrable is small (~12 validated collapse loci, ~6 free of
domain-sharer doubt) and bounded by identifiability **and** by which 'families' are real.
This is a capability StringTie lacks, not a transcript-recovery count."*

## Reproduce

```bash
PY=/home/juanfra/miniforge3/bin/python
$PY bench/copy_assign.py real               # -> copy_assign_real.out (resolvability)
$PY bench/stringtie_vs_copy_headtohead.py   # -> the tables above (incl. labeling-vs-recovery, DSFAM0 split)
```


---

## quant_benchmark

# Quantification benchmark: the value of multi-mapper resolution

Two near-identical paralog copies (A,B) at abundance theta:1-theta, differing at K PSVs; HiFi reads (P_cov=0.6,
eps=0.005). Estimate per-copy abundance via naive 1/k, abundance-EM (soft, PSV likelihood), psv_hard (per-
molecule majority vote). Metric = |theta_hat - theta|, 25 reps. (quant_benchmark.py / .tsv / .png)

## Findings
1. **K=0 (no PSVs): ALL methods fail equally** (error = |theta-0.5|). The identifiability floor — no information,
   no resolution. The same wall as detection.
2. **The value is at SKEWED abundance, not symmetric.** At 50:50, naive 1/k is already exact (0 error). The
   payoff is at asymmetric expression — which is the COMMON paralog case (one copy dominant / a pseudogene):
   value (naive_err - psv_err) reaches **+0.19 at 70:30 and +0.40 at 90:10**. (This corrects the earlier
   "gap at similar abundances" framing — for point estimates it's the opposite.)
3. **Value scales with PSV density K — the long-read edge.** 90:10: value +0.24 (K=1) -> +0.40 (K>=8). Each
   molecule covering more PSVs => exact quantification; short reads (few PSVs/read) can't.
4. **Soft EM beats hard assignment at SPARSE PSVs.** K=1, 90:10: EM err 0.007 vs psv_hard 0.160; they
   converge by K>=4-8. ACTIONABLE: for quantification, copy-assignment should be SOFT/probabilistic (use the
   likelihood), not hard labels — hard assignment wastes thin partial evidence.
5. **At symmetric abundance the value is DIFFERENTIAL POWER, not point error.** naive always reports 50:50
   (zero variance but ZERO power to detect a real copy difference); only resolution (EM/PSV) recovers the
   actual ratio and can test copy-specific change. So resolution is essential for differential analysis even
   when abundances are equal.

## Takeaway
Resolving multi-mapper ambiguity adds nothing to ASSEMBLY (structure is in the primaries) but is decisive for
QUANTIFICATION of asymmetrically-expressed paralogs, scaling with PSV density (long-read advantage) and bounded
only at the K=0 identifiability floor. The right estimator is the SOFT likelihood (Canzar's frame), not hard
labels. This is the concrete payoff of the PSV machinery: per-copy abundance, not better assembly.

## Real-reads validation (sim5x, ground truth) — sim5x_quant_validate.py

Subsampled the sim5x HiFi reads to a SKEWED abundance (40:30:20:10:10 = [0.364,0.273,0.182,0.091,0.091],
truth in read names) and swept the PSV ladder K=0..8, running the actual copy_assign soft quant:

| K | soft EM err | naive(uniform) err | hard-count err | soft abundance |
|---|---|---|---|---|
| 0 | 0.473 | 0.473 | 1.27 | [.2,.2,.2,.2,.2] (uniform: no info) |
| 1 | 0.276 | 0.473 | 1.06 | partial resolution |
| 2 | **0.007** | 0.473 | 1.02 | [.37,.27,.18,.09,.09] ≈ truth |
| 4 | 0.007 | 0.473 | 1.02 | ≈ truth |
| 8 | 0.025 | 0.473 | 1.02 | ≈ truth |

CONFIRMS the simulation benchmark on REAL reads, and proves the EM RESOLVES (not just returns the prior):
- K=0 (no PSV): soft == naive == uniform (honest identifiability floor; the soft EM does not fabricate).
- K=1: soft already pulls toward the truth (0.276 < naive 0.473) -- uses the single PSV's PARTIAL evidence.
- K>=2 (PSVs clear the HiFi error floor, psv_acc=1.0): soft = 0.007, recovering the SKEW exactly -- genuine
  per-molecule resolution, not the 0.2-each prior. Naive stays flat at 0.473 (never resolves).
- Tracks the identifiability ladder (per-read psv_acc K0=0.2 -> K2=1.0). The quant error IS the assignment error.
- (hard read-count is ~1.0 throughout: each read has 5 alignments and the wrong-frame ones default to c0 --
  itself an illustration of why naive best-alignment counting fails on multimappers; the clean contrast is
  soft-resolves vs naive-uniform.)

