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
