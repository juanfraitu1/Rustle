# Minimizer-free copy discovery via intron-chain alignment (structural axis)

**Question (user):** the family identification uses minimizers, but we're not restricted to them — a
full alignment of intron chains should also work. Try it without minimizers.

**Answer:** intron-chain alignment is a viable, discriminative, **minimizer-free** family criterion.
It recovers the flagship cross-chromosome dup, rejects domain-sharers by construction, and at high
precision is **94% concordant with the sequence pipeline** — i.e. a strong *independent confirmation*
axis. It is not a standalone recall expander in this data (it is blind to retrocopies, which were a
large fraction of the sequence finds).

## Method (bench/extract_intron_chains.py + family_intronchain_discovery.py)
1. Per gene, the **intron chain** = ordered units `(exon_len, intron_after)`, multi-exon only
   (20,502 genes; 2,565 single-exon are the explicit retrocopy blind spot).
2. **Candidate generation (no sequence minimizers)** — two modes (a documented tradeoff, below).
3. **Full Needleman-Wunsch alignment of the two intron chains.** A unit matches iff exon lengths
   agree (±max(6bp,10%)) AND the flanking intron lengths agree (±max(40bp,20%)). Gaps model intron
   gain/loss / exon skip. Structural score = matched units / shorter chain; **gate: score≥0.6 AND
   ≥4 matched units** (RABL2 matches 6; few-exon and chance matches fall below 4).

**Coupling exon+intron is essential.** Exon lengths alone are NOT discriminative (internal exons
cluster ~120–150 bp) → matching them gave 173,283 chance "copies". Adding intron lengths (which span
4–5 orders of magnitude and are preserved between copies) collapsed that to a precise set.

## Recall — the criterion works
- **RABL2A (9 exons, NC_073235.2) ↔ RABL2B (10 exons, NC_086018.1): matched 6, score 0.67** —
  recovered despite the 9-vs-10-exon intron gain/loss (the NW gap absorbs it).
- Discriminative against domain-sharers: the universe "family" LOC101144552 (19 vs 4 vs 16 exons) is
  **correctly rejected** (score 0.06–0.25) — structure is naturally hostile to the FP modes (domain
  share, shared transposon) that needed extra filtering in the sequence pipeline.

## The headline: independent confirmation (structure ⟂ sequence)
At the high-precision gate (matched≥4, length-window): **360 cross-chromosome structural copies,
340 (94%) of which the POA-sequence pipeline also found.** Two fully independent axes — intron-chain
structure vs sequence alignment — agree on the same copies. A pair confirmed by BOTH is robust;
this is a provable, minimizer-free second line of evidence (advisor-aligned).

| axis overlap (cross-chrom) | count |
|---|---|
| confirmed by BOTH structure + sequence | 340 |
| structure-only (sequence missed) | 20 |
| sequence-only (structure blind) | 7,964 |

The 7,964 sequence-only are dominated by **retrocopies / single-exon** copies (no intron chain) plus
copies whose structure diverged past tolerance — the structural axis cannot see these by construction.

## Candidate-generation tradeoff (a finding; both modes kept, `--cand`)
| mode | candidate pairs | cross-chrom confirmed | recall (univ 5) | structure-only | runtime |
|---|---|---|---|---|---|
| length-window, matched≥4 (default) | 5.98M | 360 | 1/5 | 20 | ~5 min |
| length-window, matched≥3 | 5.98M | 993 | 1/5 | 527 | ~5 min |
| 2-intron-shingle, matched≥3 | 14.5M | 38,588 | 3/5 | 37,138 | ~50 min |

- A **single intron length is not discriminative** (common sizes shared by thousands of genes → an
  intron-length index over-collapses). **2-intron shingles** (consecutive intron-size pairs) are
  specific enough to index — they recover more (partial / large-indel copies the length-window
  excludes, e.g. LOC129529434's 5-exon↔10-exon contained match) but are broad → need a tighter NW
  gate and ~10× the compute. The length-window's precision partly comes from its length constraint,
  which is also why it misses partial copies — an inherent tradeoff.

## Honest limitations
- **Retrocopies / single-exon invisible** (no intron chain) — the structural axis's hard blind spot;
  the sequence axis covers them. The two are complementary, neither subsumes the other.
- **Few-exon genes** (≤3 exons) are structurally under-determined (can't reach matched≥4); universe
  2-exon families (LOC101127159, LOC129529611) are missed for this reason.
- **Universe recall 1/5** at high precision, but the misses decompose as: 2 few-exon (under-determined),
  1 structurally-different (correctly rejected = precision win), 1 partial-copy (candidate-gen excluded);
  RABL2 (the one structurally-tractable non-partial case) is recovered.
- **Input = RefSeq gene reps**, one representative isoform per gene.

## Verdict
Minimizer-free intron-chain alignment is a **viable, discriminative, independent structural axis**,
best used as a **confirmation layer** alongside sequence (94% concordant; a both-axes copy is robust),
not a standalone recall expander here. Strengths the sequence axis lacks: domain-sharer resistance and
sequence-divergence robustness. Blind spot the sequence axis covers: retrocopies. **Best system = both
axes.** Next: the user's planned use is exactly this — structure as a second, minimizer-free signal.

## Reproduce
- `python3 bench/extract_intron_chains.py` (gene chains → /tmp/gene_chains.tsv)
- `python3 bench/family_intronchain_discovery.py --cand length` (default; `--cand shingle` for recall mode)
