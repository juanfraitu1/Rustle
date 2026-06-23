# More signal in the BAM: NM (and rl) beyond de — retained introns in the copy model

Two questions about whether the definition extracts the BAM fully.

## 1. Retained introns ARE summed into the copy model — by design, with one subtlety

`S(v)` = union of every read's `get_blocks()` (the aligned M/=/X segments). A normally-spliced
read contributes its **exons** (the intron is an `N` gap, excluded). An **intron-retention** read
has no `N` over the retained intron — it is aligned there — so that intron *is* a block and is
summed. This is the intended "full gene minus introns *except retained ones expressed in the RNA*."

**Subtlety (a real refinement).** The union takes *all* reads' blocks with no support threshold, so
a **single** intron-read-through read pads `S(v)` with the whole intron. Because `~B` scores
**reciprocal coverage**, padding `S(i)` with intronic sequence the properly-spliced homolog `S(j)`
lacks **lowers** `cov_i` and can push a real pair toward the `τ` edge — a plausible cause of the
validated families that sit at reciprocal coverage ~0.31. **Fix:** gate block inclusion by read
support (include a base only if ≥k reads align through it), so single-read read-through does not
enter `S(v)`; this should *raise* real families off the `τ` edge. (Not yet applied.)

## 2. de is not the only signal — NM discriminates bridges that de cannot

The new BAM (`-N 50 -p 0.1 --eqx`) carries `NM` (raw edit distance), `AS`/`ms` (DP scores),
`s1`/`cm` (chaining), and `rl` (length of query in repetitive seeds) — only `de` is used. Per
chromosome, classify `~R` cross-mappings as `~B`-copy vs `~B`-bridge and compare the reads:

**NC_073240.2 — 69 copy / 59 bridge edges, 4,560 read-placements**

| signal | real-copy (median) | bridge (median) | separation |
|---|---|---|---|
| `de` (gap-compressed divergence) | 0.0094 | 0.0083 | 0.9× — **none** |
| `NM`-rate (edit distance / aligned bp) | 0.0013 | 0.0035 | **2.6× — discriminates** |
| `rl`-fraction (repeat-seed bp / aligned bp) | 0.0041 | 0.0066 | 1.6× — weak |

**`de` does not separate bridges from copies at the read level (0.9×)** — confirming that `~B`
(copy-level homology) is doing the bridge-pruning, not the read-level tie test. But **`NM`-rate
separates them 2.6×**: a bridge read is forced onto a non-homologous locus with large indels, so
its *raw* edit distance is high while its *gap-compressed* `de` stays low. `de` compresses away
exactly the indel burden that flags a bridge; `NM` keeps it.

### Implication
- `de` is the correct **tie** criterion (it measures sequence confusability, which is genuinely
  similar for copy and bridge cross-mappings — the bake-off already preferred it over `AS`, which
  is biased on intronless retrocopies). `DE_max` as an absolute ceiling is weakly load-bearing.
- **`NM` adds an orthogonal, read-level bridge signal** that `de` misses. It could (a) tighten
  `~R` by down-weighting high-`NM` (forced-indel) placements before they form edges, or (b) serve
  as a cheap read-level pre-filter complementing `~B`'s copy-level pruning. `~B` remains the
  principled pruner (it directly tests "are these copies"); `NM` is a corroborating proxy — which
  is consistent with the conflict-criterion bake-off finding that NM *corroborates* de.
- `rl` (repeat-mediated-bridge flag) is directionally right but weak here; `AS`/`ms` are
  length/structure-biased. `NM` is the highest-value currently-unused tag.

Scope: one chromosome; `~B`-class via reciprocal coverage of exon-union models (`COV_MIN=0.30`).
Reproduce: `python bench/family_def_bam_signals.py NC_073240.2`.
