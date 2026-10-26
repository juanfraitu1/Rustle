# Pre-registration: cross-copy dominance prune of pool-added transcripts (FLAIR's transcriptome competition, ties abstain)

**Written 2026-10-26, before any prune run.** Branch `assembler-stringtie-gaps`. Spec: `FAMILY_LAYER_SPEC.md`. A default
flip is the user's call.

## Seen before this file (dev, human A119b genome-wide; never the verdict)
Stratified bakeoff (Amendment 1 strata): multi-copy precision 14.9 (primary-only run) vs 12.6 (pooled run, shipped
`RUSTLE_GTF_SECONDARY` 0.98); recall 17.7 vs 20.7; single-copy identical (19.8 / 22.1). Pool-added transcripts (chain
absent from the primary-only run): 12,184 = multi 4,853 (291 `=`), single 226 (0 `=`), no reference gene 7,105.

## Rule (annotation-free; `bench/stgaps/dominance_prune.py`)
1. Candidates = pool-added multi-exon transcripts (`pool_added.py`: pooled-run chain absent from the primary-only run).
2. Support of a candidate c = molecules with a record (primary, or secondary with AS >= 0.98 × the molecule's best AS,
   the pool's own admission) at c's contig/strand whose exact intron chain equals c's.
3. Every emitted multi-exon transcript of the pooled run → spliced sequence; support molecules (full sequences from
   their primary record) aligned with `minimap2 -x map-hifi -c -N 100 -p 0.8`. AS(m, t) = the `AS:i` of m's best hit on t.
4. m is **won elsewhere** for c when max_t AS(m, t) over transcripts t NOT overlapping c's span is **strictly greater**
   than AS(m, c) (no hit on c = −∞). Ties (equal AS) are NOT won elsewhere.
5. **Drop c iff every support molecule is won elsewhere** (and c has ≥ 1 support molecule found). Nothing else changes.

## Substrates
Dev: human A119b genome-wide (`rustle_figures/runs/human_A119b/human_A119b{,.primary}.gtf`, `A119b.t2t.bam`).
**Held-out: gorilla OR6737 genome-wide** (`runs/gorilla_OR6737/gorilla_OR6737{,.primary}.gtf`, `GGO_mm.bam`), strata by
the same Amendment-1 definition on GGO RefSeq (`GCF_029281585.2`) mapped with a gorilla splice index (built for this).

## Clauses (each substrate; strata metrics from `stratify.py`; PRIM / POOL / PRUNE = primary-only / pooled / pruned)
- **D1 precision:** multi-copy precision PRUNE >= PRIM's.
- **D2 keep the gain:** multi-copy `=` reference transcripts recalled: PRUNE − PRIM >= 0.8 × (POOL − PRIM).
- **D3 single-copy untouched:** single-copy recall and precision |PRUNE − POOL| <= 0.1 point each.
**ADOPT as opt-in** = D1-D3 on held-out gorilla (dev reported, not judged). Otherwise a register row.
Descriptive: pool-added no-reference-gene transcripts dropped / kept (they may be real unannotated copies: an O3 caution,
not counted as gain); class n/o/x shares per stratum.
