# DNA-derived PSV identifiability catalog (Phase 1) — design

*A reference-only (T2T DNA) catalog of exonic paralogous sequence variants (PSVs) per multi-copy-family copy
pair, yielding a genome-wide table of what RNA can vs. provably cannot resolve, validated against the RNA
census. Phase 2 (the per-read signature decoder) is designed for but deferred.*

## Goal

For every multi-copy family, compute from the T2T reference alone — per copy-pair — the number of **exonic**
PSVs **K** (positions a spliced read could use to tell the copies apart). Produce a genome-wide identifiability
table (resolvable K≥1 vs. genuine-K=0), and cross-check it against the 87%/13% RNA census so that **DNA truth
independently confirms (or refines) the RNA identifiability split**.

Why this matters: today PSVs are found *from the RNA* (coverage-biased, confoundable by error / A-to-I editing —
exactly longcallR's problem). The DNA reference has every copy's full sequence and *which allele belongs to which
copy*, so identifiability becomes a reference computation, and per-read assignment (Phase 2) becomes a supervised
nearest-signature decode rather than the NP-hard unsupervised phasing that SDA / longcallR approximate.

## Non-goals (this phase)

- No per-read decoder (Phase 2): the per-family allele matrices are persisted for it, but read→copy assignment
  and per-copy quantification are out of scope here.
- No re-alignment of reads, no Rust/pipeline integration. Pure Python bench tooling + minimap2.
- No claim about copies absent from the T2T reference (assume the T2T assembly contains the copies; flag the
  rare unaligned copy rather than infer it).

## Substrate (all reference-only)

- Families: `bench/denovo_families.tsv` (1130 families; member ids `DN_<chrom>_<start>_<n>` — same input as the
  census, so the co-located pairs line up for the cross-check).
- Reference: `/home/juanfra/winloci_scratch/GGO.fasta` (+ `.fai`).
- Annotation: `/home/juanfra/winloci_scratch/GGO_genomic.gff` (exon intervals).
- Tools: `minimap2` (2.30), `samtools` — present.

## Method (per family)

1. **Copy intervals.** For each member locus, the genomic interval = the annotated gene span (GFF gene
   overlapping the locus start) if present; else `[max(0,start−2000), start+W]` with `W=40_000` default.
   Alignment trims to the homologous extent, so the window only needs to be generous, not exact.
2. **Reference copy + projection.** Pick the longest interval as the family reference copy `ref0`. Align every
   other copy's sequence to `ref0` with `minimap2 -cx asm20 --cs=long ref0.fa copy_i.fa` (the asm20 mode the
   K=0 frontier ref-vs-ref attack used). Parse the `cs` tag to walk `ref0`/`copy_i` coordinates and record, at
   each `ref0` position, `copy_i`'s base — building a per-column **allele matrix** `allele[copy][ref0_pos]`.
   One alignment per copy (≈1130 minimap2 calls total, not per-pair ≈6362). A copy that fails to align to
   `ref0` (rare, divergent) is recorded as `unaligned` and excluded from its pairs with a flag.
3. **PSV columns.** `ref0` positions where ≥2 copies carry different bases (substitutions). Copy-private
   insertions/deletions and exon gains/losses are not PSV columns but are summed as `private_exon_bp`
   (secondary distinguishing signal, reported, not folded into K).
4. **Exonic restriction.** A PSV column is exonic iff its `ref0` position falls inside a GFF exon overlapping
   the `ref0` window. (Assumption: homologous exon structure, so `ref0`-exonic ≈ copy-exonic; Phase 2's read
   decoder uses per-copy exons precisely.) If `ref0`'s locus has **no overlapping GFF exon**, the pair is
   flagged `unannotated_exons` and reported with `psv_exonic = NA` — excluded from the resolvable/K=0 verdict
   rather than fabricating exons (no dependency on a de-novo-transcript exon file). The annotated-exon coverage
   is reported in the summary so the NA fraction is explicit.
5. **Per-pair K.** For each copy-pair `(i,j)`, `K = ` number of exonic PSV columns where `allele[i]` and
   `allele[j]` both exist and differ (Hamming over exonic PSV columns). `aln_identity(i,j)` from the projected
   matches/length. Verdict: `unannotated` if exonic restriction is NA (step 4); else `genuine_k0` if `K==0`;
   else `resolvable` (K = the identifiability margin).

## Output

- `bench/dna_psv_catalog.tsv` — one row per copy-pair:
  `family, copyA, copyB, chromA, chromB, co_located, aln_identity, psv_total, psv_exonic(=K), private_exon_bp, verdict`.
- `bench/dna_psv_catalog_summary.md` — K distribution; % resolvable vs % genuine-K=0 (co-located subset and
  all-pairs); top genuine-K=0 families; and the cross-check (below).
- Per-family allele matrices persisted to `/home/juanfra/winloci_scratch/dna_catalog/<family>.json` (scratch;
  the Phase-2 decoder's input — each copy's signature = its row of the matrix at exonic PSV columns).

## Validation — the cross-check (the headline)

Restrict to co-located copy pairs that the RNA census classified (reuse `copy_resolution_census.classify_pair`
to get `frac_same`); RNA-K0 := `frac_same ≥ 0.95`. Build the 2×2 confusion of **DNA-K=0 vs RNA-K0** and report
concordance. Interpretation of the cells:
- concordant DNA-K0 ∧ RNA-K0, and DNA-K≥1 ∧ RNA-resolvable → DNA independently confirms the RNA split.
- DNA-K≥1 ∧ RNA-K0 (RNA-tied despite a reference PSV) → the PSV sits in an unexpressed / low-coverage exon →
  "expressible but not expressed" (an honest refinement: reference identifiability ≥ read identifiability).
- DNA-K0 ∧ RNA-resolvable (RNA NM differs with 0 exonic substitution PSVs) → pseudo-K=0 via indels / splice
  shifts → caught by `private_exon_bp` (the indel signal), reconciling the two.

The discordance map is itself a result, not noise.

## Testing / verification (deterministic check)

A `check()` runnable as `python3 bench/dna_psv_catalog.py`:
- **Panel anchors (exon-identical ⇒ K=0):** MAGEA `MAGd0` (161251228~161458538), `MAGd2`
  (164381222~164442447), `MAGd3` (164397061~164426194) must have **exonic K = 0** (exons identical — matches
  the per-read finding that pair0's divergence is purely intronic, and pair2/pair3 are fully identical).
- **Resolvable anchors (K≥1):** at least one of the cross-chrom copies (`AK6`~`LOC115934278`,
  `CCDC196`~`LOC129526440`, `RABL2A`~`RABL2B`) must have **exonic K ≥ 1**.
- **Cross-check concordance:** DNA-K=0 vs RNA-K0 agree on ≥ 80% of the classified co-located pairs.
- **Internal consistency:** `psv_exonic ≤ psv_total`; every classified-census co-located pair appears in the
  catalog.

## Architecture & data flow

```
denovo_families.tsv ─┐
GGO.fasta ───────────┼─> per family: extract copies → minimap2 asm20 --cs vs ref0
GGO_genomic.gff ─────┘      → allele matrix → PSV cols → exonic restrict → per-pair K
                                   │
                                   ├─> dna_psv_catalog.tsv (+ summary)         [Phase-1 deliverable]
                                   ├─> cross-check vs copy_resolution_census   [headline validation]
                                   └─> per-family allele matrices (scratch)    [Phase-2 decoder input]
```

## Risks / open points

- **Runtime:** one minimap2 call per family (~1130) over gene-sized sequences — bounded (minutes). Mega-family
  DNFAM0 (728 members) is the outlier; cap copies-per-family processed (e.g. 200) with a logged note, since
  its members are mostly genome-scattered (few true co-located pairs) — do not silently drop.
- **Ref-projection drops copy-private insertions** → tracked as `private_exon_bp`, not lost.
- **Exon-source proxy** (ref0-exonic for both copies) is an approximation; acceptable for an identifiability
  *count* and flagged. The cross-check's discordant cells surface any case where it matters.
- **Cross-chrom pairs** (AK6/RABL2) align the same way, flagged `co_located=False`; they are mostly resolvable
  and not the identifiability bottleneck, but completeness is cheap.
