# Copy-assignment recompute on the COMPLETE multimapping BAM (GGO_mm.bam)

**Why a recompute.** The original `GGO.bam` was the multimapping UNDERCOUNT (minimap2 default `-N` cap →
few secondaries). The correct substrate is **`GGO_mm.bam`** (`minimap2 -ax splice:hq -uf --eqx -Y -N 50
-p 0.1 --secondary=yes`) — chr19 alone has **63% secondary alignments** vs nearly none before. Copy
assignment lives on those secondaries (they are the candidate copies), so it was recomputed end to end.

## Three corrections that made the recompute trustworthy

The recompute surfaced — and fixed — three problems that had been inflating the family count and the
copy-assignment numbers:

1. **PSV discovery: poasta → minimap2 asm20.** `discover_psvs` aligned every copy to copy[0] with
   poasta (exact Dijkstra DP) on the full ~10 kb spliced sequences — the dominant cost (a single rich
   family took **>5 min**). Replaced with minimap2 asm20 (`-c --eqx`), reconstructing the same 2-row
   gapped MSA from the CIGAR → **same PSVs, ~100× faster** (the genome-wide sweep went from *5 min per
   region* to **40 regions / 75 s**). Gated `RUSTLE_PSV_MINIMAP2=1`; poasta still the default.
   *Crucially, the minimap2/poasta discrepancy on real families is what exposed the over-merges below.*

2. **Same-strand only (`colocated_families` partitions by chrom AND strand).** The family path did not
   enforce strand, so **antisense pairs** (a `+` gene and a `−` gene that share sequence via an
   inverted-repeat) were merged as "copies." A read carries its own strand, so there is no copy-level
   ambiguity to resolve — and force-aligning a `+` spliced sequence to a `−` one manufactures spurious
   PSV columns. Verified on the **genome splice motifs** (authoritative; the BAM's minimap2 `-uf` flag
   makes the `ts:A` tag uninformative): the DSFAM817 "flagship" copies are **GT-AG (+) vs CT-AC (−)** —
   a genuine antisense pair, not a paralog family.

3. **Disjoint loci (`prune_same_locus`).** Copies of a multi-copy family must be at DISTINCT loci. Two
   same-strand transcripts are the SAME locus iff they **share a junction** (an intron coordinate pair is
   an exact splice-site fingerprint → isoforms of one gene) or one is a **structureless (≤1-exon) span
   contained** in the other (an unspliced read-through). Tandem paralogs have disjoint junction sets, so
   they are preserved. This removed the CAFAM0 case (a 12-exon gene grouped with its own 1-exon
   read-through). TDD'd (drop-readthrough / dedup-shared-junction / preserve-disjoint-tandems).

**Validation that the fixes are not over-aggressive** (the key check): real same-strand disjoint families
SURVIVE intact — a 2-copy `−/−` pair (silver 100/100), a 5-copy tandem recovered to 8 copies (silver
7/7) — while only the over-merges (antisense pairs, spliced/unspliced) collapse.

## Genome-wide result (281 family regions → clean families)

| | OLD (over-merged, poasta force-aligned) | **NEW (clean: minimap2 + strand + disjoint-loci)** |
|---|---|---|
| "families" | 281 regions, many antisense/spliced over-merges | **58 clean families** — a **LOWER BOUND**, see below |
| copies | inflated by spurious force-alignment | **144** (45×2-copy, 6×3, 3×4, 3×5, 1×9) |
| collapsed copies recovered | — | **978** (`--recover-copies`, from the complete multimapping) |

> ⚠ **"58" was a lower bound** (re-scan of the old 281-region `core_recip≥0.13` catalog windows). The
> **true genome-wide count is now computed — see below.**

## ⭐ The genome-wide de-tie conflict-graph catalog (O1 / L1 — first run at scale)

`gw_family_catalog` (`src/bin/gw_family_catalog.rs` → `denovo_pipeline::detect_conflict_catalog_genome_wide`)
runs the **principled, threshold-free family definition** genome-wide: a family is a connected component
of the **read-conflict graph** (loci among which reads are genuinely confused — de-tied multimappers,
NO `core_recip≥0.13` similarity bar), with the same-strand + disjoint-loci fixes applied. This is the O1
criterion the prior audit flagged as "never run genome-wide" (L1).

Genome-wide on `GGO_mm.bam`: **164,795 skeletons → 15,769 reps → 245 conflict components → 82 clean
families / 207 copies** (no OOM; per-chrom, peak <7 GB).

| catalog | criterion | result |
|---|---|---|
| OLD `denovo_families.tsv` | `core_recip≥0.13` (arbitrary threshold) | 281 "families"; DNFAM0 = **728 members chr1→chrY** (over-merge) |
| lower bound (re-scan) | per-region conflict graph on old windows | 58 |
| **NEW genome-wide** | **de-tie conflict graph, threshold-free** | **82 families** (61×2-copy, 14×3, 2×4, 1× each 5/6/7/8/9), 18 contigs, balanced strand (43+/39−) |

**82 > 58** — the genome-wide conflict graph finds *more* than the 58 lower-bound (which was a re-scan of
the old catalog windows); this is *consistent with* 58 being a lower bound (it does NOT prove strict
set-containment 82 ⊇ 58 — different scope/criterion). The catalog is now the principled Canzar-aligned
definition (connected components of read-conflict, no similarity threshold) instead of the over-merging
`0.13` bar. Real paralog arrays surface: a **9-copy** tandem (NC_073247.2), an **8-copy** (NC_073228.2),
a **7-copy** (NC_086018.1), plus 6/5/4-copy families. Output: `gw_conflict_catalog.{families,copies}.tsv`.

> ⚠ **SCOPE — the catalog is co-located same-chrom tandem arrays ONLY, by construction (not a finding).**
> The "0 mixed-strand / 82/82 single-chrom" properties are *structural guarantees* of the method (per-chrom
> edge building + `colocated_families` partitioning by `(chrom, strand)` + `min_copies≥2`), **not** evidence
> that no cross-chrom families exist. The method is triple-blind to **cross-chromosome dispersed paralog
> families** — DAZ/DAZL, RABL2A/B, dispersed retrocopies, cross-chrom segdup families are biologically real
> but **out of scope by design** (a memory-bounded per-chrom choice). So this is a *genome-wide catalog of
> co-located same-chrom tandem arrays*, not the complete multi-copy-gene-family catalog. Also excluded
> (intentional): <2-copy families, conflicts with <3 de-tied reads, antisense pairs.

**Per-read assignment — on the 82-family CONFLICT-GRAPH catalog** (the catalog and the assignment are now
the SAME substrate; 187,037 multimapping reads, canonical engine, `τ=6.9`, minimap2 PSV + strand +
disjoint-loci fixes):
- **assigned: 87,791 (46.9%)** — resolved to a copy, no 1/k
- ambiguous: 11,937 (6.4%) — decisive feature but margin < τ → **abstain**
- tied: 87,309 (46.7%) — no distinguishing feature = the identifiability floor → abstain
- **Of the decisive reads (those carrying ≥1 PSV/junction), 88.0% were assigned.**
- **collapsed copies recovered: 1,406** (`--recover-copies`).

**Silver standard** (unique-mapper agreement; a *circular* proxy): **31,443 / 32,081 = 98.0%**. Per-copy
abundance median 0.49. (78/82 regions; the 4 densest tandem arrays timed out at 150 s and are logged.)

> Note: the higher TIED fraction here (46.7%) vs the earlier **58-substrate** run (23.1% tied, 72.3%
> assigned) is honest — the conflict-graph catalog is *enriched for genuine multi-copy tandem arrays*
> (more identifiability-floor reads with no distinguishing feature), so the raw assigned-% is lower, but
> *of the decisive reads it is still 88%*, and the engine still abstains rather than guess (no 1/k). The
> earlier 58-substrate numbers were a narrower re-scan of the old catalog windows; the 82-catalog numbers
> above are the principled, same-substrate result.

### Flagship clean families (the honest replacements for the retired over-merged ones)
- **NC_073236.2 ~46.6 Mb — 4-copy, silver 1235/1235 (100%)**, 1,181 junction-resolved reads, 55 PSV cols.
- **NC_073224.2 ~129.2 Mb — 4-copy**, 9,193 reads, 62 PSV cols, 3,072 junction-resolved, silver 279/290.
- **NC_073247.2 ~59.7 Mb — 9-copy tandem** (the largest), 2,357 reads, silver 57/57.
- **NC_073228.2 / NC_073239.2 ~122–145 Mb — 5-copy**, 83–109 PSV cols, silver 418/439 & 447/448.

## ⭐⭐ Cross-chromosome paralog families + exon-sum (FLNC) validation (L3)

The same-chrom catalog above is **blind to dispersed paralogs by construction** (RABL2A/B, DAZ-class,
retrocopies). `gw_family_catalog --cross-chrom` removes the per-chrom restriction: a read confused between
copies on *different* chromosomes forms a conflict edge, so families are no longer one-chromosome. To avoid
the transitive repeat-driven mega-merge (a first run produced a 958-copy/25-chrom blob), the global conflict
graph is decomposed with Louvain communities + a Web/density guard (`family_split::decompose_families`).

Genome-wide on `GGO_mm.bam`: **164,795 skeletons → 15,769 reps → 277 communities (10 Web-dropped) → 266
conflict families, 155 cross-chromosome.**

### The exon-sum (FLNC) validation — annotation-free, de-circularising

Mapping copies to RefSeq gene names (the `gene_at` route) injects an annotation artifact: a copy's span can
overlap a *neighbouring* gene and get mislabelled (a GST paralog tagged "SEC22B"). The clean check aligns
each copy's **exon-sum** — its spliced consensus sequence (`DenovoTranscript.seq`, built from the FLNC
reads), now emitted as `<out>.copies.fa`. A family's copies are all-vs-all aligned (minimap2 asm20). This is
**independent of both** the conflict-graph family definition (read-conflict fires on a shared *fragment*;
full-length exon-sum alignment with a coverage floor cannot) **and** of annotation. `bench/validate_exon_sum.py`.

> **purity(family) = (largest mutually-aligning component) / (#copies)** — 1.0 = a pure paralog family.

| | same-chrom | cross-chrom |
|---|---|---|
| median exon-sum purity | **1.00** (51% perfectly pure) | **0.50** (only 15.5% pure) |

**The raw cross-chrom conflict graph OVER-MERGES — and it is a genuine Alu-SINE artifact, not an alignment
miss.** Adversarially verified (4-agent review): the worst families merge *unrelated* genes
(GWFAM43 = ZNF827/PDLIM1/BRCA1-class; GWFAM41 = BRCA1/HSPA4L/CLIP1) bridged by a shared **150–300 bp Alu**
covering only **4–7%** of each transcript — aligned to the AluY consensus (mapq 56–60, 83–90% identity,
canonical Alu signature, high dinucleotide entropy ruling out low-complexity). The over-merge is **robust to
the alignment threshold** (cross-chrom purity stays ~0.50 / ~15% pure across identity 0.70→0.90), so asm20 is
*not* wrongly splitting real divergent paralogs — confirmed by the negative control: GWFAM6's one genuine
divergent core (ident 0.92–0.98, full coverage) is **retained** while its Alu halo is pruned.

### The refined catalog — read-confused AND mutually homologous AND multi-locus

The principled fix makes family membership a **two-criterion combinatorial object**: copies must be
(i) in one **exon-sum-homology** component *and* (ii) at **≥2 spatially-DISTINCT loci** (distinct paralogs
occupy disjoint genomic spans → two copies are the same locus iff their spans overlap on the same
chrom+strand; threshold-free). This removes both the Alu over-merges and the same-locus isoform/fragment
contamination (a gene + its own nested low-read fragment, which exon-sum homology alone is blind to).

**266 raw conflict families → 157 refined families, 54 cross-chromosome.** Size distribution:
111×2-copy, 24×3, 9×4, 5×5, 2×6, 2×7, 3×8, 1×9.

- **⭐ RABL2 recovered (GWFAM50/REF9): a full 5-clique across 5 chromosomes** — RABL2A (NC_073235.2),
  RABL2B (NC_086018.1), + 3 retrocopies; all 10 exon-sum pairs align; Compara + mmseqs (fident 0.875–0.995)
  confirm. This is the dispersed-paralog family the same-chrom catalog was structurally blind to.
- Large refined arrays, all annotation+protein-confirmed: **GAGE** cancer-testis (9-copy tandem),
  **GSTM** mu-class (8-copy), two **cross-chrom segdup arrays spanning 7 chromosomes** (8-copy, ident
  0.98–0.997), **JMJD7** (7-copy), **NBPF**/Olduvai (6-copy).

### DNA-level corroboration (genomic-span homology) — pushing past the annotation floor

The 70.9% annotation confirmation (below) is a *floor* set by what gorilla annotation covers. To reach past
it we checked **DNA-level homology**: for each refined family, extract every copy's **genomic span (introns
included)** from the assembly and all-vs-all align (minimap2 asm20, identity≥0.90 = segdup stringency,
coverage-of-shorter≥0.50). **DNA-confirmed: 140/157 = 89.2%** — and **statistically indistinguishable
between same-chrom (89.3%) and cross-chrom (88.9%)** (Fisher p≈1.0; cross-chrom n=54), so the cross-chrom
capability adds **no** extra false families. **RABL2 (GWFAM17) = a genuine 5-clique** (all 10 genomic pairs
align directly, identity 0.98–0.99). Adversarially verified (sensitivity 85.4%→91.7% across id 0.85–0.95;
confirmed alignments are 98–99% over multi-kb extents, not fragment artifacts; of the 17 "silent" families
15 are ≥90%-identity homologs that miss only on *coverage*, so ~98% (154/157) show real homology and only
~3 are genuinely suspect).

> ⚠ **This is DNA-level CORROBORATION, not a fully-orthogonal precision, and not a real SEDEF/BISER segdup
> map.** It is *partly circular* — the genomic span contains the same exons the catalog was built on; the
> independent content (introns/flanks) dominates only for intron-rich families (RABL2, the 4-exon/19kb
> cross-chrom arrays) and is minor for compact many-exon families. Report **89.2% as a lower bound on
> segdup-grade homology.** The genuinely independent check — a genome-wide SEDEF/BISER segdup map — is
> **blocked**: BISER's precompiled aligner segfaults under WSL2 (its unrefined putative SDs are Alu-repeat
> noise that link gene loci to repeats, not to their paralogs), so the gold-standard map remains TODO
> (build SEDEF from source, or run BISER off-WSL2).
> - The 17 genomic-silent families are **partial/structurally-divergent HIGH-identity paralogs failing on
>   coverage, NOT retrocopies** (0/17 intronless; mean n_exon 9.26 ≥ confirmed 8.29 — the retrocopy
>   hypothesis is rejected).
> - **GWFAM99 is a confirmed FALSE family** the DNA check caught: its two "copies" overlap in genome
>   coordinates on OPPOSITE strands (666 reads `+` vs 3 reads `−`) — a sense/antisense mis-split that the
>   same-strand-only distinct-locus rule lets through (the antisense-overlap edge case; a minority-strand
>   collapse would remove it).

### Independent (orthogonal) confirmation of the refined catalog

Of 141 RefSeq-mappable refined families, **100 (70.9%) are confirmed by orthogonal evidence**
(gene-name root / Compara / mmseqs protein homology fident≥0.40 & qcov,tcov≥0.50) — **flat across same-chrom
(71.4%) and cross-chrom (69.8%): cross-chrom families are NOT enriched for false merges.** The 41 unconfirmed
are non-coding (lncRNA/snRNA, absent from the protein DB) or protein-coding paralogs sitting just under the
0.50 coverage gate; re-checked at the nucleotide level, **all 41 are themselves homologous — zero genuine
false merges in the sample.**

## Robustness fixes (post-review) — corrected catalog = 152 families

Three fixes from the DNA review, all in-engine + unit-tested (644 tests):

1. **Antisense edge case** (`distinct_locus_reps`). The DNA check caught GWFAM99 — two "copies" overlapping
   in genome coordinates on OPPOSITE strands (666 reads `+` vs 3 `−`), a sense/antisense mis-split. Fix:
   opposite-strand overlapping copies collapse when one is a clear read-minority (a strand artifact),
   keeping the dominant copy; a BALANCED overlapping antisense pair stays distinct. Verified: the corrected
   catalog has **0 intra-family coordinate overlaps**.

2. **Divergent introns** (`RefineParams.include_introns`, `--refine-introns`). The exon-sum is exons-only,
   so two copies with identical exons but divergent introns look identical (the K=0 identifiability
   frontier). The flag refines on the **genomic span (introns included)** so that divergence separates them
   (and adds the intron/flank signal). Tradeoff: introns diverge faster, so this is STRICTER and drops
   older paralogs whose introns diverged — hence a knob, not the default. (`validate_genomic_dna.py` is the
   intron-inclusive view at the validation layer.)

3. **Length protection from random reads** (`pass1_skeletons_robust`, `DenovoConfig.min_terminal_support=2`).
   The transcript extent was `min(start)/max(end)` over all reads, so a SINGLE runaway read (chimeric /
   intra-primed / mis-clipped terminal exon) set the boundary and inflated the length. Now the boundary
   must be reached by ≥2 reads (the k-th most-extreme position), trimming lone outliers. FLNC full-length
   reads reach the true ends in bulk, so only runaways are trimmed. Legacy 2-arg `pass1_skeletons` (k=1) is
   unchanged → no existing test/pipeline shifts.

**Corrected catalog** (`gw_xchrom_refined2`, all three fixes, default exon-sum refine): **152 refined
families, 50 cross-chrom** (was 157/54). RABL2 intact (5 copies / 5 chromosomes). **DNA-confirmed 137/152 =
90.1%** (same-chrom 89.2%, cross-chrom 92.0%) — slightly *up* from 89.2% because removing the antisense
false family and the length-inflated marginals made the catalog cleaner.

## Divergent-paralog tiers

asm20 cannot **seed** a pair below ~80% identity, so lowering its identity gate is a no-op (the catalog is
unchanged from id≥0.80 down to id≥0.70). The ~0.82 floor is asm20's **seed-sensitivity envelope** — a
structural property, not a tuned threshold. To reach more divergent paralogs the *detector* must change, so
the refinement UNIONs additional edge detectors per family onto the asm20 core. The shared coverage floor
(over most of the copy length) is the false-merge defense; tiers are **additive** (merge divergent copies
into a family, never split asm20's calls).

**Sensitive nucleotide tier — PROMOTED to default.** `minimap2 -k11 -w5` (id≥0.70) recovers
nucleotide-divergent paralogs (rapidly-evolving KZNF / SAFB families, incl a reference-unannotated 98.6% SAFB
copy) down to ~73–76% identity. Cheap (one extra minimap2/family, run time unchanged), validated to add real
families with **0 false merges**, so it is now part of the default `--refine`: **asm20-only 152 → asm20+sensitive
155 families** (the new default). `--no-sensitive` reverts to the exact asm20-only 152 baseline.

**Protein tier — opt-in `--protein-tail`.** longest-ORF 6-frame → `mmseqs` (fident≥0.50, qcov/tcov≥0.50)
recovers **synonymous-divergent CODING paralogs** — the RABL2B retrocopy family at 87–99% protein but only
~70% genomic identity, which nucleotide seeds can never anchor (the qcov/tcov guard rejects lone-domain
merges like the IFITM1 case). **BATCHED into ONE mmseqs run** across all families (within-family hits only) —
byte-identical to the per-family result (158/158) but at **~0 added wall time** (6.2 min ≈ the 6.4 min
default), vs the per-family version that was the bottleneck. Degrades gracefully if `mmseqs` is absent.

Full union (`--refine --protein-tail`): **155 → 158 families** (+3 over the default; vs asm20-only the full
gain is +6 — the union *merges* divergent copies into existing families, e.g. RABL2B's divergent retrocopies
join RABL2, not the sweep's naive +14). Additions are real: DNA-confirmed **89.9%** at a divergent-appropriate
0.65 genomic gate (vs 87.3% at the strict 0.90 gate the divergent families intentionally fail). RABL2 intact;
a new 10-copy array appears.

> Catalog tiers: `--no-sensitive` = **152** (asm20-only; the substrate of the DNA-corroboration 90.1% and
> all caveats above) · default `--refine` = **155** (asm20+sensitive, promoted) · `--protein-tail` = **158**
> (adds synonymous-divergent coding paralogs). `cargo test --lib longest_orf_aa`.

## Honest caveats
- ⚠ **The old flagship/headline copy-assignment numbers are RETIRED.** "DSFAM817 = clean 3-copy, 90%
  assigned" and "CAFAM0 = 213 assigned @ 99.1%" were on OVER-MERGED false families (antisense +
  spliced/unspliced) that poasta's force-alignment masked with thousands of spurious PSVs. Do not cite
  them; cite the clean families above.
- The **silver standard is circular** (agrees with minimap2 where minimap2 was already confident); the
  load-bearing identifiability evidence is the sim5x labeled-truth ladder, not silver.
- The 281 family **regions** come from the older de-novo catalog (coordinate windows only); `copy_assign`
  re-detects copies per region on `GGO_mm.bam`. A fully clean pass would re-derive families genome-wide
  with the conflict graph + these fixes (the remaining O1 step). 1/281 regions timed out (150 s).
- `--min-copies 2` (real paralog PAIRS are the common case; `min-3` wrongly drops them).

### Cross-chrom / exon-sum caveats (from the adversarial review)
- **"Pure by construction" is retired.** The refinement certifies *homology* + *distinct-locus*, NOT
  functional-family purity: a residual case (GWFAM18 GSTM array) still carries one NBPF10 locus that does
  not link to GSTM at the protein level. Treat 157 as "homology-validated multi-locus families," not "157
  verified gene families."
- **70.9% orthogonal confirmation is a FLOOR, not precision.** Max-overlap gene assignment can grab a
  *neighbour* (EML5/ZC3H14), and the protein DB cannot test lncRNA/snRNA paralogs, so non-coding and
  partial-coverage families are *systematically under-confirmed* rather than wrong. Do NOT claim ~100%
  precision — the nucleotide re-check that found "zero false merges" uses the same homology criterion that
  built the catalog (circular); only the protein/name/Compara legs (71%) are orthogonal.
- **The 0.50 coverage gate is itself a non-canonical threshold** (Canzar-relevant): it drops genuine
  high-identity, partial-coverage paralogs at the boundary (e.g. fident 0.98 / qcov 0.47). Mitigated by the
  measured stability of the refined count across identity 0.70–0.90 (159→170 before the distinct-locus
  filter; 157 at the 0.80/0.50 operating point), but name it explicitly.
- **Refinement is now in-engine** (`gw_family_catalog --refine`, `denovo_pipeline::refine_families_exon_sum`):
  it emits the 157-family refined catalog directly (Rust reproduces the Python `validate_exon_sum.py`
  result byte-for-byte — 157 families / 54 cross-chrom / identical size distribution / RABL2 5-clique).
  `bench/validate_exon_sum.py` is retained as the independent cross-check + purity/precision reporter.
- **The 70.9% orthogonal-confirmation precision is a FLOOR set by annotation coverage, not catalog error**
  (0 false merges in the checked sample). It is depressed by three *resource* gaps, each separately
  fixable: non-coding paralogs absent from the protein DB (~16; need an Rfam/nucleotide resource),
  partial-coverage paralogs just under the 0.50 gate (~11; need ORF/domain-aware coverage), and
  `gene_at` neighbour-mislabels (~14; require the exon-sum to align the gene CDS). **Hard ceiling:**
  de-novo copies at *unannotated* loci have no orthogonal annotation resource and are irreducibly
  circular to validate from annotation — only a **DNA segmental-duplication map (SEDEF/BISER, the Soto
  2025 parCN approach)** is a truly orthogonal confirmation, and it also covers the non-coding and
  unannotated tails. That is the principled route to a *true* precision rather than a floor (not yet run).

## Reproduce
```bash
cargo build --release --bin copy_assign gw_family_catalog
bash /home/juanfra/winloci_scratch/ca_gw/run_clean_gw.sh      # minimap2+strand+dedup, min-copies 2, resumable
python /home/juanfra/winloci_scratch/ca_gw/aggregate_clean.py # genome-wide totals
# cross-chrom catalog, exon-sum (FLNC) homology + distinct-locus refinement IN-ENGINE (needs minimap2):
./target/release/gw_family_catalog --bam GGO_mm.bam --fasta GGO.fasta --out gw_xchrom_refined \
    --cross-chrom --refine                                    # -> 157 refined families, 54 cross-chrom
python bench/validate_exon_sum.py gw_xchrom_catalog minimap2 0.80 0.50   # independent cross-check + precision
cargo test --lib prune_same_locus minimap2_psv distinct_locus_reps homology_components   # the fixes' tests
```
