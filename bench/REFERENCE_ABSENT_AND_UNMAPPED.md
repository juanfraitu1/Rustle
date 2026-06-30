# Reference Absent And Unmapped (consolidated)

> Merged from 6 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [reference_absent_catalog](#reference-absent-catalog) · [hidden_collapse_headroom](#hidden-collapse-headroom) · [unmapped_rescue_poc](#unmapped-rescue-poc) · [crosschrom_discovery](#crosschrom-discovery) · [genomic_copies](#genomic-copies) · [divergent_paralog_track_scope](#divergent-paralog-track-scope)


---

## reference_absent_catalog

# Reference-absent gene-family copy catalog (gorilla IsoSeq vs T2T)

Gene-family copies present in the gorilla IsoSeq reads but **absent from the T2T reference assembly**.
Per `unmapped_rescue_poc.md`, such copies split into two populations — and the POC settled that the
gene-family-relevant ones live in the **mapped** pile, not the unmapped pile.

## Population A — collapsed / CNV-absent (mapped reads) — the main yield

A reference copy that is actually ≥2 collapsed copies (or has an extra CNV paralog the reference
lacks) shows, among its **primary** alignments (paralog-bleed firewall), a coherent **second
haplotype**: ≥12 balanced-frequency (0.20–0.60) alt columns co-segregating on ≥5 reads
(`hidden_copy.rs::detect_hidden_copy`, mirrored in `bench/hidden_copy_scan.py`). A diploid het is
1–2 columns; sequencing error never reaches the balanced band.

**Scan:** 357 reference-copy loci across 70 multi-copy families, 19 s. **30 loci flagged.**

**Validation (`hidden_copy_validate.py`):** the main RNA confound is **RNA editing** (A→I read as
A>G), which also makes co-segregating alt columns. A genuine collapsed copy has a *diverse* ref→alt
substitution spectrum; editing is A>G/T>C-dominated. Splitting the 30 flags by spectrum:

| verdict | n | criterion |
|---|---|---|
| **COLLAPSED-COPY** | **18** | ≥6 substitution subtypes, A>G+T>C < 40%, alt-fraction ≥ 0.20 — diverse divergence |
| ambiguous | 4 | borderline A>G+T>C (0.38–0.47) |
| weak / minority | 3 | alt-fraction < 0.20 (a minor haplotype; real low-expression copy or noise) |
| RNA-editing? | 5 | A>G+T>C ≥ 50% (e.g. GSTM 0.67 / 6 subtypes; DSFAM0 0.67 / 3 subtypes) — filtered out |

**The 18 collapsed-copy candidates** concentrate exactly where CN variation is expected:
- **Zinc-finger clusters** (DSFAM5, DSFAM8 on NC_073244.2) — textbook CN-variable segmental
  duplications; DSFAM5 carries collapsed copies at several loci.
- **Novel families** (DSFAM213 at 3 loci, DSFAM58 at 3, DSFAM104 at 2, DSFAM292, DSFAM26, DSFAM205,
  DSFAM243, DSFAM377).
- Strongest (balanced alt-fraction ~0.5, high depth): DSFAM5 (NC_073244.2:65003440, 478/1002 reads),
  DSFAM205 (166/382), DSFAM26 (193/340) — co-equal second haplotypes = a full collapsed copy.

## Population B — divergent-absent (unmapped reads) — dry on T2T

`unmapped_rescue_poc.md`: 5,519 unmapped reads (0.13%), 79% already-present (99.7% id), residual =
expressed TEs, leaving **exactly one** genuine hit: **cl0** (38 reads, 775-aa ORF, ~94% absent;
novel locus vs contaminant pending an `nt` BLAST). The route is a real lever only against an
*incomplete* reference (GRCh38-style); on complete T2T it is dry.

## Discipline & honest caveats

- **Detect + flag, never place** (DAZ3): each entry is *evidence the reads imply an unmodeled copy at
  this locus*, not a placed sequence. The next step to promote a flag to a catalogued copy is to
  **assemble the alt-haplotype reads into a consensus**, confirm a coherent ORF, and show paralog
  homology **distinct from every reference copy** (and ideally a **depth-doubling** cross-check —
  a collapsed copy carries ~2× the single-copy read depth).
- The ≥12-balanced-column firewall excludes diploid hets; the spectrum filter excludes RNA editing;
  the primary-only rule excludes in-reference paralog bleed. The 4 "ambiguous" + 3 "weak" need the
  consensus/ORF step before counting.
- Scope: the 70 de-tie multi-copy families. A genome-wide single-locus scan (a 1-copy reference gene
  that is actually 2) would extend coverage — the current scan requires the family to be known.

Artifacts: `hidden_copy_scan.py` · `hidden_copy_validate.py` ·
`/home/juanfra/winloci_scratch/refabsent/hidden_copy_catalog{,_validated}.{json,tsv}`.

## Promotion of the 18 strong flags (`promote_hidden_copies.py` + `promote_evidence_plus.py`)

Each flag's **alt-haplotype reads** were isolated, **POA-assembled** (pyabpoa) into the copy's
consensus transcript, then characterised: longest **ORF**, **genome homology** (map the consensus to
the whole T2T assembly — the correct "absent" measure), and **depth-doubling** (flagged-copy reads/kb
vs family siblings). The decisive RNA signal is **divergence**: the consensus maps back to its **own
locus** but the reference there is a *divergent paralog* — too divergent to be the copy itself or a
het allele, so the actual expressed copy is **absent**. (Depth-doubling is confounded by expression;
identity to a *family-model* copy mislead — the **genome** match is what counts.)

| class | n | meaning |
|---|---|---|
| **ABSENT-DIVERGENT-PARALOG** | **4** | consensus maps to own locus at **3.9–20.4 % divergence** from the reference — the copy is absent. DSFAM26 (4.8 %, ORF 362 aa, 2.2× depth), DSFAM58 (3.9 %, 266 aa, 2.0×), DSFAM58 (8.7 %, 250 aa), DSFAM292 (20.4 %, ORF 97 aa — weakest). All clustered on **NC_073229.2 ~49–50 Mb**. |
| **CNV-CANDIDATE (needs DNA)** | 4 | present in reference (98.2–99.8 % id) but **depth-doubled 2.0–10.2×** → expressed as more copies than modelled; het-vs-CNV is unresolvable from RNA alone (DSFAM377 at 99.8 %/10.2× = a near-identical collapsed expansion). DNA parCN would settle these. |
| ambiguous | 5 | 97.7–99.4 % id, depth not clearly doubled |
| het-likely (allelic) | 2 | ≥99.5 % id, low depth |
| ARTIFACT (consensus failed) | 3 | all DSFAM213 loci — POA consensus chimeric (genome id 36–81 %, cov < 0.9); inconclusive (heterogeneous alt-reads), not promotable as-is |

**Defensible result:** of the 18 strong flags, **4 promote to reference-absent divergent paralog
copies** (clustered at NC_073229.2 49–50 Mb; DSFAM26/58 are the strongest — divergent, coherent ORF,
*and* depth-doubled). A further **4 are CNV candidates** that RNA cannot separate from het — the exact
case **DNA parCN** resolves. The rest are het/ambiguous/consensus-artifacts. The funnel
(357 loci → 30 flags → 18 spectrum-clean → 4 promoted) is the honest yield.

Promotion artifacts: `promote_hidden_copies.py` · `promote_evidence_plus.py` ·
`/home/juanfra/winloci_scratch/refabsent/promoted/{consensus.fa,promoted_final.json,final_catalog.json}`.

## Locking the 4 copies (annotation + protein BLAST; `blastx` vs the 22,614-protein gorilla proteome)

**The 4 reference-absent copy CANDIDATES are all MHC genes**, clustered in the gorilla MHC on NC_073229.2
49–50 Mb — the single most polymorphic / copy-number-variable / reference-divergent region in the
genome, exactly where reference-absent expressed copies are expected:

| copy | locus | MHC gene (annotation) | nt div | gorilla-protein hit (blastx) | ORF |
|---|---|---|---|---|---|
| DSFAM26 | NC_073229.2:49051184 | class I (Gogo-B\*0201) | 4.8 % | LOC101142134 86.0 %, e=0.0 | 362 aa |
| DSFAM58 | NC_073229.2:50166630 | DRB1 (class II) | 3.9 % | LOC101131235 95.1 %, e=2e-172 | 266 aa |
| DSFAM58 | NC_073229.2:50380202 | DQ-β (class II) | 8.7 % | LOC101133052 88.7 %, e=2e-139 | 250 aa |
| DSFAM292 | NC_073229.2:50361413 | DQ-α (class II) | 20.4 % | LOC109027175 97.9 %, e=1e-57 | 97 aa |

**Contamination ruled out** (the `cl0` check, passed): every consensus hits the gorilla proteome at
86–98 % protein identity with e ≤ 1e-57 → endogenous gorilla MHC, not exogenous mRNA. The 4 CNV
candidates likewise name to endogenous genes: DSFAM243 → **PRAME/MAGE** (cancer-testis antigen, 99.8 %
protein, 2.2× depth), DSFAM5 → **ZNF766** (98.7 %, 2.2×), DSFAM377 → **uS5m** mito-ribosomal (99.8 %,
10× depth).

**MHC caveat (honest):** the MHC is hyperpolymorphic, so "reference-absent paralog copy" vs
"hyperdivergent allele" is genuinely blurry here — DNA/haplotype data settles which. But either way the
finding holds: **expressed MHC sequence 4–20 % divergent from the T2T reference's gene at that locus**,
endogenous, coding.

**DSFAM213 (3 artifacts) — rejected after retry.** Per-isoform ava-clustering + per-cluster POA still
yields consensuses mapping to their loci at only **59–68 % identity** with weak/no protein homology
(48.8 %, 56 % over 57 aa, no-hit) — repeat/chimera, not a paralog (paralogs > 80 %). Lesson: DSFAM213
carried the **highest** raw alt-column counts (235/123/35) yet is the **least** real — raw flag
strength does not predict a copy; the assembly + genome-homology + protein gate is the filter.

**Final result:** **4 expressed endogenous divergent MHC paralog CANDIDATES** (named, protein-confirmed
*endogenous* — protein homology does NOT establish copy-vs-allele) + 4 DNA-resolvable CNV candidates
(PRAME/MAGE, ZNF766, uS5m, +1). Locking artifacts: `promote_evidence_plus.py` · `blastx.tsv` · the
`proteindb` gorilla-proteome BLAST DB.
> ⚠ **Corrected (adversarial review #4):** these are **candidates, not confirmed copies**. (i) Protein BLAST
> confirms gorilla-MHC *endogeneity*, not copy status. (ii) The "divergence ⟹ copy" rule is **circular** — the
> reads are SELECTED for ≥12 reference-divergent ALT columns, so the consensus is divergent *by construction*.
> (iii) All 4 sit in the most hyperpolymorphic region of the genome (MHC), where a 3–20% divergent second
> **haplotype** is the *expected* diploid het, not a copy. The copy-vs-allele resolver (DNA parCN) was **never
> run** (no committed output). Report as "4 expressed divergent MHC haplotype candidates pending DNA"; the same
> applies to the "15 dispersed" and "905 collapsed" figures (mechanism-detected, not copy-validated).

## Genome-wide extension (`hidden_copy_scan_genomewide.py` + `promote_genomewide.py`)

Scan ALL 12,793 expressed de-novo loci (not just the 70 families) — so a single-copy reference gene
expressed as ≥2 copies is catchable. Funnel: **12,793 loci → 1,015 flags (8 %)** → 246 RNA-editing +
~560 het (>97 % id) + 42 chimeric/repeat consensus (>20 % div, the DSFAM213 mode genome-wide) + 60
weak/partial dropped → **73 paralog-copy candidates** (3–20 % divergence, endogenous gorilla-protein
hit, full coverage; `gw_clean_candidates.json`).

**These 73 are candidates, NOT confirmed copies.** A 3–20 % divergent second haplotype of a single-copy
gene is as plausibly a processed-pseudogene/retrocopy (reads cross-map to the parent), a segdup
paralog, or an alignment artifact, as a true CNV copy — RNA cannot separate them per-candidate. The
biologically compelling ones land in the expected rapidly-evolving / CN-variable families: **PRDM9**
(13.5 %, 817 aa — recombination-hotspot ZnF array), **SRGAP3** (8.0 % — the SRGAP2C family), **ZNF208/
ZNF558/ZBTB46/KBTBD3** (ZnF segdup clusters), plus immune/membrane genes.

**Conclusion: reference-absent expressed copies are CONCENTRATED, not uniform.** The clean confirmed set
stays the 4 MHC copies (clean because hotspot + defined multi-copy family). Genome-wide single-locus
scanning trades precision for recall — a candidate-generator (73, needing per-candidate validation:
intronless?=retrocopy, maps-elsewhere?=dispersed paralog, depth/DNA), landing in duplicated/immune/
rapidly-evolving families. The **high-precision lever is family/segdup-scoped scanning**.

GW artifacts: `hidden_copy_scan_genomewide.py` · `promote_genomewide.py` ·
`/home/juanfra/winloci_scratch/refabsent/{genomewide_flags*.json,gw_promoted/gw_clean_candidates.json}`.

### Discriminating the 73: multi-mapping + intronless (`gw_discriminated.json`)

**Multi-mapping** (minimap2 `-N20 -p0.2`): a real reference-absent copy of a *multi-locus* family maps
to its own divergent locus **plus** the reference paralogs; an in-place second haplotype (het/tandem)
maps to one locus only.
- **15 DISPERSED-PARALOG** = the promotable genome-wide set (own locus 3–15 % divergent + paralog loci):
  **PRDM9** (13.5 %, maps 2 loci, only 69 % id to its best paralog — a distinct member of the
  hypervariable recombination-hotspot family), **ZNF208** (13.6 %, 4 loci), ZDHHC13, plus large
  multi-locus families (LOC129525331: 12 loci / 28 introns / 1350-aa ORF) and CMTM7/ATP6V1E1/NOM1/etc.
- **58 single-locus** = het / tandem-CNV / absent-retrocopy — RNA cannot separate (DNA parCN needed).

**Intronless test — uninformative (honest negative).** An *absent* retrocopy's reads map to the
**parent** gene and splice out the parent's introns, so they look like a normal spliced transcript; the
intronless signature only appears where the retrocopy's own locus exists in the reference (it doesn't,
by definition). Only 2/58 single-locus were intronless (single-exon genes). So the multi-mapping test is
the one that discriminates; intronless cannot promote absent retrocopies from RNA→parent mapping.

**FULL CATALOG (all CANDIDATES — copy-vs-allele pending DNA, review #4):** 4 family-scoped MHC candidates +
**15 dispersed-paralog reference-absent candidates** (genome-wide; PRDM9/ZNF208/…) + 58 single-locus & 8
family-CNV candidates needing DNA. The candidates land — as biology predicts — in immune (MHC),
recombination (PRDM9), and zinc-finger/segdup
families. Discriminator artifacts: `gw_promoted/{clean73.paf,gw_discriminated.json}`.


---

## hidden_collapse_headroom

# Task (c): hidden collapsed-copy PSV signal at single-copy recall loci

Direct test of the headroom probe's undercount caveat: do read-coherence recall isoforms at
ANNOTATED single-copy loci actually sit on HIDDEN collapsed/cross-mapped paralog copies with
copy-discriminating PSV signal (annotation-free, called from the BAM pileup)?

Scanned single-copy recall loci with a linked PSV block (n_copies>=2 dumped). HIGH_COL=8.

## Verdict totals (loci with a linked PSV block)
| verdict | loci | recall isoforms | FSM |
|---|---|---|---|
| COLLAPSED_LIKE | 306 | 895 | 59 |
| AMBIGUOUS | 0 | 0 | 0 |
| HET_LIKE | 1397 | 4393 | 405 |

## n_coseg distribution (het-vs-collapse axis) — looking for a valley
```
 coseg   HET_LIKE  COLLAPSED
     0          0          0  
     1          0          0  
     2          0          0  
     3        707         46  ########################################
     4        313         29  ########################################
     5        173         16  ########################################
     6        113         10  ########################################
     7         91         11  ########################################
     8          0         52  ########################################
     9          0         32  ################################
    10          0         15  ###############
    11          0         21  #####################
    12          0         18  ##################
    13          0         10  ##########
    14          0         10  ##########
    15          0          9  #########
    16          0          4  ####
    17          0          5  #####
    18          0          2  ##
    19          0          1  #
   20+          0         15  ###############
```

## COLLAPSED_LIKE tiers (by evidence strength)
| tier | meaning | loci | recall isoforms | FSM |
|---|---|---|---|---|
| A_ge3groups | >=3 allele groups (not diploid het) | 22 | 158 | 6 |
| B_dense2copy | >=8 linked PSVs, 2 copies, low multimap | 167 | 437 | 40 |
| C_multimap | multimapping-driven (segdup/repeat possible) | 117 | 300 | 13 |

## Raw tier counts are NOT the headroom — confounds dominate (adversarial verification)
The adversarial workflow (5-mode methodology panel + hands-on BAM re-derivation, bench/wf_hidden_collapse_verify.js) confirmed every raw tier is dominated by a false-positive mode:
- **TIER B (≥8 linked cols, 2 copies) = diploid HETEROZYGOSITY.** A real 2nd genomic copy MUST
  multimap, but the TIER-B loci are UNIQUELY mapped (frac_mq0=0). The ≥8-column bar is mis-calibrated:
  scan windows average ~94 kb (whole genes, not transcripts), so a polymorphic gene trivially phases
  8–46 het SNPs (extreme cases: 46- and 30-coseg, both frac_mq0=0, MAPQ all 60, balanced minor frac).
- **TIER C (multimap) = mostly SEGDUP SPILLOVER.** frac_sec is structurally invalid: secondary
  records carry no SEQ → contribute zero alleles; the 'copies' are built from RESOLVED MAPQ-60
  primaries then OR-overridden by evidence-free spill-in reads.
- **RNA EDITING penetrates even TIER A.** Pure A>G/T>C transition spectra co-segregate as fake
  haplotypes and mint spurious ≥3-'copy' splits (36% of TIER A editing-suspect).

## Deterministic confound-controlled headroom
frac_mq0 bands over 306 COLLAPSED_LIKE: ≥0.3 (genuine local multimap)=2; [0.10,0.3)=4; <0.10 (uniquely mapped = het/editing/spillover)=300.
**Joint gate (frac_mq0≥0.3 AND n_coseg≥8) → 0 loci, 0 recall isoforms, 0 FSM.**
The 2 genuinely-multimapping loci all have n_coseg<8 (diploid-het floor, all 0 FSM) → not confidently collapse vs het-at-a-multimap-locus.

## VERDICT
- raw COLLAPSED_LIKE: 306 loci (895 recall isoforms, 59 FSM) — ALL confounds.
- naive 'TIER A+B' (BEFORE verification): 189 loci / 595 iso / 46 FSM — refuted (het + editing).
- **confound-controlled hidden-collapse headroom = 0 loci / 0 iso / 0 FSM.**
- **GO/NO-GO: NO-GO.** The direct-BAM scan, after confound control, finds 0 PSV-resolvable hidden
  collapse at single-copy recall loci — confirming the geometric probe's 0. The undercount caveat
  does NOT open real copy-resolution headroom.

Completeness gaps (honest, the detector cannot see): identical-sequence copies emit no PSVs (invisible); copies whose reads map to a separate locus/contig (RABL2/DAZ regime) never appear in this window; strand-blind (editing detected but not auto-gated); indel/STR copy differences discarded; no independent (segdup/Compara) paralog ground truth.


---

## unmapped_rescue_poc

# Unmapped-read rescue PoC — can the unmapped pile recover gene copies missing from the reference?

**Date:** 2026-06-21 · **Data:** T2T gorilla (`GGO`) PacBio IsoSeq, mapped to the
T2T gorilla assembly · **Reproduce:** `bench/unmapped_rescue_poc.py` (every number
below is regenerated by the script from the preserved artifacts)

**Verdict (negative as a lever; one genuine n=1 hit):** On a near-complete (T2T)
reference, "rescue the unmapped reads → find expressed gene copies missing from the
genome" does **not** pan out as a productive lever. The unmapped pile is 0.13% of
reads, and once de-novo clustered and remapped it is dominated by reads whose
sequence is **already present** in the genome (single-locus, median 99.7% identity)
plus a chimera/read-through minority. Following up the residual (below) dissolves it
into **expressed transposable elements** (cl19/cl31) and one divergent single locus
(cl56), leaving **exactly one genuine hit: cl0** — a complex, reproducible (38 reads),
**protein-coding (775-aa ORF)** ~5.5 kb transcript that is ~94% **absent** from the
assembly. So the advisor's concern is *real in microcosm* (cl0 exists) but
**quantitatively negligible** on a complete reference (yield ≈ 1 in 5,519 reads). The
lever would matter against an **incomplete** reference (GRCh38-style); here it is dry.

---

## Why test this

The advisor's standing concern is **reference bias**: a copy of a gene that is
expressed but absent from the reference assembly can never align, so an
alignment-keyed assembler is blind to it. The cheapest place such copies would
surface is the **unmapped read pile**. So we tested the hypothesis directly:
extract every unmapped read, cluster the pile de-novo (no reference), and ask
whether the reproducible clusters are genuinely-absent expressed copies.

## Method

```
unmapped.bam      = samtools view -b -f 4 GGO.bam                  # primary-unmapped reads
unmapped.fa       = samtools fasta unmapped.bam
ava.paf           = minimap2 -x ava-pb -X unmapped.fa unmapped.fa  # all-vs-all overlaps
clusters          = connected components of the ava.paf overlap graph
reps.fa           = one representative read per component of size >= 5
reps.relaxed.sam  = minimap2 -ax map-pb -p 0.1 -N 50 GGO.fasta reps.fa   # NON-spliced remap
```

The relaxed (non-spliced) remap is the probe. `-p 0.1 -N 50` reports up to 50
secondary placements, so paralog homology is visible. A read's **"pieces"** are its
**primary + supplementary** alignments (the colinear split); **secondary**
alignments are alternative paralog placements and are tracked separately, *not*
counted as chimeric pieces (see the correction note). Reading the remap:

- chimera / read-through → pieces map to ≥2 loci (primary + supplementary);
- splice-rejected real transcript → maps colinearly to one existing locus;
- genuinely novel/absent sequence → most of the read aligns nowhere.

## Results

**The pile is tiny and full-length — not junk:**

| | |
|---|---|
| Unmapped reads | **5,519** = **0.125% ≈ 0.13%** of 4,409,946 primary reads (flagstat) |
| Read length | median **2,920 bp**, mean 2,934, range 50–10,644 — full-length |
| Reads forming reproducible overlap clusters (≥3 reads) | **345 clusters / 1,605 reads** |
| Clusters of ≥5 reads (representatives remapped) | **103** (largest = 38 reads) |

**Relaxed remap of the 103 cluster representatives:**

| | Count | Reading |
|---|---|---|
| Unmapped even in relaxed mode | 2 | — |
| Got a primary alignment | **101** | |
| └ primary **MAPQ 60** | 88 (87%) | primary uniquely resolved |
| └ primary MAPQ < 60 | 13 | some paralog ambiguity in the *primary* |
| └ primary MAPQ 0 | **0** | every read gets a decisive primary |
| └ primary overlaps an annotated gene | 93 (92%) | lands on existing genes |
| └ best-chain coverage ≥80% of the read | 63 (62%) | read ~fully explained colinearly |
| **reports ≥1 secondary placement** | **47 (47%)** | **paralog homology IS present** (41 cross-chromosome) |

**Classification of the 101 mapped reps** (pieces = primary + supplementary):

| Category | Count | What it is |
|---|---|---|
| **Single-locus, present** | **80 (79%)** | colinear to one existing locus; **median substitution identity 99.7%** → the sequence is in the genome. Splice-rejected (unspliced / intron-retention / genomic), not a missing copy. |
| **Multi-locus chimera** | **18 (18%)** | pieces map to ≥2 loci → read-through / fusion transcripts (3 are clean 2-piece, ≥80% covered) |
| **Novel-sequence candidate** | **3 (3%)** | <50% of the read aligns anywhere (threshold-sensitive: 1 / 3 / 10 at cov<0.4 / 0.5 / 0.6) |

## Residual-candidate follow-up (option b — done here)

We chased the only clusters where a "missing copy" could plausibly hide. The
follow-up (reproduced by `step5_candidates` in the script) resolves all of them:

| cluster | len | aligns | ORF | MAPQ | #aln/#loci/#chr | classification |
|---|---|---|---|---|---|---|
| **cl0** | 5,734 | 5% | **775 aa** | 1 | 5 / 5 / 4 | **NOVEL protein-coding transcript, absent from assembly** |
| cl19 | 3,156 | 85% | 98 aa | 1 | 51 / 22 / 11 | expressed **transposable element** (dispersed, ~50 genomic copies) |
| cl31 | 3,451 | 65% | 91 aa | 1 | 51 / 18 / 13 | expressed **transposable element** (dispersed, ~50 genomic copies) |
| cl56 | 3,028 | 63% | 130 aa | 60 | 5 / 2 / 2 | divergent/unannotated **single locus** (present; ~36 kb from LOC101128475) |
| cl34 | 4,762 | 43% | 153 aa | 60 | 8 / 1 / 1 | partial transcript anchored in ATP1A1 |
| cl102 | 2,434 | 43% | 181 aa | 30 | 1 / 1 / 1 | partial transcript anchored in ATP5F1B |

- **cl19 / cl31 are expressed transposable elements / dispersed repeats**, not gene
  copies: each read has ~50 near-equal-scoring genomic placements (primary AS ≈ best
  secondary AS) scattered across 11–13 chromosomes, MAPQ 1, weak ORF. The "~80%
  identity to intergenic" is repeat-to-repeat matching.
- **cl56** is a single-copy (MAPQ 60, AS 571 ≫ 360) intergenic locus expressed at
  ~80% identity — an unannotated/divergent locus, but the sequence *is* in the
  assembly. Minor.
- **cl0 is the one genuine find.** It is *not* low-complexity (near-maximal k-mer
  entropy, no tandem periodicity — top period match 0.286 vs 0.260 random) and *not*
  a repeat (only 4 short scattered anchors, not dispersed copies). It is a complex,
  reproducible ~5.5 kb molecule (38 reads; full-length mode + truncation spread) with
  a **775-aa stop-free ORF** (random expectation ≈ 21 aa) — i.e. a bona-fide
  **protein-coding mRNA** — of which only a 306 bp (5%) fragment anchors weakly
  (MAPQ 1, ~84% id) to PARP4, leaving **~94% with no genomic match**. (No polyA tail,
  but IsoSeq `refine` strips polyA from FLNC, so that is expected/uninformative.)

**What cl0 still needs:** its species of origin — a genuinely novel/missing gorilla
locus / assembly gap vs. an exogenous contaminant mRNA — cannot be called offline; it
requires an `nt`/`nr` database (BLAST) search. That is the single clear next step if
anyone pursues it. Either way it is one transcript, so it does not change the
"dry as a lever" verdict.

## Correction note (supersedes the in-session figures)

An earlier in-session pass reported "~74/101 multi-locus, chimera-dominated, 1 novel
(cl102)." That over-counted multi-locus reads by treating scattered **secondary**
alignments as separate loci. Re-derived here with pieces = primary + supplementary
only, the chimera fraction is **18%, not ~73%**; the bulk (79%) is single-locus
splice-rejected reads at 99.7% identity; and the biggest cluster cl0 is a single
~5%-coverage MAPQ-1 anchor, **not** a clean read-through across 5 genes. The
headline conclusion (dry for missing copies) is unchanged.

## What this means for the advisor's reference-bias concern

1. **"Copies not in the annotation"** — the real win is the population that *does*
   map but isn't annotated: our **de-novo expressed loci** (the de-tie conflict-graph
   families define copies from reads, not gene models). The unmapped pile is not
   where that signal lives.
2. **"Copies not in the genome"** — a divergent/CNV copy with residual homology
   gets **force-mapped onto its paralog**, not flagged unmapped. Consistent with that,
   the **hard, unresolvable** copy-ambiguity regime (tied paralogs, MAPQ-0 *primaries*)
   does **not** concentrate in the unmapped pile: every remapped read here gets a
   decisive primary (0 MAPQ-0 primaries, 87% MAPQ 60). Paralog homology *is* present
   (47/101 reps report secondary placements) — but it is resolvable, i.e. exactly the
   kind that minimap2's primary already picks. The genuinely tied regime therefore
   lives in the **mapped** MAPQ-0 pile (the copy-assignment problem), not here.
   *(Caveat: the MAPQ figures are from the relaxed re-map of cluster reps, not the
   original mapping; they characterise the reps, they are not a direct measurement of
   the original unmapped event.)*
3. **Scope** — the unmapped-rescue lever would matter against an **incomplete**
   reference (GRCh38-style), where SDA-style absent copies are common. The concern is
   real in general; it is just dry in this T2T gorilla data. Defensible statement:
   *"I tested it directly — in T2T it is dry (0.13% unmapped; 79% present-sequence at
   99.7% identity; 18% chimeras; 0 MAPQ-0 primaries; a few residual diverged/novel
   candidates); it is a real lever only against incomplete references."*

## Reproduce

```bash
# preserved artifacts (too large to commit) live in the scratch dir:
#   /home/juanfra/winloci_scratch/unmapped_poc/
#     unmapped.fa, ava.paf, reps.relaxed.sam, genes.bed, ggo_flagstat.txt,
#     clusters_ge5.json, low_cov_anchors.json   (+ unmapped.bam, reps.fa)
python3 bench/unmapped_rescue_poc.py [scratch_dir]
# -> prints the tables above and writes unmapped_rescue_poc_result.json
```

## Caveats / scope

- **Single dataset** (one T2T gorilla IsoSeq sample). The *quantitative* yield is
  specific to a near-complete reference; the qualitative conclusion (the unmapped pile
  is not where missing copies live, on a complete reference) is the transferable part.
- **Novel cutoff is fragile:** the "3 novel" count depends on the `COV_NOVEL=0.5`
  threshold (1 / 3 / 10 at 0.4 / 0.5 / 0.6); reported as a sensitivity, not a hard number.
- **"Present" ≠ 100% covered:** the *single-locus present* bucket includes 35 reps at
  50–80% best-chain coverage — most, not all, of the read aligns (at high identity over
  the aligned part).
- Clusters are connected components of an all-vs-all overlap graph at minimap2
  `ava-pb` defaults; a stricter/looser threshold shifts cluster counts but not the
  classification of representatives.
- Gene assignment uses largest-overlap against the `gene` features of the gorilla
  RefSeq/Gnomon annotation (`GGO_genomic.gff`); intergenic primaries (8/101) are
  reported as such. Substitution identity = 1 − (NM − inserted − deleted) / M, i.e.
  mismatches per aligned base (avoids counting indels/introns as divergence).


---

## crosschrom_discovery

# Genome-wide cross-chromosome copy discovery (RNA-level)

**Question (user):** the production family grouping gathers copies per genomic region (position-overlap bundles), so DISPERSED paralogs whose copies sit on DIFFERENT chromosomes (RABL2A/B; the headroom probe's 17 DISPERSED families) are never co-considered. Can the harness be improved to find cross-chromosome resemblant copies?

**Answer: yes.** A genome-wide discovery harness removes the chromosome restriction and finds them precisely.

## Method (bench/extract_gene_reps.py + family_crosschrom_discovery.py + crosschrom_grade.py)
1. Extract one representative spliced transcript per gene, genome-wide: **22,983 gene reps** (longest transcript per RefSeq gene, +strand-oriented).
2. **Minimizer-LSH with NO chromosome restriction** (k=15/w=10 canonical; inverted index, skip repetitive minimizers >200 genes; pairs sharing ≥4 minimizers, Jaccard≥0.03) → 18,934 candidate pairs. Minimizer-Jaccard is only a loose PREFILTER (real diverged dups like RFPL sit at Jaccard ~0.13).
3. **POA contiguous-core gate** (the validated RNA-level definition, bench/poa_family_definition.py): largest single ungapped co-aligned block ≥ T=0.13 of the shorter gene.
4. **Grade by POA core-block IDENTITY** — the discriminator core_cov alone lacks.
- The user's note (minimizers are useful but not the only option) is respected: minimizers are just the prefilter; the gate is alignment. A full **intron-chain alignment** is a planned second (structural) candidate axis — dispersed copies keep their relative intron-chain structure.

## Recall — it finds the known cross-chromosome families
- **8/8 universe cross-chromosome families recovered (RABL2A/B + 7 LOC families).** RABL2A (NC_073235.2) ↔ RABL2B (NC_086018.1): core_cov 0.17, **core_ident 0.99** (recent dup; short but near-identical core — exactly why core_cov alone is low and T=0.13 was calibrated to keep it).
- Named recent cross-chrom paralogs found (all core_ident≥0.97): CRIPTO/CRIPTO3, GK/GK3, EIF2S3/EIF2S3B, HNRNPA1/HNRNPA1L2, PGAM1/PGAM4, SLC25A51/SLC25A52, METTL2A/METTL2B — textbook dispersed paralogs.

## Precision — the per-pair signal is clean (core-identity, not Jaccard)
- **All 8304 cross-chrom pairs have POA core-block identity ≥ 0.4** (7387 ≥0.9, 5846 ≥0.95); **none sit at the ~0.25 chance baseline.** There are no chance-alignment false positives.
- Minimizer-Jaccard was a BAD precision axis: the apparent 'chance' pairs (EEF1A1↔LOC etc., Jaccard<0.1) are real retrocopies/processed pseudogenes at 0.89–0.99 core-identity — their low Jaccard is just a length-mismatch artifact (a short copy vs a long parent). Core-block identity is the right axis.
- Recent-duplicate subset (core_ident≥0.95): **5846 pairs**.

## Honest caveats / residual false-positive modes
- **Family-level transitive over-merge.** Connected components over a permissive pair gate chain distinct subfamilies through domain hubs: the largest 'families' are 145- and 114-gene components (gene SUPERFAMILIES, not single families). The 428 **size-2** families are clean copy pairs; large components need tighter clustering (mutual-best / all-pairs-above-bar) than transitive closure. (703 components total.)
- **Short high-identity shared elements.** A few pairs sit just above the core_cov=0.13 gate with high identity over a SHORT block = a shared transposon/element between otherwise-unrelated genes (e.g. IGFL3↔USP12: 16% core at 98% id). The recency filter does NOT remove these (they are high-identity); whole-gene-fraction or element-masking would. Raising core_cov would drop them but also drops real short recent dups like RABL2 (0.17) — the binding tension.
- **Recency spectrum.** Pairs span recent duplicates (core_ident≥0.95) → older paralogs/retrocopies → ancient domain-based families. 'Recent duplicate' (the advisor-defensible claim) = the high-identity subset; broader 'resemblant copies' = the full set.
- **One representative isoform per gene** (gene_rep); a different isoform could shift coverage at the margin. **Input = RefSeq gene reps** (not assembled transfrags) — swapping in rustle/StringTie output is a one-line input change (same harness).
- **Scope = the 25 large contigs' RefSeq genes** (22,983); whole-genome contigs add more.

## Verdict
The harness improvement WORKS: removing the chromosome restriction (genome-wide minimizer-LSH → POA contiguous-core → core-identity grading) recovers **8/8** known cross-chromosome families and finds hundreds more, all genuinely homologous per-pair. The headline gap it closes: production groups copies per region, so it structurally cannot see RABL2A/B; this finds them. Remaining work is family-level clustering (de-transitive-merge) and an element-sharer filter — both refinements, not blockers.

## Reproduce
- `python3 bench/extract_gene_reps.py` (genome-wide gene reps → /tmp/gene_reps_gw.fa)
- `MINIFORGE python bench/family_crosschrom_discovery.py --stage all` (LSH → POA → families)
- `MINIFORGE python bench/crosschrom_grade.py` (core-identity grading → crosschrom_graded.tsv)
- `MINIFORGE python bench/crosschrom_writeup.py` (this writeup + figure)


---

## genomic_copies

# Genomic self-alignment — completing the copy roster (pseudogene + unannotated copies)

The DNA/protein tier (mmseqs2 on translated CDS) catches ancient *coding* families but is blind to two
copy classes: **pseudogene copies** (no CDS → no protein) and **unannotated copies** (no annotation at
all). This pass finds them by genomic self-alignment: map every gene-rep transcript back to the genome
(`minimap2 -cx splice -N 20 -p 0.4 --secondary=yes`) and classify each strong homology hit (≥50% of
the transcript aligned at ≥80% identity) by the annotated feature it lands on.

## Result
PAF hit classes: self=23,123 · annotated_gene(paralog)=18,495 · **pseudogene_copy=2,825** · **unannotated=1,949**.

Merged to distinct loci: **1,978 NEW copy loci** the annotation + protein-clustering missed:
- **pseudogene copies: 1,313** (annotated as pseudogenes — the protein-clustering blind spot)
- **unannotated copies: 665** (no annotation at all — the annotation blind spot)

## RNA overlay — what transcribes (the functional readout)
- new copies are **mostly silent**: pseudogene copies 18% transcribed (239/1,313), unannotated 16%
  (104/665) — vs **72%** for coding-family copies. Exactly the expected biology: duplicated/pseudogenized
  copies are largely dark, with a transcribed minority.
- the transcribed minority is the interesting part — e.g. an **SDHA retrocopy** (335 reads, 94.8% id),
  and several high-identity LOC copies expressed at 180–390 reads. These are real transcribed copies
  invisible to both the annotation and the protein tier.

## The complete picture (three tiers + RNA overlay)
| copy class | source | copies | transcribed |
|---|---|---|---|
| annotated coding-family | mmseqs2 protein clusters | 14,545 | 72% |
| pseudogene copy | genomic self-align → annotated pseudogene | 1,313 | 18% |
| unannotated copy | genomic self-align → outside all annotation | 665 | 16% |

**DNA enumerates every copy** (annotated coding + pseudogene + unannotated); **RNA shows which are
actually transcribed**, copy-by-copy. The silent vs transcribed split *is* the genome-vs-transcriptome
result, and copy-assignment (the identifiability theorem) then says which of the transcribed copies are
individually resolvable.

## Honest scope
- Copies found by mapping ANNOTATED transcripts → genome, so a fully novel family with NO annotated
  member anywhere would still be missed (would need de-novo transcript assembly as the query).
- ORF-intactness of the new copies (true pseudogene vs intact-but-unannotated gene) is the natural next
  annotation — read the ORF off the transcribed copies' reads.
- Identity/coverage thresholds (≥0.8 id, ≥0.5 qcov) set the recent/divergence horizon for "a copy."

## Reproduce
- annotated intervals: extract `gene`/`pseudogene` features (chrom,start,end,biotype,gene) from
  `/tmp/gw/ref_*.gff3` → `annot_intervals.tsv`
- `MINIFORGE python bench/build_genomic_copies.py` (minimap2 self-align + classify + RNA overlay)
- `python3 bench/genomic_copies_fig.py`


---

## divergent_paralog_track_scope

# Scope: a divergent-paralog detection track (beyond the contiguous-core horizon)

_Status: SCOPE / not started. Author handoff doc. Builds on `compara_validation.py`, `family_detection_validation.py`, `denovo_family_pipeline.md`._

## 1. Motivation — the gap, with evidence

The shipped family detector defines a family by POA **contiguous-core coverage ≥ `t_core` (0.13)** — the longest *contiguous exact* shared run / min(len). This is tuned for **recent tandem duplicates** and is correctly conservative.

It is blind to **divergent paralogs**. Measured on DSFAM45 (`diag_prefilter_homology`, NC_073234.2:48446103-49179309):

| pair | shared 18-mers | longest contiguous run | core cov (poasta == LCS) |
|---|---|---|---|
| [3]×[4] (≈2 kb copies) | **374** | ~76 bp | **0.036** |
| unrelated pair | 2–9 | — | 0.003–0.01 |

374 shared k-mers vs 2–9 for unrelated ⇒ unmistakable common descent; but the longest contiguous exact run is ~76 bp, so contiguous-core (poasta **and** the alignment-free LCS agree — not an artifact) puts it at 0.036, far under 0.13. **These are real paralogs past the recent-duplicate horizon.**

**Why it's worth detecting (the through-line):** divergent paralogs carry MANY fixed differences → by the identifiability theorem they are the *easy, high-PSV* regime for copy-assignment (the advisor's #2 interest). The families we currently can't *detect* are the ones we could most cleanly *assign*. Detection is the gate; assignment is the payoff.

## 2. The central question (sharp, falsifiable)

> Among de-novo locus pairs the contiguous-core criterion **rejects**, can a sequence feature — validated against an **independent, LOC-covering** paralog truth set — separate TRUE divergent paralogs from FALSE domain/repeat-sharers, at a **calibrated FDR**?

Reject criterion to beat: a clean, ROC-calibrated edge predictor with reported FDR (no arbitrary threshold) that recovers DSFAM45-type families without admitting domain-sharers.

## 3. The make-or-break — the truth set — is largely already in place

The original sin (`family_detection_validation.py`) was **circularity**: a minimizer-Jaccard clustering scored against a minimap2-built universe — both reward shared subsequence. `compara_validation.py` de-circularized it with Ensembl Compara (gene-tree + species-tree, never sees our alignments) but hit a **coverage wall: 154/195 genes are LOC-only → unmappable to Ensembl → only 12 checkable pairs.** Too small to calibrate a threshold.

**The fix we can now afford — OrthoFinder on `proteins.fa`:**
- `/home/juanfra/winloci_scratch/proteins.fa` = the RefSeq gorilla proteome. RefSeq annotates **LOC genes with protein products (XP_)**, so OrthoFinder clusters them too → the truth set **covers the LOC paralogs Compara could not** (the majority, and exactly the interesting copies).
- Self-contained: no Ensembl-ID mapping, no gorGor↔T2T liftover.
- Independent of a DNA k-mer feature (DIAMOND protein all-vs-all → MCL → gene trees).
- The advisor-endorsed modality ("need protein-level (OrthoFinder) for a real claim").

**Truth plan:** OrthoFinder orthogroups = PRIMARY truth (LOC-covering, hundreds–thousands of paralog pairs). Compara (already cached, `compara_paralog_relation.json`) = independent CROSS-CHECK on the mappable named subset (guards against an OrthoFinder-specific bias). Granularity caveat carries over: orthogroups are COARSER than copy-families ⇒ score **precision** (our pairs that the truth also relates), not recall.

## 4. Feature candidates (test both; lead with the independent one)

| feature | independence vs protein truth | sensitivity to divergence | cost |
|---|---|---|---|
| **DNA k-mer containment / Jaccard** on the de-novo nt seq | FULL (cross-modality) | moderate | cheap (have the k-mers already) |
| **Protein-ORF identity** (translate longest ORF → AA → identity/k-mer) | PARTIAL (both protein; mitigate: pairwise vs OrthoFinder's tree-clustering) | high | needs ORF-finding on de-novo tx |
| shared-k-mer **count**, length-normalized | FULL | moderate (the DSFAM45 signal: 374 vs ~5) | cheap |

Lead claim = **DNA/k-mer feature vs OrthoFinder truth** (clean cross-modality). Report protein-ORF as the sensitive-but-caveated upper bound. If the cheap DNA feature already separates on independent truth, that is the strongest, simplest result.

## 5. Phased plan — gates, kill-criteria, effort

**Phase 0 — Truth set + labeled pairs (≈1–2 d). THE gate.**
- Run OrthoFinder (or a direct DIAMOND-all-vs-all + MCL if OrthoFinder is too heavy) on `proteins.fa` → orthogroups.
- Map each de-novo locus → overlapping RefSeq gene (`GGO_genomic.gff`) → protein → orthogroup.
- Build labeled pair set: POSITIVE = same orthogroup; NEGATIVE = length/GC-matched different-orthogroup pairs.
- **GATE:** ≥ a few hundred labeled paralog pairs INCLUDING LOC genes. If not (truth too sparse even with OrthoFinder) → fall back to Compara-only, smaller claim, and flag the reduced scope. **KILL if no independent truth materialises.**

**Phase 1 — Is there signal? (≈1 d).**
- Feature AUC vs the independent label, restricted to the contiguous-core-REJECT population (cov < `t_core`) — i.e. exactly the pairs the shipped detector misses.
- **KILL if best AUC < ~0.7 on independent truth** → divergent paralogs are sequence-indistinguishable from domain-sharers (an honest, publishable negative; stop).
- Decide the feature.

**Phase 2 — Calibration / principled edge (≈1 d).**
- ROC-calibrate the operating threshold on the truth; report FDR at the operating point. NO arbitrary cutoff.
- Advisor-principled option: emit the feature as a weighted edge and let the EXISTING weighted-modularity Louvain (`family_split.rs`) decompose — the threshold folds into the objective.

**Phase 3 — Integration as a SEPARATE, gated track (≈2–3 d).**
- New `divergent_family_detect` path running ALONGSIDE (not replacing) contiguous-core: same candidate_pairs prefilter, but the REJECTED-by-core pairs get the validated divergent edge criterion; reuse `decompose_families`.
- Opt-in (`RUSTLE_DIVERGENT_PARALOGS`), default-OFF → byte-identical baseline (the project's standard for every lever). TDD, adversarial review.
- Apply to DSFAM45 → does it recover the family at the reported FDR?

**Phase 4 — The payoff: does it help copy-assignment? (≈1–2 d).**
- On the recovered divergent families, run the existing copy-assignment + measure resolvable-PSV / identifiability. Hypothesis: divergent paralogs are HIGH-identifiability (many PSVs) → clean assignments. This closes the loop to the advisor's interest and is the real success bar.

## 6. Principled framing (for Canzar)

- **No arbitrary threshold:** the edge criterion is ROC-calibrated on an independent biological truth with a reported FDR; ideally folded into the existing weighted-modularity objective so the cutoff is emergent, not hand-set.
- **Precision, not recall:** truth is coarser (superfamilies) ⇒ we score "of the pairs we join, what fraction the independent truth also relates," and intentionally do NOT penalise correctly splitting ancient superfamilies into copy-clusters (the established framing).
- **Through-line:** detection feeds the identifiability story — divergent paralogs are the high-PSV, cleanly-assignable regime. The track expands the set of families we can *resolve*, which is the thesis core.

## 7. Risks / open decisions for the user

1. **Truth source priority.** OrthoFinder-primary (more setup, breaks the 12-pair ceiling, the real claim) vs Compara-only (already cached, but only 12 mappable pairs → underpowered). *Recommend OrthoFinder-primary.*
2. **Feature priority.** DNA-independent (cleanest claim) vs protein-ORF (most sensitive, needs de-novo ORF-finding, semi-circular caveat). *Recommend lead DNA, report protein.*
3. **Success bar.** Phase 3 (recover DSFAM45-type families at a calibrated FDR) — or push to Phase 4 (show the recovered families assign cleanly), which is the advisor-facing payoff. *Recommend committing to Phase 4 as the bar.*
4. **Honest-negative acceptance.** The kill-criteria can end in "no separable feature" — a real result. Agree up front that a clean negative is an acceptable outcome (it bounds the contiguous-core decision as provably optimal on sequence alone).

## 8. First concrete step (SUPERSEDED — see §9)

~~Phase 0: get orthogroups from `proteins.fa` (OrthoFinder)...~~ Replaced by the own-tool, RNA-only direction below.

## 9. CONSTRAINTS (user) — own-tool, RNA-only, Phase-4 bar

Three decisions reshape the track:

1. **No OrthoFinder dependency — extract & implement its pertinent logic ourselves.** OrthoFinder's algorithm that matters here = a **length/composition-NORMALISED all-vs-all similarity graph → graph clustering (MCL) → orthogroups**. We already have the graph (`candidate_pairs`) and a clustering engine (`family_split::decompose_families`, weighted Louvain). So the divergent-paralog DETECTOR is a small extraction, not a new tool:
   - (a) all-vs-all similarity on the REJECTED-by-core pairs, with a **normalised** score, and
   - (b) **low-complexity / tandem-repeat masking** before scoring — this is the load-bearing piece for DSFAM45: its 374 shared k-mers must be shown to come from an orthologous divergent core, NOT a shared repeat. Masking + length-normalisation is exactly how OrthoFinder-style scoring separates true paralogs from repeat/domain-sharers. We implement that principle on RNA.
   - (c) cluster with the existing weighted Louvain. No new arbitrary threshold beyond the ROC-calibrated edge.

2. **RNA-only at ship time.** The shipped feature is computed from the de-novo transcript **nucleotide (RNA)** sequence only — k-mer containment / normalised shared-k-mer count / LCS — needing NO genome FASTA and NO proteome at inference. (Any ORF translation, if ever used, is derived from the RNA itself and stays self-contained; lead with nucleotide features to honour "no protein.") The genome/protein are allowed ONLY in the offline validation, never in the shipped detector.

3. **Success bar = Phase 4** (recovered divergent families must ASSIGN cleanly — high per-read PSV identifiability).

### What this changes vs §3–5
- **Truth stays external, but as DATA not a tool.** Ensembl **Compara** paralog relations (already cached, `compara_paralog_relation.json`) are a published gene-tree oracle that never sees our RNA — using them as ground truth is NOT a tool dependency and keeps validation non-circular. We do **not** run OrthoFinder.
- **The cost of dropping OrthoFinder = the LOC-coverage wall reopens.** Compara maps only NAMED genes (154/195 universe genes are LOC → unverifiable). Mitigation: the prior 12 pairs were a curated UNIVERSE subset — map **all** de-novo loci genome-wide → named RefSeq genes (`GGO_genomic.gff`) → cached Compara paralogy for a far larger NAMED calibration set; calibrate the RNA feature there, then APPLY genome-wide incl. LOC (predict at the calibrated FDR; LOC stays unvalidated — honest, as the prior report already stated). **KILL** if the genome-wide named+Compara set is still too small to calibrate → fall back to a qualitative DSFAM45 demo + caveat.

### Revised Phase 0 (two parallel sub-tasks)
- **(i) Detector (RNA-only, our own):** implement normalised all-vs-all similarity + low-complexity masking over the core-rejected candidate pairs; feed weighted edges to the existing Louvain. Default-OFF (`RUSTLE_DIVERGENT_PARALOGS`), byte-identical baseline, TDD.
- **(ii) Truth (offline, Compara-as-data):** map de-novo loci → named RefSeq genes → cached Compara paralogy → a labeled NAMED-pair table, genome-wide. **GATE:** ≥ ~100 named paralog pairs incl. divergent ones; else underpowered → qualitative only.

### First concrete step
Build the offline Compara-based labeled NAMED-pair table genome-wide (sub-task ii) — it is cheap (reuses `compara_fetch.py`/cache + the GFF) and it is the GATE: it tells us whether there is enough independent truth to calibrate an RNA feature at all, before we invest in the detector. If the gate passes, implement sub-task (i).

## 10. GATE RESULT (Phase 0) — PASS

Run: `divergent_truth_gate.py` + `divergent_compara_label.py` (artifact: `divergent_calibration_pairs.tsv`).

- locus→gene coverage (edge loci): 67% named / 28% LOC / 5% unmapped.
- The divergent signal is in DIFFERENT-named low-core pairs (NOT same-root — those are mega-family ZNF/OR blobs).
- **GOLD calibration set: 2,081 divergent (core<0.20) pairs → 2,075 Compara-labelable → 781 (38%) paralog / 1,294 non-paralog.** Balanced, independent (Compara gene-trees), clean classes (pos: ABCA8/10, ACOT2/6, ACSL3/4, AMPD2/3, BLK/HCK; neg: ABCA13/ITSN1, ALB/AVEN — prefilter false-groupings).
- ⚠️ CAVEAT: this set is from ACCEPTED edges (core 0.13–0.20, MILDLY divergent). The DEEP regime (core<0.13, the DSFAM45=0.036 case) is the REJECTED population, not yet generated — a feature calibrated here must be re-checked on the deep tail (Phase 1.5: generate core<0.13 pairs + label).

VERDICT: enough independent truth to calibrate an RNA-only feature → proceed to Phase 1 (feature AUC on this set).

## 11. PHASE 1 RESULT — RNA-only detector is a NO-GO (best AUC 0.629)

`divergent_phase1_features.py` on the gold set (781 paralog / 1294 non-paralog divergent pairs):

| feature class | best AUC |
|---|---|
| DNA k-mer (raw/jaccard/containment/span_cov) | 0.558 |
| + low-complexity MASKING | 0.558 (masking removed ~0 k-mers → false sharing is REAL domain seq, not repeats) |
| + EXON STRUCTURE (exon count/ratio) | 0.556 |
| **protein-ORF (RNA-translated)** | **0.629** ← the only real signal |

VERDICT for the user's RNA-only / no-protein goal = **NO-GO**:
- DNA + structure features (cleanly independent of Compara) FAIL: AUC ≤ 0.56 ≈ random. Divergent paralogs and
  prefilter domain-sharers are indistinguishable by RNA nucleotide sequence.
- The ONLY signal is PROTEIN-level (0.629) — confirming divergent paralogs conserve protein under synonymous
  DNA divergence. BUT: (a) still < 0.70; (b) violates "no protein"; (c) measurement is CONFOUNDED — crude ORF
  deflates it, while Compara's protein-tree truth inflates it (semi-circular). Net: protein is the right axis
  but unresolved whether a proper impl crosses 0.70, and resolving it leaves the RNA-only frame.
- => the contiguous-core recent-duplicate boundary is the DEFENSIBLE detection limit on RNA sequence alone.
  Rigorously confirms the prior "defensible = recent-duplicate detection only" with an INDEPENDENT gold truth.

The ONE door ajar: if "no protein" is relaxed to allow RNA-derived ORF translation, a proper protein analysis
(real ORF + alignment coverage, de-confounded) is one more test — but it is a different (protein) track.

---

## MILESTONE_reference_absent_copies

# Milestone — Reference-absent gene-family copies from gorilla IsoSeq

**One line:** expressed, endogenous, divergent gene-family *second haplotypes* that are **absent from the
T2T gorilla reference** at their locus, found at the RNA level — strongest (protein-confirmed) in the MHC
and supported genome-wide in the recombination (PRDM9) and zinc-finger families. The RNA-level analog of
Eichler's novel TBC1D3 groups and Soto's paralogs missing from GRCh38. **Caveat up front:** from RNA
alone "reference-absent *copy*" vs "hyperdivergent *allele*" is not separable — these are strong
*candidates* pending DNA copy-number (parCN); see Honest caveats.

## The question (advisor's standing interest)

"Find copies of a gene that are *not in the genome*." A copy expressed in the individual but absent
from the reference assembly can never align correctly, so an alignment-keyed method is blind to it. The
question: can we recover such copies from gorilla IsoSeq against the (near-complete) T2T reference?

## Approach — two populations, and where they actually live

A copy absent from the reference surfaces in one of two read piles, and we tested both:

1. **Unmapped reads (divergent-absent).** Tested directly (`unmapped_rescue_poc.md`): the unmapped pile
   is tiny (0.13 %), 79 % already-present (99.7 % identity), residual = expressed transposable elements,
   leaving **one** genuine hit (`cl0`). On a *complete* T2T reference this route is **dry**. Its own
   conclusion told us where to look: a divergent/CNV copy with residual homology gets **force-mapped
   onto its paralog**, not flagged unmapped — so the gene-family copies live in the **mapped** pile.

2. **Mapped reads (collapsed / CNV-absent) — the yield.** A reference copy that is really ≥2 collapsed
   copies (or carries an extra CNV paralog the reference lacks) shows, among its **primary** alignments
   (the paralog-bleed firewall), a coherent **second haplotype**: ≥12 balanced-frequency (0.20–0.60)
   alt columns co-segregating on ≥5 reads (`hidden_copy.rs::detect_hidden_copy`). A diploid het is 1–2
   columns; sequencing error never reaches the balanced band.

## The pipeline (detect → validate → promote → lock)

Run directly on the copy loci (19 s for 357 family loci; the de-novo + POA path was 100× slower and was
abandoned). Each flag is then put through a funnel that each removes a distinct false-positive class:

- **RNA-editing filter** — A→I (read as A>G) makes co-segregating columns; a genuine copy has a
  *diverse* substitution spectrum, editing is A>G/T>C-dominated. (Load-bearing: removed 5/30 family
  flags, incl. a 67 %-A>G GSTM signal.)
- **Assembly + genome-divergence** — POA-assemble the alt-haplotype reads into the copy's consensus,
  map to the whole T2T genome. The decisive signal is **divergence in a band**: the consensus maps back
  to its own locus at **3–20 % divergence** — too divergent to be a het allele (<0.5 %) or the copy
  itself, too conserved to be a chimera (>20 %). (The >20 % band is consensus chimera/repeat — the
  DSFAM213 failure mode.)
- **Protein BLAST** — blastx vs the 22,614-protein gorilla proteome: names the paralog and **rules out
  contamination** (endogenous gorilla protein hit, e ≤ 1e-57 — the check `cl0` failed and these pass).
- **Multi-mapping** (genome-wide) — a real copy of a *multi-locus* family maps to its own divergent
  locus **plus** the reference paralogs; an in-place second haplotype maps to one locus only.

## Result

| tier | n | what | confidence |
|---|---|---|---|
| **MHC copy candidates (protein-confirmed)** | **4** | class I (Gogo-B\*0201), DRB1, DQ-α, DQ-β — 3.9–20.4 % divergent, ORFs 97–362 aa, endogenous (86–98 % gorilla-protein), clustered in the MHC (NC_073229.2 49–50 Mb) | high protein homology; copy-vs-allele needs DNA |
| **Dispersed-paralog copies** | **15** | genome-wide; map to own divergent locus + paralog loci — **PRDM9** (13.5 %, distinct member of the recombination-hotspot family), **ZNF208**, ZDHHC13, large multi-locus families | medium |
| single-locus / family-CNV | 66 | het / tandem-CNV / absent-retrocopy / depth-doubled — RNA cannot separate | needs DNA parCN |

**Funnel:** 12,793 expressed loci → 1,015 flags (8 %) → 73 paralog candidates → 15 dispersed-paralog
copies; plus the 4 family-scoped MHC. The copies land **exactly where biology predicts** — immune
(MHC), recombination (PRDM9), and zinc-finger / segmental-duplication families — the most polymorphic,
copy-number-variable, reference-divergent regions of the genome.

## Why it matters (thesis through-line)

This is the **identifiability thesis turned outward**: the reference cannot *model* these copies (one
gene where the reads imply more), and we detect that from the reads while **abstaining from placing the
copy** (DAZ3/Canzar discipline — detect + flag the evidence, never fabricate a sequence). It complements
the within-reference copy-assignment work: there we assign reads to the copies the reference has; here we
flag copies the reference lacks.

## Measured specificity — false-positive BOUND (`bench/o4_fp_bound.py`)

Quantified (no BAM re-scan; intersect the genome-wide flags with single-copy genes). On **17,334
single-copy genes** (unique-name protein_coding, non-LOC — definitionally cannot host a reference-absent
*copy*), the raw hidden-copy flag still fires on **828/11,206 = 7.39%** of overlapping loci — essentially
the genome-wide background rate (7.93%). **The raw flag is a non-specific SCREEN, not a specific copy
detector**: ~7.4% is an upper bound on its false-positive rate, dominated by **heterozygosity** (a
diploid gene with enough balanced SNP columns trips the ≥12-column threshold). This is the het-vs-copy
wall, measured. *Implication:* the 1,015 genome-wide flags are a candidate screen, not a copy catalog;
specificity comes from the DOWNSTREAM filters (divergence >3% ≫ het ~0.5%; RNA-editing filter; protein
BLAST), and only the **4 protein-BLAST-confirmed MHC candidates** survive an external check. Do not
assert the flagged pool as a copy set.

## Honest caveats (the loose ends)

- **Het vs copy (the central one) — now quantified at ~7.4% raw-flag FP (above).** For the
  hyperpolymorphic MHC and the 66 single-locus candidates,
  "reference-absent paralog *copy*" vs "hyperdivergent *allele*" is genuinely blurry from RNA alone. The
  resolver is **DNA copy-number (parCN)** — the concrete experimental ask. The finding holds regardless:
  expressed sequence 4–20 % divergent from the reference's gene, endogenous and coding.
- **No independent validation yet.** No DNA, no second individual, no assembled full copy locus. The 4
  MHC are protein-confirmed *endogenous* but not *cross-validated* as copies.
- **Genome-wide is a candidate-generator, not a clean catalog.** Precision drops off the
  family/segdup-scoped regions; the 15 dispersed-paralog are multi-mapping-supported candidates, the 66
  single-locus need DNA.
- **Intronless/retrocopy detection is not possible here** — an absent retrocopy's reads map to the
  parent and splice out its introns, so they look like a normal spliced transcript; the intronless
  signature only appears where the retrocopy's own locus is in the reference (absent by definition).

## Reproduce

`hidden_copy_scan.py` · `hidden_copy_validate.py` · `promote_hidden_copies.py` ·
`promote_evidence_plus.py` · `hidden_copy_scan_genomewide.py` · `promote_genomewide.py`. Consensus
sequences + evidence: `/home/juanfra/winloci_scratch/refabsent/`. Full technical writeup:
`bench/reference_absent_catalog.md`.


---

## UNANNOTATED_COPIES

# Can we find UNANNOTATED copies (not even a LOC)? — yes, at three levels of "unknown"

The advisor's question: not reference-absent copies (those mapped to *annotated* MHC/etc. loci), but
expressed multi-copy-family members at loci with **no annotation at all** — not even an uncharacterised
`LOC######`. Checkable by intersecting the de-novo expressed loci with the *entire* gorilla annotation
(gene/lnc_RNA/pseudogene/mRNA/tRNA/…, 30,191 spans, 1,510 Mb) and asking which unannotated ones are
paralog copies. **Answer: yes.**

## The funnel

- **12,793 de-novo expressed loci → 291 unannotated** (2.3%; overlap *zero* annotation).
- Of the 291: **42 are paralog copies** (transcript maps to ≥2 loci at ≥80% identity). Copy-count
  2–8 loci (31 are 2–4, 11 are 5–8) — **gene-family-like, none >8 = no dispersed-TE contamination**.
- blastx vs the gorilla proteome: **203/291 unannotated transcripts hit a protein** (most "unannotated
  loci" are coding paralogs the annotation simply missed at that spot); 88 have no protein hit.

## Three levels of "unknown"

1. **Unannotated copies of KNOWN families (the robust win).** ~35–36 paralog loci whose paralog *and*
   protein are annotated — a copy of a known gene sitting at a locus with **no annotation, not even a
   LOC**. Concrete: **ROBO2** — a 94%-identity copy family at **5 unannotated loci** (likely
   retrocopies / segmental duplications of the axon-guidance gene); ZDHHC13 paralog (53%); several
   `LOC#####` copies (75–81%). Plus **6 read-conflict-VALIDATED unannotated copies inside the de-tie
   families** — including **DSFAM817**, the difficult-case flagship win: *we assigned reads to a copy
   that isn't even annotated.*
2. **Loci whose paralogs are also unannotated, but the protein family is known.** e.g. **ERVFRD-1**
   (endogenous-retrovirus env / syncytin-2 family, 31% — a retroelement family, not a classic gene),
   and copies of annotated `LOC` genes at unannotated loci.
3. **LITERALLY unknown (the rarest).** **NC_073240.2:23.57 Mb** — a **2-copy** family with **no
   annotation, no protein homology, no LOC** — genuinely uncharacterised. The broader pool of **88
   unannotated loci with no protein hit** are the candidates for this class (novel genes / lncRNA /
   to-be-vetted).

## Honest caveats

- "Unannotated copy of a known family" is the solid, common result (40+ loci, gene-family copy-counts,
  protein-confirmed, 6 read-conflict-validated). This robustly answers "can we find unannotated copies."
- "Literally nothing known" (no annotation *and* no protein *and* no LOC) is **rare** — 1 clear
  multi-copy case here, plus the 88-no-protein-hit pool that needs ORF/coherence vetting (novel gene vs
  lncRNA vs assembly/mapping artifact). Some (ROBO2 ×5 @94%) are likely processed pseudogenes/segdups;
  ERVFRD-1 is a retroelement — both are "copies" but not novel *gene* families.
- The annotation is gorilla RefSeq (fairly complete: 34k genes incl. 7k pseudogenes, 10k lncRNAs), so
  "unannotated" here is a genuine gap, not a sparse-annotation artifact.

## Bottom line for the advisor

**Yes — we find expressed paralog copies at loci the gorilla annotation does not cover at all (no LOC):
42 genome-wide (ROBO2 ×5, ZDHHC13, LOC-copies…) + 6 read-conflict-validated (incl. the DSFAM817 win
copy).** Most are unannotated copies of families whose protein is known elsewhere; a small residue
(1 clear multi-copy + an 88-locus no-protein pool) is the "literally unknown" frontier worth vetting.

## Vetting the 88 no-protein-hit pool (`vet_unannotated_novel.py`)

Separating novel gene families from lncRNA/noise (ORF, exon count, paralog count, 4-mer complexity):
- **4 NOVEL multi-copy gene-FAMILY candidates** (multi-copy + ORF ≥100 aa + spliced + complex, no
  annotation, no gorilla-protein hit): **NC_073236.2:139.89 Mb — a 5-copy family, 183-aa ORF, spliced,
  complexity 0.96** (the flagship); NC_086018.1:20.93 Mb (2 cp, 157 aa); NC_086018.1:30.31 Mb (2 cp,
  142 aa, 3 ex); NC_073240.2:23.57 Mb (2 cp, 117 aa).
- 50 novel **single-copy** gene candidates (ORF + spliced, 1 locus — novel genes, not families).
- 15 **lncRNA** candidates (spliced, no/low ORF). 19 ambiguous. 0 low-complexity/noise (complexity
  filter clean).

**Caveat — "truly novel" needs a cross-species check.** "No protein hit" is vs the *gorilla* proteome
only; these 4 could still be homologs of a protein in another species (= gorilla-unannotated, not
universally novel). The airtight confirmation is blastp/diamond vs SwissProt/nr — the next step. Reads
are modest (3–11/locus); the 183-aa ORF is well above the ~21-aa random expectation, so coding potential
is real; complexity 0.96–1.0 + spliced + multi-copy argue against artifact.

Artifacts: `unannotated_loci.json` · `unannotated_copies.json` · `unannotated_novel_vetted.json` ·
`unannot_blastx.tsv` · `unannot_tx.{fa,paf}` (in `/home/juanfra/winloci_scratch/refabsent/`).

