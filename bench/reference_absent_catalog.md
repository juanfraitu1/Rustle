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

**The 4 reference-absent copies are all MHC genes**, clustered in the gorilla MHC on NC_073229.2
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

**Final locked result:** **4 endogenous reference-absent divergent MHC paralog copies** (named,
protein-confirmed) + 4 DNA-resolvable CNV candidates (PRAME/MAGE, ZNF766, uS5m, +1). Locking artifacts:
`promote_evidence_plus.py` · `blastx.tsv` · the `proteindb` gorilla-proteome BLAST DB.

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

**FULL CATALOG:** 4 confirmed MHC copies (family-scoped) + **15 dispersed-paralog reference-absent
copies** (genome-wide; PRDM9/ZNF208/…) + 58 single-locus & 8 family-CNV candidates needing DNA. The
copies land — as biology predicts — in immune (MHC), recombination (PRDM9), and zinc-finger/segdup
families. Discriminator artifacts: `gw_promoted/{clean73.paf,gw_discriminated.json}`.
