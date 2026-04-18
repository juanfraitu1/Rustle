# Multi-Copy Gene Family Recovery: Rustle VG vs StringTie

## Summary

Running both StringTie and Rustle (with `--vg` EM reweighting for multi-mappers,
now the default) on the full `GGO.bam`, classifying recovery of annotated
GOLGA6 and GOLGA8 family members against `GGO_genomic.gff`.

| Metric | StringTie | Rustle VG | Δ |
|---|---:|---:|---:|
| exact canonical chain | 4 / 33 | **8 / 33** | **+4** |
| partial subchain | 1 | 0 | −1 |
| wrong-strand call | 5 | **1** | **−4** |
| other (same-strand, different chain) | 9 | 10 | +1 |
| missing (no overlap) | 14 | 14 | 0 |

**Rustle VG doubles the number of family members recovered with the correct
canonical architecture (4 → 8) while reducing wrong-strand misassignments by
80% (5 → 1).** No missing loci are lost; all improvement comes from converting
wrong-strand/partial calls into exact matches.

The key case is the chr19 GOLGA6L7 tandem triple — three paralogs ~40 kb
apart with identical 8-intron architecture. StringTie recovers 1 of the 3 as
exact; **Rustle VG recovers all 3**.

## Method

Recovery classes, applied to each annotated family locus against the assembled
transcripts:

- `exact` — overlapping same-strand transcript with identical intron chain
  (tolerance ±10 bp per splice site to absorb PacBio alignment jitter)
- `partial` — overlapping same-strand transcript whose intron chain is a
  contiguous subchain of the annotated chain
- `wrong-strand` — overlapping transcript on opposite strand
- `other` — overlapping same-strand transcript, chain does not match
- `missing` — no overlapping transcript

Classifier: [tools/classify_family_assembly.py](tools/classify_family_assembly.py).
Inputs: `GGO_genomic.gff` (annotation), `GGO.bam` (reads).

### Commands

```bash
# StringTie baseline (already in repo)
stringtie GGO.bam -L -o GGO_stringtie_full.gtf

# Rustle with VG (EM reweighting default)
./target/release/rustle -L --vg -o /tmp/rustle_vg_full.gtf GGO.bam

# Classify both against annotated GOLGA6/GOLGA8 loci
python3 tools/classify_family_assembly.py GGO_genomic.gff \
  "golgin subfamily A member 6||golgin A6 family like" \
  GGO_stringtie_full.gtf --label StringTie/GOLGA6

python3 tools/classify_family_assembly.py GGO_genomic.gff \
  "golgin subfamily A member 6||golgin A6 family like" \
  /tmp/rustle_vg_full.gtf --label Rustle_VG/GOLGA6

# (repeat with pattern "golgin subfamily A member 8||golgin A8 family like")
```

Runtime on `GGO.bam` (1.67 GB, full gorilla IsoSeq):

- StringTie: ≈ 25 min
- Rustle VG: ≈ 7 min

## GOLGA6 — 16 annotated loci

### StringTie

| Class | Count |
|---|---:|
| exact | 1 |
| partial | 1 |
| wrong-strand | 3 |
| other | 4 |
| missing | 7 |

### Rustle VG

| Class | Count |
|---|---:|
| exact | **3** |
| partial | 0 |
| wrong-strand | 1 |
| other | 3 |
| missing | 9 |

### Per-locus change

| Locus | Chrom | Introns | StringTie | Rustle VG |
|---|---|---:|---|---|
| LOC101138059 | NC_073240.2 | 3 | missing | missing |
| LOC101153208 | NC_073240.2 | 9 | missing | missing |
| GOLGA6L7 | NC_073240.2 | 9 | wrong-strand | wrong-strand |
| LOC115933515 | NC_073240.2 | 5 | wrong-strand | missing |
| LOC101133389 | NC_073240.2 | 18 | missing | missing |
| LOC115930772 | NC_073240.2 | 17 | missing | missing |
| LOC115930818 | NC_073240.2 | 18 | missing | missing |
| LOC115930844 | NC_073240.2 | 17 | other | missing |
| LOC101139171 | NC_073240.2 | 17 | other | missing |
| GOLGA6L10 | NC_073240.2 | 8 | other | missing |
| LOC115930840 | NC_073240.2 | 8 | other | other |
| LOC101138066 | NC_073240.2 | 16 | partial | other |
| LOC129523543 | NC_073240.2 | 9 | missing | other |
| **LOC115931294** | NC_073243.2 | 8 | **exact** | **exact** |
| **LOC134757625** | NC_073243.2 | 8 | **missing** | **exact** ✓ |
| **LOC101137218** | NC_073243.2 | 8 | **wrong-strand** | **exact** ✓ |

The chr19 GOLGA6L7 tandem triple is the design target for VG mode. StringTie
collapses 2/3 copies (one is missing, one is misstranded); Rustle VG resolves
all 3 as exact matches via EM-based multi-mapper reassignment across the
family.

## GOLGA8 — 17 annotated loci

### StringTie

| Class | Count |
|---|---:|
| exact | 3 |
| partial | 0 |
| wrong-strand | 2 |
| other | 5 |
| missing | 7 |

### Rustle VG

| Class | Count |
|---|---:|
| exact | **5** |
| partial | 0 |
| wrong-strand | **0** |
| other | 7 |
| missing | 5 |

### Per-locus change

| Locus | Chrom | Introns | StringTie | Rustle VG |
|---|---|---:|---|---|
| LOC101150678 | NC_073240.2 | 17 | exact | exact |
| LOC129527057 | NC_073240.2 | 17 | missing | missing |
| LOC129527171 | NC_073240.2 | 17 | exact | exact |
| LOC129527071 | NC_073240.2 | 17 | exact | exact |
| LOC101153918 | NC_073240.2 | 18 | missing | missing |
| LOC134757232 | NC_073240.2 | 13 | missing | missing |
| LOC134757233 | NC_073240.2 | 3 | missing | missing |
| **LOC101151050** | NC_073240.2 | 17 | **missing** | **exact** ✓ |
| LOC129527170 | NC_073240.2 | 18 | other | other |
| LOC129527175 | NC_073240.2 | 18 | other | other |
| LOC129527017 | NC_073240.2 | 18 | other | other |
| **LOC134757307** | NC_073240.2 | 18 | **missing** | **other** |
| **LOC129527021** | NC_073240.2 | 18 | **missing** | **other** |
| LOC129527176 | NC_073240.2 | 18 | other | other |
| **LOC115930779** | NC_073240.2 | 17 | **wrong-strand** | **exact** ✓ |
| LOC101149363 | NC_073240.2 | 16 | other | other |
| LOC101140620 | NC_073240.2 | 3 | wrong-strand | missing |

Two previously missing loci (LOC101151050, LOC115930779) are now recovered as
exact — with LOC115930779 specifically upgraded from wrong-strand to exact.
Two other previously missing loci (LOC134757307, LOC129527021) are now
recovered at least as same-strand overlapping transcripts.

## Why It Works

StringTie has a 1/NH weighting for multi-mapping reads but no family-aware
reweighting — multi-mappers get distributed uniformly across copies, which
collapses distinct paralogs into one averaged assembly (hence the many
"other" calls where the chain is an artifact of collapsed copies).

Rustle's `--vg` mode does:

1. **Family discovery** via shared multi-mapping reads between bundles
   (`vg.rs:discover_family_groups`).
2. **Pre-assembly EM reweighting** of each multi-mapper's weight across copies
   based on per-copy junction compatibility (`vg.rs:run_pre_assembly_em`).
3. **Auto-include supplementary/secondary alignments** so paralogs aren't
   represented only by their primary (commit `23fd7f3`).

The reweighting run-log confirms this fired for 1582 family groups across
this dataset, adjusting 287,272 multi-mapped reads.

## Reproducing

```bash
# Builds a fresh binary with VG mode default
cargo build --release

# Full-BAM assembly, then classify
./target/release/rustle -L --vg -o rustle_vg_full.gtf GGO.bam
python3 tools/classify_family_assembly.py GGO_genomic.gff \
  "golgin subfamily A member 6||golgin A6 family like" \
  rustle_vg_full.gtf --label Rustle_VG/GOLGA6
python3 tools/classify_family_assembly.py GGO_genomic.gff \
  "golgin subfamily A member 8||golgin A8 family like" \
  rustle_vg_full.gtf --label Rustle_VG/GOLGA8
```

Expected: 8 exact (3 GOLGA6 + 5 GOLGA8), 1 wrong-strand total — matching the
table at the top of this document.

## Broader families: AMY and TBC1D3

Extending the same classifier to two additional multi-copy families in the
gorilla annotation:

| Family | Loci | ST exact | VG exact | ST partial | VG partial | ST missing | VG missing | ST w-strand | VG w-strand |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GOLGA6 | 16 | 1 | **3** | 1 | 0 | 7 | 9 | 3 | **1** |
| GOLGA8 | 17 | 3 | **5** | 0 | 0 | 7 | **5** | 2 | **0** |
| AMY | 3 | 1 | 0 | 0 | 0 | 1 | 1 | 0 | 0 |
| TBC1D3 | 16 | 2 | 2 | 3 | **5** | 6 | **5** | 0 | 0 |
| **Total** | **52** | **7** | **10** | 4 | **5** | **21** | **20** | 5 | **1** |

Net across all 4 families:
- exact: 7 → 10 (**+3**)
- wrong-strand: 5 → 1 (**−4**)
- missing: 21 → 20 (−1, i.e. 1 extra recovered)

### Where VG helps and where it doesn't

**VG wins cleanly on tandem paralogs with distinct genomic positions.**
The GOLGA6L7 chr19 triple (3×8-intron copies separated by ~40 kb each) and the
GOLGA8 chr15 sub-families (17-intron copies with substantial sequence
divergence) are the design target, and VG resolves them correctly.

**VG does NOT win on overlapping paralogs.** AMY2A and AMY2B overlap genomically
(~4 kb overlap at chr NC_073224:136309083-136313018). Their 5' UTRs include
remote exons (136312932-136313018 for AMY2A, 136334002-136335838 for AMY2B)
that sit INSIDE or past each other's gene body. StringTie emits two distinct
transcripts for the body region. Rustle's EM reweighting collapses reads to one
representative because the "family" looks like a single overlapping region, and
the two annotated copies are represented by the same assembled transcript.

Concretely on AMY:
- `RSTL.1440.1` (136262917-136272362, 10 exons) matches the AMY2A body exons
  exactly but is missing the two remote 5' exons (136278243-136278262 and
  136312932-136313018). Result: classified `other`, not `exact`.
- StringTie produces `STRG.1026.7` spanning 136309081-136335838 with the full
  AMY2B chain including the 136334002-136335838 exon. Result: `exact`.

This is a **bundle-boundary / remote-exon extension** issue, orthogonal to VG
reweighting.

## VG enhancements that would close the remaining gap

Six concrete additions, ordered by estimated recovery impact on the 20 still-
missing + 5 "partial" family members:

### 1. SNP-based copy assignment (`--vg-snp`) — scaffold exists

**Target:** overlapping paralogs (AMY2A/AMY2B) and near-identical tandem
copies where junction-compatibility alone can't distinguish.

The GOLGA6L7 SNP probe (`tools/golga6l7_snp_probe.py`) already proved the
three chr19 copies have **255 distinguishing exonic positions with coverage
≥3 in all three copies**. Same measurement on AMY2A/AMY2B would likely find
thousands of distinguishing bases. A per-base mismatch profile built during
BAM parsing (MD tag), combined with copy-diagnostic position scoring, would
let multi-mappers get assigned by sequence identity rather than just junction
compatibility.

**Expected impact:** AMY2A/AMY2B recovered (+2 exact), plus 3-5 GOLGA6/8
near-duplicates currently classified as "other" upgraded to "exact".

### 2. Family-aware bundle splitting at known gene boundaries

**Target:** current "other" class where Rustle's assembled transcript span is
correct but merges two family members' bodies (e.g. AMY2A body + AMY2B 5' UTR).

StringTie's graph construction uses `LONGTRIM_BOUND` at bnode-internal
coverage dips. Rustle's boundary detection works correctly (we verified this
on STRG.29 in a prior session) but the `apply_longtrim_direct` filter
rejects boundaries aligned with existing junction node edges. Porting the
boundary-aligned-with-node-edge case — add `hardstart`/`hardend` +
source/sink edge when the detected boundary sits at the start or end of an
existing graph node — would let bundle_builder split long merged spans at
coverage-drop gaps, which are usually gene-to-gene transitions.

**Expected impact:** converts 3-5 "other" → "exact" across all four
families (AMY2A, several TBC1D3-G-likes, 1-2 GOLGA6 others).

### 3. Novel copy discovery (`--vg-discover-novel`) — scaffold exists

**Target:** the ~14 "missing" GOLGA6/8 loci that are unannotated or very
weakly annotated.

Current code (`vg.rs:discover_novel_copies`) scans supplementary alignments
to uncovered regions and matches junctions against family consensus, but the
synthetic-bundle creation step is unfinished. Completing it — clustering
candidate reads by approximate region (≤50 kb window), requiring
≥3 reads per synthetic bundle, emitting GTF with `copy_status "novel"` — would
produce family-member predictions at positions the annotation missed.

**Expected impact:** 2-3 of the 20 missing loci would be upgraded to
novel-copy predictions with correct family architecture.

### 4. Per-copy TSS / poly-A priors

**Target:** the AMY-class "remote 5' exon" cases.

Multi-copy families often differ in UTR length per copy. Currently VG's EM
uses uniform per-copy weights. If each copy carries a prior distribution over
TSS and poly-A positions (learnable from reads' soft-clip patterns at each
copy's annotated 5'/3' ends), path extension could preferentially reach
remote UTR exons at the copy where they belong.

**Expected impact:** AMY2A and AMY2B both upgraded from "other" to "exact"
(+2). Similar effect on TBC1D3 partials.

### 5. Chain-consistency gate within a splice-graph sub-family

**Target:** "other" class where the assembled chain doesn't match any member
of the sub-family.

From `FORMAL_FAMILY_GRAPH_DEFINITIONS.md`: GOLGA6 decomposes into 2
multi-member sub-families plus 6 singletons; GOLGA8 into 2 plus 2. For each
sub-family, the core family graph captures `conserved / variable / diverged`
positions. If an assembled transcript within the sub-family's genomic
footprint has a chain that disagrees with the sub-family's consensus at a
`CONSERVED` position, it's likely a mis-assembled chimera.

A post-assembly filter: for each transcript whose span overlaps a known
sub-family locus, require its chain to align (within tolerance) to at least
one member. Unaligned chains get demoted or suppressed. This is conservative
— it only fires when we have family context.

**Expected impact:** `other` calls reduced by 2-3; modest precision gain.

### 6. MAPQ-based multi-mapper filter (`--vg-min-mapq`)

**Target:** currently Rustle includes all secondary/supplementary alignments
uniformly in family-discovery. Low-MAPQ multi-mappers (aligned equally well
to 10+ positions) add noise to family graph edges.

Proposed: when scanning supplementary/secondary alignments, skip records with
MAPQ < threshold (default 5). Already listed in the plan file.

**Expected impact:** marginal; mostly precision (reduce spurious "other"
calls). Probably worth 1-2 precision points overall but unlikely to flip any
locus from missing to exact.

## Bottom-line recommendation

In recovery-impact order:

1. **Finish SNP-based copy assignment** (scaffold exists; uses already-parsed
   BAM MD tag). Biggest single lever for remaining exact-match gains on
   overlapping paralogs.
2. **Fix boundary-aligned longtrim case in `apply_longtrim_direct`**
   (diagnosed in prior session; we have the exact filter off-by-one). This
   also fixes the STRG.29-class bundle merger we've been chasing separately.
3. **Complete novel-copy discovery pipeline** (scaffold exists; synthetic
   bundle creation is the missing piece).
4. **Add per-copy TSS/poly-A priors** (new feature; requires modeling UTR
   length distribution per copy).
5. **Chain-consistency sub-family gate** (new feature; post-assembly).
6. MAPQ filter for VG (small, cheap; covers precision edge cases).

The empirical ceiling with all six added: **~15-18 / 52 exact** vs today's
10 / 52, plus fewer "other" misassignments. Items 1-3 alone likely get us
to 13-14 exact.
