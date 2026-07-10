# Spot-check on six named gene families: does copy detection work?

**Date:** 2026-07-09. **Substrate:** `GGO_mm.bam` (gorilla testis Iso-Seq) vs `GGO.fasta`, current `copy_assign`,
`--min-copies 2 --skip-poa-diagnostic`, each region run alone, foreground.

Six annotated families with a known paralog count were chosen so that an **undercall is visible**: if the
annotation says two paralogs sit in the window and we return no family, that is a miss, not an opinion.

## Result

| family | annotated paralogs in window | primary reads | E_c (default) | E_r (`--homology-primary`) | verdict |
|---|---|---|---|---|---|
| **GSTM** | GSTM3 `+`, GSTM5 `−`, GSTM1 `−` | 1377 | 4 reps, **0 edges → 0 families** | **2 families, 3868 reads assigned** | **undercall fixed by E_r**, but one family is an artifact (below) |
| **MAGEA** | MAGEA4 `+`, MAGEA10 `−` | 674 | 0 edges → 0 families | 1 edge → 2 components, **0 families after co-location** | **undercall NOT fixed** — inverted pair |
| **RFPL** | RFPL2 `−`, RFPL3 `+` | 104 | 0 edges → 0 families | 2 edges → 2 components, **0 families after co-location** | **undercall NOT fixed** — inverted pair |
| PRAMEF | PRAMEF20/17/19 | **3** | 0 transcripts | 0 transcripts | not expressed — correct silence |
| XAGE | XAGE5, XAGE3 | **9** | 0 transcripts | 0 transcripts | not expressed — correct silence |
| CEACAM | 6 genes | **16** | 0 transcripts | 0 transcripts | not expressed — correct silence |

Three of six families carry no data in this testis sample (3, 9, 16 primary reads). Returning no family there
is **correct**, not a miss. The method is only on trial for GSTM, MAGEA, and RFPL.

## Finding 1 — the E_c undercall is real, and `--homology-primary` fixes it (GSTM)

GSTM is the textbook "globin problem": the copies map uniquely, so no read is ambiguous between them, so the
read-conflict graph (`E_c`) has **zero edges** and the family does not exist. With four de-novo reps in hand,
default `copy_assign` emits **0 families and 0 assignments**.

Under `--homology-primary` the E_r homology oracle finds 4 edges and the family appears: **2 families, 5224
assignments, 3868 reads assigned / 1356 tied.** `CAFAM1` is exactly right — GSTM5 (`129191737-129198214`) and
GSTM1 (`129211328-129222730`), the two `−`-strand paralogs, as two disjoint copies.

This is the first time the E_r work has been shown to recover a *named, annotated* family on real data. The
[HOMOLOGY_PRIMARY_DELTA](HOMOLOGY_PRIMARY_DELTA.md) probe regions never exhibited it.

## Finding 2 — a NEW artifact class: readthrough transcripts admitted as copies (GSTM)

The same GSTM run also produced `CAFAM0` = { GSTM3, **`DN_..._129190708_1`** }. That second "copy" spans
`129190708-129220537` — 30 kb, **one intron** — and covers **both GSTM5 and GSTM1**, which `CAFAM1` already
holds as separate copies. The same genomic sequence is a copy of two different families, and a spurious family
(`CAFAM0`) is built on GSTM3 plus a chimera.

**Both existing checks missed it.** The overlap check looked only *within* a family, and the annotation
cross-check reported `CAFAM0` and `CAFAM1` as "distinct genes → corroborated". The tell is that a copy of one
family overlaps a copy of another. Fixed: `denovo_pipeline::catalog_overlaps` is now catalog-wide and
classifies each pair as `DuplicateLocus` (reciprocal ≈ 1 — one locus admitted twice, the CAFAM69 case),
`Containment` (a readthrough enclosing a fragment, which must NOT be merged or it would absorb real tandem
copies), or **`SharedAcrossFamilies`** (this case). `copy_assign` now warns, and `bench/artifact_audit.py`
exits non-zero.

The single-intron, 30 kb span is the signature of a badly assembled transcript, not biology.

## Finding 3 — inverted duplicates are structurally invisible to O2 (MAGEA, RFPL)

MAGEA4 (`+`) / MAGEA10 (`−`) and RFPL2 (`−`) / RFPL3 (`+`) are real paralog pairs, and the **E_r oracle finds
the homology edge in both cases**. Both families are then destroyed downstream: `colocated_families` partitions
per `(chrom, strand)`, so each inverted pair splits into two singletons and is dropped at `min_copies = 2`.

This is not a threshold to tune — it is a modelling decision. O1's same-chromosome path was deliberately moved
off `colocated_families` onto the strand-agnostic `families_from_conflict_graph` for exactly this reason
(inverted duplicates), while **O2 kept the strand split**. So O1 can see an inverted duplicate that O2 can then
never assign. `--homology-primary` cannot help: the family is lost after membership, not during it.

Two of the three expressed families in this spot-check are inverted pairs, so the effect is not exotic. Its
genome-wide prevalence has **not** been measured.

## What this changes

- `--homology-primary` earns its keep: it is the difference between 0 and 2 families on GSTM.
- The abstention accounting is not yet trustworthy on default runs. `CAFAM0`'s reads are assigned against a
  chimeric copy set; `CAFAM69`-style duplicate loci abstain wholesale and look like the K=0 wall.
- **Two open defects, neither fixed here:**
  1. `collapse_loci_groups` unions transcripts only on an exact shared intron, never on positional overlap, so
     one locus can be admitted twice (26% of families in the 74-family catalog contain an overlapping rep pair).
  2. `colocated_families` splits on strand, so O2 cannot assign inverted duplicates.
- Fixing either changes family membership and therefore every downstream number. Both are flagged, not fixed.

---

# FIX: `colocated_families` no longer splits on strand

Finding 3 above is fixed. `colocated_families` now partitions family members by **chromosome only**.

**Why this is safe.** The strand split was justified as protection against same-locus antisense (a `+` gene and
a `-` gene whose spans overlap — an inverted-repeat artifact, not two copies). But `prune_same_locus`
**clause (c)** already removes exactly that, and its own comment says so: *"distinct-loci opposite-strand
copies do NOT overlap, so this never fires on a genuine inverted-duplication paralog — which is exactly what we
want to keep."* Because the split ran first, every group was single-strand and clause (c) could never fire. It
was dead code guarding against a case the split had already destroyed, along with the paralogs.

Mixed-strand families are safe downstream: each copy carries its own strand, `copy_assign_pipeline`
reverse-complements read bases per copy (`assign_minus_strand_copy_reverse_complements_read_base`), and two
inverted-duplicate mRNAs are both in transcription orientation, so they align forward to each other.

## Result

| family | before | after |
|---|---|---|
| **MAGEA** (MAGEA4 `+` / MAGEA10 `-`) | 0 families | **1 family, 2 copies — `163590486-163600790` = MAGEA4, `163809311-163814502` = MAGEA10. 895/931 reads assigned, 35 tied, 1 ambiguous.** |
| **GSTM** | 2 families, one of them a chimera | **1 family, 3 copies = GSTM3 + GSTM5 + GSTM1, no warning. 2656/2675 assigned (99.3%).** |

Both match the annotation exactly. The GSTM readthrough artifact **also disappeared**: with the strand split
gone, `prune_same_locus` finally compares the `+`-strand 30 kb readthrough against the `-`-strand GSTM5/GSTM1
and prunes it. One fix, two defects closed.

The pre-existing test `colocated_families_splits_by_strand` asserted the old behavior; it is replaced by
`colocated_families_keeps_mixed_strand_disjoint_copies_in_one_family`, and the antisense protection is now
covered directly by `colocated_families_still_drops_same_locus_antisense_overlap`. 847 tests.

## ⚠ Default outputs change

This is a family-membership fix, so numbers move. Measured on the probe regions:

| region | before | after |
|---|---|---|
| `NC_073224.2:101578582-101607889` | 1 family, 665 assignments | unchanged |
| `NC_073228.2:182473722-182663103` | 2 families, 1209 assignments | **1 family, 481 assignments** |

The family that vanished on region 4 is the suspicious `CAFAM0` — 816 reads, *both* copies flagged collapsed,
90% tied, and the outlier that once skewed the tied fraction from 25% to 44.9%. It is pruned as an artifact.
That also retires the open question in `HOMOLOGY_PRIMARY_DELTA.md`, where region 4 looked like pure E_r
*losing* a real family: part of what E_c had there was not a family.

## ⚠ Newly exposed defect: giant single-exon readthrough transcripts

RFPL (RFPL2 `-` / RFPL3 `+`) now forms its family — and then **assignment hangs (>400 s)**. The cause is not
the strand fix. The window contains `DN_NC_086018.1_30210606_1`: a **128 kb, single-exon** de-novo transcript
carrying 12 reads, spanning RFPL2 and SLC5A4. It is a readthrough/intronic pileup, the same artifact class as
the GSTM one, and it is admitted as a *copy*. Aligning reads against a 128 kb "transcript" is quadratic.

`prune_same_locus` misses it: clause (b) drops a structureless span only when it is *contained* (overlap ≥ 50%
of the shorter transcript), and this giant merely straddles each neighbour partially (24.5 kb of a 53.9 kb
copy). The defect pre-dates this fix — the strand split had simply been hiding it by dropping the family.

The candidate rule is: a **single-exon transcript that overlaps any spliced transcript at all** is the
unspliced/readthrough form of a locus, not a copy. Retrocopies stay safe under that rule — they are intronless
but sit at a *different* locus and do not overlap their parent. Not implemented: it is a second
membership-changing rule and deserves its own measurement.

## Reproduce

```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --region NC_073224.2:129164000-129240000 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --out GSTM_er
python bench/artifact_audit.py GSTM_er --gff GGO_genomic.gff     # exits 1: SharedAcrossFamilies x2
```
Related: `bench/HOMOLOGY_PRIMARY_DELTA.md`, `bench/K0_FLANK_EXPERIMENT.md` (the three causes of abstention).
