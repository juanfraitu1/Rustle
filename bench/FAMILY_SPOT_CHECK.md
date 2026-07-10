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

## Reproduce

```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --region NC_073224.2:129164000-129240000 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --out GSTM_er
python bench/artifact_audit.py GSTM_er --gff GGO_genomic.gff     # exits 1: SharedAcrossFamilies x2
```
Related: `bench/HOMOLOGY_PRIMARY_DELTA.md`, `bench/K0_FLANK_EXPERIMENT.md` (the three causes of abstention).
