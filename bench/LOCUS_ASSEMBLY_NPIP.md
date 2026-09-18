# Locus assembly at NPIP — measured against `docs/PREREG_locus_assembly_2026-09-18.md` (md5 02237f9d6da31ec694a8e53cd8cbdf2b)

**Verdict: the bottleneck is NOT ambiguity and NOT depth — it is EXACT-CHAIN COLLAPSE in the assembler.**
O2 ran, assigned 206 contested molecules, and changed the reconstruction by **exactly zero junctions**.
The assembler emits 20 points FEWER junctions than simply taking the primary reads' junctions.

Data `/mnt/linuxdisk/home/juanfraitu/locus_asm/`. No new algorithm was written; `copy_assign` was run as
shipped (rebuilt `--release`).

## Arms

| arm | command | families found |
|---|---|---|
| **A** annotation-free | `--regions npip.regions --gtf --read-provenance` (27 per-copy windows) | **0** |
| **C** copies supplied (O2 isolated) | `--regions npipC.regions --families npipC.copies.tsv --copies-fa npipC.copies.fa --gtf --read-provenance` (chr16:11.9–19.0 Mb, 9 co-located copies) | 1 |

**Arm A found 0 families.** NPIP is a DISPERSED family — 26 chr16 copies from 11.96 to 80.44 Mb, median gap
between consecutive copies **265,515 bp**, max 45.2 Mb. `copy_assign` detects CO-LOCATED copies, so with one
copy per window it had nothing to pair. This is a scope fact about the binary, not a failure of assembly:
it still assembled 750 transcripts (653 spliced).

## AS-1 — integrity: **PASSED** (nothing exceeds the §6m4 ceiling)

## The 9-copy co-located cluster, 106 annotated junctions

| read set / arm | junctions | complete chains | exact chains |
|---|---|---|---|
| primary-read union (baseline to beat) | 85 (**80.2%**) | 2/9 | — |
| ALL-alignment union (**ceiling**, perfect O2) | 91 (**85.8%**) | 4/9 | — |
| **arm A — assembler, no O2** | **64 (60.4%)** | **0/9** | **0/9** |
| **arm C — assembler + O2** | **64 (60.4%)** | **0/9** | **0/9** |

**AS-2 FAILS on both bars**: 60.4% < 75.5% required (a); 0/9 < 7/26-equivalent (b).
**Arm C is byte-identical to arm A on every copy.** O2 is not the constraint.

Whole-family arm A (26 copies, 249 junctions): union 149 (59.8%), exact chains 3/26 — same picture.

## AS-5 — O2 did run, and abstained properly

| | |
|---|---|
| molecules in the swept region | 23,884 |
| AS-tied, entered the certificate | 1,456 (9,646 records) |
| evaluated | 955 — origin-rejected 294 (O3's), **CONTESTED 661** |
| **assigned** | **206 (31.2%)** |
| tied | 176 (26.6%) |
| ambiguous (abstained) | 279 (42.2%) |

Assign-or-abstain works: 69% of contested molecules were refused. But 206 assignments out of 23,884
molecules cannot add an isoform, and the chains they support were already carried by unique reads.

## Why the assembler loses 20 junctions it can see

Of the 85 annotated junctions visible in primary reads, the assembler emits 65 and **drops 20**. It is not
a depth floor — 10 of the 20 have depth ≥ 10, five of them **255–287 reads**:

| | n | median depth | median #distinct intron chains | **median largest single chain** |
|---|---|---|---|---|
| emitted | 65 | 126 | 63 | **32 reads** |
| **dropped** | 20 | 12 | 11 | **2 reads** |

The five high-depth casualties:

| junction depth | distinct chains carrying it | largest single chain |
|---|---|---|
| 287 | 131 | 42 |
| 275 | 134 | 42 |
| 274 | 132 | 42 |
| 270 | 129 | 42 |
| 255 | 117 | 42 |

**The mechanism.** The assembler collapses reads by EXACT intron chain and then applies a per-chain floor.
A junction carried by 287 reads spread over 131 different chains produces no chain big enough to survive,
so the junction disappears — while a junction on one dominant chain sails through. Dropped junctions have a
median largest-chain of **2 reads** against **32** for kept ones: a 16× difference, and it is the whole story.

This is §6m0's fragmentation finding in the assembler rather than in `shared_definition`: 5′-truncated
Iso-Seq reads shatter one transcript into many distinct chains. `shared_definition` already ships the fix —
**read-isoform widening** (`widen_with_read_isoforms`, k = 5), which admits a chain when every junction has
≥ k support instead of demanding one exact chain clear a floor. `copy_assign`'s assembler does not have it.

## Conclusions

1. **O2 is not the bottleneck for reconstruction.** It contributed 0 junctions here. Its value is copy
   ASSIGNMENT (206 resolved, 455 correctly refused), which is a different deliverable.
2. **The bottleneck is exact-chain collapse**, and the fix already exists in the codebase.
3. **NPIP is dispersed**, so annotation-free `copy_assign` finds no family at all on per-copy windows.
   Reconstructing a dispersed family needs the copy set supplied (O1's job) or a non-co-located detector.
4. Even with all three fixed, §6m4's ceiling stands: **≤ 4/9 complete chains in this cluster, ≤ 10/26
   family-wide**, because 16.1% of annotated junctions have no read in this library at all.

**Next lever, in order of measured value:** port read-isoform widening into the assembler (worth up to
+21 junctions, 60.4% → 80.2%), then the all-alignment pool (+6, → 85.8%). O2 changes neither.
