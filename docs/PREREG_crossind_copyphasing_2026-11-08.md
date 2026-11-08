# PREREG: does copy_phasing certify the cross-individual strong loci? (register 1582, 2026-11-08)

Follows `PREREG_crossind_missing_copies_2026-11-07.md` (rows 1578-1579). That test left 12 "strong" loci in OR6737 (an OR-only
three-haplotype flag plus a 1 % divergence pile; `/anvme/workspace/v138da13-rustle/crossind/strong.tsv`). The three-haplotype
test is a crude Strainy step (no refinement, no statistical column test, no reference-copy slots). copy_phasing
(`src/copy_phasing/`, Strainy's phasing + IsoCon's test + the allele-slot certificate) is the full version. Question: does it
certify an extra copy at these loci in OR6737, and not at the same loci in Jim (KB3781), whose own genome the reference is?

## Disclosure (written before the lookup)
- **No new run.** The products already exist: `~/o2ties_real/{or,kb}.A.*`. That is job 4403643 (2026-11-01, register 1483, arm A),
  `copy_phasing.port run --semi --guard` over all OR (`or.copies.tsv`, 446 families) and KB (`kb.copies.tsv`, 280 families)
  copies on GGO.fa. The BAMs are the same as rows 1578-1579. `--semi`/`--guard` act only in the O2 loop after certification
  (`src/copy_phasing/family.rs:410-426`), so the certificate is the default one.
- **What was already seen.** The first five family lines of `or.A.log` and `kb.A.log` (MCL0/1/10/100/101 certified counts), and
  register 1483's O2 numbers. Nothing at the strong loci.
- **Depth confound.** The KB run used the full BAM (13.5 M primaries), not the 0.35 subsample matched to OR (4.4 M). More depth
  gives more power to certify, so the KB null is not depth-matched. It is conservative for "OR > KB".
- **Catalogs differ.** The OR and KB catalogs are separate de novo catalogs with unrelated MCL ids. Loci are matched by coordinates
  only.
- **Slots.** Every copy gets 2 slots. The one strong locus on chrX (NC_086017.1, APOOL/SATL1) has no family copy expected, so it
  should be untestable anyway. If it is testable, it is reported but excluded from the verdict.

## Definitions
- **Locus window.** A `strong.tsv` row's contig:start-end.
- **Testable in a sample.** At least one copy in that sample's `copies.tsv` overlaps the window. A locus with no copy in KB is
  "untestable on the null". That is listed, not counted as a pass.
- **Hit.** A `flagged=1` row of `candidates.tsv` or any row of `putative.tsv` whose home copy (`nearest_locus` / `copy`, resolved
  through `copies.tsv`) overlaps the window. Counting putative copies follows Theorem A: a real extra copy can be held only as
  putative when a homozygous survivor has a free slot. Certified-only counts are reported beside the primary counts.
- **Background.** The per-copy rate over the whole sample's catalog: the fraction of `copies.tsv` copies that are home to at least
  one hit. Reported, not barred.

## Bars (verdict on the OR-testable loci among the 12)
- **PASS (copy_phasing corroborates OR-specific extra copies)** requires all three:
  1. OR hits at >= 4 testable loci, or >= half of them if fewer than 8 are testable.
  2. KB hits <= 1 at the same loci (among KB-testable ones).
  3. Per copy, the OR hit rate inside the windows is >= 2x the OR background.
- **FAMILY COMPLEXITY.** Bar 1 passes and KB hits >= 2. The certificate fires on these families in Jim too, so a certification there
  does not separate a missing copy from CNV-family structure.
- **NOT CORROBORATED.** Bar 1 fails.

## Scorer
`bench/crossind/copyphasing_lookup.py` reads the two candidates/putative/copies tables and strong.tsv and prints one row per
locus plus the verdict. Report only. No default changes.

## Amendment (after the lookup, disclosed)
The Slots note above names NC_086017.1 as chrX. It is chromosome 2 (the ZNF717 / MCL197 locus). ChrX is NC_073247.2 (the APOOL/SATL1 locus,
as the note said), which has no family copy. The first scorer run excluded chr2 by mistake: OR testable 6, hits 5, KB 0, per-copy
10/18. The corrected run excludes NC_073247.2, as the note intended. The verdict is the same.

## Result (2026-11-08): PASS as registered
Scorer output: `/anvme/workspace/v138da13-rustle/crossind/copyphasing_lookup.txt`.

| locus | OR copies (families) | OR certified / putative | KB copies (families) | KB certified / putative |
|---|---|---|---|---|
| NBPF, NC_073224.2:117.29-117.53 Mb | 2 (MCL5) | 1 / 0 | 3 | 0 / 0 |
| IGHV-like, NC_073242.2:18.35-18.60 Mb | 6 (MCL3, MCL67) | 10 / 2 | 4 | 0 / 0 |
| NACA-like, NC_073242.2:22.13 Mb | 1 (MCL7) | 1 / 0 | 1 | 0 / 0 |
| INTS4-like/CCT6A, NC_073230.2:74.82-74.95 Mb | 2 (MCL35, MCL271) | 2 / 0 | 1 | 0 / 0 |
| SMPD4-like, NC_073235.2:31.68 Mb | 1 (MCL310) | 0 / 0 | 2 | 0 / 0 |
| NBPF, NC_073224.2:223.16-223.94 Mb | 6 (MCL5, MCL224) | 0 / 1 | 3 | 0 / 0 |
| ZNF717, NC_086017.1:85.70-85.82 Mb | 1 (MCL197) | 0 / 0 | 1 | 0 / 0 |
| 5 loci with no family copy (incl. chrX NC_073247.2) | 0 | untestable | 0 | untestable |

- **Bar 1 PASS.** OR hits at 5 of 7 testable loci (certified at 4), against a need of 4.
- **Bar 2 PASS.** KB has 0 hits at the same 7 loci, even at about 3x OR's depth.
- **Bar 3 PASS.** The OR per-copy hit rate in the windows is 10/19 = 0.53, against a background of 76/1,468 = 0.052 (10x).
- **Background.** KB's genome-wide per-copy rate is 53/812 = 0.065. These are certifications on Jim's own genome, so they are false by
  construction, or collapsed assembly copies. The certificate does fire on Jim elsewhere; it just does not fire at these loci.

**Reading (caveats that limit the PASS):**
1. **The same reads chose the loci and certify them.** The loci were picked as OR-only three-haplotype loci, and copy_phasing
   phases the same OR reads. So this corroborates the method, not the data. Bar 2 is partly built into the selection too: a KB
   three-haplotype flag at a locus excluded it. Independent evidence would have to come from OR6737 DNA, which we do not have.
2. **The IGHV-like locus (10 of the certified copies) has a non-copy explanation.** B-cell transcripts carry V(D)J rearrangement
   and somatic hypermutation, and both make extra haplotypes. Without that locus: 4 of 6 hits, need 3, and the result still passes.
3. What survives is NBPF, NACA-like and INTS4-like/CCT6A: certified extra haplotypes in OR6737, none in Jim. NBPF is the textbook
   copy-number-variable family; NACA and INTS4 have known paralogs or pseudogenes. This fits row 1579: a different individual's missing copies are few and sit in CNV families.
