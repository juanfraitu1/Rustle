# minimap2 `-N50` vs `-N1000` — the secondary-alignment cap is airtight for the shipped pipeline

**Question.** The shipped alignment is `minimap2 -ax splice:hq -uf --eqx -Y -N50 -p0.1 --secondary=yes`
(`GGO_mm.bam`, minimap2 2.31, confirmed by `@PG`). `-N` caps the number of secondary alignments per
read; `-p0.1` keeps a secondary only if it scores ≥ 0.1× the primary. Would bumping `-N50 → -N1000`
reveal more copies for the >50-copy **mega gene families** (KRAB-ZNF / satellite / high-copy segdup),
or would the within-`p0.1` secondary census plateau below 50 anyway (so `-N1000` = no gain)?

**Method.** Read the `-N50` census directly from the shipped BAM (it *is* the `-N50` result), for
13,968 primary reads over 5 mega-families spanning divergent→near-identical (ZNF92/ZNF578 KRAB-ZNF,
LOC115930164 tandem, XRCC2 & FRG1 near-identical segdups). Re-aligned the **cap-saturated** reads with
`-N1000` (and a matched `-N50` cross-check: N50 ⊆ N1000 always, N50 = min(N1000, 51); the only exceptions
are off-by-one at the cap boundary). Deterministic; intermediates in `winloci_scratch/`.

## Result — two clean regimes

| family | copies | per-read census (`-N50`) | `-N1000` effect |
|---|---|---|---|
| **XRCC2** (near-identical, id 0.991) | 11 | median 11, **max 49, 0 saturating** | **nothing** — plateaus at true CN, `-N`-insensitive |
| FRG1 / LOC115930164 (segdup/tandem) | 11–15 | ~19–22 | negligible |
| **ZNF92 / ZNF578** (KRAB-ZNF) | 41+ | **56% of ZNF92 reads saturate at 51** | census → median 89, **max 1001** |

- **Cap-hit fraction = 1,665 / 13,968 = 11.9%** saturate at exactly 51 (= 1 primary + 50 secondaries);
  0 records exceed 51 → the `-N50` cap is respected. Cap-hits are **overwhelmingly KRAB-ZNF**; clean
  high-copy families plateau *well under 50* and are completely `-N`-insensitive.

## The KRAB-ZNF "gain" is repeat fragmentation, not new copies

Re-aligning the cap-hit reads with `-N1000` (513-read sample, all at the `-N50` cap): 96.1% have a true
census > 51 (median 89, p99 745, max 1001). **But:**
- **records per distinct 100 kb locus = 2.75** — ≈3 fragmentary sub-alignments per gene (the C2H2
  zinc-finger tandem array producing partial hits), *not* distinct copy-loci.
- **distinct-locus census rises only +7 (median 33 → 43)** — diffuse KRAB-domain homology, **all
  MAPQ-0** (the permissive-`p0.1` tail; the primaries are often confident MAPQ-60).

## Cost

- Output inflation: **records 2.68×, BAM 2.32×** (on the cap-hit sample).
- Compute, same input: **~1.8–1.9×** (wall 1.81×, CPU 1.89×). *(An earlier ">3×" claim was
  unsupported — corrected.)* The absolute slowness of these reads (`mid_occ≈1136`) is intrinsic chaining
  cost shared by both `-N` values.

## Does it change any thesis quantity? No.

- **χ(H) / conflict census:** the +7 distinct loci are pure KRAB-ZNF over-enumeration — already the
  **K-frontier / abstain** regime, not counted as clean copies.
- **Copy assignment (O2):** abstains on exactly these MAPQ-0 mega-family reads → extra secondaries unused.
- **Copy-number depth probe:** does not consume secondaries at all.

## Verdict — keep `-N50`

It captures the true census for every clean/near-identical family (all plateau under 50) and truncates
only the KRAB-ZNF repeat-fragmentation tail, which the pipeline abstains on anyway. `-N1000` buys
2.3–2.7× larger BAMs and ~1.8× slower alignment for diffuse, MAPQ-0, abstained multimappings. If a real
KRAB-ZNF copy census is ever wanted, the fix is a **distinct-locus / CIGAR-collapsed metric** (collapsing
the ~2.75 fragments/locus), **not** raising `-N`.
