# Simulating TSPY: is the 0/5 tie-break invariance honest, or a method miss?

**Date:** 2026-07-11. **Substrate:** the 6-copy TSPY tandem array on gorilla chrY (`NC_073248.2`,
`34731504-34847734`), copies LOC129530275/277/276/278/279/280. **Question:** the real-data certificate reports
TSPY `0/5` copies tie-break-invariant. Is that honest (the copies are genuinely RNA-unresolvable) or a failure of
the method? A ground-truth simulation — reads with known copy-of-origin — is the only way to tell.

## The copies are 100% identical — the answer is information-theoretic

All-vs-all `minimap2 asm5` over the 6 genomic copy sequences:

| pair | identity |
|---|---|
| c275 = c277 = c278 = c279 | **100.000%** (2782/2782 bp) |
| c280 vs the group | 99.964% (1 bp) |
| c276 vs the group | 99.782% (~6 bp) |

Four of the six copies are **byte-identical over their entire 2782 bp**. A read from one is indistinguishable
from a read from another by *any* method — there is zero distinguishing information. TSPY's near-total
unresolvability is therefore not a method limitation; it is the definition of the K=0 frontier.

## Ground-truth simulation confirms it end-to-end

`tspy_sim.py`: build each copy's spliced mRNA, simulate 40 HiFi reads/copy (0.3% error) tagged with true
copy-of-origin, align (`minimap2 -ax splice:hq -N 50`), run `copy_assign`, evaluate against the tags.

**A subtlety that had to be fixed first (methodological note).** The annotation assigns these *identical* copies
**inconsistent exon boundaries** (annotated spliced lengths 1147/1140/1108). A naive simulation using per-copy
annotation exons injects fake copy-specific *junctions* and makes identical copies look resolvable. The faithful
simulation applies **one shared splice structure** to every copy's genomic sequence, so the identical copies
produce byte-identical transcripts (true splice sites are identical when the sequence is).

**Result (40 reads/copy):**

| group | minimap2 primaries | `copy_assign` outcome |
|---|---|---|
| 5 near-identical copies (200 reads) | **spread arbitrarily** across the identical positions, all MAPQ 0 | **100% tied / abstain, 0 misassignment** |
| c276 (40 reads, ~6 bp private) | maps **uniquely** (MAPQ > 0) | **100% resolved to c276** |

The certificate flags exactly this: **1/6 invariant = c276** (anchored 40, `tie_invariant=TRUE`); the five
identical copies are `anchored=0, tie=false, junction=false`. The pipeline **does not fabricate** assignments for
indistinguishable copies — it ties them — and it **does** resolve the one copy carrying real sequence divergence.

This is the honest outcome, proven with ground truth: where the copies carry information (c276), assignment is
100% correct; where they carry none (the identical four), the method abstains rather than guessing. The `0/5` on
real data is the certificate correctly reporting that the *expressed* TSPY copies are the identical ones (c276,
the one resolvable copy, had 0 reads in the real sample — it is exactly the copy the real run missed).

## Why injected exonic SNVs do NOT rescue it (and why DAZ2 needed junctions)

A second arm injected 5 copy-specific exonic SNVs into every copy. It **did not** increase resolution
(`psv_cols=1` in both arms, identical group still 100% tied). Cause: when reads spread arbitrarily across the
identical genomic positions, `copy_assign` assembles each position's reads into a per-position consensus that
**mixes all copies**, averaging the injected SNVs away. Only divergence large enough to make reads **map uniquely**
(c276's ~6 bp) escapes the mixing. This is precisely why the real DAZ2 was recovered by **copy-specific junctions**
(structural, position-anchored) rather than exonic PSVs — and why exonically-identical tandem arrays are the
genuine K=0 wall.

## Would a better aligner (tuned minimap2 / winnowmap) help? No — measured

The barrier is information-theoretic, so no aligner can cross it. Confirmed on the ground-truth sim (arm A):

| aligner | uniquely-mappable reads (MAPQ>0) | copies tie-break-invariant |
|---|---|---|
| minimap2 `splice:hq` (default) | 40 / 240 | 1/6 (c276) |
| minimap2 `-k11 -w5 -p0.1 -N200` (sensitive) | 40 / 240 | 1/6 |
| **winnowmap `splice:hq`** (repeat-specialist, meryl-weighted) | **40 / 240** | **1/6** |

All three place c275's reads *arbitrarily* across the identical positions and resolve only c276. Winnowmap — the
tool built specifically for repetitive mapping — gives byte-for-byte the same answer, because its advantage is
mapping reads to the right *region* and calibrating MAPQ in repeats, not manufacturing sequence differences
between 100%-identical copies. There are none to find.

There is also **no flank lever**: the whole ~8 kb tandem unit (gene + intergenic) is **99.005% identical**, so
even reads extending into the flanks (which Iso-Seq mRNA reads do not) would mostly tie. The only escapes remain
copy-specific divergence large enough to map uniquely (c276), DNA spanning the ~1% divergent flank positions, or
aggregate quantification — none of them an alignment-parameter change.

## Conclusion

TSPY `0/5` is **honest and correct**. The unresolvable copies are literally identical (0 distinguishing bases);
the tie-break invariance certificate flags exactly the copies ground truth confirms are unresolvable, and passes
exactly the one (c276) that carries real signal. No method — ours, minimap2, IsoCon, phylogeny — can resolve
100%-identical copies from RNA; the correct answer is to abstain, which is what the certificate reports.

## Reproduce

```
python3 tspy_sim.py 40                       # -> tspy_simA.fastq (real) + tspy_simB.fastq (SNV-injected)
minimap2 -ax splice:hq -N 50 GGO.fasta tspy_simA.fastq | samtools sort -o tspy_simA.bam - ; samtools index tspy_simA.bam
copy_assign --bam tspy_simA.bam --fasta GGO.fasta --region NC_073248.2:34731504-34847734 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --out tspy_resA
# quant.tsv: c276 tie_invariant=TRUE (anchored 40); the 5 identical copies anchored=0, all tied.
```

Related: `bench/PRIMARY_SECONDARY_INVARIANCE.md`, `bench/KNOWN_FAMILY_REGRESSION.md` (TSPY 5/6),
`project_k0_frontier_unresolvable`, `project_daz2_locus_support` (junction-defined recovery).
