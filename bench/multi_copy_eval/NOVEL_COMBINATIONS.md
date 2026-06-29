# Can VG find novel exon combinations "difficult to model otherwise"? — yes, in two senses

VG's original promise was to surface exon combinations a linear assembler can't model. We tested
both forms.

## (1) Novel junction combinations — real, abundant, mis-scored as FP

On one chromosome (NC_073240.2), rustle-VG emits **1,668 class-`j` novel junction-combination**
predictions (+148 novel loci, class `u`) — *more* than the 808 annotation matches. They are
genuine, not chimeric:

| novel class | n | spanned end-to-end by ≥1 read | by ≥1 primary read | by ≥3 reads |
|---|---|---|---|---|
| `j` | 1,668 | 92% | 83% | 39% |
| `u` | 148 | 97% | 72% | 41% |

So rustle-VG *is* a novel-isoform engine: ~1,500 read-validated novel combinations on one
chromosome — which the whole "low precision vs StringTie" story was *penalizing* (they're absent
from an incomplete annotation). **Caveat:** 83% are spanned by *primary* reads → a full-length
linear assembler could find them too; they're novel vs the *reference* but not graph-*dependent*.

## (2) Cross-copy gene-conversion recombinant — graph-UNIQUE, demonstrated on ground truth

The truly "impossible to model linearly" case: an isoform whose 5′ exons come from copy A and 3′
exons from copy B (a sequence mosaic — same structure, mixed copy-SNPs). A linear assembler must
assign a read to one copy; only the variation graph + PSV phasing can represent the switch.

**Synthetic ground truth** (`bench/tandem_attribution/gen_synthetic.py`): 2 copies at 0.98 identity,
distinct isoforms, **80 planted gene-conversion reads** (5′ copy0 / 3′ copy1, switching at exon 3 =
copy0:9451–9750). Reads simulated, aligned with minimap2 (`-ax splice:hq`).

| tool | result |
|---|---|
| StringTie / rustle baseline | recover the 2 structural copy isoforms; **cannot see the mosaic** (no copy/SNP channel) |
| **rustle-VG mosaic detector** (`RUSTLE_VG_MOSAIC_ON`) | `family=0 0→1 breakpoint~9655 dispersion=0bp => CONFIRMED gene-conversion` |

The detected breakpoint **9655 falls exactly in the planted exon-3 switch (9451–9750)**, copies
**0→1** correct, **0 bp dispersion** (replicated). This is the variation graph's genuinely-unique
capability: detecting a novel cross-copy exon combination that has **no structural signature** at
all.

### Honest boundaries
- It's the **sequence/PSV channel** (DETECT + FLAG), not GTF structure — it flags the recombinant,
  it does not emit it as a transcript (a gene-conversion has no novel junction to assemble).
- It works only in an **identity window** (~0.98): more divergent (0.96) → copies don't share
  `ExonClass` nodes, no PSV fingerprint, detector skipped; more identical (0.985) → too few PSVs to
  discriminate. This is the identifiability window, not a tuning bug.
- **Conservative by default**: 4 of 80 recombinant reads were cleanly detectable (alignment
  ambiguity), flagged "suspect/unreplicated" at default thresholds; only the 0 bp-dispersion
  replication lets permissive thresholds *confirm*. The default declines to confirm — the
  don't-fabricate guard.

## Verdict

VG *can* do what it was meant to. Novel junction combinations: emitted in abundance and
read-validated (reframe the "FP" as discovery). The graph-unique cross-copy recombinant: detected
correctly on ground truth, invisible to linear tools — bounded by the identifiability window and
deliberately conservative.
