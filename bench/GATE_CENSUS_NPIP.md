# Where pass-1 skeletons die, and why the "ceiling" was not a depth limit

**Headline: §6m4's "16.1% of NPIP junctions have no read in this library" was the WRONG explanation.
Those junctions are NON-CANONICAL, and the reads are right to be absent.** 209 of 249 annotated NPIP
junctions carry a canonical motif — **exactly the 209 that §6m4 called the all-alignment ceiling.**

Lib suite **881 passed / 0 failed**. Instrumentation: `denovo_assemble::assemble_gate_census` + `GateCensus`
(`kept + the four rejection buckets == skeletons.len()` is an asserted invariant), printed by
`copy_assign` under `RUSTLE_GATE_CENSUS=1`. Print-only — the transcripts are unchanged.

## 1. The census: the gate's only real filter is the sequence/motif test

| arm | skeletons | kept | rej reads | rej span | **rej seq (motif/coords)** | rej len |
|---|---|---|---|---|---|---|
| k=0 | 1,093 | 750 | 0 | 0 | **343 (31.4%)** | 0 |
| k=3 | 2,978 | 1,124 | **1,369** | 0 | 485 | 0 |

At the default, **nothing** is lost to the read floor, the span cap or the length window — the entire
pass-1 → GTF loss is `build_spliced_seq` returning `None`. At k=3 the widening adds 1,885 skeletons and the
read floor immediately kills 1,369 of them, which is why §6m6's W-2 failed.

**Neither published knob reaches this path.** `RUSTLE_JUNCTION_NC_MAX_BP=100000000` and
`RUSTLE_GATE_MIN_READS=2` both left the census **bit-for-bit identical** (343/1,369 rejections unchanged).
The non-canonical *size* tolerance is not what binds; the "at least one canonical junction" requirement is.

## 2. Why: the annotation's missing junctions are the non-canonical ones

Motifs of all 249 annotated NPIP junctions, strand-aware, against whether any alignment (primary OR
secondary) carries them:

| | observed in reads | NOT observed |
|---|---|---|
| **canonical** (GT–AG / GC–AG / AT–AC) | **205** | 4 |
| **non-canonical** | 4 | **36** |

- P(observed \| canonical) = **0.981**
- P(observed \| non-canonical) = **0.100**
- 209 canonical total = **the exact §6m4 "ceiling" of 209/249**.

209 of 249 are `GT..AG`; the rest are one-sided variants (`GT..AA` ×6, `GT..TG` ×5, `GT..GA` ×3 …) —
the signature of annotation, not biology. No copy has zero canonical junctions; 17 of 26 have at least one
non-canonical one.

**So the reads and the assembler agree with each other and disagree with the annotation.** The assembler
refusing to build a transcript through a non-canonical junction is correct behaviour, not a defect.

## 3. Rescored against a canonical-only truth

| arm | all junctions (9-copy / 26) | **canonical-only (9-copy / 26)** | **complete chains** |
|---|---|---|---|
| k=0 baseline | 65/106, 151/249 | **65/94, 151/209** | **3/26 → 8/26** |
| k=3 widening | 67/106, 161/249 | **67/94, 161/209** | **3/26 → 9/26** |

Removing 40 junctions that are almost certainly annotation artifacts nearly **triples** complete-chain
recovery, from 3/26 to 8–9/26, with no change to the assembler at all.

## Corrections this forces

- **§6m4** framed the 16.1% as "absent information, not misassignment", implying depth. It is neither:
  the junctions are very likely not real. The claim "a second library would be needed" is withdrawn —
  a second library would not create GT..AA sites.
- **§6m5/§6m6**'s gap of 60.4% vs an 85.8% "ceiling" is partly an artifact of scoring against those 40.

## What is still genuinely open

Against the canonical-only truth the assembler still reaches only 151/209 (72%) at k=0 and 161/209 (77%)
at k=3, so **~48 canonical junctions are still lost** — and the census says they die at
`build_spliced_seq` too. Next: separate that 343-skeleton bucket into *fetch failure*, *no canonical
junction anywhere in the chain*, and *one bad junction in an otherwise canonical chain*. Only the third is
recoverable, and it is the one worth a rule.
