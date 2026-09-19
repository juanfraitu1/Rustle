# The chr16 arm `build_spliced_seq_with` demands before `RUSTLE_JUNCTION_MAJORITY` can be a default

`denovo_assemble.rs` documents a harm measured on chr16 and states: **"REQUIRED BEFORE FLIPPING: a chr16
arm."** §6m8 ran only 27 NPIP windows and explicitly disclaimed being that arm. This is it: a full
chr16 RNA catalog, both ways, nothing else changed.

Command (twice, identical but for the env var):
`gw_family_catalog --bam chr16.bam --fasta chm13v2.0.fa --threads 3 --out {off,on}`
`chr16.bam` = 682,958 primary reads. ~57 min per arm.

## Result

| metric | OFF | ON | delta | direction the warning feared |
|---|---|---|---|---|
| **families** | 282 | **290** | **+8** | FEWER (fusion) — **did not happen** |
| copies | 1,404 | 1,418 | +14 (+1.0%) | more |
| **strictly-engulfed copies** | 79 | **86** | **+7 (+8.9%)** | more |
| **max family size** | 71 | **71** | **0** | bigger (hub fusion) — **did not happen** |
| families of size 2 | 127 | 134 | +7 | — |
| median family size | 3.0 | 3.0 | 0 | — |

The documented harm was **families 121 → 117, copies 678 → 700, engulfed 60 → 63**.

## Reading

**The failure mode the flag was held back for — family FUSION — does not occur at chr16 scale.** Families go
UP (+8), the largest family is unchanged at 71 copies, and the median family size is unchanged. The +8
families are almost entirely new 2-copy families (+7), i.e. the flag finds small families that the strict
motif rule was deleting outright, rather than gluing existing ones together.

**The residual cost is engulfment: +7 strictly-engulfed copies (+8.9%)**, the same direction as the
documented harm (+5%) and of comparable size. That is the real price, and it is a precision cost on copy
boundaries, not a family-structure cost.

## Recommendation

The blocking condition in the code comment is now satisfied and its specific fear is refuted. What remains
is a judgement call about 8.9% more engulfed copies against the §6m8 gains (+11 canonical junctions,
complete NPIP chains 9 → 10 on that substrate, 96% of added junctions read-supported).

**I am not flipping the default.** Two reasons: the decision is the user's, and one chromosome is one
substrate — the standing "hold a substrate back" rule means a second chromosome (or the gorilla contigs)
should agree before a shipped default moves. The arm is recorded so that decision can be made on numbers.
