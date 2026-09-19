# `RUSTLE_GTF_SECONDARY` — admitting secondary alignments into the `--gtf` read pool

**Verdict: the largest single gain measured in this line of work — complete NPIP chains 10/26 → 20/26 and
canonical junctions 170/209 → 197/209 (81% → 94%) — at a 9x transcript cost. Shipped OPT-IN, default off,
byte-identical when unset.** Lib suite **882 passed / 0 failed**.

## What the diagnostic found

Hunting the "selection defect" of §6n1 — 18 of 26 copies have a read carrying their whole canonical chain,
yet only 10 are emitted — the cause turned out not to be selection at all:

| read pool | copies with a complete canonical chain available |
|---|---|
| **PRIMARY only** | **10/26** |
| primary + secondary | **18/26** |

**10/26 is exactly what the assembler emits.** It was never dropping anything: it was already at the
primary-only ceiling. The 8 copies in between have their complete chain present at that locus *only as a
secondary record* — reads minimap2 placed primarily at a different NPIP paralogue.

⚠ So §6n1's framing ("the assembler drops chains that exist whole in single reads") was wrong in mechanism.
It is a **read-pool scope** decision, not a selection defect. Corrected here.

## The change

`reads_in_region` now builds its site pool with `alignment_read_from_record(.., gtf_secondary_enabled())`
instead of `primary_read_from_record`. `RUSTLE_GTF_SECONDARY` defaults OFF, so the pool is primary-only and
output is byte-identical (verified: the OFF arm is byte-identical to the §6m8 `mj3.gtf`). Supplementary
records stay excluded in both modes.

This is the O1⊥O2 abstention the module already documents for site construction: *a read may support
several sites at once, because "which copies exist" (O1) must not require first answering "where did this
molecule come from" (O2)*. It is explicitly **not** a copy assignment — that is O2's job, and the certificate
path is untouched.

## Result (26 spliced NPIP copies, §6m7 canonical truth, `--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1`)

| arm | transcripts | canonical / 209 | **complete chains** | chains containing the full set |
|---|---|---|---|---|
| primary only (shipped default) | 1,456 | 170 (81%) | **10/26** | 7 |
| **+ secondary** | **13,132** | **197 (94%)** | **20/26** | **15** |
| — oracle: union of all alignments | — | 205 (98%) | 22/26 | — |

**+10 copies and +27 junctions**, reaching 90% of the way to the all-alignment oracle.

**The cost is real: 9x more transcripts (1,456 → 13,132).** Most of the additions are multimapper echoes of
the same molecule at several paralogues — by construction, since that is what a secondary record is. Any
downstream consumer that counts transcripts, or treats each as an independent observation, will be wrong
under this flag.

## Recommendation

Use it when the question is **"which junctions/copies exist"** (O1, reconstruction, this ladder). Do NOT use
it when the question is **"how much of what"** — quantification, expression, or anything per-molecule — until
the echoes are collapsed. Default stays off for that reason.

**Not done:** collapsing the echoes (one representative per molecule across its placements) would keep the
+10 copies and undo most of the 9x. That is the obvious follow-up and it is untested.
