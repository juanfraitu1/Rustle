# PREREG row 1667: which single-exon rule removes the no-model family members (2026-11-13)

Assembly-layer diagnosis (step 0 only, descriptive). User request after row 1666.

A119b genome-wide `assemble`, same binary and tree as row 1666's frozen copy (workspace `r1666_mono_root/src3`,
switch unset), one rule off per arm:
- `p`: priming filter off (`RUSTLE_POLISH_MONO_PRIMED=0`);
- `s`: shadow rule off (`--polish-mono-shadow` removed);
- `q`: read floor off (`--polish-mono-quantile 0.0`);
- `psq`: all three off.
Fixed truth set: the 291 no-model members of the baseline funnel (`r1666_mono_root/g_off/funnel.tsv`), split by the
earlier read classes. Count per arm: members that gain a model (famtruth.expressed_from_gtf on `PREFIX.families.gtf`,
the funnel's own test), and the cost: single-exon models added and multi-exon transcripts changed.
Reading rule, frozen: a rule is a candidate for a registered fix only if its arm alone gives >= 10 members a model;
if only `psq` does, the loss is shared and no single rule is the lever; if none does, the members are lost before
the polish (seeding or clustering) and the next step is there. No rule is changed by this row.

## Result (job 4468872)
Baseline: 13,508 single-exon and 218,258 multi-exon transcripts; multi-exon count identical in every arm.

| arm | members gaining a model (of 291) | single-exon genes among them | single-exon models added | added models per member |
|---|---|---|---|---|
| `p` priming filter off | 11 | 3 | 13,884 | 1,262 |
| `s` shadow rule off | 22 | 16 | 42,648 | 1,939 |
| `q` read floor off | 25 | 3 | 29,809 | 1,192 |
| `psq` all off | 103 | 30 | 236,544 | 2,297 |

Overlaps between single arms: p&s 0, p&q 1, s&q 2; 48 members appear only with all three off.
By the frozen reading rule each rule alone is a candidate (>= 10), but none is a lever as a blunt switch: every arm
pays more than a thousand added single-exon models per member. The loss is spread over the three rules and nearly
half of the recoverable members need two of them lifted together. 188 of the 291 gain no model even with all three
off (132 of the 291 have fewer than 3 own reads). Ceiling of any single-exon polish change: 103 members.
No rule changed.
