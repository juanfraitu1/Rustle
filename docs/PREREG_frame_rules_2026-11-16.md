# PREREG, register 1766: reading-frame rules for chains every tool misses (2026-11-16)

User goal: rules that raise sensitivity and precision against RefSeq, from chains StringTie, FLAIR, IsoQuant, bambu and the
reference-free arm (`isoseq cluster2`, standing in for IsoCon, which needs a family panel) do not output although reads carry them.
Samples: chimp PTR testis (development), bonobo PPA testis (held out; same tissue, sister species: a weak transfer test).
Workspace `r1766_allmiss`; scripts `bench/famcmp/allmiss.py` and the `bench/famcmp/frame_rules/` (`cand.py`, `orf.py`, `join.py`, `rules.py`; copies of the workspace scripts).

## Frozen on PTR

No PPA candidate, feature or precision table was opened before the freeze. Seen before it: the PPA counts of RefSeq chains by
evidence tier and by what ours has at the locus (RefSeq-conditioned, no rule precision).

Features of a chain: `ptc` = the longest ATG..stop ORF of the model, on its strand, stops > 50 nt before the last junction;
`orf` = that ORF's length in nt; `longer` = primary reads whose chain contains this chain as a contiguous part and is longer;
`locus` = no default model shares an intron (`none`), a default model is a contiguous part of the chain (`extends`), the chain is
inside a default model (`inside`), else `alt`; `sister` = primary reads of the sister species' testis sample with exactly the lifted
chain (`xspecies_chains.py`, PTR <-> PPA maps).

- **P1 (precision):** drop every default chain with `ptc`.
- **S1 (sensitivity, one sample):** admit a chain of the 1-read-floor run that the default lacks when: no `ptc`, `orf` >= 900,
  one primary read has exactly the chain, `longer` = 0, `locus` is `none` or `extends`.
- **S2 (sensitivity, two samples):** as S1 with `orf` >= 300 and any `locus` but `inside`, plus `sister` >= 1.

PTR values: default 57,039 chains / 21,904 correct / 0.384; 2-read chains 0.226. P1 drops 9,599 chains at 0.111 (4.8% of the
correct chains); kept set 0.439 against 0.398 for the read cut at the same correct count (final.gtf 0.513 against 0.478).
S1 3,443 new chains, 818 correct, 0.238. S2 3,287 / 1,275 / 0.388. P1 + S1 + S2 on `ours.gtf`: 53,586 / 22,613 / 0.422.

## Bars on PPA (all computed on PPA itself)

- P1 PASS: kept-set precision beats the read cut at the same correct count by >= 0.02 on `ours.gtf` and on `ours.final.gtf`, and
  loses <= 6% of the correct chains.
- S1 PASS, S2 PASS (each): precision of the new chains >= precision of PPA's default 2-read chains, and new correct chains >= 2% of
  the default's correct chains.
- Combined PASS: P1 + S1 + S2 on `ours.gtf` has more correct chains and higher precision than the default.

Not tested here: another tissue, a non-sister pair, NMD-annotated or non-coding RefSeq records lost by P1 (counted after the verdict).

## Results (2026-11-16)

Who outputs each RefSeq multi-exon chain, by read evidence (`PTR.t1.tsv`, `PPA.t1.tsv`; tools = StringTie, FLAIR, IsoQuant, bambu at
its default NDR and at NDR 1, `isoseq cluster2`; ours = `ours.gtf`). `full` = a primary read has the whole chain; `linked` = no such
read, every intron and every adjacent intron pair in a read; `nodata` = an intron in no primary read (left out of everything below).

| | ours only | tools only | both | nobody |
|---|---|---|---|---|
| PTR full | 352 | 5,254 | 20,927 | 6,783 |
| PTR linked | 366 | 461 | 207 | 18,917 |
| PTR nodata | 46 | 0 | 2 | 71,925 |
| PPA full | 321 | 3,098 | 17,835 | 4,931 |
| PPA linked | 176 | 224 | 103 | 12,538 |
| PPA nodata | 105 | 1 | 3 | 50,111 |

RefSeq chains per arm, PTR / PPA: ours 21,904 / 18,545; StringTie 21,532 / 17,288; FLAIR 19,730 / 15,934; cluster2 17,417 / 15,339;
IsoQuant 15,761 / 13,109; bambu NDR 1 13,835 / 12,429; bambu default 10,164 / 9,631 (bambu could not train on either sample and used
its pretrained model).

The `nobody, full` class is one read per chain (PTR 6,700 of 6,812 before bambu was added). What ours has at those loci (PTR, 6,812):
another combination of the same introns 3,905, no model sharing an intron 1,683, a model that is a part of the chain 955, the chain
inside a longer model 269.

Verdicts on PPA, held out:

| Rule | PPA | Bar | Verdict |
|---|---|---|---|
| P1 on `ours.gtf` | 0.370 vs read cut 0.323 at 18,149 correct; 2.1% of correct lost | +0.02, <= 6% | **PASS** |
| P1 on `ours.final.gtf` | 0.412 vs 0.367; 2.0% lost | +0.02, <= 6% | **PASS** |
| S1 | 3,563 new chains, 681 correct, 0.191 | >= 0.156 (2-read chains), >= 371 | **PASS** |
| S2 | 3,341 new, 920 correct, 0.275 | same | **PASS** |
| P1 + S1 + S2 | 55,445 / 19,533 / 0.352 vs default 58,953 / 18,545 / 0.315 | more correct, higher precision | **PASS** |

`ours.final.gtf` with all three: PTR 45,114 / 21,767 / 0.482 (was 45,418 / 20,993 / 0.462); PPA 48,052 / 18,529 / 0.386 (was 48,736 /
17,501 / 0.359). gffcompare on the emitted GTFs (`TAG.rules.gtf` from `emit.py`), intron-chain sensitivity / precision: PTR 16.8 / 38.5
-> 17.4 / 42.3; PPA 20.2 / 31.5 -> 21.3 / 35.3.

Of the chains nobody outputs that a read carries whole, S1 or S2 recover 1,014 of 6,783 (PTR) and 892 of 4,931 (PPA); the 1-read-floor
run holds 2,961 and 2,191 of them. The rules recover none of the `linked` class.

What ours has at the loci of the `nobody` chains (full / linked): another combination of the same introns PTR 3,905 / 10,937,
PPA 2,988 / 6,985; a model that is a part of the chain 955 / 5,036, 502 / 2,158; no model sharing an intron 1,683 / 1,975,
1,093 / 1,626; the chain inside a longer model 269 / 975, 373 / 1,775 (counts before bambu was added).

Post hoc, after the verdict:

- **S1 does not beat the margin once P1 is applied.** With P1 on, the weakest accepted tier is the 2-read chains without the stop:
  PTR 0.271, PPA 0.196 (ORF >= 900 only: 0.326, 0.238). S1 is 0.238 and 0.191, under it on both; S2 is 0.388 and 0.275, over it on
  both. S1 is a floor trade inside the frame-clean class; S2 is the one admit rule above the margin, and it needs a second sample.
  Against P1 alone on PPA (49,021 / 18,149 / 0.370) the admits give 55,445 / 19,533 / 0.352.
- `orf()` passes a toy self-check (three strings); strand handling of `orf.py` is not tested directly, the 1% loss among coding
  chains is the indirect evidence.

- **P1 on five more samples** (default assemblies of `apetools/`; kept precision vs read cut at the same correct count, `ours.gtf`):
  PPY testis 0.234 vs 0.197, PAB testis 0.263 vs 0.235, PTR lymphoblastoid 0.388 vs 0.355, PPA fibroblast 0.280 vs 0.224, PAB
  fibroblast 0.386 vs 0.354; correct chains lost 1.9% to 6.1% (PAB testis is over the 6% bar by 0.1 point).
- **P1 is a coding-transcript rule.** Correct chains lost by RefSeq class: coding (NM / XM) 1.3% (PTR) and 0.8% (PPA); non-coding
  (NR / XR) 46% (804 of 1,748) and 53% (253 of 481). RefSeq's XM models are themselves ORF-based predictions, so part of P1's gain is
  agreement with the annotator's own criterion.
- **Read-only generators stay closed** (PTR, bar 0.226): one-read chains by locus class, intron support and novelty, best cell 0.227;
  joins of two read chains sharing introns 155,183 candidates at 0.014; a read's local event set into a default model 192,441 at 0.031.
  The `linked` class has no rule. Chains with a non-canonical intron: ours emits none; StringTie 1,724 (38 correct), FLAIR 9,018 (43),
  IsoQuant 747 (26).

Not done: port to the assembler (S1 needs the longest-ORF scan and the `longer` / `locus` features at emission; P1 is a filter on the
finished GTF); S1 / S2 on the orangutans or cell lines; a P1 variant that spares non-coding loci; IsoCon proper (panel-only).

## User decisions (2026-11-16)

- **P1 is an attribute, not a filter.** `tools/gtf_ptc.py` (self-test `--selftest`) adds `orf_len` and `ptc` to the transcript lines
  of `PREFIX.final.gtf` as the last step of assemble (`RUSTLE_PTC=0` = off; nothing dropped). PPA `ours.final.gtf`: 7,013 of 48,736
  multi-exon transcripts carry `ptc "1"` (the test above counted 7,108 chains with the ORF of the `ours.gtf` model; layered ends differ).
- **S2 is not ported:** no sister-species sample can be assumed at analysis time. S1 stays unported (under the margin once P1 applies).

## Addendum: one read plus agreeing reads (user 2026-11-16), frozen on PTR

User: a 1-read rule where more than half of the other reads agree; look at the reads on each junction and their colinearity.
`frame_rules/agree.py`: for a new chain C with one full read, every other primary read sharing an intron with C is `agree` (its chain
is a contiguous part of C), `longer` (C is a contiguous part of it) or `conflict`. Per junction, `jsolo` counts junctions only the
full read has and `jmin` is the lowest share of colinear reads among the reads on a junction. Rules (all: one full read, `longer` = 0):

- **M1, the rule as asked:** agree / (agree + conflict) > 0.5. PTR 3,264 new chains, 0.142.
- **S3:** `jsolo` = 0 (every junction also in another read), no `ptc`, `orf` >= 900. PTR 3,017 / 926 / 0.307.
- **S3b:** S3 with `jmin` > 0.05. PTR 1,417 / 522 / 0.368.

PTR, seen while choosing: the agreeing share separates little on its own (<= 0.2: 0.093; > 0.5: 0.14), because at an expressed gene
most reads on a shared junction belong to the main isoform; `jsolo` = 0 is 0.175 against 0.08 with one solo junction. PPA's
`agree.tsv` was written by the same command and not opened. PPA is not blind for S3's frame half (no stop, ORF >= 900 is S1's, already scored
on PPA); only the junction features were unseen there. With this addendum six rules have been judged on PPA.

Bars on PPA (each rule): new-chain precision >= PPA's 2-read chains without the stop (0.196, the margin once `ptc` is in use) and
new correct chains >= 371.

### Addendum verdicts on PPA, held out

| Rule | PPA new chains / correct / precision | Bar (0.196, 371) | Verdict |
|---|---|---|---|
| M1 majority agreement | 2,817 / 216 / 0.077 | | **FAIL** both bars |
| S3 every junction in another read + frame | 3,719 / 792 / 0.213 | | **PASS** |
| S3b S3 + `jmin` > 0.05 | 1,497 / 444 / 0.297 | | **PASS** |

Added to the default (chains / correct / precision):

| | default | + S3 | + S3b |
|---|---|---|---|
| PTR `ours.gtf` | 57,039 / 21,904 / 0.384 | 60,056 / 22,830 / 0.380 | 58,456 / 22,426 / 0.384 |
| PTR `ours.final.gtf` | 45,418 / 20,993 / 0.462 | 48,435 / 21,919 / 0.453 | 46,835 / 21,515 / 0.459 |
| PPA `ours.gtf` | 58,953 / 18,545 / 0.315 | 62,672 / 19,337 / 0.309 | 60,450 / 18,989 / 0.314 |
| PPA `ours.final.gtf` | 48,736 / 17,501 / 0.359 | 52,455 / 18,293 / 0.349 | 50,233 / 17,945 / 0.357 |

Of S3's correct chains, 346 (PTR) and 331 (PPA) are in no other tool (S3b 144 and 128): 5% and 7% of the whole-read chains nobody
outputs. The rest of S3's gain is chains another tool has and ours drops at the 2-read floor.

Reading: a majority of agreeing reads is the wrong test for a minor isoform (its junctions are shared with the main one, whose reads
conflict by construction). What carries is that no junction rests on the one read alone, plus the frame. S3 is one-sample and
replaces S1 (0.191 on PPA). S3 raises correct chains at a precision cost of 0.004-0.010; S3b is close to neutral on precision.
Not done: port; S3 on another tissue or genus; S3 with `ptc` filtering of the default.

### Like-for-like tier, and what the junction feature adds (post hoc)

- The accepted chains that match S3 on frame are the 2-read chains with no stop and ORF >= 900: PTR 0.326, PPA 0.238. **S3 is under
  that tier on both (0.307, 0.213): a floor trade, as S1 was. S3b is over it on both (0.368, 0.297).** S3b is the rule to port if one is.
- One read, no longer read, no stop, ORF >= 900, split on the junction feature: a junction resting on the one read alone 0.159 (PTR)
  and 0.109 (PPA); none 0.307 and 0.213.
- "More than half agree" read per junction (`jmin` > 0.5): PTR 2,861 chains at 0.133, PPA 2,572 at 0.066. A junction only the full
  read has scores 1.0 there; with every junction also in another read only 83 and 73 chains remain (19 and 17 correct).
- Chains whose every adjacent junction pair is in another colinear read, no stop, ORF >= 300: PTR 222 at 0.374. Cleanest signal, tiny yield.

### The reads on each junction

`frame_rules/junc.py` -> `r1766_allmiss/{ptr,ppa}.s3.junctions.tsv`: one row per junction of every S3 chain (reads on the junction,
the full read, agree, longer, conflict, colinear share). Two PPA chains on NC_073249.2:

| chain | junction | reads | full | agree | conflict | colinear share |
|---|---|---|---|---|---|---|
| correct, 5 introns, from 21,200,406 | 1 | 38 | 1 | 0 | 37 | 0.026 |
| | 2 | 2 | 1 | 0 | 1 | 0.500 |
| | 3 | 34 | 1 | 0 | 33 | 0.029 |
| | 4 | 43 | 1 | 1 | 41 | 0.047 |
| | 5 | 42 | 1 | 2 | 39 | 0.071 |
| wrong, 4 introns, from 1,132,377 | 1 | 2 | 1 | 0 | 1 | 0.500 |
| | 2 | 21 | 1 | 0 | 20 | 0.048 |
| | 3 | 21 | 1 | 0 | 20 | 0.048 |
| | 4 | 11 | 1 | 0 | 10 | 0.091 |

Both look alike at this level: the junctions are shared with a main isoform whose reads conflict, and the one rare junction (2 reads)
is what makes the chain new. The correct one has a few agreeing partial reads at its 3' junctions; that is the S3b signal.

Limits: colinearity compares intron chains only. A read whose end exon runs into one of the chain's introns counts as agreeing, and
unspliced reads are not counted. The "+ S3 on `final.gtf`" rows assume admitted one-read chains pass the chain-score filter, which is
untested: a porting question. The `final_ptc` driver step was run end to end on a 20 Mb PPA region (290 multi-exon transcripts
tagged, 42 with `ptc "1"`; `PREFIX.gtf` untouched).

## S3b ported as opt-in; orangutan test for the default (user 2026-11-16), bars set before the run

Port: `tools/oneread.py` (self-test; reproduces the prototype's PPA set, 1,497 of 1,497 chains), driver stage `oneread` ->
`PREFIX.oneread.gtf`; `RUSTLE_ONEREAD=1` makes assemble append those models to `PREFIX.final.gtf` (after the chain-score filter,
before the `ptc` tag). User: test on an orangutan and make it the default if promising.

Samples, unseen for this rule: PPY testis-library `PPY_mm`, PAB `PAB_mm`, and PAB fibroblast `PAB_AG06213_fib_mm` (another tissue).
Fresh `assemble` with `RUSTLE_ONEREAD=1` from a `tools/` snapshot (`r1766_allmiss/orang/`). Per sample, on that sample:

1. precision of the admitted chains >= precision of the default's 2-read chains with no early stop and ORF >= 900;
2. admitted correct chains >= 2% of the default's correct chains;
3. `final.gtf` precision with the models appended falls by <= 0.005.

**Default if all three samples pass all three bars.** Anything else: stays opt-in, reported as is.

### Outcome (job 4477777, 12 min; scored with `frame_rules/orang_score.py` against `pongo_asm/{PPY,PAB}.refseq.gtf`)

Chains / correct / precision. Tier = the default's 2-read chains with no early stop and ORF >= 900.

| sample | default `o.gtf` | admitted | tier | bar 1 | yield bar (2%) | bar 2 | `final.gtf` without -> with | bar 3 |
|---|---|---|---|---|---|---|---|---|
| PPY_mm | 113,109 / 21,347 / 0.189 | 1,482 / 299 / 0.202 | 0.186 | PASS | 427 | FAIL | 86,535 / 20,620 / 0.238 -> 87,964 / 20,907 / 0.238 | PASS |
| PAB_mm | 75,338 / 16,931 / 0.225 | 1,026 / 268 / 0.261 | 0.268 | FAIL | 339 | FAIL | 55,876 / 16,418 / 0.294 -> 56,860 / 16,671 / 0.293 | PASS |
| PAB_AG06213_fib_mm | 73,064 / 25,253 / 0.346 | 577 / 169 / 0.293 | 0.244 | PASS | 505 | FAIL | 59,495 / 23,648 / 0.397 -> 60,052 / 23,810 / 0.396 | PASS |

**FAIL as registered: not the default.** The yield bar fails on all three samples (1.0-1.6% of the default's correct chains
against 2%), and PAB testis misses bar 1 by 0.007. `final.gtf` precision moves by at most 0.001. `oneread` stays an opt-in stage
(`RUSTLE_ONEREAD=1` to append to `final.gtf`): +160 to +290 correct chains per sample at about the precision of the file it joins.
