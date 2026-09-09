# PREREG — three-tool bakeoff on one multi-copy family, scored by ONE derived rule (2026-09-08)

**Written before any comparison output was produced.** md5 recorded in
`/mnt/linuxdisk/home/juanfraitu/bakeoff/PREREG.md5`.

## The question the user asked
*"Can we make this more of an aligner-adjacent pipeline so it is comparable to StringTie and flair, but
emphasising the pieces they miss or conflate because of secondaries?"*

## Why a new instrument is needed at all
⚠⚠ **The obvious metric is vacuous and must not be used.** Scoring on "copy attribution accuracy" gives
StringTie, flair and isoseq a zero **by construction** — none of them emits a copy attribute (§6gj, F1–F5 of
§6gs). A metric only one tool can produce proves nothing; that is the failure class already tracked in
`project_vacuous_instruments`.

⭐ **The fix: derive the copy call from every tool's output by ONE identical rule.** All four tools emit a
GTF whose transcripts have genomic coordinates, and a read's intron chain either matches a transcript or does
not. So for *every* tool: molecule → matching transcript(s) → the copy interval those transcripts sit in.
Nothing tool-specific enters the scorer. `bench/tool_bakeoff.py` implements exactly this and is run with the
same three inputs for each arm.

**Derived states (identical for all arms):** `derived_one` (all matching transcripts sit in one copy),
`derived_multi` (matching transcripts span ≥2 copies — the tool CONFLATED copies), `derived_none` (no
transcript carries this molecule's chain). ⚠ `derived_multi` is a property of the emitted transcript set, not
an accusation of error; the truth arms below are what make it an error.

## Substrate (fixed now, not chosen after seeing results)
- family `MCL1_073242`, **80 copies**, `/mnt/linuxdisk/home/juanfraitu/mcl_ann/sweep_v19/fam_MCL1_073242/copies.tsv`
- reads `/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam` (fibroblast IsoSeq, 15,961 primaries in-family)
- reference `/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa` (3 contigs)
- flair **3.0.0**, installed at `/home/juanfra/miniforge3/envs/flair` — ⭐ all six §6gs findings were
  re-verified line-for-line in the INSTALLED version before this file was written (the source read was
  3.0.1+master; `flair_align.py:150/156`, `-N 4` in quantify/collapse, the `tlen`/`tname` sort keys, the
  `mapping_quality = 60` overwrite are byte-identical). The findings therefore describe what will actually run.
- StringTie: the build already in `tools/stringtie`.

## Predictions
| # | prediction | what refutes it |
|---|---|---|
| **P1** | flair places a **smaller** fraction of molecules in `derived_one` than we do, because `flair align --secondary=no` gives its BED one placement per read | flair ≥ ours |
| **P2** | flair emits **exactly zero** molecules in any tied/abstaining state — it has no such state (§6gs F3) — while ours is non-zero | any flair abstention |
| **P3** | flair's `derived_multi` count **exceeds** ours: its transcripts absorb reads whose origins we separate | flair ≤ ours |
| **P4** | **excision arm** — with one copy removed from the reference, flair and StringTie assign **≥90 %** of that copy's reads to a surviving copy, while ours abstains (§6ff measured 100 % abstention, 0 silent) | either tool abstaining on ≥10 %, or ours failing to abstain |
| **P5** ⚠ | **guard against our own framing** — flair recovers **more** read-supported junctions than we do, as StringTie did (58.6 % vs 42.1 %, §6gj) | flair recovers fewer, which would be a surprise worth reporting rather than a win to claim |

## ⛔ The falsifier that kills the framing
**If flair's `derived_one` calls agree with our assigned calls on ≥ 95 % of the molecules we assign, the
secondary-alignment difference does not matter on this substrate** and the "they cannot see the ambiguity"
line must be dropped from the thesis framing, however true it is of the source code. Reading code is not
evidence of a measurable difference; this line is what converts F1–F5 into a result or retires them.

## Rules held
- ⚠ **Never pool human and gorilla numbers** — this is gorilla only.
- ⚠ Each arm is scored by the SAME `bench/tool_bakeoff.py` invocation shape; no per-tool tuning.
- ⚠ Tools run at their **documented defaults**. flair's `--allow_paralogs` is OFF (its default, and `quantify`
  cannot pass it at all); a second flair arm WITH it is a separate, labelled arm, not a substitute.
- ⚠ Heavy runs strictly one at a time, foreground; outputs to `/mnt/linuxdisk/home/juanfraitu/bakeoff`.
- ⚠ This substrate developed the framing. Per trap 15 (`feedback_hold_a_substrate_back`), **a second family is
  held back** (`fam_MCL2_073244`, 64 copies) and is NOT to be looked at until the MCL1 arms are complete and
  written down. A conclusion that does not survive it is not a conclusion.

---
## Amendment 1 (2026-09-08, written BEFORE the flair arm was run)
⚠ **P3 as written cannot fire, and that is a defect in the instrument, not a result.** `derived_multi` asked
whether the same intron chain is asserted in ≥2 copies — but a chain is a tuple of GENOMIC coordinates, so two
copies at different coordinates can never share one. Measured on the two arms already run: `derived_multi` = 0
for **both** ours and StringTie, which is arithmetic, not a finding.

**P3 is restated, on measures that can actually distinguish the arms:**
- **P3a** — `transcripts spanning ≥2 copies`: a transcript whose span overlaps more than one copy interval has
  merged copies into one model. Prediction unchanged in direction: **flair and StringTie exceed ours.**
- **P3b** — `copies with ≥1 transcript`: how many of the 80 copies receive any model at all. Prediction:
  **ours covers more copies**, because a tool that never sees a read's secondary placements cannot emit a model
  at a copy whose reads were all given to a paralog.

⚠ Both are computed by the same tool-agnostic code for every arm. `derived_multi` is kept in the output for
honesty — it will read 0 nearly everywhere and that is expected, not suppressed.
