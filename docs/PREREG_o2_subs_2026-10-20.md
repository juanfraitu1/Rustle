# Pre-registration: O2 by substitutions only, on fresh substrates (2026-10-20)

**User (2026-10-20):** "how about the analysis of CIGAR and bam tags such as md, nm, de, ms etc" — and, on the proposal to
test the rule it revealed on fresh data: "yes please".

## The rule

For a read tied within 30 intron-free points on two or more copies of one family (`o2window.py` window 30, stitched
alignments excluded, each copy by its highest-`ms` alignment): count each alignment's **substitutions** over the read bases
every tied alignment aligns (`alignment_events` "S" events in the common window: what the MD tag counts, less indels). The
copy with the **unique fewest** takes the read; a shared minimum keeps minimap2's primary (abstain). Indels, clips, AS,
ms, NM and de are not used.

**Why (seen on HG002, descriptive, `PREREG_o2_smn_showcase` "Tags"):** substitutions pick the true copy 2,603 / 1 on SMN,
5,482 / 5 on THOC3, and decide 354 / 0 of the reads the Pareto set abstained on; AS, ms, NM and de are misled by indels.
Register 1343's preamble said the same rule was a coin flip outside SMN/THOC3 on the 10-AS domain, with the edit answer
key that over-penalized clips; this test settles it.

## Substrates, never used for this rule

- **S1, one genome:** simulation round 27, `families27.json` (sha1 ebedbf33389e, 30 of the 100 qualifying non-PAR pool
  families no earlier round drew, seed 2027, 83 copies), the same build as round 26 (`heldout26/run_cell.sh` with
  families27, seeds 1701-1703, profiles A119b / public, INTACT / CROSS / SILENT). Truth = the simulated copy. Note: 3 drawn
  families (DUSP22, GTF2H2C, HSPA6) appeared in the HG002 tag table above; the rule was not tuned on them.
- **S2, real reads, another donor:** A119b testis fragments, seed 2 (`o2frag_a119b2/frag.bam`, 33,840 fragments of 11,280
  reads, `--eqx`); truth = the source read's copy (`truth.tsv`), as in register 1343. Known bias: a systematic full-read
  error counts against the rule. Qualifying families only (`famident/qualify_nonpar.json`).
- **S3, replication, seen families:** gorilla KB3781's 143 held-out families (`personal/window/g`), intron-free key from
  his assemblies (answer key only). Reported, no verdict.

## Scoring (`o2window.py pick-bars`, column `subs_pick`)

Scored read: determined truth among the window's copies. Per substrate: reads, families, the rule's decisions, right,
minimap2's primary right, fixes (rule right, primary wrong), breaks (the reverse), families fixing more than breaking and
the reverse; by band (first tied within 10 / 20 / 30).

## Bars (S1 and S2 each; as register 1343's)

- **P:** >= 300 scored reads from >= 11 families, minimap2 wrong on >= 50.
- **G:** accuracy >= minimap2's + 2.0 points.
- **S:** fixes >= 5 × breaks.
- **F:** >= 11 families fixing more than breaking, no more than a fifth as many (rounded up) the reverse.

**Verdict per substrate:** PASS = all four; FAIL = P holds and another fails; UNDERPOWERED = P fails. S1 and S2 are not
pooled; a claim needs both to PASS.

## Rules

Code committed with this file; each substrate prepared and scored once; bugs → amendment + rerun; foreground, one heavy
job at a time; audit on a frozen copy before reporting; human and gorilla never pooled.

## Amendments

(none yet)

## Results (code a45c91bc; each substrate scored once; frozen copy `personal/window/audit_frozen_subs_a45c91bc/`)

| substrate | scored / families | minimap2 right | substitution rule right | fixes / breaks | families better / worse | verdict |
|---|---|---|---|---|---|---|
| **S1** one genome, round 27 (18 worlds) | 23,145 / 27 | 21,125 | 21,125 | 0 / 0 | 0 / 0 | **FAIL** (G, F) |
| **S2** A119b fragments, seed 2 | 3,646 / 107 | 3,259 | 3,251 | 4 / 12 | 2 / 11 | **FAIL** (G, S, F) |
| S3 gorilla held-out (replication) | 2,711 / 29 | 2,464 | 2,439 | 1 / 26 | 1 / 6 | (no verdict; bars fail) |

**S1, read by read.** Of 23,145 scored reads the substitution count decides 18,713, and on **every one** it picks minimap2's
own copy (right on 18,710; both wrong on 3). The other 4,432 are substitution ties — 3,086 with equal counts on every copy —
and there minimap2 is right on 2,415 and wrong on 2,017: a coin toss. Every minimap2 error in S1 is a read whose bases
carry no substitution that tells its copies apart. Bands: first tied within 10, 11,501 scored, minimap2 right 9,483;
10-20, 7,394 / 7,392; 20-30, 4,250 / 4,250.

**S2.** Fixes 4, breaks 12, in 107 families (2 better, 11 worse). The truth is each fragment's source read's own call
(register 1343's known bias: a systematic full-read error counts against the rule).

**S3.** 1 fix, 26 breaks: 23 on reads carrying his other parent's alleles (SM233 10, SM502 8, SM556 5), 3 on reads whose
copy equals the reference. SM392's stitched reads abstain (one plain candidate, nothing to compare).

**Reading (pending the audit).** The HG002 SMN tag result (seen data) does not generalize. In one genome the substitution
count adds no correction — minimap2's primary already is the fewest-substitution copy whenever one exists — but it does
say exactly which primaries are coin tosses. In real donors the count is misled by the donor's own alleles, as every
allele-blind statistic is.
