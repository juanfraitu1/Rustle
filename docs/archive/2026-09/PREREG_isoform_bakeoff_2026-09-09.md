# PREREG — the ISOFORM-level (lift-aware) hard-locus bakeoff (2026-09-09, before the metric runs)

**Why.** `tool_bakeoff.py` matches a read's intron chain to a transcript at the read's primary locus. The
evidence-placed GTF (`ours_copyset.gtf`) emits an isoform once, at its evidence copy or with a copy SET, so the
exact-locus rule turns every read whose primary sits at a dropped phantom into `derived_none` (§6hn P5).
**Rule (`bench/isoform_bakeoff.py`).** For every hard molecule (gate row with a primary inside a copy span) with
≥ 1 intron: the read's chain at its primary copy P, and its lifts to every other copy c (± 5 bp at each boundary).
A tool CARRIES the molecule if any transcript's chain equals the read's chain at P or its lift at some c. The
tool's ADDRESS SET for that isoform = the copies where it emits the chain, plus, for our GTF, the `copies`
attribute (evidence copies, or the undecided set). Tools: ours (copy-set GTF), ours-raw (`ours_final2.gtf`),
flair, StringTie, isoseq (fuzz 0 for all; fuzz 5 as a secondary table).

| # | prediction | refuted by |
|---|---|---|
| P1 | lift-aware carried on the hard set: **ours ≥ 0.85** (the copy-set GTF loses nothing once matching is lift-aware), ours-raw within 1 pt of ours; flair, StringTie **≤ 0.60**; isoseq **≤ 0.80** | ours < 0.80, or any competitor ≥ ours |
| P2 | the 262 O2-assigned molecules that a tool carries: the tool's address set contains O2's copy — **ours ≥ 95 %**; ours-raw ≥ 90 %; competitors reported (their address is the aligner's primary) | ours < 90 % |
| P3 | for the molecules O2 leaves undecided that a tool carries: **competitors commit to exactly ONE address in ≥ 95 %** (a transcript has one locus); ours reports a SET (≥ 2 addresses incl. `outside`) in **≥ 95 %** | either < 85 % |
| P4 | per copy: ours-only vs each competitor ≥ 1, competitor-only ≤ 1 (§6hh P4, now lift-aware) | as §6hh |
Human only; support-matched (chain multiplicity ≥ 2) as a secondary table.

## Outcome (2026-09-09) — `bench/isoform_bakeoff.py` (`bakeoff/human/hard/isoform_bakeoff.txt`); per-tool phantoms via `bench/isoform_copy_lift.py`
Hard spliced molecules with a primary in a copy span: 3,633 (contested 696 = 207 assigned / 489 undecided).
| # | verdict |
|---|---|
| P1 | ⚠ ours **0.848** (boundary miss of ≥ 0.85; ours-raw 0.849 within 1 pt ✓); flair 0.596 ✓, StringTie 0.524 ✓, isoseq 0.756 ✓. On the contested stratum isoseq 0.889 > ours 0.661 (support-matched 0.910 vs 0.877): **lift-aware "carried" rewards a tool for emitting the same isoform at several copies** — the very duplication the next table measures |
| P2 | ✓ of the O2-assigned molecules a tool carries, the address set contains O2's copy: **ours 128/128 = 100 %**, ours-raw 95 %; flair 78 %, StringTie 82 %, isoseq 86 % — the competitors put 14–22 % of the certified isoforms at another copy |
| P3 | ⛔ refuted in letter: competitors report ONE address for only 35 % (flair) / 14 % (StringTie) / 26 % (isoseq) of the undecided molecules — not because they represent uncertainty but because they emit SEPARATE transcripts at several copies; ours reports a set for **95 %** (18 single-address cases are undecided reads carried by an evidence-backed shared isoform) |
| P4 | ✓ flair (ours-only copies 15, 18), ✓ StringTie (0, 15); ⛔ isoseq (no ours-only copy; isoseq emits at all 25) |

### The table that answers "other assemblers struggle here" (same lift rule on every tool's GTF; evidence = unique mapper or certified read)
| tool | multi-intron transcripts in copies | multi-copy isoform groups | groups with a PHANTOM copy | phantom transcripts | isoforms on abstaining reads only |
|---|---|---|---|---|---|
| **ours, copy-set GTF** | 534 | 45 | **2 (4 %)** | **2** (both lifted placements, definitional) | **40 — each emitted ONCE with `copies "A,B[,outside]"`** |
| ours, raw | 542 | 53 | 13 (25 %) | 14 | 39 (aligner's address) |
| flair | 733 | 14 | 4 (29 %) | 4 | 76 (aligner's address) |
| StringTie | 433 | 34 | 13 (38 %) | 15 | 20 (aligner's address) |
| isoseq collapse | 2,699 | 60 | **34 (57 %)** | **45** | **200** (single-read transcripts at whichever copy won the tie-break) |
Phantoms sit beside a sibling ≥ 0.985 in 87–100 % of cases for every tool. The claim that survives, stated
exactly: **every assembler emits phantom transcripts at near-identical copies and gives coin-toss isoforms a
single arbitrary address; ours is the only GTF that (a) drops the phantoms it can identify, (b) places certified
isoforms at O2's copy (100 % of carried certified molecules), and (c) says "A or B" when the reads cannot tell.**
