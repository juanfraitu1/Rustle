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
