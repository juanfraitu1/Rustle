# PREREG — emit the GTF O2 believes: evidence-placed transcripts and copy SETS (2026-09-09, before the post-processor runs)

**Rule (post-processor `bench/gtf_copy_set.py`, human MCL0, on `ours_final2.gtf`; the binary is untouched in
this step).** Member reads of a transcript = primaries with the same intron chain (the collapse rule).
Evidence at copy c for an isoform GROUP (§6hl lift groups) = ≥ 1 member read that is a UNIQUE mapper at c
(not admitted by the AS-tied gate) or read-level `assigned` to c. Then, per group:
1. **Evidence at ≥ 1 copy:** emit one transcript per evidence copy — the existing transcript where one sits
   there; where a certified read's copy has no transcript in the group, LIFT the chain to that copy through the
   copy-to-copy alignment (every exon boundary within ± 5 bp) and emit it there (`placed_by "assigned_read"`,
   `lifted_from`). Transcripts of the group at copies WITHOUT evidence are dropped (the phantoms); their
   abstaining reads join the group's `reads_undecided` and `copies_undecided` (the union of those reads'
   candidate copies from the read-star dump).
2. **No evidence anywhere (abstaining reads only):** emit ONE transcript, at the aligner's majority primary
   locus, `copy_status "undecidable"`, `copies "A,B,…"` = the union of its reads' candidate copies.
3. Attributes on every family transcript: `evidence_unique "c:n,…"`, `evidence_assigned "c:n,…"`,
   `copies` (evidence copies or the undecided set), `pooled_copy`/`pooled_margin` from §6hm (report only,
   never a coordinate). Non-family transcripts pass through unchanged.

| # | prediction | refuted by |
|---|---|---|
| P1 | the 14 phantom copies (§6hl) are gone: re-running `isoform_copy_lift.py` on the new GTF finds **0** multi-copy groups with an evidence-less copy | ≥ 1 |
| P2 | the 39 single-copy abstaining-only isoforms are emitted once with a `copies` set of **≥ 2** members in ≥ 95 % | < 80 % |
| P3 | the 43 shared groups keep every evidence copy (no shared isoform loses an address) | any loss |
| P4 | family transcripts **670–697** (684 − 14 phantoms + ≤ 13 lifted placements for certified reads whose transcript sat at another copy) | outside |
| P5 | hard-locus bakeoff on the new GTF (fuzz 0): ours carried on the hard set within **± 1 point** of 0.852; copy attribution of the 262 O2-assigned rises from 166/182 to **≥ 95 %** of the carried ones; transcripts spanning ≥ 2 copies still **0** | carried drops > 2 points, attribution < 90 %, or any spanning transcript |
Human only. If the object holds, the same rule is ported into the binary's `--gtf` emitter (a later, separate step).

## Outcome (2026-09-09) — `bench/gtf_copy_set.py` → `bakeoff/human/ours_copyset.gtf` (684 → 676 family transcripts: 488 kept, 13 phantoms dropped, 1 duplicate undecided dropped, 6 lifted placements added (3 lifts failed), 40 undecided isoforms emitted once)
| # | verdict |
|---|---|
| P1 | ⚠ 2, not 0: both are LIFTED placements (copy 12) that the lift script cannot see evidence for — its evidence rule wants chain-matching reads AT the locus, and a certified read's primary is at the copy it was lifted from. Definitional; the 13 original phantoms are gone |
| P2 | ⛔ as pre-registered (read-star candidate union): 13/40 had a set ≥ 2 — 27 undecided isoforms are made of single-candidate reads whose tie partner lies OUTSIDE the catalog. **Amended before re-running:** the copy set = the units at the read's AS-TIED placements (the gate's tie set, from the BAM) plus `outside` when a tie partner overlaps no unit. Result: **40/40 have a set ≥ 2** (27 include `outside`, 11 are catalog pairs, 2 triples) |
| P3 | ✓ the 43 shared groups keep every evidence copy |
| P4 | ✓ 676 (670–697) |
| P5 | ⚠ mixed, and the metric is the problem: derived carried on the hard set 0.835 (−1.7 pts), on the contested stratum 0.673 (−5.6 pts); attribution 166/181 = 92 % (unchanged); spanning transcripts 0. The derived rule matches a read's chain to a transcript AT THE READ'S PRIMARY LOCUS, so a phantom dropped at copy 7 turns its reads into `derived_none` even though the isoform is emitted at copy 8 with `copies_undecided "7"`, and a lifted placement is never matched. The exact-locus metric cannot see copy sets or O2 placements — hence the lift-aware metric of the next PREREG |
