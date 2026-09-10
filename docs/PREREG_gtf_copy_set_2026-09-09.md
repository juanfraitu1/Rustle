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
