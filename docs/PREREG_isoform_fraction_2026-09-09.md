# PREREG — `--min-isoform-fraction`: demote low-relative-depth transcripts, StringTie-style (2026-09-09, before the code change)

**Finding it follows from** (§6hq, the width-bias chase): our own `--gtf` emitted a 2-read, 10-exon
transcript at copy 4 spanning 15 kb PAST the copy's true 17.8 kb boundary. Verified NOT a mis-assembly:
flair's raw per-read transcripts (uncollapsed) show individual CCS reads genuinely aligning across the same
coordinates, and StringTie's own output at the identical locus reports the SAME readthrough — as a MINOR
alternative isoform (`STRG.55764.4/.10`, coverage 6.1/1.7) beside the dominant, correctly-bounded isoform
(`STRG.55764.1`, coverage 50.5). StringTie's defense is not "bundling" (we already group by locus via
`collapse_loci_groups`/`gene_tid`) — it is a per-locus RELATIVE-DEPTH classification from flow decomposition:
a chain far below the locus's dominant chain's coverage is reported as a minor isoform, not an equal call.
We collapse by exact chain with one flat `min_reads` floor — a 2-read and a 50-read chain at the same locus
are both just "a transcript."

**Rule (behind `--min-isoform-fraction <f>`, default 0.0 = off, byte-identical GTF).** `TranscriptRec.gene_tid`
already carries the locus-grouping key (from `collapse_loci_groups`, threaded through unconditionally, not
only under `--gtf-copy-set`). Before emission, per `gene_tid`: `group_max = max(n_reads)` over its
transcripts. Per transcript: `isoform_fraction = n_reads / group_max`; `low_confidence = isoform_fraction <
min_isoform_fraction` (only evaluated when the flag is > 0). When the flag is > 0, EVERY family transcript
gains `isoform_fraction "X.XXX"; low_confidence "true"/"false";` (report-only, non-family transcripts
untouched). When `--gtf-copy-set` is ALSO active: a `low_confidence` transcript is excluded from the
evidence-placement machinery (falls to the plain pass-through branch, tagged but never grouped/lifted/used
as evidence for a copy's span) — the direct fix for the width-bias mechanism (§6hq row 805/806).

| # | prediction | refuted by |
|---|---|---|
| P0 | flag at 0.0: GTF byte-identical to `ours_final3.gtf`/`ours_final3_nocs.gtf` on both species; assignments untouched | any byte differs |
| P1 | at 0.10 (StringTie's own `-f` default is 0.01, but HiFi depth here is orders of magnitude lower, so a coarser floor is expected to matter): `DN_chr16_15368428_10` (copy 4, isoform_fraction ≈ 2/max) is flagged `low_confidence "true"` | not flagged |
| P2 | at 0.10, **≤ 15 %** of human MCL0 family transcripts are newly low-confidence (this is meant to catch rare outliers, not a large fraction of real alternative isoforms) | > 30 % |
| P3 | re-running `bench/width_deficit.py`'s containment sweep with `--gtf-copy-set --min-isoform-fraction 0.10`'s GTF (excluding low-confidence transcripts from the span) at containment floor 0.5 (the LOOSE floor that let the readthrough through) recovers **most** of the improvement the 0.8/0.9 floor gave: median \|rel deficit\| for the previously-outlier copies (0, 2, 4, 8, 9, 12) drops by **≥ 70 %** | < 30 % recovered |
| P4 | the phantom table (§6hn/§6ho rule) on the new GTF: phantom transcripts **do not increase** (excluding low-confidence transcripts from evidence should not manufacture NEW phantoms, since a low-confidence transcript was already contributing spurious width, not spurious evidence) | phantoms increase |
| P5 | suite passes; excision sweep (§6hg) numbers on `ours_final2`/`ours_final3`-equivalent assignments are untouched (this flag only touches the `--gtf` emitter, never `copy_assign`'s own O2 verdicts) | any assignment-side number changes |
Human MCL0 primary; gorilla MCL1 secondary (report only, competitor comparison not repeated — testis
substrate rule stands). Default stays the user's decision after the measurement.
