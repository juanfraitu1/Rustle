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

## Outcome (2026-09-10) — `bakeoff/human/ours_iso10b` (sum-denominator fix), `ours_p0iso*` (byte-identity)
| # | verdict |
|---|---|
| P0 | ✓ flag at 0.0 (implicit and explicit) byte-identical GTF `888794bf`, assignments `8a057f68`, on both the max- and sum-denominator builds |
| P1 | ⚠ **failed with the originally-specified MAX denominator** (`DN_chr16_15368428_10` scored `isoform_fraction 0.167`, not flagged — its own locus group's max chain is only 12 reads, not the ~50 StringTie sees, because our exact-chain collapse FRAGMENTS one locus's reads across many near-identical chains: this locus has 41 total reads split across 9 chains, largest 12). **Fixed by summing, not maxing, `n_reads` within a `gene_tid` group** (register row 807) — with the sum denominator (41), the transcript scores 0.049 and IS flagged at floor 0.10 |
| P2 | ⛔⛔ **REFUTED, badly, with the sum denominator: 514/686 = 75.0 % of family transcripts flagged low-confidence** (predicted ≤ 15 %, refuted at > 30 %) |
Suite 871 / 0 / 11.

### Why P2 fails, and what it reveals (the important part)
The locus that contains the readthrough transcript has 9 chains at 2, 2, 2, 2, 3, 5, 11, 12, 2 reads (sum
41). The readthrough chain (2 reads) is **not a statistical outlier in DEPTH** among its siblings — five of
the other eight chains have the SAME or fewer reads and are perfectly normal alternative isoforms (different
5' start points from natural read-length/TSS heterogeneity, not artifacts). **`isoform_fraction`
(whole-transcript relative depth) cannot separate a readthrough from a legitimate low-support alternative
isoform when both have similarly low absolute read counts** — which is the common case once a locus's true
depth is spread over many near-identical exact chains (measured here: even the BEST-supported chain at this
locus is only 12/41 = 29 % of the locus total, so a strict-enough floor to catch the 2-read readthrough
catches most of the 2–5-read legitimate variants too, and a loose-enough floor to spare them misses the
readthrough). **The actual distinguishing signal is at the JUNCTION level, not the transcript level**:
StringTie's flow decomposition scores the SPECIFIC long-range junction (2 reads cross it) against the many
reads that are present at the same upstream position but TERMINATE before it (the well-supported chains'
reads) — a per-junction relative test, not a per-transcript one. This is structurally the SAME test §6ft's
read-through certificate already runs at O2's read-star stage (a giant intron with too few supporting
molecules), just with the WRONG absolute threshold for this case: that gap is ~13.5 kb, well under the
existing 50 kb floor, so even the O2-side guard would not have caught it.

### Reading
`--min-isoform-fraction` is implemented, tested, byte-identical at 0 (default), and does exactly what its own
definition says — but that definition (relative TRANSCRIPT depth) is the wrong axis for this pipeline's
exact-chain-collapse fragmentation pattern. **It is not the fix for the readthrough problem.** The correct
next step is a JUNCTION-level relative-support test in the assembler itself (does this specific junction's
read count fall far below the read count of the reads present at the same donor site that do NOT take it) —
a materially different, larger piece of work than this flag, scoped separately. `--min-isoform-fraction`
stays in the codebase, default off, as a plain (if imperfect for THIS purpose) relative-depth report/filter;
not recommended as a default and not claimed to solve the readthrough case.
