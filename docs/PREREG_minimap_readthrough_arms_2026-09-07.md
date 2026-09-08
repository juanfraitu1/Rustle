# PREREG — can a minimap2 setting remove the read-throughs without destroying the real ones? (2026-09-07)

**Written before the arms.** md5 in `mcl_ann/adj/readthrough/PREREG_arms.md5`.

## Why this is now testable rather than arguable
§6gc gave a **truth set**: **30 of the 42** guarded gorilla read-through junctions are used by chimpanzee
reads (negative control 0/42). An aligner setting can therefore be scored by **how many CONSERVED junctions it
destroys**, not by whether it reduces the read-through count.

## Arms
Same 1,009,396 primary FLNC reads (those the shipped run placed on the three contigs), same reference index
(`GGO.splice.mmi`, full genome), only the mapping options differ:
- **A0 control** — `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (the shipped BAM's own line).
- **A1 no long join** — the same plus `-r 500,0`: the second value is minimap2's **long-join bandwidth**, the
  heuristic that patches distant anchors into one chain. Zero disables it. This is the only untested
  aligner-side lever that targets chaining across a gap.
⛔ `-G 50k` is not re-run: register 749 already measured it hitting 1 of 46 junctions while damaging 8.2 % of
gorilla transcripts. ⛔ `--junc-bed` is refused on principle (register 729).

## Metric
Per arm, per junction, the number of primary reads whose CIGAR carries an intron at exactly those coordinates
(± 5 bp), over the **30 conserved**, the **12 non-replicating** and the **6 the §6fw guard already removes**.

| # | prediction |
|---|---|
| **M1** | A0 reproduces the shipped junction support within 10 % on the 30 conserved junctions — if not, the re-alignment is not comparable and no arm may be read |
| **M2** | A1 **destroys conserved junctions**: ≥ 5 of the 30 lose all support |
| **M3** | A1 does **not** selectively remove the 6 guard-rejected ones — its effect on them is no larger, in proportion, than on the 30 |

## Interpretation fixed in advance
- M2 + M3 holding ⟹ **no aligner setting is the right place for this**; the guard stays downstream, and this
  closes the aligner line for good.
- M2 failing (A1 leaves the conserved set intact) ⟹ `-r 500,0` deserves a full evaluation as a default,
  which would be a genuine surprise and would be reported as one.

---
## AMENDMENT 1 (2026-09-07, before A1 ran)
`-r 500,0` is **rejected by minimap2**: `NUM1 (500) can't be larger than NUM2 (0)` — the long-join bandwidth
must be at least the chaining bandwidth. A1 therefore uses **`-r 500,500`**, the smallest legal setting, which
leaves the long-join heuristic **no reach beyond ordinary chaining**. That is the same intervention the
prediction was written about; no prediction changes.
