# The fair comparison: `isoseq collapse` vs ours — cluster recipe

**Why this is the right comparison, and StringTie was not.** StringTie's long-read mode is a short-read
assembler adapted to long reads: it carries a flow model because short-read coverage is fragmentary, and it
will assert junction combinations no single read shows. Comparing against it invites the obvious objection
that we picked the wrong tool. **`isoseq collapse` is PacBio's own route for exactly this data** and does what
full-length reads actually need — group observed structures — so it isolates the only thing in dispute: **it
cannot say which copy an isoform came from, and we can.**

## What to run on the cluster — TWO arms, because they answer different questions

⚠ **Corrected 2026-09-08 (user: "we don't need to cluster first with cluster2?").** The canonical IsoSeq
route is `refine → cluster2 → pbmm2 align → collapse`. Skipping `cluster2` is a real, supported choice
("cluster-free"), not an oversight — but the two are **not the same experiment** and both are worth having.

### Arm 1 — cluster-free: `collapse` on the aligned FLNC (the apples-to-apples arm)
This is the one that matches what we do: group observed structures by genomic mapping, no consensus building,
no polishing. Our isoform set is exactly that, so this arm asks *is our grouping equivalent to PacBio's?*

```bash
samtools sort -@4 -o flnc.sorted.bam npip3.bam && samtools index flnc.sorted.bam
isoseq collapse --do-not-collapse-extra-5exons flnc.sorted.bam flnc.collapsed.gff
```
⚠ If the installed version refuses aligned FLNC as input, that route is unavailable in that build — use arm 2
and say so; do not work around it by pre-clustering and calling it cluster-free.

### Arm 2 — canonical: `cluster2` then align then `collapse` (the arm a reviewer expects)
`cluster2` clusters FLNC by similarity and builds a **polished consensus** per cluster. That is an extra
inference and error-correction step we do **not** perform, so this arm is not apples-to-apples — it asks the
different and equally fair question *is the standard PacBio pipeline better than ours?*

```bash
isoseq cluster2 flnc.bam clustered.bam                       # unaligned FLNC in
pbmm2 align --preset ISOSEQ --sort ref.fa clustered.bam mapped.bam
isoseq collapse --do-not-collapse-extra-5exons mapped.bam clustered.collapsed.gff
```

⚠ `--do-not-collapse-extra-5exons` in both arms: without it, 5′-truncated molecules merge into longer models
and the comparison tilts in OUR favour by collapsing distinct structures we keep.
⭐ Both arms write `*.read_stat.txt` mapping each read to its isoform — keep it, it is what lets the copy
comparison run on the same molecules.

## What to send back
`flnc.collapsed.gff` and/or `clustered.collapsed.gff`, plus the matching `*.read_stat.txt` and
`*.abundance.txt`. Say which arm each came from — the two are scored the same way but read differently.

## Scoring it here — one command
```bash
python3 bench/isoform_set_compare.py  <ours>.gtf  flnc.collapsed.gff  <bam>  copies.tsv \
        --labels ours,isoseq
```

It reports, on the same molecules and the same family copies: transcripts per set, distinct intron chains and
how many are shared, the fraction of spliced molecules whose exact chain is present, read-supported junction
recovery, junctions asserted below 2 molecules, transcripts over-extending their copy, and how many
transcripts in each set carry a copy attribute at all.

## The same table for StringTie, so the shapes are already known
| | ours | StringTie |
|---|---|---|
| transcripts overlapping a copy | 403 | 175 |
| distinct spliced intron chains | 305 | 164 (64 shared) |
| spliced molecules whose exact chain is present | **0.892** | 0.794 |
| read-supported junctions recovered | 0.435 | **0.586** |
| junctions asserted with < 2 molecules | **0** | 11 |
| transcripts > 2× their best copy | **47** | 59 |
| transcripts carrying a COPY attribute | **404** | **0** |

**Predictions to write down before the run.**
- **Arm 1 (cluster-free)** should land **closer to us than StringTie did** on chains and molecule agreement —
  it groups rather than infers — and should assert **few or no** junctions below 2 molecules.
- **Arm 2 (cluster2)** should emit **fewer** isoforms than arm 1 and than us, because consensus building
  merges near-identical structures, and should show **lower** exact-chain agreement with the raw molecules for
  the same reason. If arm 2 instead agrees with the molecules better than arm 1, our reading of what
  clustering does is wrong and should be corrected.
- If **either** arm matches or beats us on molecule agreement, that is the honest outcome and the right
  response is to **use it for the grouping and keep only the attribution**, which is the contribution anyway.
⚠ Neither arm can carry a copy attribute, so that column stays the point of the exercise in both.
