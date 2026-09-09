# The fair comparison: `isoseq collapse` vs ours — cluster recipe

**Why this is the right comparison, and StringTie was not.** StringTie's long-read mode is a short-read
assembler adapted to long reads: it carries a flow model because short-read coverage is fragmentary, and it
will assert junction combinations no single read shows. Comparing against it invites the obvious objection
that we picked the wrong tool. **`isoseq collapse` is PacBio's own route for exactly this data** and does what
full-length reads actually need — group observed structures — so it isolates the only thing in dispute: **it
cannot say which copy an isoform came from, and we can.**

## What to run on the cluster

The aligned FLNC BAM is all that is needed; `isoseq collapse` takes an aligned, position-sorted BAM.

```bash
# input: the same aligned FLNC used everywhere here
#   gorilla NPIP three contigs : npip3.bam   (subset of GCA_029281585.2_flnc_mm.bam)
#   human chrY (DAZ)           : chrY.bam    (subset of A119b.t2t.bam)

samtools sort -@4 -o flnc.sorted.bam  npip3.bam
samtools index flnc.sorted.bam

# one command; writes flnc.collapsed.gff plus .abundance.txt and .read_stat.txt
isoseq collapse --do-not-collapse-extra-5exons flnc.sorted.bam flnc.collapsed.gff
```

⚠ Use `--do-not-collapse-extra-5exons` so 5′-truncated molecules are not merged into longer models. Without
it the comparison is unfair in OUR favour: it would collapse distinct observed structures that we keep.

If the unaligned FLNC is easier to reach, the full route is
`isoseq cluster2 flnc.bam clustered.bam` → `pbmm2 align --preset ISOSEQ` → `isoseq collapse`. The collapse-only
route above is closer to what we do and is the cleaner comparison.

## What to send back
`flnc.collapsed.gff` (and `flnc.collapsed.abundance.txt` if it is written).

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

**Prediction to write down before the run.** `isoseq collapse` should land **closer to us than StringTie did**
on chains and molecule agreement — it groups rather than infers — and should assert **few or no** junctions
below 2 molecules. If it matches or beats us on molecule agreement, that is the honest outcome and the right
response is to **use it for the grouping and keep only the attribution**, which is the contribution anyway.
