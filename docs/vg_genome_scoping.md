# VG genome scoping (fast `--vg` on region-scoped BAMs)

`rustle --vg --genome-fasta <genome>` loads a reference FASTA to phase reads to
gene-family copies. On a region-scoped BAM (e.g. one chromosome, or a single
locus), loading the *whole* genome is the bottleneck — for the gorilla genome it
is ~6 minutes / several GB of RAM just to phase a 60 kb region.

There are two ways this is handled; the first is automatic.

## 1. Automatic scoping (default, in-tool)

In `--vg` mode rustle now loads **only the contigs that actually have mapped
reads in the BAM**, not the whole genome. It does this without scanning the BAM:

- it reads the BAM index (`.bai`) to find which reference sequences have
  `mapped_record_count > 0`, then
- seeks straight to those contigs in the FASTA via the FASTA index (`.fai`).

A region slice of `chrY` therefore loads ~66 MB instead of 3.6 GB.

**Result (DAZ region):** ~6 min → **~1.2 s**, identical output. (`[VG] Loading
genome FASTA scoped to 1 BAM contig(s) ...` appears in stderr.)

**Requirements (both are standard and cheap to create):**

```bash
samtools index  reads.bam     # -> reads.bam.bai   (used to find contigs with reads)
samtools faidx  genome.fa     # -> genome.fa.fai   (used to seek to those contigs)
```

**Safe fallback:** if either index is missing/unreadable, or the BAM header lists
contigs not in the FASTA, rustle transparently loads the whole genome exactly as
before — scoping never makes a run *worse*, it only speeds up the common case.

**Note:** a genome-wide BAM has reads on every contig, so scoping naturally loads
everything (no change). The win is specifically for region/locus-scoped BAMs.

## 2. Manual pre-subset helper (optional)

`bench/subset_genome_for_vg.sh` writes a small FASTA containing only the BAM's
mapped contigs, for pipelines that want to materialize a reduced reference (or
for other tools):

```bash
bench/subset_genome_for_vg.sh  genome.fa  reads.bam  subset.fa
rustle --vg --genome-fasta subset.fa -L reads.bam -o out.gtf
```

It uses `samtools idxstats` + `samtools faidx` (no full-genome scan); on DAZ it
produces a 66 MB `subset.fa` from the 3.6 GB genome in ~1 s.

## Capacity-confidence + abstain (related)

With the genome loaded, `--vg` apportions multi-mapping reads across copies and
emits, per transcript:

- `capacity_confidence "X"` — fraction of the copy's coverage that comes from
  reads decisively anchored *to this copy* (1.0 = fully anchored, 0.0 = phantom).
- `low_confidence "true"` — set when `capacity_confidence < RUSTLE_VG_ABSTAIN_FLOOR`
  (default 0.05). The transcript is kept (not dropped), just flagged.

Example (DAZ): the real DAZ1 copy reads `capacity_confidence "1.000"`; the
unexpressed DAZ3 copy reads `capacity_confidence "0.000"` + `low_confidence
"true"`.
