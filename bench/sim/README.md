# Simulated IsoSeq Paralog Fixture

Ground-truth fixture for testing Rustle's **Layer-2 paralog recovery**: the ability to
recover a lowly-expressed ("starved") paralog copy whose reads multi-map primarily to a
well-expressed homologous copy.

---

## Design

### Synthetic genome (`bench/fixtures/sim_paralog.fa`)

One contig `chrSIM`, 50,000 bp, generated from a deterministic pseudo-random sequence
(Python `random.Random(seed=42)`).

### Gene copies

| Copy | Region (0-based, half-open) | Exons                                  | Expression |
|------|-----------------------------|----------------------------------------|------------|
| A    | ~10,000–12,000              | [10000,10300), [10800,11100), [11700,12000) | 40 reads (full 3-exon) + 10 reads (2-exon skip) |
| B    | ~30,000–32,000              | [30000,30300), [30800,31100), [31700,32000) | 5 reads (full 3-exon only) |

Both copies have the same exon lengths (3 × 300 bp) and the same intron lengths
(500 bp intron 1, 600 bp intron 2), on the `+` strand.

### Homology and PSVs

Copy B's exon sequences are derived from copy A's by applying ~1.5% base substitutions
at deterministic positions (fixed seed = 43), yielding:

- **14 PSVs** in 900 exonic bp
- **98.44% exon identity** between copies

The 14 PSV positions (0-based) in the genome:

```
30074, 30149, 30808, 30832, 30893, 30970, 31033, 31058, 31083, 31785, 31837, 31847, 31917, 31991
```

These PSVs are what Layer-2's PSV machinery uses to distinguish the two copies during
constrained flow decomposition.

### Reads (`bench/fixtures/sim_paralog.reads.fa`)

55 reads total, generated as **exact spliced transcript sequences (0 sequencing errors)**,
forward-oriented (IsoSeq-refine convention).  HiFi accuracy is ~99.9%, so error-free reads
are a valid approximation that avoids a heavy simulator dependency.

| Read set      | Count | Isoform                    |
|---------------|-------|----------------------------|
| copyA_full_*  | 40    | A: all 3 exons (900 bp)    |
| copyA_skip_*  | 10    | A: exon 1 + exon 3 (600 bp, skips exon 2) |
| copyB_full_*  | 5     | B: all 3 exons (900 bp)    |

No read simulator was installed.  The Python script `gen_paralog_fixture.py` emits
reads directly.

### Alignment (`bench/fixtures/sim_paralog.bam`)

Command:
```
minimap2 -ax splice:hq -uf --secondary=yes -N 5 --MD -Y \
    bench/fixtures/sim_paralog.fa \
    bench/fixtures/sim_paralog.reads.fa \
  | samtools sort -o bench/fixtures/sim_paralog.bam
samtools index bench/fixtures/sim_paralog.bam
```

---

## Verified multimapping structure

```
samtools idxstats bench/fixtures/sim_paralog.bam
#  chrSIM  50000  110  0
#  *       0      0    0

samtools view -c -f 0x100 bench/fixtures/sim_paralog.bam
#  55
```

Breakdown by origin and landing region:

| Count | Type | Origin copy | Lands in region |
|-------|------|-------------|-----------------|
| 50    | PRIMARY   | A | copy A (pos ~10001) |
| 50    | SECONDARY | A | copy B (pos ~30001) |
| 5     | PRIMARY   | B | copy B (pos ~30001) |
| 5     | SECONDARY | B | copy A (pos ~10001) |

**Interpretation:** The 5 copy B reads (starved copy) get their PRIMARY alignment to copy B
and a SECONDARY alignment to copy A — exactly the cross-copy secondary signal that Layer-2's
side-index ingests.  All copy A reads reciprocally produce secondary alignments into the copy B
region.  The 14 PSVs in the exons allow the PSV machinery to distinguish origin.

Spliced alignments confirmed (CIGAR contains `N` operators):

```
copyA_full_000  FLAG=0    POS=10001  CIGAR=301M500N299M600N300M
copyB_full_000  FLAG=256  POS=10001  CIGAR=301M500N299M600N300M  (secondary; primary at 30001)
```

---

## Expected Layer-2 recovery

Layer-2 should recover **copy B's full 3-exon isoform** from the secondary side-index.
Without Layer-2, copy B's bundle receives only its 5 primary reads (coverage too low for
assembly alone under default thresholds).  With Layer-2, the secondary signal from copy A's
50 reads is borrowed into the copy B graph, and the PSVs confirm the copy-B-specific path,
enabling emission of the 3-exon isoform.

---

## Regenerating

```bash
bash bench/sim/make_paralog_fixture.sh
```

All outputs are deterministic (fixed Python `random.Random` seeds).

---

## Files

| File | Description |
|------|-------------|
| `bench/sim/make_paralog_fixture.sh` | Shell driver (re-runnable) |
| `bench/sim/gen_paralog_fixture.py`  | Python genome + read generator + alignment |
| `bench/fixtures/sim_paralog.fa`     | Synthetic genome (50 kb) |
| `bench/fixtures/sim_paralog.fa.fai` | FASTA index |
| `bench/fixtures/sim_paralog.reads.fa` | 55 HiFi-like reads (FASTA) |
| `bench/fixtures/sim_paralog.bam`    | Sorted, indexed alignments |
| `bench/fixtures/sim_paralog.bam.bai` | BAM index |
