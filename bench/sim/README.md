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

---

# Genuine Starvation Fixture

**`bench/fixtures/sim_starved.*`** — controlled test for copy B DROPPED from default
`--vg` output while its bundle is thin-but-present (1 primary read survives so the
side-index can attach B's secondaries to its locus).  This is the target case for
Layer-2 homology coordinate-transfer recovery.

---

## Design

### Synthetic genome (`bench/fixtures/sim_starved.fa`)

One contig `chrSIM_STARVED`, 50,000 bp, generated from a deterministic pseudo-random
sequence (Python `random.Random(seed=42)`).

### Gene copies

| Copy | Region (0-based, half-open) | Exons                                  | Role |
|------|-----------------------------|----------------------------------------|------|
| A    | ~10,000–12,000              | [10000,10300), [10800,11100), [11700,12000) | Dominant — well assembled |
| B    | ~30,000–32,000              | [30000,30300), [30800,31100), [31700,32000) | Starved — absent from default --vg output |

Both copies have the same exon/intron structure: 3 exons × 300 bp, introns 500 bp +
600 bp, `+` strand.

### Homology and PSVs

Copy B's exon sequences are derived from copy A by applying 6 deterministic base
substitutions (2 per exon, placed at 1/3 and 1/2 through each exon):

- **6 PSVs** in 900 exonic bp
- **99.33% exon identity** between copies

PSV positions (0-based genome coordinates in copy B exons):

```
30100, 30150, 30900, 30950, 31800, 31850
```

The 6 PSVs are what makes `copyB_primary_000` prefer copy B as its primary alignment
over copy A (PSVs give it exact matches at B vs mismatches at A).  The remaining 70
secondaries at copy B carry the B-specific splice coordinates that Layer-2 needs.

### Read mix (genuine starvation)

| Read set | Count | Sequence | Primary lands at | Purpose |
|----------|-------|----------|------------------|---------|
| copyA_full_* | 40 | exact copy A 3-exon | copy A | dominant expression |
| copyA_skip_* | 10 | copy A exon1+exon3 skip | copy A | alt isoform evidence |
| copyB_primary_* | 1 | exact copy B 3-exon (has 6 PSVs) | copy B | bundle anchor |
| copyB_xmap_* | 20 | copy A sequence | copy A | B expression cross-mapped to A; secondary at B |

Total: 71 reads.

### Why this produces starvation

- `copyB_primary_000` uses copy B's transcript sequence (including PSVs). minimap2
  assigns it PRIMARY at copy B (MAPQ=60) and secondary at copy A.  This 1 primary
  read causes copy B's bundle to materialise internally — but 1 read is below rustle's
  assembly threshold, so no transcript is emitted.
- The 60 copy-A-sequence reads (`copyA_full_*` + `copyB_xmap_*`) and 10 skip reads
  all primary at copy A. They produce SECONDARY alignments at copy B (FLAG=256), each
  with a spliced CIGAR spanning copy B's exon coordinates.  These 70 secondary reads
  are the side-index signal that Layer-2 can use to reconstruct copy B's isoform.

## Verified alignment structure

```
samtools idxstats bench/fixtures/sim_starved.bam
#  chrSIM_STARVED  50000  142  0

samtools view -c -F 0x100 bench/fixtures/sim_starved.bam 'chrSIM_STARVED:10000-12000'
#  70   (copy A primaries)

samtools view -c -F 0x100 bench/fixtures/sim_starved.bam 'chrSIM_STARVED:30000-32000'
#  1    (copy B primary — bundle anchor)

samtools view -c -f 0x100 bench/fixtures/sim_starved.bam 'chrSIM_STARVED:30000-32000'
#  70   (secondary reads carrying copy B's splice coords)
```

Copy B primary read:
```
copyB_primary_000  FLAG=0  POS=30001  CIGAR=298M500N304M600N298M  MAPQ=60
```

Sample secondary reads at copy B (first 3):
```
copyA_full_000  FLAG=256  POS=30001  CIGAR=298M500N304M600N298M
copyA_full_001  FLAG=256  POS=30001  CIGAR=298M500N304M600N298M
copyA_full_002  FLAG=256  POS=30001  CIGAR=298M500N304M600N298M
```

## Empirical starvation verification

```bash
RAYON_NUM_THREADS=1 ./target/release/rustle -L bench/fixtures/sim_starved.bam \
  --vg --genome-fasta bench/fixtures/sim_starved.fa -o /tmp/starved_default.gtf 2>/dev/null
```

Expected output: **copy B is ABSENT** (no transcript between 30000-32000).
Copy A is assembled into 2 isoforms (full 3-exon + exon-skip).

```
awk '$4 > 20000' /tmp/starved_default.gtf
# (no output — copy B absent)
```

## Expected Layer-2 recovery

When Layer-2 (`--vg-layer2`) is implemented:

1. The side-index detects 70 secondary reads at copy B's locus (30000-32000).
2. These secondaries carry the splice CIGAR `298M500N304M600N298M`, reconstructing
   copy B's 3-exon structure.
3. The 6 PSVs at positions 30100/30150/30900/30950/31800/31850 confirm copy B origin.
4. Layer-2 emits copy B's full 3-exon isoform (previously absent from default --vg).

## Regenerating

```bash
bash bench/sim/make_starved_fixture.sh
```

All outputs are deterministic (fixed Python `random.Random` seeds: 42, 141).

## Files

| File | Description |
|------|-------------|
| `bench/sim/make_starved_fixture.sh`  | Shell driver (re-runnable) |
| `bench/sim/gen_starved_fixture.py`   | Python genome + read generator + alignment |
| `bench/fixtures/sim_starved.fa`      | Synthetic genome (chrSIM_STARVED, 50 kb) |
| `bench/fixtures/sim_starved.fa.fai`  | FASTA index |
| `bench/fixtures/sim_starved.reads.fa`| 71 HiFi-like reads (FASTA) |
| `bench/fixtures/sim_starved.bam`     | Sorted, indexed alignments |
| `bench/fixtures/sim_starved.bam.bai` | BAM index |
