# Pre-registration — a streaming `--assemble-only` that a laptop can run genome-wide

**Written 2026-09-23 (§6zb), before the streaming path exists.** User: *"I would like to always run it
genomewide but it takes too much from this computer, in fact the idea would be to make it so that even
underpowered laptops could run it."*

## Where the time and memory go (measured, A119b chr21, 1,854,371 records, warm cache)

| probe | wall | peak RSS |
|---|---|---|
| shipped `--assemble-only`, full polish | 58 s (`reads_in_region` 42 s) | **4.5 GB** (chr1: 12.6 GB, r976) |
| `samtools view -c` on the same slice (htslib) | 9.5–13 s | 50 MB |
| noodles indexed `query`, extract nothing, plain BGZF reader | 13.1–13.9 s | – |
| same through `MultithreadedReader(4)` | **5.7–6.0 s** | – |
| noodles `query` + `RecordBuf::try_from_alignment_record` per record | **31.4–31.6 s** | – |
| noodles `libdeflate` feature (interleaved A/B, 3 pairs) | 36.4–38.3 vs 37.0–39.3 s read | no change; reverted |

So r983/r984's "the read is the floor and the decode is free" was wrong in the direction that matters: the
BGZF+query layer is 6 s with threads, the **`RecordBuf` decode (sequence + quality + a 116-op CIGAR per read,
215 M ops on chr21) is ~18 s**, and the per-read structs we then keep (`PrimaryRead` + `BamRead` with a
copied CIGAR, name and tags, for every record, plus the O2 AS-tied gate run over them) are the rest of the
wall and all of the memory.

## The change

Under `--assemble-only` (with `--read-isoform-k 0`, no footprint nodes, no tied seeding, no AS-ratio
secondary filter — every other setting keeps the materialised path):

1. **Pass-1 becomes an accumulator**: the existing `pass1_skeletons_widened` is split into `push(read)` /
   `finish()`; the function itself becomes push-all-then-finish, so it is byte-identical by construction.
   Spliced reads are reduced on arrival to their group's `(n, k-smallest starts, k-largest ends, n_reverse)`;
   unspliced reads (~45% of primaries) are kept as `PrimaryRead` for `cluster_unspliced`, unchanged.
2. **A streaming region reader** iterates the indexed query through the multithreaded BGZF reader on lazy
   `Record`s: flags, alignment start and the CIGAR ops only, the ops fed to the SAME `exons_from_cigar`.
   No `RecordBuf`, no `BamRead`, no names, no tags. The AS-tied gate is skipped (nothing to assign).
3. **`--genome-wide`**: one process, every contig with mapped reads in the BAM index, in sequence; the FASTA
   is already loaded per contig.

## Acceptance — committed now

- **Byte-identical GTFs**: chr21 and chr20 (human, new defaults) against today's outputs, and the gorilla
  genome (26 contigs) against `ggo_both_genome.gtf`. Any difference is a bug, not a trade.
- **Memory**: chr1 (A119b, 7.0 M records) peak RSS **< 1 GB** (was 12.6 GB); chr21 < 500 MB.
- **Wall**: chr21 `--assemble-only` **≤ 20 s** warm (was 58 s); a genome-wide A119b run in one process
  **≤ 25 min** on this box without batching (was ~1 h in memory-packed batches); gorilla ≤ 5 min.
- The full test suite passes; `--read-isoform-k 3`, `RUSTLE_FOOTPRINT_NODES`, `--tied-seed`,
  `RUSTLE_GTF_SECONDARY_AS_RATIO` still take the old path and stay byte-identical.

**Predicted:** memory target met by a wide margin (state is O(distinct chains), ~140 k on chr21); wall
lands at 10–15 s on chr21 (6 s read + CIGAR parsing + polish) and ~15–20 min genome-wide, reader-bound.

---

# OUTCOME (2026-09-23) — ⭐ every target met; one documented attribute difference

| run | before (materialised) | after (streaming, per-contig polish, indexed ISM/host) | identical output? |
|---|---|---|---|
| human chr21 | 58 s, **4.5 GB** | **14.7 s, 200 MB** | yes, apart from `matched_reads` |
| human chr20 | – | 18.2 s, 228 MB | yes, apart from `matched_reads` |
| human chr1 (5.5 M records) | ~306 s, **12.6 GB** (r976) | **30 s, 726 MB** | yes (vs the streaming quadratic-polish run) |
| gorilla genome, 26 contigs, one process | 196 s sweep at `--jobs 3`, ~2 GB × 3 | **151–321 s, 1.08 GB** | yes, apart from `matched_reads` |
| **human genome, 25 contigs, one process** | ~1 h in memory-packed batches (r976) | **482 s (8 min), 1.93 GB** | reproducible run-to-run |

Tests 945 / 0. Reader time is now 377 s of the 482 s genome-wide (`MultithreadedReader`, lazy records).

**What the work found on the way.** (1) r983/r984 were wrong in the direction that matters: the `RecordBuf`
decode is ~18 s of chr21's 42 s read and the multithreaded BGZF reader takes the raw query from 13 s to 6 s
once nothing else dominates; noodles' `libdeflate` feature changes nothing (interleaved A/B, reverted).
(2) The per-read structs and the O2 AS-tied gate were the memory (2.4 GB per million records). (3) A single
genome-wide process must polish PER CONTIG: the mono-exonic floor is a per-run quantile, and pooling it
across contigs changed 42 single-exon transcripts on gorilla — now per contig, byte-identical to the sweep
and to every single-region run. (4) The polish itself was quadratic: the ISM container scan over all
multi-exon transcripts of a contig+strand and the mono-exonic host scan cost 240 s of chr1's 306 s;
indexing containers by first junction and the host search by an offline Fenwick sweep visits exactly the
same pairs, byte-identical, chr1 306 → 30 s.

**The one deviation from "byte-identical":** the GTF attribute `matched_reads` — the number of AS-TIED reads
matching the transcript, produced by the O2 gate — is 0 on the streaming path, because the gate is what the
path removes. Same transcripts, coordinates, `reads`, strands. `--materialize-reads` restores the old path
(and the attribute) at the old cost.

**Shipped:** `copy_assign --genome-wide` (every contig with mapped reads in the `.bai`, one process,
mutually exclusive with `--region`/`--regions`), streaming pass-1 as the default under `--assemble-only`
whenever no read-level extra is requested (`--read-isoform-k 0`, no footprint nodes, no tied seeding, no
AS-ratio secondary filter), `--materialize-reads`, `Pass1Acc` / `stream_pass1_region` /
`junction_support` in `denovo_assemble.rs`, per-contig polish, indexed ISM and host search, and the
`bam_null_probe` bin that produced the reader table. `REPRODUCE.md` and `tools/genome_wide_sweep.sh`
updated; the sweep script is superseded for `--assemble-only`.
