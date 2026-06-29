# Input formats (SAM/BAM/CRAM) + aligner ties

## Native SAM / BAM / CRAM input (SHIPPED)

rustle now accepts **SAM, BAM, and CRAM** input natively (pure noodles — no samtools/htslib needed).

* **Detection** is by magic bytes (`bam.rs::detect_format`): `CRAM` container magic → CRAM; BGZF/gzip
  (`0x1f 0x8b`) with first decompressed bytes `BAM\1` → BAM, else bgzipped SAM; plain text → SAM.
* **BAM** uses the existing fast multithreaded BGZF reader, **completely unchanged** (the hot path —
  BAM decode is ~40 % of wall-clock — is untouched).
* **SAM / CRAM** are transcoded **once** to a temporary BAM at startup (`bam.rs::transcode_to_temp_bam`)
  and the unchanged multi-pass pipeline runs on it. This is far less invasive than re-plumbing every
  read loop, and the pipeline reads the alignments several times, so one up-front transcode is also
  cheaper than repeatedly streaming a slow SAM/CRAM reader. The temp BAM is deleted on exit
  (`TempBam` drop guard).
* **CRAM** needs an indexed reference FASTA: `--cram-ref <ref.fa>` (defaults to `--genome-fasta`); the
  FASTA must have a sibling `.fai`. Missing reference → a clear, actionable error.

Why this is safe: BAM is a lossless container of the same alignment records, so the assembly is
**identical** to running on an equivalent BAM.

**Verified (chr19 / NC_073243.2):** running rustle on `GGO_19.bam`, a SAM of it, and a CRAM of it
(`--cram-ref GGO.fasta`) produced **byte-identical GTFs** (BAM-vs-SAM and BAM-vs-CRAM = 0 differing
lines; 2006 transcripts / 23565 exons each). Temp files auto-cleaned; 618-test suite green.

**Operational caveats (SAM/CRAM only — BAM is untouched):**
- The transcode writes a temporary BAM **the same size class as the input** (a multi-GB SAM → a
  multi-GB temp BAM); its size is printed to stderr. It goes to `$TMPDIR` (falls back to `/tmp`) — point
  `$TMPDIR` at a roomy disk for whole-genome SAM/CRAM.
- The transcode is currently **single-threaded** (the multi-pass assembly that follows is parallel; only
  the one-time up-front transcode is serial). For a large whole-genome CRAM this adds a few minutes.
- The `TempBam` drop guard deletes the temp on normal exit and on errors, but **not on `SIGKILL`/OOM** —
  a hard kill can leave a `rustle_input_<pid>_*.bam` behind; clean `$TMPDIR` if a run is force-killed.
- Byte-identity is verified on chr19 only; it holds by construction (BAM losslessly carries the same
  records), but a whole-genome equality check has not been run.

```bash
rustle -L reads.bam  -o out.gtf                       # BAM   (fast path)
rustle -L reads.sam  -o out.gtf                       # SAM   (native)
rustle -L reads.cram --cram-ref ref.fa -o out.gtf     # CRAM  (native; ref must be .fai-indexed)
```

## Can we stop minimap2 / winnowmap from making ties?

**Short answer: no — and that's a feature, not a gap.** A genuine tie means the read is
*information-theoretically identical* across the copies over its mapped span (no PSV in the read → no
distinguishing column). That is the identifiability floor (the `n_decisive = 0` / K-frontier case): no
aligner score, `--eqx`, scoring matrix, or winnowmap weighting can break it, because there is no
sequence difference to score. Forcing the aligner to "pick one" doesn't remove the tie — it **hides**
it behind an arbitrary primary flag, which is exactly the failure mode to avoid.

What you **can** (and should) do — expose ties and resolve them, abstain only at the floor:

1. **Expose all tied/near-tied placements** instead of hiding them:
   `minimap2 -ax map-hifi --secondary=yes -N <K> -p <low> --eqx` (Eichler used `-p0.5 -N20`; for full
   tie exposure `-N50 -p0.1`). Then every co-/near-optimal hit is in the BAM with its `AS`, and rustle's
   PSV layer adjudicates.
2. **Don't trust the arbitrary primary among MAPQ-0 paralogs** — rustle already doesn't. Copy
   assignment uses the **AS / decisive-margin gate** (`psv_linkage::assign_read_to_copy`), which
   *abstains on a tie* ("no copy beats the runner-up by the margin → `None`"), the PSV-space restatement
   of Eichler's AS ≥ 10, with **no 1/k** splitting.
3. **Reduce the *breakable* (spurious) ties**: `--eqx` for base-level PSVs; soft-mask repeats to cut
   repeat-driven multimaps; a merged minimap2(`-N50 -p0.1`)+winnowmap BAM maximizes exposed candidates
   in segdups (winnowmap changes *which* copy is primary / reduces mismapping, but identical-copy ties
   remain). Longer reads / less 5′ degradation span more PSVs → fewer ties (an input property, not an
   aligner knob).

So the lever is never "make the aligner not tie"; it is "expose the ties, resolve with PSVs, and
**abstain** exactly when the molecule carries no distinguishing information" — which is precisely the
thesis's identifiability boundary (and Canzar-aligned: resolve ambiguity, no arbitrary 1/k).
