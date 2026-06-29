# Simulated HiFi reads + extended-CIGAR: is *information* the limit, or *identifiability*?

Two follow-up questions, both answered empirically:

1. Does the gene-conversion / copy-assignment result survive a **realistic** long-read error
   model, or was it an artifact of the toy uniform-error generator?
2. Would re-running minimap2 with the **extended CIGAR** (`--eqx`: `=` matches, `X` mismatches)
   give rustle more discriminating information than the aggregate `NM` tag?

**Bottom line:** information is *not* the binding constraint. The error model and the CIGAR
encoding both leave the result unchanged; what moves it is **identifiability** — the density of
copy-diagnostic positions (PSVs), set by how recently the copies duplicated.

---

## 1. Realistic HiFi reads (badread `pacbio2021`) — result holds, and gets *stronger*

Installed via mamba (`badread` 0.4.2 + `pbsim3` 3.0.5, env `readsim`). Replaced
`gen_synthetic.py`'s uniform per-base error with badread's homopolymer-aware `pacbio2021`
model at ~99 % read identity, 4× depth, on the same 2-copy / 0.98-identity / distinct-isoform
genome with 80 planted gene-conversion molecules (5′ copy0 / 3′ copy1, switch at exon 3 =
9451–9750).

| reads | mosaic detector verdict |
|---|---|
| toy uniform-error (prior) | `breakpoint~9655 molecules=4 dispersion=0bp => CONFIRMED` |
| **badread realistic HiFi** | `breakpoint~9655 molecules=11 dispersion=0bp => CONFIRMED` |

Breakpoint **9655 ∈ the planted exon-3 switch**, **0 bp dispersion**, 36 diagnostic sites, all
37 family reads decisive. Realistic error did **not** obscure the PSV channel — the
copy-diagnostic SNPs (18/copy at 0.98) stand ~18× above the ~1 % random error, so they separate
cleanly. Higher depth gave *more* confirming molecules (11 vs 4).

**Methodology gotcha (real, worth recording):** badread randomly reverse-complements reads
(genomic dsDNA convention). Fed straight to the splice aligner this makes a **mixed-strand
family** → no joint graph → *"0 diagnostic sites — pileup fallback"* → PSV channel disabled and
no mosaic call. Real IsoSeq/MAS-seq cDNA reads are **strand-oriented** (the `isoseq refine`
step). Re-orienting badread's `-strand` reads to forward (one `revcomp` pass) restores the
single-strand family and the detection. Any long-read **RNA** simulation driven by a **DNA**
simulator needs this orientation step.

## 2. The controlled identifiability sweep — vary *only* PSV density on the same realistic reads

| identity | PSVs/copy | realistic-HiFi result |
|---:|---:|---|
| 0.980 | 18 | CONFIRMED · bp 9655 · 11 mol · 0 bp dispersion |
| 0.985 | 13 | CONFIRMED · bp 9655 · 11 mol |
| 0.990 |  9 | CONFIRMED · bp 9625 · 13 mol |
| 0.995 |  5 | CONFIRMED but bp drifts to 6802 · 3 mol · 10 diag-sites · precision degrades |
| 0.997 |  3 | **detector skipped — no diagnostic sites — floor reached** |

As the copies approach identical, detection **degrades then fails** — not from missing
information but because there are too few copy-distinguishing positions for any read to phase.
At 3 PSVs/copy no read can decisively own a copy; no error model and no CIGAR encoding
manufactures discrimination the genome doesn't contain.

Note the window is **coverage-dependent**: at 4×/120-reads the realistic run confirms down to
0.99, where the earlier low-coverage toy run failed at 0.985. The two knobs are **PSV density**
(identifiability, biological) and **coverage** (replication, buys margin) — *neither* is
"information per read."

## 3. Extended CIGAR (`--eqx`) — inert for rustle

Re-aligned the synthetic reads with and without `--eqx`:

- **NM is byte-identical** — 0 of 189 alignments differ. `--eqx` does not change the alignment
  or its edit distance; the strict-NM "ownership" signal (which copy a read fits best) is
  unchanged.
- The CIGAR *does* change — `400M5200N350M…` → `284=1X20=1X94=…` — i.e. `--eqx` localizes each
  mismatch as an `X` op. **But rustle never reads that distinction.**
  `bam.rs::extract_mismatches_vs_fasta` (line 526) treats `M | = | X` identically: for every
  aligned position it compares the read base against the in-memory reference FASTA and records
  `(ref_pos, query_base)` itself. rustle therefore already has full **per-position** mismatch
  localization — the exact information `--eqx` would surface — derived independently of the
  CIGAR.

So `--eqx` would feed the SNP/PSV/mosaic channel a byte-identical mismatch list. It helps tools
that lack the reference and must read mismatch positions off the CIGAR/MD tag; rustle keeps the
genome resident and reconstructs them directly.

## Verdict

> **Information is not the issue.** A realistic homopolymer-aware HiFi error model leaves the
> gene-conversion detection intact (in fact stronger), and the extended CIGAR adds nothing
> rustle doesn't already extract. The binding constraint is **identifiability** — copy-PSV
> density (with coverage as the second knob). Detection holds 0.98–0.99, degrades at 0.995,
> and vanishes at 0.997 (3 PSVs/copy), regardless of error model or CIGAR encoding. That is the
> measured ceiling, not a tooling gap.

Reusable harness: `gen_synthetic.py --error-rate 0` (clean molecules) → badread
`--error_model pacbio2021` → forward-reorient `-strand` reads → minimap2 `-ax splice:hq` →
`RUSTLE_VG_MOSAIC_ON`. See [NOVEL_COMBINATIONS.md](NOVEL_COMBINATIONS.md).
