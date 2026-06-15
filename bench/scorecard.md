# Genome-wide three-way scorecard — rustle-VG vs StringTie vs NCBI RefSeq

Generated 2026-06-14 on the chimpanzee (GGO) long-read dataset. Quantifies, genome-wide,
whether rustle in VG mode (`--vg --vg-layer2 -G stringtie.gtf`, the guided deployment mode)
(1) reproduces StringTie's transcripts (parity) and (2) recovers additional NCBI-RefSeq-
annotated transcripts that StringTie misses — with an illustrative "yield vs StringTie" that
can exceed 100%.

## How it was produced

Per-chromosome **serial** (whole-genome `-L` OOMs at ~18 GB); resumable via `.done` sentinels.

- `bench/gw_threeway.sh` — per contig (>1 Mb, with reads): slice BAM from `GGO.bam`, run
  StringTie (`-L`), run rustle-VG-guided (`--vg --vg-layer2 -G st.gtf --genome-fasta`), then
  three `gffcompare` runs: StringTie vs NCBI, rustle-VG vs NCBI, and parity (VG vs StringTie).
- `bench/gw_aggregate.py` — parses the per-chrom `.tmap`/`.stats` and prints the scorecard.

Reproduce:

```
bash bench/gw_threeway.sh      # ~27 min, writes /tmp/gw/
python3 bench/gw_aggregate.py  # prints the scorecard below
```

Inputs (paths hardcoded in the harness): `GGO.bam`, `GGO.fasta`, `GGO_genomic.gff` (NCBI
RefSeq annotation), `tools/stringtie/stringtie` (v3.0.1), `gffcompare` v0.12.10.

## Result (25 chromosomes, full coverage, ~27 min, 0 failures)

|                              | StringTie | rustle-VG |
|------------------------------|-----------|-----------|
| Total transcripts            | 68,157    | 71,015    |
| **Parity** (% of StringTie reproduced) | — | **100.00 %** (68,153/68,153), **0 regressions on every chromosome** |

**NCBI RefSeq recovery (unique annotated transcripts matched):**

| Tier            | ST-matched | VG-matched | VG-only (net-new) | ST-only (regr.) | illustrative yield |
|-----------------|-----------|-----------|-------------------|-----------------|--------------------|
| strict `=`/`c`  | 25,148    | 25,567    | **+419**          | **0**           | **101.7 %**        |
| broad `=`/`c`/`j` | 31,567  | 32,018    | **+451**          | **0**           | **101.4 %**        |

rustle-VG is a **strict superset of StringTie's annotation-corroborated output** (ST-only = 0
at both tiers) and recovers **419 NCBI-corroborated transcripts StringTie misses entirely**,
across **401 genes**, enriched in paralog / multi-copy families (the X chromosome alone
contributes 29; multi-copy genes RBMS3, MGA, PABIR2 and several LOC paralog families
contribute 2–3 each). Full list: `/tmp/gw/vgonly_strict.tsv`.

Speed: total ~27 min; no chromosome exceeded the 6-min rustle-VG flag (slowest = chr1
NC_073224.2 at 110 s for 481k reads).

## Caveats

- The **>100 % "yield" is illustrative** — true sensitivity caps at 100 %; the ratio exceeds
  it only because the denominator is StringTie's matched set, not the full annotation.
- The raw absolute Sn vs the full per-chrom annotation is low for *every* tool because the
  annotation contains far more loci than any single dataset expresses; the meaningful
  comparison is the **cross-tool delta** on the identical reference (denominator cancels).
- **Attribution:** the 419 VG-only transcripts combine *both* rustle's assembler beating
  StringTie *and* the VG secondary-alignment layer. The 4-way attribution
  (`bench/gw_attribution.sh` + `bench/gw_attribute.py`) adds a rustle guided-baseline
  (no `--vg`) per chromosome and splits the 419 — see below.

## 4-way attribution (isolating the VG layer's own contribution)

Adds rustle **guided-baseline** (`rustle -G st.gtf`, NO `--vg`) per chromosome, then splits the
419 strict (`=`/`c`) VG-only-vs-StringTie recoveries by whether the guided baseline already
finds them. Run after the three-way:

```
bash bench/gw_attribution.sh   # ~9 min, writes /tmp/gw/gd_*.gtf + gc_gd_*
python3 bench/gw_attribute.py   # prints the split below
```

Genome-wide (strict `=`/`c`, unique NCBI ref transcripts):

| Output                                   | NCBI-matched |
|------------------------------------------|--------------|
| StringTie                                | 25,148       |
| rustle guided-baseline (no `--vg`)       | 25,505       |
| rustle-VG (`--vg --vg-layer2`)           | 25,567       |

Splitting the **419** VG-only-vs-StringTie net-new recoveries:

| Source | count | meaning |
|--------|-------|---------|
| **assembler** | **357 (85 %)** | also found by the guided baseline — rustle's *assembler* beats StringTie; the VG layer was not needed |
| **VG-layer**  | **62 (15 %)**  | found **only** with `--vg --vg-layer2` — the secondary-alignment mechanism's own contribution |

So genome-wide, most of rustle-VG's advantage over StringTie is **rustle's base assembler**
(357 net-new), and the **VG secondary-alignment layer adds 62** NCBI-corroborated transcripts
that *both* StringTie *and* rustle's own assembler miss — across **54 genes, enriched in
paralog / multi-copy families** (LOC134756437 ×3; LOC101146541, LOC129527020/24, and other LOC
paralog families ×2), concentrated on chr1 (14), the X chromosome (8), and other read-dense
chromosomes. This is the VG mechanism's defensible genome-wide contribution; it is modest in
absolute count but real, strictly additive (0 regressions), and lands where the mechanism is
designed to help. Full lists: `/tmp/gw/vglayer_genes.tsv`, `/tmp/gw/assembler_genes.tsv`.

## Intron-chain (multi-exon) recompute + baseline-parity finding

The transcript-level `=`/`c` counts above include single-exon transcripts, which is where
rustle is *permissive* and the headline is weakest. Re-cut at **intron-chain level** (multi-exon
only — the rigorous gffcompare metric):

**Genome-wide intron-chain Sn/Pr vs NCBI** (from the per-chrom gffcompare stats):

| Tool | matching intron-chains | Sn | Pr |
|------|------------------------|------|------|
| StringTie | 23,304 | 24.39 % | 34.58 % |
| rustle-VG | 23,438 | 24.52 % | 34.27 % |

Near-identical genome-wide. **Multi-exon decomposition of the net-new recoveries:**
- Assembler net-new (357 transcript-level): **113 are multi-exon** (244 = 68 % single-exon —
  the permissiveness-driven part, the weak headline).
- VG-layer (62 transcript-level): **62 = 100 % multi-exon, intron-chain-corroborated.** The
  *entire* `--vg-layer2` secondary-alignment contribution is multi-exon real-structure recovery,
  no single-exon inflation — the pushback-resistant headline.

**Baseline-parity finding (de-novo rustle vs StringTie, no flags — the "baseline must equal
StringTie" question).** On chrY / chr19 / chr1, rustle de-novo already reproduces **92.8–95.6 %
of StringTie's intron chains** (Pr 84–90 %). The ~5–7 % divergence (1,384 rustle-only chains over
the 3 chroms) is **83 % NCBI-corroborated** (real annotated transcripts StringTie misses, not FPs).
**It is architectural, not threshold-gateable:** rustle's scalar thresholds already match StringTie's
`-L` defaults (`-c -f -j -m -a -g`); the only difference is `-s` (single-exon), where rustle is
*stricter* (4.75 vs 1.5) — matching StringTie there *explodes* single-exon output (wrong direction).
The only knob that shrinks the divergence (strict junction acceptance) removes just 17 % of it while
destroying ~5× as many *shared* real chains. So the divergence is the documented strict-early /
lean-downstream architecture (the inverse of StringTie), not a tunable — and it is mostly *real
recall*, so suppressing it to match StringTie would discard true positives. Achieving exact
baseline parity would require an architectural re-port (gating rustle's junction-acceptance and
downstream-retention behavior), not a threshold flag.

## PSV-linkage channel (`--vg-layer2-psv-linkage`) — genome-wide result

The PSV→junction linkage channel (`bench/gw_psvlink.sh` + `bench/gw_psvlink_aggregate.py`)
assigns a junction to the correct paralog copy when a single read spans both a PSV and the
junction. Genome-wide (25 chroms, ~23 min), comparing PSV-linkage ON (`pl`) vs the rest of
the VG layer (`vg`, Part A on):

| Tier | vg-matched | pl-matched | **PSV net-new** | dropped (pl<vg) |
|------|-----------|-----------|-----------------|-----------------|
| strict `=`/`c` | 25,567 | 25,567 | **0** | 0 |
| broad `=`/`c`/`j` | 32,018 | 32,018 | **0** | 0 |

**Correct and safe, but inert as an additive channel.** No chimeras (0 transcripts > 1.5 Mb —
the per-locus-frame fix holds genome-wide), 0 regressions, 100% StringTie parity preserved.
But **PSV net-new = 0**: the channel adds no NCBI-corroborated transcripts over the VG layer.

Why: PSV-linkage emits a PSV-*validated subset* of what Part A already emits — both build the
same per-copy chains from the same reads (`build_exons_from_chain`); PSV-linkage merely
restricts to reads whose alleles confirm the copy. So as an *additive* channel alongside Part A
(the design choice), it can only produce chains Part A already produced → union-by-chain dedups
them → net contribution ≈ 0. The mechanism *is* real (in isolation, with Part A off via
`RUSTLE_VG_LAYER2_NO_MULTI_ISOFORM=1`, it recovers a copy-specific skip isoform resolvable only
via the linked PSV — harness leg 9, fixture `sim_psvlink`). Its value (precision-safe per-copy
assignment, no phantoms) would only materialize if it **replaced or filtered** Part A rather
than augmenting it. Kept **default-off**, like Part B; the engine + `PsvCertificate` are reusable.
