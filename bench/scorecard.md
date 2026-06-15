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

## Genome-wide StringTie-EXACT floor (`-e -G st.gtf`) — the baseline anchor

The "baseline must equal StringTie exactly, then VG adds on top" question, answered
genome-wide. `bench/gw_eonly_parity.sh` + `bench/gw_eonly_aggregate.py` run rustle in
**eonly-guided** mode (`-e -G st.gtf`: emit ONLY transcripts matching the guide) per
chromosome and gffcompare the result against StringTie's own GTF.

```
bash bench/gw_eonly_parity.sh      # ~9.5 min (568s), reuses st_$C.gtf
python3 bench/gw_eonly_aggregate.py
```

**Result (25 chroms, after the two fixes below): BYTE-EXACT on every chromosome.**

| metric | value |
|--------|-------|
| Transcript **Sensitivity** vs StringTie | **100.0 on every chromosome** |
| Transcript **Precision** vs StringTie | **100.0 on every chromosome** |
| floor transcripts vs StringTie | **68,157 == 68,157** (difference 0) |
| chromosomes byte-exact (ee_tx==st_tx, Sn=Pr=100) | **25 / 25** |
| NCBI-matched (strict =/c) vs StringTie | **25,148 == 25,148** (same annotated set) |

`rustle -e -G st.gtf` reproduces StringTie's *exact transcript set on all 25 chromosomes* —
the StringTie-exact floor the VG layer builds on. Reaching it took two fixes (below): the
journey was **+50 → +6 → 0**.

### The duplicate-emission bug (found here, fixed)

The first floor run showed **+50** extras (68,207 vs 68,157), all class `=`. Tracing them:
every one was a **duplicate guide emission**. Root cause:

> A guide transcript spanning a coverage gap is split across bundles; guided/eonly mode
> **force-emits the guide once per overlapping bundle** with bundle-local coverage,
> producing identical-coordinate duplicates. StringTie emits each guide exactly once. No
> existing dedup caught them — the intron-chain / same-span passes all **exempt guides and
> skip single-exon**, which is exactly the duplicated population. The split copies' coverage
> sums to the true single-bundle value (e.g. `1.004 + 43.433 = 44.437`).

It was a **general rustle bug, present in every mode** (StringTie never duplicates):

| output | exact-duplicate-transcript extras | single-exon of those |
|--------|-----------------------------------|----------------------|
| StringTie       | **0**  | 0  |
| rustle-eonly    | **44** | 41 |
| rustle-VG       | **20** | 5  |
| rustle guided   | **19** | 5  |

**Fix:** `transcript_filter::dedup_identical_transcripts` collapses byte-identical
transcripts (same chrom/strand/exons), keeps a guide-preferred representative, and **sums
coverage** onto it (feeding the correct value to TPM/FPKM). Wired at the global stage before
quantification, gated `!precise_mode()` (RUSTLE_PRECISE stays byte-identical to 4705ab1),
opt-out `RUSTLE_DEDUP_IDENTICAL_OFF=1`. Paralog-safe: VG copies live at distinct coordinates
and are never byte-identical. Verified: 4 unit tests; `NC_073247.2` eonly 2396→2391 (==
StringTie); chr19 **de-novo default unchanged** (2013==2013, no-op where no dups exist);
`RUSTLE_PRECISE` dup preserved. Genome-wide eonly floor **+50 → +6**, 20/25 chroms byte-exact.

### The residual +6 (`_micro5merge`) — gated out of eonly

After the dedup fix, 6 extras remained: **8 `_micro5merge` variants** genome-wide (net +6).
rustle deliberately emits a micro-exon-**merged** alternative of a guide alongside the
faithful copy (e.g. fusing two ≤73 bp micro-exons), tagged `source
"guide:STRG.x.y_micro5merge"`, class `m`. This is a *separate, intentional* micro-exon
transformation — a rustle "different decision," not a duplicate — but in eonly mode (strict
reference reproduction) it leaks one extra isoform per affected guide.

**Fix:** the 5' micro-exon merge block (pipeline.rs ~18146) is now suppressed in eonly mode
(`config.eonly && !precise_mode()`); the faithful unmerged copy always survives, so
Sensitivity is unaffected. It remains default-on in guided / de-novo mode (where the extra
isoform is a legitimate prediction) and under RUSTLE_PRECISE. With both fixes the eonly floor
is **byte-exact on all 25 chromosomes** (68,157 == 68,157, Sn=Pr=100, NCBI 25,148 == 25,148).

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

## PSV / multimapping — final verdict (2026-06-15): would PSVs beat StringTie? Not meaningfully.

Rigorously tested both PSV routes to beating StringTie. Conclusion: **PSVs are a real, principled
mechanism but the defensible margin over StringTie is small**; the broad primary-read route is
alignment noise. Scripts: `bench/psv_sizing/`.

1. **Secondary-driven copy recovery** (`RUSTLE_VG_RECOVER_COPIES`, commit e6a51f7): **+77
   NCBI-corroborated additive isoforms** over the primary-only/StringTie-level baseline on 2
   paralog-dense chroms (21 solid `=`/`c`; 56 `j`; **only 3 exact `=`**; 1/27 starved-copy targets
   cleanly improved). Real but small and mostly partial; exactness is coverage/identifiability-limited.

2. **Primary-read PSV-phasing** (sizing only, not built): sizer flagged 84 loci with PSV-linked
   structural divergence — but validation collapsed it: 14 already in StringTie, **2 RefSeq-
   corroborated**, 68 novel. The 68 are **not** cross-mapping artifacts (reads align strictly best
   at their own locus, 0/69) — but junction-QC of the 89 novel divergent junctions found **88/89
   non-canonical, 0 credible-real** (canonical + support≥5 + non-RT): they are **systematic paralog
   alignment artifacts** (spurious non-canonical `N`-gaps), not real splicing. Defensible headroom ≈ **2**.

**Bigger lever for "beat StringTie" is non-PSV:** read-coherence (`--read-chain`) adds **+1,857
strict-real (FSM/NIC) isoforms** genome-wide over the flow baseline (`bench/readcoherence_finding.md`)
— ~20× the PSV margin. Pivoting effort there.
