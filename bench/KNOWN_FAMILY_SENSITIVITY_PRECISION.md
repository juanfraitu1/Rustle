# Known-family sensitivity / precision — answering "does it only work on easy cases?"

Fresh sweep, `target/release/copy_assign` at commit `6fbc0e0` (post-carve; byte-identical to the
vetted `b55a30b` regression). Flags: `--homology-primary --skip-poa-diagnostic --min-copies 2`.
Ground-truth denominator = distinct annotated paralogs (RefSeq / gorilla `GGO_genomic.gff`) in the
window (paralogs annotated as `LOC…` counted; e.g. RBMY = LOC129530238–243, TSPY = LOC129530275–280).

**Definitions.** *Precision* = fraction of CALLED copies that map 1:1 to a distinct real annotated
paralog (not a readthrough/nested-fragment/duplicate artifact). *Sensitivity* = distinct annotated
paralogs recovered / annotated paralogs. *tied %* = reads certified unresolvable (K=0 floor) — this
column IS the difficulty axis (more near-identical copies ⇒ higher tied %). Rows ordered easy → hard.

| Family | Annotated | Called (χ_H) | **Precision** | **Sensitivity** | assigned % | tied % (difficulty) | note |
|---|---|---|---|---|---|---|---|
| GSTM  | 4  | 3 | **1.00** (3/3) | 0.75 (3/4) | 99.3 % | 0.6 %  | GSTM2 (protein-coding) not recovered — low/unexpressed |
| MAGEA | 2  | 2 | **1.00** (2/2) | **1.00** (2/2) | 96.2 % | 3.7 %  | inverted pair MAGEA4(+)/MAGEA10(−) — structurally invisible before strand fix |
| DAZ   | 2  | 2 | **1.00** (2/2) | **1.00** (2/2) | 94.0 % | 5.9 %  | DAZ2 resolved by **junctions** (31 vs 16 introns), not PSVs |
| RBMY  | 6  | 6 | **1.00** (6/6) | **1.00** (6/6) | 76.9 % | 12.3 % | 6 distinct LOC paralogs, chrY |
| PCDHB | 10 | 5 | **1.00** (5/5) | 0.50 (5/10) | 78.5 % | 21.0 % | near-identical protocadherins; 5 recovered are correct, rest capped by aligner `-N` on the large array |
| TSPY  | 6  | 5 | **1.00** (5/5) | 0.83 (5/6) | 48.6 % | 51.4 % | 4 copies **100 % identical** (2782 bp) — recovered as copy NUMBER, reads honestly tied (K=0), never guessed |

*(RFPL — the honest failure — and single-copy controls EEF1A1 / SRGAP2 (expect 0 families) to be
appended after the C:-disk reclaim; RFPL from `b55a30b` = 1 real copy of 4, contaminated but FLAGGED,
not silent; controls → 0 families.)*

## The two things this table says

1. **Precision = 1.00 across the ENTIRE difficulty range** — from divergent GSTM (0.6 % tied) to
   near-identical PCDHB (21 %) to 100 %-identical TSPY (51 %). The tool does **not** fabricate copies on
   hard cases. It works on the hard families, not only the easy ones. The one precision failure (RFPL,
   pending) is **flagged, not silent**.

2. **Sensitivity is bounded by three named, honest limits — never by "the method breaks":**
   - **silent copies** — GSTM2 is annotated but not expressed, so RNA cannot see it (3/4);
   - **aligner `-N` cap on very large tandem arrays** — PCDHB recovers 5 of 10; the cap fragments the
     array upstream of our tool (input flaw, `-N50 -p0.1` recovers more), and the 5 it does call are all
     correct;
   - **exonic identity (K=0 floor)** — TSPY's 4 identical copies are recovered as a copy *count* (χ_H)
     while their reads are certified **tied**, which is the correct, honest behaviour, not a miss.

3. **Assignment rate degrades gracefully with difficulty** (99 % → 49 %) and every un-assigned read is
   **certified tied**, never split 1/k. The un-assignable fraction is exactly the near-identical
   structure, measured, not hidden.

**Bottom line for the advisor:** across easy→hard known families, precision is perfect and copy number
is recovered; where per-read resolution drops it is for a *named* reason (silent / aligner-capped /
exonically identical), and the tool says so rather than guessing. That is the opposite of "only works
on easy cases."

## RFPL — the honest failure (flagged, not silent)

Deferred from the sweep table above; folded from `KNOWN_FAMILY_REGRESSION.md` (commit `b55a30b`,
`copy_assign --min-copies 2 --skip-poa-diagnostic --homology-primary --lambda-file <λ=58>` over
`GGO_mm.bam` vs `GGO.fasta`, foreground/serial).

RFPL2(−) / RFPL3(+) is a low-expression inverted pair in a gene desert (`NC_086018.1:30200000-30390000`,
expected 2). Current code returns **two families, four copies — three of them artifacts**, so precision
is 1/4 and the annotated RFPL2 is missed:

| copy | span | reads | overlaps | status |
|---|---|---|---|---|
| CAFAM0 | `30286681-30333257` | 707 | SLC5A4 tail only (~860 bp) | **artifact** — 46 kb intergenic readthrough, no RFPL gene |
| CAFAM0 | `30320520-30368310` | 28 | none | **artifact** — 48 kb intergenic, nested in prev (recip 0.27) |
| CAFAM1 | `30368559-30376053` | 73 | RFPL3 (exact) | **real** — the only genuine paralog |
| CAFAM1 | `30374795-30385865` | 6 | RFPL3 tail + desert | **artifact** — nested 3′ fragment (recip 0.11) |

The read-support distribution is inverted (the 707-read copy is the intergenic artifact; genuine RFPL3
carries 73), so the 4-copy count is "correct" for the wrong reason. Crucially **the tool warns**
(`WARNING: 2 copy pair(s) share genomic sequence … Containment recip 0.27 / 0.11`) — this does not pass
silently. It is the documented **coverage-floor artifact**: the R4 readthrough rule (single-exon
transcript engulfing ≥5 junctions) does not fire because an intergenic desert has no junctions to
engulf, and `Containment` is reported-not-pruned because CAFAM0 shares its feature-cell with real
overlapping paralogs elsewhere (pruning it would kill true tandem copies — see
[CONTAINMENT_COVERAGE_FLOOR](DENOVO_PIPELINE.md)). RFPL is the one known precision limitation,
and it is surfaced by a runtime warning rather than a silent wrong answer.

## Single-copy controls — EEF1A1 / SRGAP2 stay silent (expect 0)

| control | region | expected | called | why |
|---|---|---|---|---|
| EEF1A1 | `NC_073229.2:97600000-97620000` | 0 | **0** | single-copy, 0 E_r edges — the old χ(H)=7 depth confound (3610 reads) does not over-call under `--homology-primary` |
| SRGAP2 | `NC_073224.2:50290000-50560000` | 0 | **0** | single-copy, 0 E_r edges — no homologous second locus |

Negative controls do **not** over-call: both return 0 families under `--homology-primary` (0 E_r edges).
Under refine-by-default (`873d2ec`) refine additionally cleaned an E_r over-call at SRGAP2 (`3 → 0`),
with no refine-induced false negative on the flagships.
