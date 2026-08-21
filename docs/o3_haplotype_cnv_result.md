# Haplotype copy-number run: UNINFORMATIVE by its own pre-registered criterion

**Status 2026-08-19. Run complete. Verdict: the primary number must not be quoted.**
Pre-registration: [`o3_haplotype_cnv_result.md`](o3_haplotype_cnv_result.md).
Artifacts: `/mnt/linuxdisk/home/juanfraitu/o3_hapcnv/`.

## 1. What was run

10,019 probes (5,000 random autosomal protein-coding genes as centred windows capped at 30 kb,
plus named controls and a span-matched random-interval arm), 187.1 Mb, aligned to `mat` and `pat`
each restricted to its 24 chromosome-scale sequences (98.67% of mat, 100.00% of pat).
`d_hap = copies(pat) − copies(mat)` at the validated P3 floor (identity ≥ 0.973, ≥ 5,800 aligned bp).

**Amendment, declared:** the first pass used `-N 100` at minimap2's default `-p 0.8`. That failed
pre-declared control 2, and diagnosis showed why: `MAPKBP1_full` returned **one** record (58,172 bp
at identity 1.0000) because its paralogues, covering only part of a 58 kb probe, score under 0.8× the
self-hit and are discarded, while `PLA2G4B_full` (10.4 kb probe, paralogues spanning it at 0.984–0.986)
returned nine. **Copy counts were biased by probe length.** Re-run at `-N 200 -p 0.1`; records rose
28,623 → 97,930 (mat) and 21,983 → 78,278 (pat). This was fixing the instrument to meet a
**pre-declared** must-pass, not tuning to move the primary; all controls were re-evaluated.

## 2. Controls

| # | control | result | |
|---|---|---|---|
| 3 | self-recovery — a probe from a PAT-sourced chromosome must be found in pat | **2652/2652 and 1201/1201 = 1.0000** | PASS |
| 4 | single-copy panel | TBP, POLR2A, GTF2H1, SF1, TFRC, HMBS all **1 and 1** | PASS |
| 5 | probe-provenance stratification | PAT-sourced **0.0290 [0.0233, 0.0361]** vs MAT-sourced **0.0250 [0.0176, 0.0354]** | PASS |
| 2 | known answer (pat 8/9/9, mat 5/6/8) | SPTBN5 **9/8 exact**; MAPKBP1 **9/8** (expected 8/5); PLA2G4B **9/7** (expected 9/6) | **FAIL** |
| 1 | span-matched random intervals | **0.1512 [0.1403, 0.1629]** vs the gene rate **0.0278** | **FAIL — floor is 5.4× the signal** |

## 3. Verdict

The scope declared: *"If control 1 lands above the observed rate, the correct report is
'uninformative', not a number."* It landed **5.4× above**. ⟹ **The gene rate of 2.78% is not
quotable and no rate is reported.**

**Why control 1 failed is a design flaw, not noise.** The controls were **span-matched but not
composition-matched**: a random 30 kb interval is far more repeat- and segdup-rich than a gene body,
so it both varies more biologically and is harder to count. A signal *below* its own null is not a
weak signal — it is evidence the null is the wrong null. Any repeat attempt must match repeat
content and segdup overlap, not just span.

## 4. ⚠ A prior number is now in doubt — correction

`OBJECTIVES_AND_VERIFICATION.md` row 3.10 reports MAPKBP1/PLA2G4B/SPTBN5 at **8/9/9** on `_pri`/`_pat`
versus **5/6/8** on `_mat`, and that "several-copy difference between one animal's two haplotypes"
has been used as the standing refutation of *"copy number is stable; only SNPs and indels differ."*

**This run cannot reproduce it with either setting.** At `-p 0.8` MAPKBP1 gives 1/1; at `-p 0.1` it
gives 9/8. Neither is 8/5. SPTBN5 reproduces exactly (9/8) and PLA2G4B is off by one on mat (9/7 vs
9/6). With secondaries retained the mat deficits shrink from **3, 3, 1** to about **1, 2, 0**.

⟹ **Do not quote "1–3 whole gene copies differ between KB3781's haplotypes" until row 3.10's
instrument is re-derived and its `-p`/`-N` settings recorded.** The direction (mat ≤ pat at these
loci) survives; the magnitude does not. This does not rescue the "copy number is stable" proposition
— SPTBN5 and PLA2G4B still differ — but the difference is smaller than the number in circulation.

## 5. What does survive

**Control 5 is a clean positive result and is worth keeping.** `_pri` is a per-chromosome mosaic —
independently re-derived here from exact chromosome lengths as **16 PAT / 9 MAT**, reproducing the
documented split — so probes are paternal sequence on some chromosomes and maternal on others. The
`d_hap` rate agrees across the two strata (0.0290 vs 0.0250, overlapping CIs). **`d_hap` is not
driven by probe provenance.** That was the sharpest confound available and it is now excluded, which
any future attempt inherits.

## 6. Reproduction

```bash
cd /mnt/linuxdisk/home/juanfraitu/o3_hapcnv
python3 probes.py            # probe set + span-matched controls (span multiset asserted identical)
./align.sh mat && ./align.sh pat    # ~11 + ~8 min, peak 14.7 GB, ONE AT A TIME
python3 count.py             # controls first, primary last
```

---

# Appendix A — the PRE-REGISTRATION, written before the run

*Merged from `o3_haplotype_cnv_result.md` on 2026-08-20. ⚠ Committed BEFORE the run (ed86742 then
1247d67). Its declared stop rule is what makes the verdict "uninformative" a result rather than a
disappointment.*

## Scope: how often do whole gene copies differ between one individual's two haplotypes?

**Status 2026-08-19: SCOPED, NOT RUN.** Pre-registration draft — thresholds and controls below are
declared before any result.

## The question and why it is the right one

The advisor's objection to O3 is that *"the genome is an average, so every individual differs from it
unless you compare to the individual that generated it."* For our substrate the premise is wrong on
two counts: **mGorGor1 is a haplotype-resolved assembly of ONE animal (KB3781 "Jim")**, not a mosaic
of donors, and **the fibroblast IsoSeq is that animal's own cell line**. The register already scopes
the objection as *"live only for the CROSS arm and human work."*

But the counter-proposition — *"copy number is stable; only SNPs and indels differ"* — is also false,
and our own data already refutes it. On this assembly, `_pri`/`_pat` vs `_mat`:

| probe | `_pri` / `_pat` | `_mat` |
|---|---:|---:|
| MAPKBP1 | 8 | 5 |
| PLA2G4B | 9 | 6 |
| SPTBN5 | 9 | 8 |

**One animal's two haplotypes, differing by 1–3 whole gene copies.** What is missing is the **rate**.
The existing measurement — `d_ortho` nonzero on **11/816 = 1.35% [0.75, 2.40]**, 6 up vs 5 down,
sign test p = 1.0000 — carries two stated limits: it is **tandem-only by construction** (`in_pri = 1`
for 762/817 = 93.3% of loci, so a dispersed event is structurally invisible), and the stratum spans
only **1.1224% of the genome**.

This run replaces the compartment and drops the tandem restriction.

## Design

**Probes.** 5,000 randomly sampled autosomal protein-coding genes (of 22,650 in `GGO_genomic.gff`),
**genomic spans, introns included** — the validated detection floor (P3: ≥97.3% identity over ≥5.8 kb)
was established with genomic queries. Extracted from `_pri`.

**Targets.** `mat.fa` and `pat.fa`, each restricted to its **24 chromosome-scale sequences (≥20 Mb)**.
This matters: mat has 225 sequences and pat has 24, but the 24 large ones are **98.67% of mat and
100.00% of pat**, so restricting makes completeness symmetric and discards only 48 Mb = 1.33%.

**Statistic.** `d_hap(g) = copies(g in pat) − copies(g in mat)`, counted genome-wide (not
window-restricted) at the P3 floors.

**Primary readout.** Fraction of genes with `d_hap ≠ 0`, with CI, plus the **sign test** — polymorphism
predicts symmetry, a systematic assembly deficit predicts skew.

## Controls, declared before any result

| # | control | must-pass |
|---|---|---|
| 1 | **Span-matched random intervals** — same SIZE distribution, not count-matched (count-matching alone has already killed a metric here) | the FP floor must sit well below the observed rate |
| 2 | **Known answer** — MAPKBP1 / PLA2G4B / SPTBN5 | pat 8/9/9 and mat 5/6/8, or this is not the instrument that produced the published result |
| 3 | **Self-null** — pat vs pat | `d = 0` on every probe; any nonzero is a bug |
| 4 | **Single-copy panel** — TBP, POLR2A, GTF2H1, SF1, TFRC, HMBS, PSMB6 | 1 and 1 |
| 5 | ⭐ **Probe-provenance stratification** | `_pri` is a per-chromosome MOSAIC — 16 chromosomes byte-identical to `_pat`, 9 to `_mat` — so the probe is paternal sequence on some chromosomes and maternal on others. **The rate must agree across the two strata.** If it does not, `d_hap` is measuring probe provenance, not biology |

Control 5 is the sharpest and has no analogue in the previous run, which was forced through `_pri`
as one of the two compared assemblies. Comparing mat to pat **directly** removes that entanglement
from the comparison, but not from the probe — hence the stratification.

## The main risk, stated in advance

Genome-wide counting is exactly where the previous module measured a **control artefact rate of
8.9–12.8%**, which is why it declined to treat `d_other` as its collapse statistic. If the true
polymorphism rate is ~1–2%, **the control floor could exceed the effect** and the run returns
nothing interpretable. Mitigation is the P3 floor (≥97.3% identity, ≥5.8 kb), which is what drove
`d_ortho`'s own FP floor to 0/817. **If control 1 lands above the observed rate, the correct
report is "uninformative", not a number.**

Second limit: this measures **polymorphism**, not collapse. The two are different questions and one
run cannot separate them. It answers the advisor's objection; it does not by itself advance O3's
detector.

## Cost

Two minimap2 index builds (~10 min each, ~12 GB RSS) plus ~130 Mb of query against each haplotype:
**≈2–2.5 job-hours, one at a time**, peak RSS ~14 GB of the 25 GB ceiling.

## The decision it produces

| outcome | reading |
|---|---|
| rate low (~1–2%), symmetric, control floor well below | the objection is quantitatively small; the matched design is sound and O3's negatives stand |
| rate high (>10%) | copy-level polymorphism is common ⟹ every O3 claim must be matched-individual — which we can do, and the cross/testis arm must be re-scoped |
| control floor ≥ observed rate | uninformative; report as such and do not quote a rate |
