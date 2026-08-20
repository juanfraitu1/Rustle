# Haplotype copy-number run: UNINFORMATIVE by its own pre-registered criterion

**Status 2026-08-19. Run complete. Verdict: the primary number must not be quoted.**
Pre-registration: [`o3_haplotype_cnv_scope.md`](o3_haplotype_cnv_scope.md).
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
