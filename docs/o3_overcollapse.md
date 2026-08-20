# Over-collapse: is it happening, and could we see it?

**Status 2026-08-19.** Analysis of the original O3 premise, plus a scoped simulation that is
**specified, not run**. Companion to [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md),
[`o3_unmapped_route.md`](o3_unmapped_route.md) and the O3 restatement in
[`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md).

## 0. The premise

*"Very similar copies may be collapsed into one in the assembly. Align RNA against that assembly and
the over-collapsed copies are missed."* That is O3's original motivation. It decomposes into two
questions with **different answers**, and conflating them is the main error to avoid.

## 1. Is it happening? — UNKNOWN, and our zero does not say no

| | |
|---|---|
| DNA-side S1, three assemblies of one gorilla | **0 collapses at 816/817** |
| transcript side | **0/915** — concordant |
| literature (Yoo/Rhie) | ~**1–2 Mbp of collapse per haplotype** |
| what that predicts over the **1.1224%** of the genome those windows span | **0.47–0.94 events** |

**The zero was predicted before it was observed.** An instrument with sub-one expected yield returns
zero whether or not collapse exists. This is underpowered *by construction*, not a negative result.

⚠⚠ **And the screen's validation now needs re-checking.** It rested on recovering a published
expansion exactly — MAPKBP1/PLA2G4B/SPTBN5 at 8/9/9. On 2026-08-19 that count proved
**`-p`-sensitive**: at minimap2's default `-p 0.8`, MAPKBP1 returns **1** copy, because its paralogues
cover only part of a 58 kb probe and score under 0.8× the perfect self-hit; at `-p 0.1` it returns
**9**. If the collapse screen ran at default `-p`, its copy counts were suppressed on long probes —
**biasing it toward zero**. See [`o3_haplotype_cnv_result.md`](o3_haplotype_cnv_result.md).
**Do not re-quote "0/817 with a validated instrument" until that screen's `-p`/`-N` are recorded.**

## 2. Could we see it? — the signature is NOT what O3 originally assumed

**A collapsed copy does not orphan its reads. It absorbs them.** Measured in the whole-genome excision
control: absorbed copies show **1.75× depth** on the surviving paralogue (FixItFelix reports 1.5×),
concordance 0.967. So the signature is **depth + PSVs** — the field-standard **S2** — and *not*
unmapped reads, and *not* clipping (no published collapse detector uses clipping).

### 2a. In RNA, the depth half is dead

Depth conflates copy number with **expression level**, which varies over orders of magnitude, so 2×
coverage at a locus means nothing without a per-gene expectation. Already on the register as
human-review-flag-only. **That leaves only the PSV half.**

### 2b. PSVs need divergence — and collapse is *caused by* the lack of it

| | TPR |
|---|---:|
| S2 detector, copies **≥ 0.01** diverged | **0.4500** |
| S2 detector, copies **< 0.01** diverged | **0.0588** |
| **fraction of true positives below 0.01** | **45.78%** |

held-out overall: **TPR 0.2703 / FPR 0.0200**.

**Copies collapse *because* they are near-identical.** The mechanism that creates the target is the
mechanism that destroys the evidence. This is **adverse selection, not bad luck**, and no threshold
choice escapes it.

At the limit it is a genuine impossibility rather than a power problem: **K = 0 — perfectly identical
copies — is unidentifiable, and that is entailed, not measured.** Between K=0 and 0.01 divergence it
is power; above 0.01 it is 0.45.

## 3. So it is not impossible — but not from RNA alone

| half | source | status |
|---|---|---|
| PSV | RNA | available, bounded, adverse-selected |
| depth | **RNA** | ⚠ **structurally dead** — expression confound |
| depth | **DNA** | ⭐ **available and unblocked** |

The gorilla WGS was confirmed usable — the pre-registered k-mer test rejected the "Y flow-sorted"
label at **17–1600×** (ENA's `sample_title` is simply wrong for every SAMN04003007 run). **The open
route is RNA PSVs + DNA depth on the matched individual, over a compartment large enough to expect
more than one event.** Both halves exist; neither has been run against the other.

## 4. Scoped simulation — a synthetic-collapse calibration ladder

**SPECIFIED, NOT RUN.** This is the calibration the collapse screen explicitly lacks: *"`d_ortho` has
never been calibrated by a synthetic collapse. The one such test on record is arithmetic, not
calibration — it proves the statistic is wired without an off-by-one, not that minimap2 plus the
clustering rule resolve two real adjacent tandem copies as two."*

### 4a. Collapse is NOT deletion — simulate the right operation

The excision control **deleted** one copy: the survivor is copy A's sequence, and copy B's reads must
find a home. **Collapse merges**: the assembly carries **one** sequence where biology has **two**, and
reads from *both* copies land on it. The operations differ, and the existing panel only tested
deletion. **Both arms must be run on the same families** so the difference is measured, not assumed.

### 4b. Design

Hybrid — **real reads, modified genome**. Simulating reads would forfeit the error and expression
structure that make the excision control credible.

1. Take the existing **162 two-copy family** panel and its matched IsoSeq.
2. Bin each family by its pair's **measured divergence** (0, ≤0.001, ≤0.002, ≤0.005, ≤0.01, ≤0.02, >0.02).
3. **Collapse arm:** in ONE modified genome, replace both copies of every family with a single
   sequence — as the excision arm masked all 162 in one genome, so it is **one alignment, not one per bin**.
4. **Deletion arm:** the existing excision genome, unchanged, for the contrast in §4a.
5. Realign the real matched IsoSeq; run the existing `detector.py`.
6. **Primary output: TPR as a function of divergence** — converting the two-bin summary
   (0.4500 / 0.0588) into a **located detection floor**.

Reuses `o3_excise/{mask_genome.py, align.sh, detector.py, panel.py}` and the frozen panel.

### 4c. Controls, to declare before running

| # | control | must-pass |
|---|---|---|
| 1 | unmodified genome, same reads, same detector | the FP floor; must be far below the collapse-arm rate |
| 2 | single-copy loci | must not fire |
| 3 | deletion arm on the same families | establishes that collapse ≠ deletion rather than assuming it |
| 4 | `-p` / `-N` recorded and swept | ⚠ the 2026-08-19 lesson — a default `-p 0.8` silently discarded 8 of MAPKBP1's 9 copies |

⚠ **Known traps.** Identical copies **fake junctions** (the TSPY simulation). Simulations need `-N 50`
— `--secondary=no` yields **0 families**. And a planted-divergence experiment scored on planted
divergence is circular: the detector must be blind to the bin.

### 4d. What it can and cannot establish

**Can:** where the detection floor actually sits, whether K=0's impossibility begins at 0.000 or bites
well above it, the FP floor under a known negative, and whether collapse and deletion really differ.

**Cannot:** ⚠ **whether over-collapse is happening in mGorGor1.** That is a *prevalence* question and
needs a larger compartment, not a calibration curve. A ladder tells you what the instrument can see;
it never tells you what is there.

### 4e. Cost

One genome modification, one index, one IsoSeq realignment against the existing panel: **≈2–3
job-hours**, one at a time, comparable to the excision run.
