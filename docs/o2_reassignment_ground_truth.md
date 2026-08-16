# O2 — non-circular ground truth for REASSIGNMENT: design + feasibility (2026-08-15)

**Status: DESIGNED, feasibility MEASURED, not yet run.** This is the last Tier-1 gap in O2. Everything
below is executable with data on disk.

## 1. The gap this closes

O2's abstention capability now has a non-circular validation (AUC 0.7995 against MAPQ's 0.4944, on
excision reads whose true copy is known by design). **Reassignment does not.** Its only real-data
evidence is agreement with minimap2's primary flag — which *is* the "98.4% restatement" complaint, so
the validation and the criticism are the same measurement.

⚠ The existing non-circular validations are **simulated**: `project_sim_ground_truth` (planted
2-chromosome genome) and `project_k0_flank_experiment` (planted 4-copy locus, 60 reads/copy). Both are
airtight and both are synthetic. **The gap is real reads.**

## 2. The idea: STRUCTURAL anchors, not PSVs

For a read at a multi-copy locus, ask what could label its copy of origin *without* an aligner's
contested decision.

* ⚠ **A PSV is a scoring decision** — and O2 scores PSVs. Using PSVs as truth would be circular.
* ⭐ **A structural anchor is presence/absence**: sequence present in copy A and simply **absent** from
  copy B. A read containing it is from A. No threshold, no score, nothing to tune.

**The protocol.** For a read overlapping an anchor: (i) truth = the copy carrying the anchor;
(ii) **trim the anchor off the read**; (iii) hand the trimmed read to the assignment method, which now
sees only shared sequence; (iv) score against the anchor truth. **The method never sees the evidence
that produced the label.**

The same trimmed reads are scored under minimap2's primary flag, giving the head-to-head comparison
that has never been possible on real data.

## 3. Feasibility — measured 2026-08-15, on the 162 two-copy excision panel

⚠ **Substrate matters and was got wrong once already.** On GENOMIC spans, 38/60 copy pairs do not align
at all — `E_r` nodes are **exon-sum spliced representatives**, not genomic intervals (rep-vs-rep
reproduces the shipped rule 103/104; genomic-vs-genomic 1/3). And since IsoSeq reads are spliced, the
anchor must be **exonic** regardless. All numbers below are on exon-sum sequence (`q915_exon.fa`).

| | families |
|---|---|
| **structural anchor ≥100 bp** | **34 / 162 = 21.0%** |
| anchor 50–99 bp | 1 |
| no anchor | 18 |
| no alignment under `asm20` | 109 |

Anchor size: min 104 bp, **median 779 bp**, max 10,724 bp.
⚠ The 109 no-alignment cases are `asm20` being far stricter than E_r's sensitive tier (`-k 11 -w 5` at
identity ≥ 0.60) — **the true anchor count is a lower bound**; re-running anchor discovery at the
shipped tier should recover more.

Read depth at ten anchored families (primary only, `-F 2308`): median **72** reads on the *lower*-depth
copy, ~1,404 reads across those ten counting only the lower copy. Across all 34, several thousand.

## 4. Protocol

1. **Anchor discovery** at the SHIPPED sensitive tier (not `asm20`): align each family's exon-sum reps
   pairwise with `--cs`, extract indels ≥100 bp. Record anchor coordinates in both rep and genomic space.
2. **Truth set**: reads whose alignment overlaps an anchor interval by ≥90% of the anchor. Label = the
   copy carrying it. ⚠ Require the overlap to be *within* the read, not at a soft-clipped edge.
3. **Trim**: remove the anchor-overlapping portion plus a 50 bp margin. Keep reads with ≥300 bp
   remaining; report how many are lost to the length floor (they are not failures, they are out of scope).
4. **Score**, per trimmed read: O2's decision (assign A / assign B / abstain), and minimap2's primary
   flag on the same trimmed read, both against the anchor truth.
5. **Report**: accuracy-when-committed, abstention rate, and the head-to-head against the primary flag,
   each with a 95% CI, split by whether the trimmed read still carries any PSV.

## 5. What it will and will not establish

**Will**: the first non-circular reassignment accuracy on real reads; the fraction of genuinely
contested reads that are resolvable at all; and whether O2 beats the primary flag where truth is known.

**Will not**: generalise beyond families that *have* a structural anchor — 21% of two-copy families,
and anchored families may be systematically more divergent than the rest. **State that as a selection
effect, do not hide it.** ⚠ Expect a large abstention/unresolvable fraction: the flank experiment showed
exon-confined reads of exonically identical copies are **100% tied — provably nothing in the BAM**.
That is the expected result for the K=0 stratum and is informative, not a failure.

## 6. Traps

* **T8** — an offline re-derivation is a hypothesis generator, never a test. Score through the actual
  pipeline, not a proxy.
* **T1** — do not condition the denominator on the method's own output. The denominator is
  *anchor-labelled reads*, fixed before scoring.
* **T12** — state the unit. This is per-READ; do not mix with per-copy or per-family rates.
* **T7** — trimming changes what the read is; do not compare trimmed-read results to untrimmed baselines.
* Pre-register thresholds (length floor, overlap fraction) **before** scoring.

## 7. Artifacts

Feasibility probe: `/home/juanfra/winloci_scratch/o2_truth_real/`. Panel: `o3_excise/panel.json`.
Exon-sum sequence: `o3_collapse/method/intervals/data/q915_exon.fa`. Reads:
`/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam` (matched individual).
