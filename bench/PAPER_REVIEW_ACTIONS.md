# Paper-grounded review — actions executed (2026-06-28)

Results of the items selected from `bench/PAPER_GROUNDED_REVIEW.md` (F1, F2, F5, F6). Each closes a
literature-grounded gap with a concrete, verified artifact.

---

## F6 — Eichler AS≥10 rule vs our significance gate (`bench/as_decisive_vs_gate.py`)
**Paper:** Guitart/Eichler 2024 (TBC1D3) — the rule the advisor keeps citing: map a read to all haplotypes,
assign iff best minimap2 AS beats the 2nd by ≥10, else ambiguous. We had only *calibrated* τ=6.9 to its
operating point; now benchmarked literally, on the labeled sim5x ladder.

| K | AS≥10 recall | gate recall | gate acc | AS margin med/p90/max |
|---|---|---|---|---|
| 0 | 0.0% | 0.0% (all Tied, min_p=1) | — | 0/0/0 |
| 1 | 0.0% | 60.0% | 1.000 | 5/5/5 |
| 2 | 0.0% | 100.0% | 1.000 | 0/5/5 |
| 4 | 0.0% | 100.0% | 1.000 | 0/5/5 |
| 8 | 0.0% | 100.0% | 0.990 | 0/0/0 |

**Divergent control** (copy0 vs a 3%-diverged copy, 30 reads): AS≥10 assigns **30/30 correctly, median margin
295 (≫10)** — so the rule is NOT broken.

**Finding (the exhibit):** on **near-identical / collapsed paralogs** — the regime our thesis targets — the raw
whole-read AS margin is **≤5 and never reaches 10**, so AS≥10 resolves **0%**: the 1–2 distinguishing PSVs are
drowned in the full-read score. Our PSV gate, which scores only the **decisive columns**, resolves **60–100% at
≥0.99 accuracy**, and at K=0 (identical copies) certifies every read **Tied (min_p=1)** instead of silently
discarding. This reproduces the literature's regime split — **AS for divergent copies (TBC1D3); PSV-level
methods (IsoCon/SDA/ours) for the collapsed regime** — and positions our gate as the principled generalization of
the advisor's favourite rule, with an identifiability certificate AS≥10 lacks.

---

## F5 — StrandOddsRatio (SOR) strand-bias filter for ASJ (`bench/asj_strand_bias.py`, wired into `asj_aggregate.py`)
**Paper:** longcallR 2026 filters allele-specific junctions with SOR ≥ 2 (reject) to kill single-strand-driven
calls; our ASJ pipeline had **no** strand-bias check and didn't even tally read strand. Implemented the GATK SOR
(`ln(symmetricRatio)+ln(refRatio)−ln(altRatio)`, +1 pseudocounts) over the 2×2 [used/not] × [fwd/rev] from
`read.is_reverse`, annotated all 475 calls, and added a `--max-sor` gate (off-by-default-equivalent: annotates
`sor`+`sor_pass`, hard-drops only when `--max-sor` is passed).

**Results (475 calls):**
- KEEP at SOR<3.0 (GATK default): **401/475 (84.4%)**; at SOR<2.0 (longcallR-strict): **337/475 (70.9%)**.
- Transversion genetic core (120): **94 pass <3.0, 80 pass <2.0**.
- ⚠ **Both flagship calls FAIL strand bias: PSMD2 SOR=10.45 (all 14 junction reads forward-strand, R−=0), DAXX
  SOR=7.08.** Combined with the M1 motif retraction, this is a third independent reason **PSMD2/DAXX should not be
  headlined** as the clean ASJ exemplars — they carry a single-strand-usage artifact pattern. The genetic *core*
  as a set is robust (most survive); the two poster-children are not.
- Caveat surfaced: 36/74 failures have 0 junction-using reads in `GGO_mm.bam` (calls were generated on the
  single-mapping `GGO.bam`); their SOR reflects spanning-read strand skew, not junction-usage skew. PSMD2/DAXX
  are NOT in that degenerate set — they fail on real junction-usage imbalance.
- Follow-up: the Rust production path (`src/bin/asj.rs`) needs `is_reverse` threaded into `AlignedRead` for a
  native gate (noted in the script).

**Action item it raises:** retire PSMD2/DAXX as the flagship ASJ pair; lead O3 with the genetic-core *set*
(transversion, SOR-passing, editing-controlled) and the DAZ1-vs-DAZL copy-specific example instead.

---

## F2 — Soto-2025 DNA/CN family validation, no SEDEF (`bench/soto_family_validate.py`)
**Paper:** Soto 2025 (Cell) defines families at the DNA level: map segmental-duplication region DNA back to the
genome, keep mappings covering >99% of each exon (shared exons), group genes sharing exons with consistent
family copy-number (famCN, MAD<1). Field standard for recent paralogs; needs no SEDEF. We mirror it as an
ORTHOGONAL DNA witness for the RNA-conflict families: extract each copy's **genomic span** from `GGO.fasta`, map
back to the genome (`minimap2 -cx asm20 --eqx -N50 -p0.5`), and test mutual exon-sharing (≥200 bp, ≥90%-id
homologous block onto another copy's locus) + famCN consistency.

**Result (82 families, 207 copies):**
- **DNA shared-exon CONFIRMED: 25/82 (30.5%)** — every copy carries a ≥200 bp ≥90%-id homologous block onto
  another copy's locus (a genuine segmental-duplication signature, computed without the read-conflict graph).
- **famCN-CONSISTENT (MAD<1): 67/82 (81.7%)**; **BOTH: 23/82 (28.0%)**.
- The famCN multiplicities expose real SD structure the RNA missed: GWFAM0 famCN=[9,9], GWFAM28=[11,11] — the
  RNA family captured 2 *expressed* copies of a larger 9–11-copy DNA duplication (expression-restricted, exactly
  the TBC1D3/green-opsin pattern Eichler reports). This is a feature: RNA sees the transcribed subset.

**Interpretation:** ~30% of the read-conflict families are independently DNA-confirmed as segmental-duplication
families by the field-standard (Soto) criterion — a non-circular witness available *without* the cluster-scale
SEDEF run. The other ~70% are either bridged by shorter shared segments (<200 bp domains/repeats — the known
over-merge hazard), more-divergent paralogs (<90% id), or single-locus signals; this is a conservative FLOOR
(the criterion demands mutual ≥200 bp/≥90% homology of every member). Honest caveat: the genomic span shares
exon sequence with the catalog, so it is not fully orthogonal — but the famCN multiplicity and mutual cross-
mapping onto *other* copies' loci are duplication signals the read-conflict graph does not encode.

## F1 — O4 mask-a-copy positive control (`bench/o4_mask_readmit.py`)
**Paper:** Soto 2025 — reference-absent copies are real (37% of human paralogs miss in GRCh38). Our bounded O4
run admitted 0 on the T2T-complete gorilla; this resolves "works-but-complete vs inert" by deleting (N-masking)
known copies from the reference, realigning the locus's reads, and checking O4 re-admits them. Control = the
identical pipeline on the unmasked reference (copies present → must admit 0).

**Attempt 1 — GWFAM71 (9-copy tandem, copies >98% identical):** masked 3 copies (25.5 kb), realigned 263 reads.
**0 admitted; routed to DNA-needs** (4× "<3 clusters", 2× "not min_p-distinct", 2× "≥98% remap identity"). This
is the gate failing **safe and correctly**: at >98% identity a copy is genuinely indistinguishable from a het
without DNA, so O4 refuses to invent it. But it is the wrong test family — near-identical copies are outside the
admission regime by design.

**Attempt 2 — GWFAM47 (6-copy family, copies ~95% identical, the admission regime):** masked copies 0/1/5
(48.5 kb), realigned 1,371 reads — they collapsed onto copy2's locus as expected (1,002 reads at ~122.74 Mb).
**Still 0 admitted, and 0 DNA-needs candidates.** Root cause (from the detection log): with the paralogs masked
out, the collapsed reads assemble into **1 transcript → 1 rep → conflict-graph 0 edges → 0 families**, so the
collapsed-copy discovery — which is **gated behind family detection** (the de-tie conflict graph needs ≥2
*conflicting reference loci*) — **never runs**.

**The real finding (an architectural scope boundary, empirically confirmed):** masking a *whole locus* tests the
**DIVERGENT** reference-absent route (a copy whose locus is entirely missing). But the implementation wires only
the **COLLAPSED** route — a copy hiding *inside* a detected multi-copy family's collapsed locus (the spec marks
the divergent route, "assemble unmapped reads → contig → realign," as DEFERRED). Deleting a locus removes the
very family structure the collapsed route needs, so the reads pile into one consensus and no discovery fires.
This is exactly consistent with the spec's stated scope; the experiment **pins it down on real data**:
- O4 as built detects **collapsed/CNV copies within a detected family**, NOT **de-novo absent loci**.
- The collapsed route's admission IS proven (synthetic `absent_copy_sim.py`: plants a copy *collapsed onto a
  co-located host within a family* → admitted as `AC_`). The real-data analog would require a family where a
  divergent copy is naturally collapsed onto a co-located host — and the H4 bounded run showed GGO's real such
  candidates fail at ≥98% remap (DNA-unresolvable) — so "0 admitted on GGO" = complete reference + conservative
  gate + the divergent route unbuilt, NOT an inert detector.

**Net F1 value:** the positive control did not produce a green "re-admitted" because it exercised the unbuilt
divergent route; but it **root-caused and confirmed O4's exact scope** (collapsed-within-family only) on real
data — a sharper, more honest O4 boundary than before, and a concrete spec for the next step (the divergent
route: assemble the collapsed pile and call it a distinct copy even with one reference locus).

## F3 — DNA-supervised copy decoding = non-circular O2 accuracy (`bench/dna_supervised_decode.py`)
**Paper:** Vollger/SDA 2019 — DNA pre-phases PSVs→copies; our identifiability theory gives the conditions under
which that recovery is exact. The consequence (and the answer to the defense audit's "what's your real-data
accuracy?", since silver is circular and sim5x is synthetic): build per-copy PSV **signatures from the T2T DNA
catalog** (copies AND their distinguishing columns defined independently of the RNA reads), decode RNA reads
against them, and validate with a **held-out-DNA-column cross-validation** — split each read's DNA-defined PSV
columns into a TRAIN half (drives the call) and a TEST half (confirms it). Non-circular on every axis: the
copies are not RNA-assembled, the columns are not RNA-discovered, and there is no minimap2-primary silver label.

**Method:** per co-located DNA-catalog family (≤3 Mb span, 2–8 members), load `{fid}.json` (ref0 + the per-copy
substitution matrix), take the **exonic** PSV columns (genomic pos = ref0 start + local pos), build per-copy
allele signatures, **realign the family's RNA reads to ref0** (one common frame — sidesteps cross-copy
projection), read each read's allele at each PSV position, decode with the production significance gate
(`copy_assign.assign_read`), and run the held-out CV.

**Result (genome-wide, 87 families, 23,504 reads):**
- **Reads decoded (margin-pass): 21,468 / 23,504 (91.3%).**
- **Held-out-DNA-column confirmation: 19,986 / 20,568 = 97.2%** (reads with ≥2 DNA PSV columns), vs a weighted
  1/K chance of 43.2% → **2.2× chance** (enrichment is modest only because the genome-wide set is K=2-dominated,
  where chance is 50%; on the K≥3-richer validation subset it was **95.6% = 3.9× chance**).
- Per-family: DNFAM483 (K=2) 99.8% on 3,769 reads, DNFAM92 (K=3) 100%, DNFAM21 (K=6) 99.9%, DNFAM55 (K=4) 100%;
  the hardest, DNFAM19 (K=8) 77.2% (more copies → more competitors → lower per-column confirmation, as expected).

**Why it matters:** this is the **first non-circular, real-data O2 accuracy**. The silver standard (99.9%) only
shows agreement with minimap2 where minimap2 was already confident, using RNA-assembled profiles — circular on
both the profile and the label. Here the reference is built entirely from DNA, and held-out DNA columns confirm
the per-read call at **97.2%**. It also realises the SDA/Vollger "supervised nearest-signature decode" the DNA-
catalog spec deferred as Phase 2 — turning NP-hard unsupervised phasing into a supervised classification when a
DNA reference exists. Honest scope: co-located DNA-catalog families with ≥2 exonic PSVs; the DNA signatures and
the RNA reads share exon sequence (not orthogonal to *sequence*), but the held-out split + DNA-defined copies/
columns remove the RNA-self-assembly and silver-label circularity that the audit flagged.
