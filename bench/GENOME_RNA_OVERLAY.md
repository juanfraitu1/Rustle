# Genome-grounded vs RNA-conflict families: cross-modal overlay

**What this is.** A test of whether grouping de-novo RNA loci into "copy families" by
read/alignment evidence agrees with a *genome-grounded* family definition built from
segmental-duplication structure — and, critically, whether that agreement is more than
genomic proximity would produce by chance.

**Artifacts**
- `bench/genome_rna_overlay.py` (overlay + null + lift)
- `bench/genome_rna_overlay.json` (committed numbers; all figures below trace to it)
- `bench/genome_rna_overlay_per_family.tsv` (per-family breakdown)
- Genome side: `bench/genome_family_def.py` → GFAM (nodes = NCBI genes; edge = one SEDEF
  segdup pair covers gene *u* on side A and gene *v* on side B, each ≥50%; family =
  connected component, |C|≥2). Reads `final.bed` (SEDEF self-alignment) + `GGO_genomic.gff`.
  **Never opens a BAM.**
- RNA side: `denovo_families.tsv` (de-novo locus catalog, members `DN_<contig>_<pos>_<idx>`),
  joined to `denovo_transcripts.meta.tsv` for `[start,end)` intervals.

---

## 1. The definition mismatch — and why it is the point, not a bug

The two catalogs do **not** mean the same thing:

| | Genome side (GFAM / SEDEF) | RNA side (denovo_families) |
|---|---|---|
| Unit | annotated gene / genomic interval | expressed de-novo locus |
| Edge | a still-detectable **segmental duplication** block links two loci | the two loci are **grouped by sequence-similarity** of their de-novo transcripts (18-mer share + POA `core_recip` ≥ 0.13) |
| Modality | genome self-alignment (homology of *DNA blocks*) | RNA-discovery grouping (loci that the read evidence cannot keep apart) |
| Bias | recency / high-identity SDs only | expression + multi-mapping ambiguity |

So a within-family RNA pair is the hypothesis "these two expressed loci are copies of one
another"; a SEDEF link is the independent genomic statement "a duplicated block physically
joins these two regions." Asking *what fraction of within-RNA-family pairs are SEDEF-linked*
is therefore a **precision-style cross-modal check**, not a like-for-like reproduction.

**Why precision + distance-matched null + lift is the right test.** SEDEF density is
strongly distance-dependent (tandem-rich at short range, near-zero across chromosomes).
A bare "X% of family pairs are segdups" number is uninterpretable because some of that is
just proximity. The airtight version fixes the question to *excess over chance*: for every
within-family pair, draw a random non-co-member pair **matched on contig and distance bin**,
under the **identical** overlap rule, and compare. Lift = observed / matched-null. The null
collapses to ~1 under random family relabeling, so a large lift is the signal.

---

## 2. Headline numbers

Within-family pairs are **distinct mapped loci** (overlapping isoform-split intervals merged
to one locus); overlap rule = ≥50% of the locus covered by a SEDEF block (pre-registered).

- **SEDEF-linkage precision (excl. DNFAM0) = 0.270 = 698 / 2585** — primary headline.
  Computed **exhaustively (no sampling)**; sampling-stable.
- All-family precision = 0.125 = 698 / 5585 (mixes in DNFAM0).
- **No-merge lower bound = 0.188 = 812 / 4328** (isoform-split contamination *included*).
- DNFAM0 = the documented Alu/repeat over-merge: 636 distinct loci, **0 / 3000 sampled
  pairs linked** — the clean over-merge signature, correctly excluded from the headline.
- **Distance- & chromosome-matched null background = 0.00061** (61 / 100 000).
- **Lift** (all-family) = 0.125 / 0.00061 = **204.9×**. (The excl-DNFAM0 lift is corrected
  below — see Limitations FLAW 1.)

**Per-stratum precision vs matched null (the mechanistically informative cut):**

| Stratum | Precision | Null | Lift |
|---|---|---|---|
| cis <10kb (genuine tandems) | 0.73 (11/15) | 0.0034 | 216× |
| cis 10–100kb | 0.46 (96/211) | 0.0115 | 39.6× |
| cis 0.1–1Mb | 0.43 (92/214) | 0.0016 | 270× |
| cis 1–10Mb | 0.36 (137/377) | 0.0016 | 223× |
| cis >10Mb | 0.117 (63/539) | 0/9771 | (null≈0) |
| trans (cross-chrom) | 0.071 (299/4229) | 0/75678 | (null≈0) |

The short-range tandem stratum is the cleanest: 73% of adjacent non-overlapping family
tandems are SEDEF-confirmed against a 0.3% background.

**Significance (honest form — see FLAW 2):** permutation floor **p < 0.001** (1/1001;
max null replicate rate 0.003 vs observed 0.125), plus a **family-level** test: 171 of 538
families are SEDEF-coherent vs ≈2.4 expected under the null (p≈0), with a family-clustered
bootstrap precision CI of [0.184, 0.377] (excl-DNFAM0).

---

## 3. Conditioned recall (and its load-bearing caveat)

**Conditioned recall = 0.608 = 113 / 186.** Conditioning chain: of 538 protein-coding GFAM
families, **186** have ≥2 member genes that are RNA-expressed (a de-novo locus overlaps the
gene); of those, **113** have ≥2 expressed members landing in **one** RNA family.

This is **not** a raw recall and is **not** comparable to 1. Uniquely-mappable genome
paralogs produce no read conflict and are *correctly* not grouped by the RNA side, so they
are excluded by design. Reported honestly as a denominator chain:
- raw vs all 538 PC families = 113/538 = **0.21**;
- raw vs all 253k SEDEF pairs = **~0.10** (most genome paralogs are unexpressed or
  uniquely-mappable). This is the structural floor, not a defect.

---

## 4. Reconciliation with the prior ~0.64 / 0.72 / 0.10

Two corrections up front:

1. **The task's premise is wrong on provenance.** The prior 0.64/0.72/0.10 did **not** come
   from a SEDEF or DFAM catalog. It came from `bench/family_def_dna_pr.py`, whose DNA truth
   is **all-vs-all cDNA homology** (minimap2 asm20, id ≥ 0.90). The "DFAM" expression overlay
   (`dna_rna_overlay.py`) is mmseqs **protein clustering** and reported *no* P/R at all.
2. **The cDNA-homology truth is partly circular.** It scores the same transcript sequence
   lineage the reads derive from; SEDEF self-alignment never sees a read.

Given that, the numbers reconcile by design:
- **Precision is *expected* to drop, 0.64 → 0.270.** SEDEF is a *stricter, narrower* oracle:
  it only catches recent high-identity SDs and structurally misses ancient / WGD /
  single-gene paralogs (RAB/RHO, keratins) and any duplication whose segdup block has
  diverged away. A pair can be a true cDNA-homolog (prior TP) yet have **no segdup**
  (SEDEF "miss"). So a lower precision against SEDEF is the correct, not the worse, number.
- **The genuine upgrade is the LIFT**, which the prior never computed (it had no null).
- **Recall reproduces.** Conditioned recall 0.608 closely matches the prior's conditional
  recall on the cross-mapping subset (0.658); the prior's raw ~0.10 reappears here as the
  same structural floor.

Net: this **replaces** the partly-circular 0.64/0.72/0.10 with a same-shape metric
(precision-as-lift + recall-conditioned-on-expression) — *but see the circularity finding
below, which qualifies how "non-circular" the replacement actually is.*

---

## 5. Limitations (every refuter finding, surfaced not buried)

### L1 — MAJOR: shared-substrate leak; "two independent estimators" is overstated
The adversarial provenance trace found a **positive leak**. The "de-novo transcript"
sequences that the RNA grouping scores (`core_recip` POA + 18-mers) are **reference-genome
substrings**: `denovo_assemble_gate.py` fetches `GGO.fasta` at read-discovered exon
boundaries (`seq = "".join(fa.fetch(c, xs, xe) ...)`). Reads contribute exon coordinates,
count, and strand — **never sequence content** (there is no read-consensus / POA-of-reads
step). SEDEF *also* scores reference self-homology. So **both sides score the same reference
bytes**; the grouping is a reference-self-homology detector and SEDEF is a reference-self-
homology detector, and segdup blocks largely subsume exons. The distance/chrom-matched null
draws reference-*dissimilar* pairs, so it **cannot** absorb this shared-substrate inflation.

Consequence: the 200–440× lift is **not** independent corroboration *via sequence content*.
The framing "two independent estimators of the same biology / genuinely non-circular oracle"
must be **retracted**. The corroboration is **not fatal** — residual non-circular signal
survives: the **reads define which loci exist** and the **exon architecture** (spliced-exonic
vs unspliced-genomic differ), and DNFAM0's 0/3000 shows the two detectors do diverge — so it
is not a pure tautology. But the load-bearing claim moves from "independent lift" to
"expressed, read-defined, exon-architected loci coincide with reference segdups, over a
shared reference substrate."

### L2 — MAJOR (statistics): two real defects in the reported numbers
- **FLAW 1 (~2×):** the excl-DNFAM0 lift of **442.7×** divides excl-DNFAM0 precision (0.270)
  by the **all-family-mix** background (0.00061). DNFAM0 is 86.7% trans (a zero-null
  stratum); removing it shifts the observed mix toward nonzero-null cis bins, so the
  correctly stratum-matched excl-DNFAM0 background is **0.00120**. The honest excl-DNFAM0
  lift is **~225×, not 442.7×.** (The all-family 204.9× is correctly matched.)
- **FLAW 2 (pseudoreplication):** Fisher / binomial **p < 1e-300 are invalid.** The 5585
  pairs are k(k−1)/2 overlapping within-family pairs over ~538 families (DNFAM0 alone = 3000
  of them), not independent Bernoulli trials. Effective N is the **families**. Honest
  significance = permutation floor **p < 0.001** + cluster-level test (171 coherent vs ~2.4
  expected) + family-clustered bootstrap CI. **Retire the p<1e-300 claims.**
- Robustness that *did* survive: cis-only lift (drop the zero-null trans stratum) = **117×**
  with a nonzero matched null; per-cis-bin lifts 40–270×; not family-dominated (top family
  = 11.2% of links; leave-one-out moves precision 0.125→0.113).

### L3 — MINOR (mapping): one disclosure spin to correct
- Coordinate integrity is **clean** (0 mismatches over 102,455 meta rows; meta and
  `final.bed` both 0-based half-open; contig names identical; mito excluded). 10 unmapped
  members dropped (negligible).
- The `merge_loci` step is the one real magnitude lever (0.270 → **0.188** no-merge). But it
  is a **correct same-locus isoform dedup**: removed pairs link at 6.5% vs 27% for retained,
  so it *de-contaminates*, it does not inflate (union-span over-merge actually *deflates*).
  **Fix the wording:** caveat #4's "conservative (lowers n_pairs)" spins it — report
  **0.188 (no-merge) as the contamination-inclusive lower bound** alongside 0.270.
- The `.loci` files cover only `NC_073243.2`; they are **not** a genome-wide coordinate
  oracle (the genome-wide locus→interval map is validated by id==meta identity instead).

### L4 — Catalog identity (pre-existing caveat, restated)
The object tested is the **POA-similarity O1 catalog** (`core_recip ≥ 0.13`,
`denovo_families.tsv`), **not** the shipped per-region de-tie read-conflict graph
(`read_conflict.rs`, which runs per-region and never genome-wide). The POA catalog still
contains the over-merges (DNFAM0). Claims here are about the O1 similarity catalog.

### L5 — Definitional and oracle ceilings
- Precision conflates two failure modes: (a) genuine RNA over-merge, (b) SEDEF
  incompleteness (a correct copy pair too diverged for SEDEF). So absolute precision is a
  **lower bound on agreement**, not an RNA error rate.
- SEDEF can never link ancient/WGD paralogs, structurally capping precision below 1.
- The annotation variant (0.1255) is **not independent** of the headline (GFAM is
  SEDEF-derived) and is symbol-keyed; a granularity cross-check only.
- Null co-member exclusion is approximate (only same-RNA-family draws excluded); a random
  null pair could be a real ungrouped paralog, which would *understate* the lift —
  conservative.

---

## 6. Verdict

**VERDICT: DEFENSIBLE ONLY AS A REFRAMED, SHARED-REFERENCE-SUBSTRATE corroboration of O1 —
NOT as the headlined "two independent estimators / 200–440× independent lift"; the RNA
grouping scores reference-genome substrings (the same bytes SEDEF scores), so the lift is
substrate-inflated and must be retracted as "independent." The real, surviving non-circular
signal is that read-defined, expression-conditioned, spliced-exon-architected loci coincide
with reference segdups far above a distance/chrom-matched null — corroborative but partial.**

### Must-fix items
1. **(L1, MAJOR)** Retract "two independent estimators / genuinely non-circular oracle."
   Either (A) rebuild the RNA catalog from **read-consensus** sequences (POA over the reads,
   carrying copy-distinguishing PSVs) and re-run the overlay — the only route to a genuinely
   cross-modal lift — or (B) re-characterize as "expressed, exon-defined, **reference-
   homologous** loci coincide with reference segdups," with corroboration carried by
   expression-conditioning + spliced-vs-genomic structure over a **shared reference
   substrate**, and stop citing 200–440× as independent-of-genome. At minimum, state plainly
   that the RNA "transcript sequences" are reference substrings.
2. **(L2/FLAW 1)** Recompute the excl-DNFAM0 lift against its **own** stratum-matched
   background (0.00120) → report **~225×, not 442.7×**. The all-family 204.9× is fine.
3. **(L2/FLAW 2)** Drop Fisher/binomial **p < 1e-300** (pseudoreplication). Report
   permutation floor **p < 0.001** + cluster-level test (171 coherent vs ~2.4 expected) +
   family-clustered bootstrap CI; lead with the fair cis-only 117× / per-cis-bin 40–270×.
4. **(L3, minor)** Report the no-merge **0.188** as the contamination-inclusive lower bound;
   fix caveat #4's "conservative" wording (the merge *raises* the headline); note `.loci`
   is not a genome-wide oracle.
5. **(L4)** State that the tested object is the POA-similarity O1 catalog, not the shipped
   per-region `read_conflict.rs` graph.
