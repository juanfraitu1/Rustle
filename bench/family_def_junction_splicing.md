# Can the VG / junction structure define gene families? (precision, not just resolution)

Synthesis of 4 angles (junction-architecture, genome-scale, vg-cothreading, splicing-unique),
each adversarially verified (confound/base-rate + scale/generalization lenses).
All headline numbers below were re-reproduced (this session) from the live scripts/data.

## Bottom line

**Q1 (VG for family DEFINITION): PARTIAL / mostly NO at scale.** The VG/junction-backbone is the
*correct structural object* for stating the criterion (Canzar-shaped, threshold-light), and it gives
one clean by-construction precision bit. But as a genome-wide family *caller* it does NOT separate the
real over-merge problem: no single VG/architecture operating point both cuts the 389-gene domain-bridge
mega-family AND keeps the validated real families connected.

**Q2 (a UNIQUE junction signal orthogonal to sequence id): YES, but narrow.** Exactly one clean,
threshold-free, sequence-orthogonal signal survives verification: **retrocopy intronlessness**. It flags
processed-pseudogene copies that cDNA identity provably cannot (id ~0.95-1.0). Everything broader
(continuous intron-count / intron-length / exon-count "architecture concordance") is either an arbitrary
panel-fit threshold, a gene-size proxy, or actively destroys real diverged paralogs.

cDNA id is *anti-discriminative* on this panel (AUC real-vs-homologous-counter = 0.029): the false
counters have id ~1.0, higher than the diverged real families (0.90-0.99). So an orthogonal signal is
genuinely needed — junctions provide one, but only for the retrocopy sub-problem.

## The single cleanest unique signal: retrocopy intronlessness (S1)

Definition: an edge (A,B) with cDNA id >= 0.9 where **min(n_intron(A), n_intron(B)) == 0 while the other
is spliced**. Exact from the GFF exon count; no threshold.

What it buys (that sequence id cannot):
- Flags processed pseudogenes / retrocopies. Panel: EEF1A1 (7 introns) ~ LOC retro (0 introns) id=0.98;
  CNN2 (6) ~ LOC retro (0) id=0.99. 2/2 retrocopies flagged, **0/5 false flags on real families**.
- The edges it cuts have median cDNA id ~0.947 (range 0.90-1.00) — sequence id is blind to them.
- Mechanistically clean & VG-native: a 0-intron copy **provably cannot thread a multi-exon spliced
  backbone** (no bubbles to traverse) → structural, by-construction exclusion, not a tuned filter.
- NOT a coverage proxy: retrocopies have high cDNA coverage (0.43-1.00).

Where it FAILS / its scope limits:
- It does **NOT** touch the actual over-merge. The 389-gene mega-family is built from **multi-exon ↔
  multi-exon DOMAIN BRIDGES**. The strict intronless-only filter (Cspliced) cuts only **196/1522 = 13%**
  of mega edges and shrinks the mega merely 389 → a **238 + 201 + 140-gene** triple blob. The domain
  bridge survives essentially intact. (Reproduced this session.)
- **Critical correction (verified):** "intronless ⇒ retrocopy ⇒ not a family" is biologically FALSE.
  Of the 5198 either-intronless edges only **1931 are the genuine one-side-intronless (retrocopy)
  pattern**; **3267 (63%) are BOTH-intronless** and contain real single-exon tandem families:
  **PCDHB (33 edges), IFNA (26), HNRNPCL, ELOA...** Dropping every intronless-incident edge therefore
  destroys real intronless tandem families ~10:1 over true retrocopies. **The valid signal must be
  restricted to one-side-intronless-vs-spliced-parent (1931 edges), NOT all-intronless.** The panel hid
  this because all 3 panel real families are multi-exon (zero intronless real families in the panel).

So S1 (restricted to one-side-intronless) is a real, clean retrocopy filter — orthogonal to id — but it
is a side problem, not the domain-bridge over-merge.

## The VG combinatorial criterion (Canzar-shaped) and what the empirics say

Criterion (threshold-light, provable-shaped):
> A candidate set S is a **family** iff its members are **parallel s-t paths through one variation graph
> that share its junction backbone** — i.e. there is a single bubble-chain spine whose branch nodes
> (junction boundaries) are co-threaded by all members in the sequence-aligned frame. A domain-sharer
> threads only a sub-path / lands in a separate component; a retrocopy is a single straight-through path
> with no bubbles and cannot thread the spine at all.

Operationalized (A3) as `shared_backbone_fraction` = (# reference-spine junctions co-threaded by all
members, in the cDNA-aligned frame, position-matched within TOL bp) / (# reference junctions).

What the empirics support:
- **Panel (n=9): clean.** REAL = {RABL2 0.89, APOBEC3 0.67, RFPL 0.50}; ALL domain-sharers and retros
  = 0.00. Margin REAL_min − COUNTER_max = **+0.50**, and it is **tolerance-insensitive** (margin stays
  +0.50 for TOL in {0,4,8,...,50} bp, even exact TOL=0). It beats intron-count-ratio (margin −0.167)
  and it is NOT a coverage proxy (APOBEC3-real at alnCov 0.39 beats CASP8-domain at alnCov 0.46 — the
  discriminator is junction *position*, not coverage). The mechanism is the right one: **position is the
  correct invariant; intron/exon LENGTH is the wrong invariant** (real paralog introns drift in length).

What the empirics do NOT support (verifier-downgraded from "strong" to "promising on n=9"):
- **The fractional position metric was never deployed at scale.** The genome-wide numbers substitute the
  exact intronless axis + a cDNA-cov<0.5 proxy, not the validated junction-position projection.
- **The strict genome-scale architecture form (Carch = equal intron count + ordered length concordance)
  DESTROYS all three real families** (verified this session): under Carch the mega breaks 389 → 32, but
  RABL2 keeps only RABL2A, and **APOBEC3 and RFPL have ZERO members left in the graph**. Reason: real
  paralog intron lengths/counts diverge post-duplication (RABL2A/B introns 455 vs 214, 4213 vs 2858;
  RFPL counts 3/4/2).
- **The continuous count_conc "family gate" (>=0.6) is panel-overfit and breaks a held-out positive.**
  The single Compara-validated paralog held out of the panel, **GGT1~GGTLC2 (cDNA id 0.98), has
  count_conc = 0.267** (15 introns vs 4) → its only edge is dropped → family destroyed (false negative).
  RFPL2~RFPL3 = 0.50. A real partial-duplication paralog has the *same* architecture signature
  (mismatched intron counts, asymmetric coverage) as a domain-sharer (count_conc up to 0.67), so the gate
  provably cannot separate them. The 0.6 threshold is non-eliminable and chosen only to keep RFPL
  connected — exactly the magic number Canzar distrusts.
- **exon-count concordance (S5, panel AUC 0.929) is a gene-SIZE proxy** (corr |dExon| vs max-exon-count
  = 0.71) and degrades to AUC ~0.755 genome-wide; it fails on size-matched domain-sharers (ASDURF dExon=1,
  GPR39 dExon=2 fall inside the real range).
- **Two goals are in direct tension:** the real-family-safe gate leaves a 204-238-gene blob and
  collaterally damages ~65% of comparable size-2/3 families; the blob-breaking gate annihilates the real
  families. No single operating point does both.

## Decisive vs suggestive

DECISIVE (survives both verifier lenses, reproduced):
- cDNA id is anti-discriminative for this over-merge (AUC 0.029); an orthogonal signal is needed.
- **One-side-intronless retrocopy flag** is a clean, threshold-free, sequence-orthogonal precision bit:
  1931 such edges genome-wide, median id ~0.95, 2/2 panel retros, 0 panel real. VG-native exclusion.
- The retrocopy flag does **not** solve the domain-bridge mega (only 13% of mega edges; 238+201+140 blob).
- Strict architecture concordance (Carch) destroys all 3 real families — rejected.
- The continuous count_conc gate (>=0.6) is panel-fit and produces a false negative on the held-out
  Compara paralog GGT1~GGTLC2 (cc=0.267).
- "all-intronless ⇒ retrocopy" is false; the both-intronless half (3267 edges) holds real tandem
  families (PCDHB, IFNA, ...). The signal must be the one-side-intronless restriction only.

SUGGESTIVE (clean on n=9, unproven at scale):
- The **junction-POSITION co-threading fraction** (shared_backbone_fraction, margin +0.50, tolerance-free)
  is the most promising VG-native criterion and the correct invariant (position not length), but it was
  never run at genome scale; the deployed proxy contradicts it. This is the one thing worth building.

## Honest answer

- VG helps DEFINE families only partially: it supplies the right *structure* and one clean by-construction
  bit (retrocopies can't thread a spliced spine), but it does not, at any tested operating point, cleanly
  cut the multi-exon domain-bridge over-merge while preserving real diverged paralogs.
- The unique junction signal orthogonal to sequence id is real but narrow: **retrocopy intronlessness**
  (restricted to one-side-intronless), which sequence id cannot see. The broader "architecture
  concordance" story does not survive scale.

## Recommended next step

Build and genome-wide-deploy the **junction-POSITION co-threading metric** (A3's
`shared_backbone_fraction`, position-matched in the cDNA-aligned frame, tolerance-free) as the actual
VG criterion — it is the only candidate that is both Canzar-shaped and panel-clean — and test whether a
*position*-based (not length/count) backbone-overlap fraction can cut multi-exon domain bridges
(partial-colinear / sub-path alignments) while keeping length-drifted real paralogs (RABL2/APOBEC3/RFPL,
GGT1~GGTLC2) connected. Ship the one-side-intronless retrocopy flag now as a standalone, by-construction
precision filter (1931 edges) — independent of whatever the position metric shows.

Artifacts: family_def_backbone_vgspine.py (A3, position metric, panel),
family_def_junction_genomewide.py (Cspliced/Ccount/Carch genome-wide),
family_def_intron_index.py (builds /home/juanfra/winloci_scratch/gene_intron_index.json),
family_def_splicing_signals_a4.py (S1-S5), family_def_junction_lever.py (count_conc gate),
compara_paralog_relation.json (held-out positives incl. GGT1~GGTLC2).
