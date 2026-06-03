# Gene-conversion mosaic-read detector (VG "unusual exon combinations")

**Status (2026-06-03):** detector BUILT, correct, opt-in (`RUSTLE_VG_MOSAIC_ON`, default OFF,
additive — EM weights untouched). Per-read core + family aggregation unit-tested (13 tests).
End-to-end on synthetics: correct detection + specificity proven; full per-locus coverage is a
scoped next step (see Limitation).

## Idea
A gene-conversion recombinant read's per-diagnostic-site copy pattern SWITCHES contiguously
(copy A for a run of sites, then copy B) — an "unusual combination" directly observed in ONE
read, not enumerated. Unlike the session's other VG levers (which were bounded by the
identifiability limit), this DISCOVERS real biology that is literally in the reads, and abstains
where there is no diagnostic signal.

## Design (synthesized from a 3-lens design panel: statistical / algorithmic / biological)
- Pure `detect_mosaic(obs, n_copies, eps, params)` (src/rustle/vg_hmm/mosaic.rs). `obs` = per-site
  cross-copy match vector (built from `fp.per_copy_site_refs` + the read's mismatches, mirroring
  `score_read_exon_fingerprint` but RETAINING the per-site pattern the scorer discards).
- Tokenize: Copy(c)=exactly one copy matches; Ambig=≥2 (non-identifiable here); Novel=0 (wildcard).
- Gates (abstain first): NonIdentifiable (ambig fraction > 0.5), LowPower (< 6 decisive sites).
- Switch: exhaustive single-changepoint scan over decisive sites maximizing explained sites for
  (copy A before k, copy B after k). With decisive-site agreements the LR reduces to
  `2·margin·ln((1−eps)/eps)`; the χ²(2df) null quantile is exactly `−2·ln(α)`, so the threshold
  `bic·ln(S) − 2·ln(α_corrected)` is closed-form (no bootstrap).
- DUAL GATE (call Mosaic): LR > threshold AND hard backstops — each tract ≥ 3 decisive sites,
  integer margin ≥ 3, per-tract purity ≥ 0.85. The LR calibrates; the integer/run-length
  backstops hold even under correlated/homopolymer error bursts (the i.i.d. null's weakness).
- Family aggregation (`aggregate_family`): a genuine conversion RECURS at a fixed breakpoint
  across independent molecules; clusters of < `family_min_supporting_reads` (3) or with dispersed
  breakpoints are ChimeraSuspect, not Confirmed. Only Confirmed events would ever seed a
  recombinant transcript (the no-fabrication guard — one unreplicated switch is never promoted).

## Validation (synthetic ground truth: gen_synthetic.py --recomb-reads/pair/exon)
- Per-read: planted recombinants are detected with the CORRECT oriented pair and a breakpoint
  bracket containing the planted crossover (e.g. exon2→exon3: bracket 6679–9615); margin/LR strong.
- Specificity: pure-copy reads → SingleCopy, **0 false mosaics**; error-flip / low-power /
  non-identifiable → abstain (unit tests).
- Family CONFIRM logic: unit-tested (3 reproducible reads → Confirmed; singleton → ChimeraSuspect).

## Limitation (scoped next step, NOT a correctness issue)
The pass is hooked in `run_fingerprint_em`, where the family bundles are reduced to the
discovery-LINKED (multi-mapping) reads — so it sees only those reads, not all primary reads at a
copy's locus. On the synthetic, recombinants that don't multimap (their 3′ tract makes the
sister-copy secondary fall below minimap2 `-p0.5`) never enter the family, so only ~1 reaches
detection — too few for the in-binary ≥3 CONFIRM (the CONFIRM logic itself is unit-proven).
**Next step:** run the pass at the assembly-bundle stage (full per-locus reads) with the
family's `fp`, so unique recombinants are also seen and real conversion events can Confirm
end-to-end. Then audit on real DAZ/RBMY/NBPF (default-OFF must stay byte-identical).
