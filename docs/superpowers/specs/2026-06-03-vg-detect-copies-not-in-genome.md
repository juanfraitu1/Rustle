# Detecting copies NOT in the reference genome (VG)

**How can a copy not be in the genome?** A gene-family copy can be present in the sequenced
individual yet absent from the reference assembly: a collapsed segmental duplication, an
assembly gap, copy-number variation (the individual carries more copies than the reference),
or a private/recent duplication. Segdup-rich families (GOLGA8, DAZ, RBMY) are exactly where
references are most incomplete.

**The detectable signal (proven).** A hidden copy's reads have no correct home, so they mismap
to the closest sibling carrying their PRIVATE SNPs — a COHERENT second haplotype the reference
(one copy at that locus) cannot explain. Ground truth: gen_synthetic --hidden-copies C (reads
from copy C present; its genome locus replaced by random). Probe (bench/hidden_copy/probe.sh):

  hidden:  copy0-locus depth=60  alt-haplotype-positions=52
  visible: copy0-locus depth=30  alt-haplotype-positions=0

So a hidden copy leaves a clean signature: ~2x depth + dozens of coherent 20-60% alt-allele
positions (the hidden copy's gene-wide divergence), vs ~0 when every copy is in the reference.
The signal cleanly separates a hidden COPY (dozens of co-occurring alt positions) from a single
heterozygous SNP (1-2) or sequencing error (0).

**Discipline (DAZ3).** The honest goal is to DETECT and FLAG the discrepancy — "the reads imply
≥N+1 copies; the reference models N" — NOT to fabricate the missing copy's sequence/location
(the phantom trap the gated-off discover_novel_copies rescue path risks). Detect, abstain on
placement.

**Status:** infrastructure built (synthetic + probe), signal PROVEN. Honest detector next.

## Honest detector BUILT (2026-06-03)
vg_hmm/hidden_copy.rs: pure detect_hidden_copy over PRIMARY reads at a reference-copy locus.
Design = independent design panel (statistical/algorithmic/honesty lenses) synthesized.
- PRIMARY-reads-only matrix = the paralog-bleed firewall (an in-reference paralog's reads are
  primary at ITS locus, secondary here, so excluded). Candidate columns = alt-allele frac in
  [0.20,0.60] (≫ 0.5% error; below fixed/ref-error). Co-segregation: a hidden copy's alt
  columns co-occur on ONE read subset (the alt haplotype H). FLAG iff ≥12 candidate columns
  (het firewall — a het is 1-2) AND ≥5 alt-haplotype reads.
- 5 unit tests: hidden-copy flagged; the 3 FP negatives REJECTED (sequencing error → 0
  candidates; heterozygous SNP → <12 positions; fixed difference → above the band); low-depth
  abstains.
- Wired per-copy-locus (opt-in RUSTLE_VG_HIDDEN_COPY) → [VG-HIDDEN] trace + GTF attrs
  hidden_copy_evidence / hidden_copy_alt_positions / hidden_copy_read_frac. Default-off
  byte-identical.
- DAZ3 discipline: DETECT + FLAG only ("reference models 1, reads imply ≥2") — NEVER places or
  fabricates the copy.
- Validated end-to-end: HIDDEN (copy2 absent from reference) → flagged, 49 alt-haplotype
  positions, read_frac 0.38; ALL-VISIBLE → 0 positions, not flagged (siblings bleed as
  secondary, filtered). Suite green.
- Limitation: the pass sees the family-REDUCED reads (8-13 vs ~60 at the locus), enough to
  carry the signal here but lower-power; full-read coverage (BAM re-fetch, like the mosaic
  pass) is the robustness follow-up. The synthesized spec also has further rigor not yet
  implemented (per-read Poisson error-null, LLR model-comparison, span/concentration checks).
