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
