# PREREG — roster admission by read origin (O2 → O1 feedback, §6fc), 2026-09-05
Input: the sweep on `rna_units_v9` (read-star default). For every family, the unit reads the origin certificate
rejects whose PRIMARY lies outside every unit of the family, clustered into loci (≥ 3 primaries within 5 kb).
Admission: each such locus becomes an ADDITIONAL COPY of the family's O2 candidate set — its exon chain built
from those primaries by the unit rule (a base is exonic if ≥ 3 reads cover it and more cover it than splice over
it), its sequence from the genome, its strand the reads' majority — flagged `read_admitted`. This is a candidate
for O2 (the read's origin must be in the set for the certificate to mean "no copy explains the read"), NOT a
family member: membership stays with the core rule (`--sedef`). Non-circular: the locus is proposed by where the
ALIGNER put the reads (their primary), the family by sequence homology.
Predictions on the families that receive ≥ 1 admitted locus (paired, same reads):
- P1 origin-rejected unit reads fall by ≥ 50 % in those families.
- P2 MAPQ-60 placement agreement stays ≥ 99.5 % and the reads at an admitted locus are assigned TO it (≥ 90 %).
- P3 NPIP's 62 audited anchors: 0 wrong.
- P4 total assigned across the sweep rises (the admitted loci's reads were tied or rejected before).
Fail rule: P3 fails ⟹ admitted loci must not enter O2's set without the core rule; P2 fails ⟹ the chain rule on
primaries alone is not a unit and the locus needs the annotation.
