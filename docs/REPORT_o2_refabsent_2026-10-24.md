# Report-only: O3 -> O2 on KB3781's real reference-absent copies (refabsent panel), 2026-10-24

**Status:** written before the port is run on these families. The refabsent panel (34 families with haplotype-only loci, truth
from KB3781's own haplotypes, registers 1205-1207) has been used before (O3 stage, 2026-10-01), so this is a REPORT, not a test:
no verdict, nothing tuned on it.

**Method.** `o3port.py --semi --guard` as of 11779a31 (P14's method), on `refabsent/R0.bam` (fibroblast reads on `_pri`), copies =
`refabsent/copies_genes.tsv`, read set = `refabsent/fam/<family>.fa`. Script: `bench/o3_graph/o2refabsent.py`.

**Truth.** Each read's best record on the haplotype the family's copies lift to: on a lifted copy -> that copy; on a B-only locus
-> absent; a tie between two loci -> unknown (dropped). Within-delta B-only loci (114 of 127) are near-identical to a `_pri` locus,
so their truth is close to allelic and is reported apart from beyond-delta loci.

**Reported.** (1) Per family: expressed B-only loci vs the port's certified / putative flags. (2) O2 on tied reads from reference
copies: placed, wrong, vs minimap2 MAPQ >= 1. (3) O2 on tied reads from B-only loci: placed anywhere but a certified new copy =
wrong; vs minimap2 MAPQ >= 1 (always wrong there).
