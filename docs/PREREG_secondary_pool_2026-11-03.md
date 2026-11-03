# PREREG SECONDARY-POOL CANDIDATE COPIES (register 1513, 2026-11-03)

Why (user, 2026-11-03): IsoCon's reference-free strength is its significance test on a pool of reads from one locus. Our port
(`candidate_copies`, register 1401: 4 reference-absent expressed copies found only with the IsoCon side) pools every primary
AND every secondary record on a family's copies, and the BAMs carry `-p 0.1 -N 50` secondaries, so the pools hold mostly
uninformative records (r1511 measures how many). This row runs the same merge on pools restricted to the secondaries that
carry information: the reads minimap2 could not place, pooled by the rival locus they land on.

Tool: `candidate_copies --molecules P.molecules.tsv --min-ms-ratio R` (new; unset = unchanged): a secondary record names its
read for a net only when its ms >= R x the read's primary ms (as_table column `primary_ms`). R = the r1511 bound (the lowest
ms-ratio bin whose useful fraction is >= 5%), fixed by that row's reading rule before this one runs; if r1511 gives no bound,
R = 0.8 (minimap2's default `-p`).

Run: KB3781, the r1509 40-family panel (`~/o2contested/panel/families.txt`, `--families`), wall limit 2 h per arm, both arms
the same binary and delta: A = default pool, B = `--min-ms-ratio R`.

Scoring (fixed now): every flagged contig (`candidates.tsv` flagged = 1, sequence in `contigs.fa`) aligned with minimap2
`-x splice:hq -uf --secondary=no` to Jim's haplotypes (`o2contested/haps/KB3781.{mat,pat}.mmi`), best hit by matches:
- `absent_from_ref`: identity >= 0.998, read coverage >= 0.9, locus outside every catalog copy's haplotype locus
  (`haplotype_truth` hits `C|family|copy`, file `truth_all10/hits.tsv`, g = 0) -> a real copy the reference lacks;
- `allele_of_known_copy`: identity >= 0.998 on a catalog copy's locus -> a false new-copy flag;
- `not_in_genome`: best identity < 0.995 -> junk.
Decision rule: arm B keeps >= 80% of arm A's `absent_from_ref` contigs (by haplotype locus, 1 kb tolerance) with a false-flag
count (allele + junk) no higher than A's, and runs in <= 50% of A's wall: then the informative-secondary pool is the pool,
and the flag becomes the driver's default for the o3 stage. Fewer true contigs or more false ones: FAIL, flag stays opt-in.
Fewer than 3 `absent_from_ref` contigs in arm A: UNDERPOWERED, report only. Report beside it: net sizes per family in both
arms, contigs per class, wall.
