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

## RESULT (2026-11-03, job 4421549, R = 0.5; workspace o2pool/, `scoreA/B.contigs.tsv`)
Arm A (default pool): pass A 188,284 primary + 434,429 secondary records; 25 flagged contigs in 13 families: absent_from_ref 3
(two loci: CM054600.2:6.73-6.76 and 6.77-6.80 Mb, families MCL239 / MCL240), allele_of_known_copy 1, not_in_genome 20,
between 1; wall 416 s. Arm B (ms ratio >= 0.5): 327,226 secondary records kept, 107,203 dropped; 28 flagged in 15
families: absent_from_ref 6 (both A loci kept = 100%, plus CM054564.2:48.57-48.67 Mb in MCL172 and MCL84, and
CM054582.2:127.83 Mb in MCL1), allele 1, not_in_genome 18, between 3; wall 447 s.
Bars: (1) B keeps >= 80% of A's absent_from_ref loci: PASS (2 of 2, and 2 more loci); (2) false flags 19 <= 21: PASS;
(3) wall <= 50% of A: FAIL (447 vs 416 s: the nets step is primaries-bound, the dropped secondaries cost nothing).
Verdict as registered: ❌ FAIL on the wall bar; the pool itself is better on both quality bars (twice the reference-absent
contigs, fewer junk). The flag stays opt-in; a follow-up may re-register without the speed bar, on more families, with the
wall reported, not barred. Caveat: 14 of the 38 not_in_genome contigs sit on two recurring loci (CM054591.2:75.61 Mb,
CM054595.2:100.75 Mb) shared by unrelated families: repeat-driven clusters, a filter candidate for `candidate_copies`.
