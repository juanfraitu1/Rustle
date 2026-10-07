# Which alignments seed the loci, and what a primary-first representative adds: NPIP and TBC1D3 on real reads (2026-10-07)

**Status: DESCRIPTIVE, human A119b DEV (the development block of both families), gorilla OR6737 second-species check (section 5). Protocol: `docs/PREREG_seed_pool_real_reads_2026-10-07.md` (committed 94e69328 before any arm existed; Amendment 1, a27786a6, before any gorilla arm). Switch: `tools/rustle_pipeline.sh --seed-pool primary|good|all`, `--seed-as-ratio`, `--contig`, `--as-table` (commits c5cd2524, dae62624, 77c3f9ef; 23 driver tests; the full `cargo test --release`, 1,193 tests, and the 40 bench unit tests pass). Runner and scorers: `bench/seed_pool/{run.sh, composition.py (12 tests), gates.py, table.py}`. Products: `/mnt/linuxdisk/tmp/seed_pool_2026-10-07/{human_chr16, human_chr17, gorilla_npip, gorilla_tbc1d3}/` (`table.txt`, `table.md`, `table.tsv`, `matrix.tsv`, `gates.json`, one directory per arm with the driver's own products and `run.log` stamps).** Human and gorilla are never pooled. User direction (2026-10-07): 'lets also make it easy to set the good secondaries and all in order to prove to my advisor how does this decision affect the outcome'.

## Answer

1. **The pool is one axis, and every count that matters moves along it.** On NPIP (chr16, 25 copies, 22 expressed) from primary-only (P) through ties (G100), the default band (G98) to every secondary (A): the number of copies whose locus representative is exact falls 12, 10, 8 (G995 and G98), 6, 3, 2 of 22; node precision (nodes of the family's clusters that lie on a truth copy) is .72 at P, .59 at the default and .41 at A; the copies that sit in the ONE largest family cluster fall 23 to 17 of 25 (A splits the family over 15 clusters); the aligner's work rises 23 times (72,118 to 1,648,588 PAF records; about 30 minutes of bounded, sharded all-vs-all against a 36-second families stage). What widening buys: one more copy with a node (24 to 25, PKD1P6-NPIPP1, at rho <= 0.95), and family recall against the human family truths (Soto sensitivity .423 to .732, Soto F .588 to .832; Compara F .667 to .674 and U2 F .644 to .638 do not move).
2. **TBC1D3 (chr17, 16 copies, 11 expressed) is a different shape.** The default beats primary-only on every count except copies with a node, which ties at 13 (exact representatives 7 to 8, node precision .81 to 1.00, copies in the largest cluster 9 to 11, Compara F .396 to .475, Soto F .434 to .487) and `all` collapses (83 nodes, 68 of them on no truth copy, node precision .18). So the default is a good point for TBC1D3. For NPIP it costs four exact representatives (12 to 8) and .13 of node precision (.72 to .59) against primary-only, and gains .07 of Soto F (.588 to .655).
3. **0.98 is a point on a continuum, not a knee.** On NPIP the representative count changes by about two copies per step of the width (S8 failed: G95 has 6 exact representatives against 8). The first loss already happens at ties (G100: NPIPB5 and LOC124907807 lose theirs), i.e. exactly-tied secondaries, the coin-toss multimappers, are enough. On TBC1D3 the band 0.90 to 0.995 stays within one copy of the default (S8 held).
4. **A primary-first representative (Rule 1) repairs the representative without touching membership, at a small cost in family F.** NPIP: G98 8 to 11 exact (NPIPA8, NPIPB15 and LOC124907807 return; NPIPA6 and NPIPB5 do not), A 2 to 9, own-node and largest-cluster counts unchanged, P+R1 identical to P. The U2 and Compara F of G98+R1 move by .022 and .017 and Soto F by .013, so the registered clause |dF| < .01 fails (S4 failed). TBC1D3: inert at the default; at A it lifts exact representatives 5 to 8 and copies with a node 12 to 15 (S5 held).
5. **The non-dominated arms on (copies with a node, exact representatives, node precision)** are P, A and A+R1 among the five primary NPIP arms (G98 and G98+R1 are dominated by P), and G98, G98+R1 and A+R1 among the TBC1D3 ones. Over the whole sweep G95+R1 (25 / 11 / .73) sits beside P (24 / 12 / .72) on NPIP: a descriptive finding of a sweep, not a registered claim.
6. **The family-level scores and node precision understate what `all` adds:** `family_score` deletes predicted members that are outside the truth universe before scoring, so the 46 nodes of A on no truth copy at chr16 (68 of 83 at chr17) are invisible to it. Its precision for A falls less than node precision does (chr16: Compara 1.000 to .857, U2 .760 to .629, Soto .968 to .963, against node precision .72 to .41). Node precision counts those nodes as false; some may be real unannotated paralogs (NPIP-core sequence outside annotated genes), so it is a floor on precision, not an error rate. It also counts every locus on a truth copy as correct, so loci piled on one copy raise it (A has 2.8 on-copy nodes per copy with a node at chr16, 2.3 at gorilla NPIP); the copy-level precision CP* (copies with a node over nodes) is .60, .44 and .15 at chr16 P, G98 and A.

7. **Gorilla NPIP reverses the sign (held-out for Rule 1, partly spent for the pool).** Outside the `all` pool no locus carries an annotated intron of any of the 9 expressed copies; the `all` pool gives 2-4 of them at 7 of the 9 and a node to 22 of 24 copies against 11 (the registered S1 and S3 fail in the direction of `all`). The expressed gorilla NPIP copies have a median of 2 primary alignments carrying their introns against 156 in human NPIP, so the primary pool has nothing to build from. Gorilla TBC1D3 (E = 7, underpowered) behaves like human: widening adds copies with a node (4, 7, 8 of 14) and the exact-representative count is flat. Rule 1 changes nothing on gorilla (floor). The `all` arms were deferred at the user's request (Amendment 2); the finished gorilla NPIP `all` arm is the best arm there, so the deferred runs (A on TBC1D3, A+R1 on NPIP) are the informative ones.

## NPIP, human chr16 (N = 25 copies, E = 22 expressed under Amendment E; arms in order of pool width)

| arm | transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M3 | M6 | compara sens / prec / F | u2 sens / prec / F | soto sens / prec / F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 8,673 | 2,579 | 72,118 | 40 (29 / 1 / 2 / 8) | 0.72 | 1.2 | 0.60 | 24 | 12 | 21 | 23 | 0.500 / 1.000 / 0.667 | 0.559 / 0.760 / 0.644 | 0.423 / 0.968 / 0.588 |
| G100 | 8,914 | 2,648 | 76,572 | 44 (31 / 2 / 1 / 10) | 0.70 | 1.3 | 0.55 | 24 | 10 | 21 | 23 | 0.481 / 1.000 / 0.650 | 0.559 / 0.731 / 0.633 | 0.451 / 0.970 / 0.615 |
| G995 | 9,108 | 2,705 | 81,388 | 55 (34 / 3 / 5 / 13) | 0.62 | 1.4 | 0.44 | 24 | 8 | 21 | 23 | 0.481 / 1.000 / 0.650 | 0.559 / 0.704 / 0.623 | 0.479 / 0.944 / 0.636 |
| G98 | 9,473 | 2,832 | 94,433 | 54 (32 / 3 / 8 / 11) | 0.59 | 1.3 | 0.44 | 24 | 8 | 21 | 23 | 0.500 / 1.000 / 0.667 | 0.588 / 0.714 / 0.645 | 0.507 / 0.923 / 0.655 |
| G95 | 9,767 | 2,953 | 112,419 | 48 (34 / 2 / 3 / 9) | 0.71 | 1.4 | 0.52 | 25 | 6 | 21 | 24 | 0.500 / 1.000 / 0.667 | 0.559 / 0.704 / 0.623 | 0.535 / 0.974 / 0.691 |
| G90 | 10,216 | 3,089 | 127,721 | 59 (39 / 3 / 9 / 8) | 0.66 | 1.6 | 0.42 | 25 | 3 | 21 | 24 | 0.481 / 1.000 / 0.650 | 0.618 / 0.656 / 0.636 | 0.634 / 0.957 / 0.763 |
| A | 14,183 | 5,951 | 1,648,588 | 168 (69 / 21 / 32 / 46) | 0.41 | 2.8 | 0.15 | 25 | 2 | 19 | 17 | 0.556 / 0.857 / 0.674 | 0.647 / 0.629 / 0.638 | 0.732 / 0.963 / 0.832 |
| P+R1 | 8,673 | 2,579 | 72,118 | 40 (29 / 1 / 2 / 8) | 0.72 | 1.2 | 0.60 | 24 | 12 | 21 | 23 | 0.500 / 1.000 / 0.667 | 0.559 / 0.760 / 0.644 | 0.423 / 0.968 / 0.588 |
| G100+R1 | 8,914 | 2,648 | 76,572 | 46 (31 / 2 / 1 / 12) | 0.67 | 1.3 | 0.52 | 24 | 11 | 21 | 23 | 0.481 / 1.000 / 0.650 | 0.559 / 0.731 / 0.633 | 0.451 / 0.970 / 0.615 |
| G995+R1 | 9,108 | 2,705 | 81,388 | 59 (34 / 3 / 5 / 17) | 0.58 | 1.4 | 0.41 | 24 | 9 | 21 | 23 | 0.481 / 1.000 / 0.650 | 0.559 / 0.704 / 0.623 | 0.479 / 0.944 / 0.636 |
| G98+R1 | 9,473 | 2,832 | 94,433 | 54 (33 / 3 / 8 / 10) | 0.61 | 1.4 | 0.44 | 24 | 11 | 21 | 23 | 0.481 / 1.000 / 0.650 | 0.559 / 0.704 / 0.623 | 0.493 / 0.921 / 0.642 |
| G95+R1 | 9,767 | 2,953 | 112,419 | 45 (33 / 2 / 3 / 7) | 0.73 | 1.3 | 0.56 | 25 | 11 | 21 | 24 | 0.500 / 1.000 / 0.667 | 0.559 / 0.704 / 0.623 | 0.535 / 0.974 / 0.691 |
| G90+R1 | 10,216 | 3,089 | 127,721 | 60 (37 / 4 / 10 / 9) | 0.62 | 1.5 | 0.40 | 24 | 10 | 21 | 22 | 0.500 / 1.000 / 0.667 | 0.618 / 0.656 / 0.636 | 0.606 / 0.956 / 0.741 |
| A+R1 | 14,183 | 5,951 | 1,648,588 | 176 (72 / 21 / 36 / 47) | 0.41 | 2.9 | 0.14 | 25 | 9 | 20 | 17 | 0.574 / 0.838 / 0.681 | 0.647 / 0.595 / 0.620 | 0.732 / 0.945 / 0.825 |

M2 is the strict Amendment E criterion (the representative starts at a capped 5' end of the copy's reads and carries its first three introns, within an own node); M3 asks the same of any transcript of the locus. Copies that change between arms (from `matrix.tsv`): P to G98 loses the exact representative at NPIPB5 and LOC124907807 (already at ties), NPIPA6 and NPIPA8 (at 0.995) and NPIPB15 (at 0.98), and gains NPIPB10P; Rule 1 at G98 returns NPIPA8, NPIPB15 and LOC124907807.

## TBC1D3, human chr17 (N = 16 copies, E = 11)

| arm | transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M3 | M6 | compara sens / prec / F | soto sens / prec / F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 10,439 | 2,794 | 70,710 | 16 (13 / 0 / 1 / 2) | 0.81 | 1.0 | 0.81 | 13 | 7 | 9 | 9 | 0.253 / 0.905 / 0.396 | 0.291 / 0.852 / 0.434 |
| G100 | 10,502 | 2,837 | 72,775 | 16 (13 / 0 / 3 / 0) | 0.81 | 1.0 | 0.81 | 13 | 7 | 9 | 9 | 0.267 / 0.909 / 0.412 | 0.304 / 0.857 / 0.449 |
| G995 | 10,565 | 2,853 | 73,616 | 19 (19 / 0 / 0 / 0) | 1.00 | 1.5 | 0.68 | 13 | 7 | 9 | 11 | 0.280 / 0.913 / 0.429 | 0.329 / 0.839 / 0.473 |
| G98 | 10,668 | 2,880 | 75,896 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 8 | 9 | 11 | 0.320 / 0.923 / 0.475 | 0.367 / 0.725 / 0.487 |
| G95 | 10,791 | 2,935 | 81,664 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 7 | 9 | 11 | 0.320 / 0.889 / 0.471 | 0.392 / 0.721 / 0.508 |
| G90 | 10,923 | 3,035 | 90,308 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 | 0.418 / 0.733 / 0.532 |
| A | 14,472 | 5,726 | 1,061,646 | 83 (15 / 0 / 0 / 68) | 0.18 | 1.2 | 0.14 | 12 | 5 | 9 | 10 | 0.373 / 0.933 / 0.533 | 0.519 / 0.707 / 0.599 |
| P+R1 | 10,439 | 2,794 | 70,710 | 16 (13 / 0 / 1 / 2) | 0.81 | 1.0 | 0.81 | 13 | 7 | 9 | 9 | 0.253 / 0.905 / 0.396 | 0.291 / 0.852 / 0.434 |
| G100+R1 | 10,502 | 2,837 | 72,775 | 16 (12 / 0 / 2 / 2) | 0.75 | 1.0 | 0.75 | 12 | 6 | 9 | 8 | 0.280 / 0.913 / 0.429 | 0.316 / 0.862 / 0.463 |
| G995+R1 | 10,565 | 2,853 | 73,616 | 13 (13 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 13 | 7 | 9 | 11 | 0.307 / 0.920 / 0.460 | 0.354 / 0.875 / 0.505 |
| G98+R1 | 10,668 | 2,880 | 75,896 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 8 | 9 | 11 | 0.320 / 0.923 / 0.475 | 0.367 / 0.725 / 0.487 |
| G95+R1 | 10,791 | 2,935 | 81,664 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 | 0.392 / 0.721 / 0.508 |
| G90+R1 | 10,923 | 3,035 | 90,308 | 22 (22 / 0 / 0 / 0) | 1.00 | 1.7 | 0.59 | 13 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 | 0.418 / 0.733 / 0.532 |
| A+R1 | 14,472 | 5,726 | 1,061,646 | 94 (27 / 0 / 0 / 67) | 0.29 | 1.8 | 0.16 | 15 | 8 | 9 | 11 | 0.373 / 0.933 / 0.533 | 0.519 / 0.804 / 0.631 |

(The other Rule-1 sweep arms are in `human_chr17/table.md`.) The `TBC1D3` truth has no U2 family, so U2 is not scored.

## The registered predictions (reported as the runner prints them; no combined verdict)

| | statement | chr16 NPIP | chr17 TBC1D3 |
|---|---|---|---|
| S1 | exact representatives P >= G98 >= A, one strict | HELD (12 >= 8 >= 2) | FAILED (7, 8, 5) |
| S2 | copies with a node G98 >= P | HELD (24, 24) | HELD (13, 13) |
| S3 | node precision P >= G98 > A | HELD (.72, .59, .41) | FAILED (.81, 1.00, .18) |
| S4 | G98+R1 >= G98 + 1 exact, no fewer nodes, \|dF\| < .01 | FAILED (11 against 8, but \|dF\| .022) | FAILED (8 against 8) |
| S5 | A+R1 >= A + 2 exact, no fewer nodes | HELD (9 against 2) | HELD (8 against 5; nodes 15 against 12) |
| S6 | P+R1 identical to P | HELD | HELD |
| S8 | G95 and G995 within one copy of G98 on M1 and M2 | FAILED (G95 has 6 against 8) | HELD |

Gates: G0 (the default through the new driver, `--contig chr16 --as-table`, is byte-identical to the 10-06 default re-score: assembled GTF, families GTF, clusters, loci GFF3), G1 (own node 24, exact 8, locus 21, U2 F .645), G1b (`composition.py` equals `nodes.py` on P and G98), G2 (P+R1 = P), G3 (`copy_support` identical under `PYTHONHASHSEED` 0 and 1) pass on chr16; G2, G3 and G4 (the sharded all-vs-all equals a single process; 6 shards on the P arm of chr17) pass on chr17, where G1b does not apply (`nodes.py` stops on three loci that share a span; `composition.py` gives such loci the cluster ids of their span, and reports 5 shared spans in A at chr16 and 13 at chr17). The headline counts were recomputed from the raw per-copy tables with separate code and agree.

## Reading

- The pool decision is a choice of where to stand on a monotone trade-off, not a free improvement: wider pools connect more paralogs (family recall rises on the Soto truth) and cost exact representatives, node precision and cohesion; at `all` the family is cut into 15 clusters on NPIP although every copy has a node. This is the same connectivity-for-precision trade the 09-22 read-pool study found on chr20 and gorilla, now measured per copy.
- The default is not wrong on either family, but it is not best on both: at NPIP primary-only is ahead on exact representatives and node precision and ties on membership and cohesion (the default is ahead only on Soto F), while at TBC1D3 the default is ahead or tied on every count. A rule that picks the pool per locus (rather than genome-wide) is not tested here.
- Rule 1 is a cheap, local repair of the representative; it is not free at the family level (up to .02 of F on the NPIP truths), and at TBC1D3 it only matters when the pool is `all`. Its native implementation (a `primary_reads` attribute, `mcl_families --representative primary-first`) is not justified by one dev library.

## Limits

One library per species; the human block is the development block of both families (nothing here is held out for them); E-found is the strict, cap-based criterion, so absolute exact-representative counts are conservative; counts are of copies (25, 16), a difference of one is one copy and the sweep (G100 to G90, their Rule-1 arms) is descriptive, so single-step differences are not claims; node precision treats unannotated nodes as false; the family-level scorer's universe intersection hides the cost of `all` (item 6); the primaries of Rule 1 come from the same reads and polish; M6 was added after the human arms (post hoc for them, Amendment 1); the redundancy and CP* columns and `secondary_support.py` were added after seeing the gorilla arms (post hoc); the aligner intermediates were deleted after each arm (the PAF record counts are kept).

## 5. Gorilla OR6737 (Amendment 1; held out for Rule 1; the all-secondaries arms were deferred by Amendment 2)

Reads OR6737 against the KB3781 assembly and RefSeq annotation (cross-individual); copy truth = the annotation-only sets of the 09-29 copy-recovery study; the criterion is Amendment A (the gorilla reads carry no cap signal), E from reads and annotation only. NPIP: contigs NC_073241.2 + NC_073242.2, **N = 24 copies, E = 9** (the 25th copy is on NC_073244.2, not assembled). TBC1D3: NC_073228.2 + NC_073224.2, **N = 14, E = 7, so UNDERPOWERED** for every prediction that reads the exact-representative counts. No family truth exists on these contigs, so M4 is not scored. Arms run: P, P+R1, G98, G98+R1, G995, G95 on both; A on NPIP only (it had finished when the all-secondaries arms were deferred; its Rule-1 arm was stopped and deleted); nothing for A on TBC1D3.

### Gorilla NPIP (N = 24, E = 9)

| arm | transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M3 | M6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 4,347 | 1,170 | 56,165 | 17 (9 / 0 / 3 / 5) | 0.53 | 1.0 | 0.53 | 9 | 0 | 0 | 8 |
| G995 | 4,449 | 1,221 | 69,802 | 20 (11 / 0 / 3 / 6) | 0.55 | 1.0 | 0.55 | 11 | 0 | 0 | 10 |
| G98 | 4,522 | 1,252 | 71,759 | 22 (11 / 0 / 3 / 8) | 0.50 | 1.0 | 0.50 | 11 | 0 | 0 | 10 |
| G95 | 4,705 | 1,336 | 75,679 | 22 (11 / 0 / 3 / 8) | 0.50 | 1.0 | 0.50 | 11 | 0 | 0 | 10 |
| A | 6,863 | 2,448 | 219,687 | 57 (51 / 0 / 0 / 6) | 0.89 | 2.3 | 0.39 | 22 | 7 | 7 | 13 |
| P+R1 | 4,347 | 1,170 | 56,165 | 17 (9 / 0 / 3 / 5) | 0.53 | 1.0 | 0.53 | 9 | 0 | 0 | 8 |
| G98+R1 | 4,522 | 1,252 | 71,759 | 22 (11 / 0 / 3 / 8) | 0.50 | 1.0 | 0.50 | 11 | 0 | 0 | 10 |

### Gorilla TBC1D3 (N = 14, E = 7)

| arm | transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M3 | M6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 12,233 | 3,127 | 238,328 | 4 (4 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 4 | 4 | 4 | 4 |
| G995 | 12,343 | 3,182 | 239,833 | 6 (6 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 6 | 4 | 4 | 6 |
| G98 | 12,462 | 3,197 | 241,374 | 7 (7 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 7 | 4 | 4 | 7 |
| G95 | 12,747 | 3,232 | 240,590 | 8 (8 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 8 | 5 | 5 | 8 |
| P+R1 | 12,233 | 3,127 | 238,328 | 4 (4 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 4 | 4 | 4 | 4 |
| G98+R1 | 12,462 | 3,197 | 241,374 | 7 (7 / 0 / 0 / 0) | 1.00 | 1.0 | 1.00 | 7 | 4 | 4 | 7 |

(Rule 1 changes neither family's clusters at P or at G98: P+R1 and G98+R1 equal P and G98 in every column.)

| | statement | gorilla NPIP | gorilla TBC1D3 |
|---|---|---|---|
| S1 | exact representatives P >= G98 >= A, one strict | FAILED (0, 0, 7: reversed) | NO VERDICT (E = 7, no A) |
| S2 | copies with a node G98 >= P | HELD (11 against 9) | HELD (7 against 4) |
| S3 | node precision P >= G98 > A | FAILED as registered (.53, .50, .89); holds with copy-level precision CP* (.53, .50, .39), see below | NO VERDICT (no A) |
| S4g | G98+R1 >= G98 + 1 exact, no fewer nodes, no smaller largest cluster | FAILED (0 against 0: nothing to repair) | NO VERDICT (E = 7) |
| S5 | A+R1 >= A + 2 exact | NO VERDICT (A+R1 deferred) | NO VERDICT |
| S6 | P+R1 identical to P | HELD | HELD |
| S8 | G95 and G995 within one copy of G98 | HELD (no difference) | NO VERDICT (E = 7; counts within one) |

Gates: G1b (composition equals `nodes.py` on P and G98), G2 and G3 pass on NPIP; G2 and G3 pass on TBC1D3 (`nodes.py` stops on shared spans there); G0, G1 and G4 are human-only.

**What gorilla NPIP shows.** The sign of the pool effect reverses. Under P, G995, G98 and G95 no locus carries any annotated intron of any of the 9 expressed copies (0, by the representative and by any transcript of the locus), while the `all` pool gives 2-4 annotated introns at 7 of the 9 (NPIPA5 and NPIPB14P stay at 0) and a node to 22 of 24 copies against 11. Each expressed copy has 2-9 primary alignments carrying its annotated introns (median 2) against 18-36 secondary ones (median 29; `bench/seed_pool/secondary_support.py`), so the primary pool has almost nothing to build the chain from. The same count in human is a median of 156 primary against 3,427 secondary alignments per expressed NPIP copy (TBC1D3: 150 and 2,037; gorilla TBC1D3: 6 and 77). The secondary-to-primary ratio is 13 to 22 on all four, so what differs is the absolute primary evidence, and with it whether the primary pool suffices (human NPIP: P finds 12 of 22 expressed copies exactly; gorilla TBC1D3, median 6 primaries: 4 of 7; gorilla NPIP, median 2: 0 of 9). This locates the earlier 'lost before the assembler' reading of gorilla NPIP (`docs/PREREG_copy_recovery_tools_2026-09-29.md`, P5): the loss is in the pool. Caveat: at these copies most of the alignments the `all` pool adds are secondary, i.e. reads whose best placement is elsewhere, so a chain 'found' at the copy may be paralog evidence placed there; the pool cannot tell the two apart (O2's question), and E requires only two primary alignments at the copy.
**NP is inflated by redundancy.** The registered node precision counts every locus on a truth copy as correct, so several loci piled on one copy raise it: A has 51 on-copy nodes for 22 copies (2.3 per copy) and an NP of .89 against .50 at the default. Copy-level precision (copies with a node over nodes, reported beside as CP*) is .39 for A against .50 at G98, and on the human contigs CP* is .15 (chr16) and .14 (chr17) for A. The registered S3 verdict on gorilla NPIP stands as printed; CP* was added after seeing it (post hoc).
**TBC1D3.** Widening the pool adds copies with a node (4, 6, 7, 8 of 14 at P, G995, G98, G95) and one more exact representative at G95 (5), with node precision 1.00 throughout and all nodes in one family cluster; the counts are too small for a verdict.

