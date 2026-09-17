# NPIP read-group nodes vs shipped de novo nodes (human CHM13, 2026-09-17)

DEVELOPMENT family (NPIP, 27 copies), substrate = the ladder's copy-centred windows. Not validation. No src change, nothing committed.
Declarations: `/mnt/linuxdisk/home/juanfraitu/sd_readgroup/DECLARATIONS.txt`, written 05:57:14 (md5 `f160935a`), before any node list (arms.pkl 05:59), mapping (06:03) or metric (results.json 06:04). **Addendum 1** (06:04:56, after the first scoring run, report-only, decision unchanged) adds "copies sharing an M1 node" and M1 node span; the part above it is unchanged (md5 `f160935a` recomputed).

## Read first

FAMILY level, triangle leaders. F = bipartite F with NPIP-only precision. Composition = distinct nodes in the matched family: M1 copy nodes + NPIP-exon fragments + other (size_strict in brackets where copies share nodes). "shared" = copies whose M1 node is also another copy's M1 node (Addendum 1).

| arm | nodes | fams | largest | R | P | F | P strict | **F strict** | comp F / F strict | pair sens / prec strict | NPIP family | full-length | opp-strand | shared |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 target | 167 | 23 | 40 | 1.000 | 1.000 | 1.000 | 0.675 | **0.806** | 1.000 / 0.540 | 1.000 / 0.450 | 40 = 27+7+6 | 27/27 | 0 | 0 |
| G0 shipped | 560 | 68 | 93 | 1.000 | 1.000 | 1.000 | 0.290 | **0.450** | 1.000 / 0.196 | 1.000 / 0.082 | 93 = 27+23+43 | 5/27 | 4 | 0 |
| G3 strand-fix | 483 | 60 | 84 | 1.000 | 1.000 | 1.000 | 0.321 | **0.486** | 1.000 / 0.234 | 1.000 / 0.101 | 84 = 27+18+39 | 5/27 | 1 | 0 |
| G1a read-group | 259 | 22 | 56 | 1.000 | 1.000 | 1.000 | 0.409 | **0.581** | 1.000 / 0.302 | 1.000 / 0.164 | 56 = 17+12+27 [66] | 26/27 | 0 | 19 |
| G1b + spliced strand | 191 | 16 | 29 | 1.000 | 1.000 | 1.000 | 0.692 | **0.818** | 1.000 / 0.450 | 1.000 / 0.474 | 29 = 17+2+10 [39] | 25/27 | 0 | 19 |
| G2a G1a + split | 379 | 34 | 75 | 1.000 | 1.000 | 1.000 | 0.318 | **0.482** | 1.000 / 0.206 | 1.000 / 0.098 | 75 = 17+15+43 [85] | 26/27 | 0 | 19 |
| G2b G1b + split | 297 | 25 | 45 | 1.000 | 1.000 | 1.000 | 0.491 | **0.659** | 1.000 / 0.277 | 1.000 / 0.236 | 45 = 17+3+25 [55] | 25/27 | 0 | 19 |

**Decision-rule outcome: NO CANDIDATE — no arm is eligible (G3 loses 11 G0 families; G1a/G1b/G2a/G2b each lose 59-66 G0 families, 5 G0 families lose every exonic base, and the largest G0 family's best-match Jaccard is 0.13-0.17), and G3, the arm closest to passing, beats G0 by only 0.036 F strict (< 0.05) anyway.**

Beside the rule: the A1 target row itself would fail the same footprint gates (65 G0 families "lost", 25 lose every base, largest best J 0.185), so the Jaccard >= 0.5 / >= 0.8 gates reject any node set whose footprints differ from G0's, not only harmful ones.

## SUBFAMILY levels (triangle; no arm splits NPIP — matched group is always the NPIP family)

| level | arm | R | P | F | P strict | F strict | pair sens | pair prec strict | comp F / F strict |
|---|---|---|---|---|---|---|---|---|---|
| SUB-1 (NPIPB→fam, J 0.704) | A1 | 0.704 | 0.704 | 0.704 | 0.475 | 0.567 | 1.000 | 0.255 | 0.704 / 0.380 |
| | G0 | 0.704 | 0.704 | 0.704 | 0.204 | 0.317 | 1.000 | 0.047 | 0.704 / 0.138 |
| | G3 | 0.704 | 0.704 | 0.704 | 0.226 | 0.342 | 1.000 | 0.057 | 0.704 / 0.165 |
| | G1a | 0.704 | 0.704 | 0.704 | 0.288 | 0.409 | 1.000 | 0.093 | 0.704 / 0.212 |
| | G1b | 0.704 | 0.704 | 0.704 | 0.487 | 0.576 | 1.000 | 0.269 | 0.704 / 0.317 |
| | G2a | 0.704 | 0.704 | 0.704 | 0.224 | 0.339 | 1.000 | 0.056 | 0.704 / 0.145 |
| | G2b | 0.704 | 0.704 | 0.704 | 0.345 | 0.463 | 1.000 | 0.134 | 0.704 / 0.195 |
| SUB-2 (A6-9→fam, J 0.148) | A1 | 0.148 | 0.148 | 0.148 | 0.100 | 0.119 | 1.000 | 0.031 | 0.148 / 0.080 |
| | G0 | 0.148 | 0.148 | 0.148 | 0.043 | 0.067 | 1.000 | 0.006 | 0.148 / 0.029 |
| | G3 | 0.148 | 0.148 | 0.148 | 0.048 | 0.072 | 1.000 | 0.007 | 0.148 / 0.035 |
| | G1a | 0.148 | 0.148 | 0.148 | 0.061 | 0.086 | 1.000 | 0.011 | 0.148 / 0.045 |
| | G1b | 0.148 | 0.148 | 0.148 | 0.103 | 0.121 | 1.000 | 0.032 | 0.148 / 0.067 |
| | G2a | 0.148 | 0.148 | 0.148 | 0.047 | 0.071 | 1.000 | 0.007 | 0.148 / 0.031 |
| | G2b | 0.148 | 0.148 | 0.148 | 0.073 | 0.098 | 1.000 | 0.016 | 0.148 / 0.041 |

## Blast radius vs G0 (whole substrate, triangle families; footprint = merged member exons, strand-blind)

| arm | families | G0 lost (best J < 0.5) | G0 losing every exonic base | largest G0 family best J | G0 changed (J < 1) | new (best J < 0.5; 0-bp) | pairs (exon / body edges) | eligible |
|---|---|---|---|---|---|---|---|---|
| A1 (report-only) | 23 | 65 | 25 | 0.185 | 68 | 20 (1) | 762 (461 / 711) | – |
| G0 | 68 | 0 | 0 | 1.000 | 0 | 0 (0) | 1919 (1349 / 1885) | – |
| G3 | 60 | 11 | 0 | 0.965 | 27 | 3 (0) | 1622 (1085 / 1583) | no: 11 lost |
| G1a | 22 | 64 | 5 | 0.134 | 68 | 18 (0) | 811 (607 / 730) | no |
| G1b | 16 | 66 | 5 | 0.142 | 68 | 14 (0) | 377 (261 / 313) | no |
| G2a | 34 | 59 | 5 | 0.170 | 68 | 25 (2) | 1441 (1076 / 1236) | no |
| G2b | 25 | 62 | 5 | 0.165 | 68 | 19 (3) | 708 (466 / 575) | no |

G3's 11 lost G0 families have 12, 6, 5, 4, 3, 3, 3, 3, 2, 2, 2 members (G0 ids 1, 12, 19, 26, 36, 37, 39, 41, 48, 49, 53). The 5 G0 families with no exonic base covered in every read-group arm are G0 ids 13, 45, 50, 51, 58.

## Per copy: G0 and G1b (highest F strict non-G0 arm; INELIGIBLE — no arm is eligible)

covered = M1 node exon bp on copy exons / copy exon bp · #ov = nodes overlapping copy exons · in = in matched family · reads = same-strand primary MAPQ>=1 reads on copy exons (ladder).

| copy | str | reads | G0 node | G0 covered | #ov | in | G1b node | span kb | G1b covered | #ov | in | G1b node shared with |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NPIPB2 | - | 422 | chr16:11963703-11968998 + 2ex | 416/1702 (0.24) | 2 | yes | chr16:11786638-12015041 - 45ex | 228.4 | 1598/1702 (0.94) | 2 | yes | - |
| NPIPA2 | + | 350 | chr16:14749588-14763868 + 7ex | 2827/4517 (0.63) | 3 | yes | chr16:14740775-15167445 + 85ex | 426.7 | 4472/4517 (0.99) | 2 | yes | NPIPA1 |
| NPIPA1 | + | 741 | chr16:14938814-14953123 + 7ex | 1021/1084 (0.94) | 1 | yes | chr16:14740775-15167445 + 85ex | 426.7 | 1084/1084 (1.00) | 2 | yes | NPIPA2 |
| PKD1P6-NPIPP1 | - | 456 | chr16:15115422-15132313 - 10ex | 1898/5391 (0.35) | 5 | yes | chr16:15029880-15262008 - 58ex | 232.1 | 5379/5391 (1.00) | 3 | yes | - |
| NPIPA5 | - | 146 | chr16:15368407-15382704 - 7ex | 2915/3741 (0.78) | 2 | yes | chr16:15368406-15469419 - 31ex | 101.0 | 3653/3741 (0.98) | 2 | yes | - |
| NPIPA6 | + | 127 | chr16:16344720-16359014 + 8ex | 1088/1570 (0.69) | 2 | yes | chr16:16328099-16438449 + 50ex | 110.3 | 1569/1570 (1.00) | 1 | yes | NPIPA7 |
| NPIPA7 | + | 81 | chr16:16391852-16406187 + 7ex | 1093/1264 (0.86) | 1 | yes | chr16:16328099-16438449 + 50ex | 110.3 | 1233/1264 (0.98) | 1 | yes | NPIPA6 |
| NPIPA8 | - | 38 | chr16:18325177-18339582 - 7ex | 1080/1588 (0.68) | 2 | yes | chr16:18147539-18493395 - 68ex | 345.9 | 1300/1588 (0.82) | 1 | yes | NPIPA9 |
| NPIPA9 | - | 849 | chr16:18372354-18376104 - 2ex | 704/2869 (0.25) | 1 | yes | chr16:18147539-18493395 - 68ex | 345.9 | 2869/2869 (1.00) | 1 | yes | NPIPA8 |
| NPIPB3 | - | 1298 | chr16:21337422-21344055 - 3ex | 2989/3674 (0.81) | 2 | yes | chr16:21314626-21842398 - 134ex | 527.8 | 3668/3674 (1.00) | 1 | yes | LOC128966608 |
| LOC128966608 | - | 485 | chr16:21692009-21693757 + 1ex | 129/3557 (0.04) | 2 | yes | chr16:21314626-21842398 - 134ex | 527.8 | 3550/3557 (1.00) | 1 | yes | NPIPB3 |
| NPIPB4 | + | 648 | chr16:22382178-22387825 + 2ex | 184/3634 (0.05) | 4 | yes | chr16:22216076-22958833 + 126ex | 742.8 | 3634/3634 (1.00) | 2 | yes | NPIPB5 |
| NPIPB5 | + | 280 | chr16:22788631-22792526 + 1ex | 3031/9237 (0.33) | 4 | yes | chr16:22216076-22958833 + 126ex | 742.8 | 9229/9237 (1.00) | 3 | yes | NPIPB4 |
| NPIPB6 | - | 659 | chr16:28625753-28628690 + 1ex | 2635/7025 (0.38) | 5 | yes | chr16:28623009-28772328 - 43ex | 149.3 | 6325/7025 (0.90) | 2 | yes | NPIPB7 |
| NPIPB7 | - | 124 | chr16:28736995-28771999 - 8ex | 1274/2952 (0.43) | 2 | yes | chr16:28623009-28772328 - 43ex | 149.3 | 2866/2952 (0.97) | 2 | yes | NPIPB6 |
| NPIPB8 | + | 44 | chr16:28935124-28939424 + 3ex | 781/1293 (0.60) | 1 | yes | chr16:28772052-29800624 + 210ex | 1028.6 | 1266/1293 (0.98) | 1 | yes | NPIPB9, NPIPB10P |
| NPIPB9 | + | 614 | chr16:29048645-29050801 + 1ex | 520/4193 (0.12) | 3 | yes | chr16:28772052-29800624 + 210ex | 1028.6 | 4105/4193 (0.98) | 1 | yes | NPIPB8, NPIPB10P |
| NPIPB10P | + | 66 | chr16:29328121-29333923 + 1ex | 441/880 (0.50) | 1 | yes | chr16:28772052-29800624 + 210ex | 1028.6 | 817/880 (0.93) | 1 | yes | NPIPB8, NPIPB9 |
| NPIPB11 | - | 155 | chr16:29663364-29679826 - 7ex | 3236/6117 (0.53) | 2 | yes | chr16:29663336-29949744 - 85ex | 286.4 | 4550/6117 (0.74) | 2 | yes | NPIPB12 |
| NPIPB12 | - | 27 | chr16:29765341-29767376 + 1ex | 1783/3737 (0.48) | 3 | yes | chr16:29663336-29949744 - 85ex | 286.4 | 3046/3737 (0.82) | 2 | yes | NPIPB11 |
| LOC124907834 | - | 487 | chr16:30507437-30523997 - 7ex | 3120/3244 (0.96) | 1 | yes | chr16:30437930-30741477 - 111ex | 303.5 | 3244/3244 (1.00) | 1 | yes | NPIPB13 |
| NPIPB13 | - | 151 | chr16:30609406-30625966 - 7ex | 3108/6144 (0.51) | 2 | yes | chr16:30437930-30741477 - 111ex | 303.5 | 6144/6144 (1.00) | 2 | yes | LOC124907834 |
| NPIPB14P | - | 1242 | chr16:75785725-75790348 - 2ex | 1574/2418 (0.65) | 2 | yes | chr16:75785693-76018542 - 25ex | 232.8 | 2418/2418 (1.00) | 1 | yes | - |
| NPIPB15 | + | 168 | chr16:80195301-80209898 + 7ex | 1605/4340 (0.37) | 1 | yes | chr16:80153261-80313480 + 15ex | 160.2 | 3581/4340 (0.83) | 3 | yes | - |
| LOC124907808 | + | 22 | chr16:80319739-80324246 + 2ex | 2777/4981 (0.56) | 2 | yes | chr16:80307233-80324271 + 6ex | 17.0 | 4839/4981 (0.97) | 1 | yes | - |
| LOC124907807 | + | 77 | chr16:80424045-80438604 + 7ex | 1805/4571 (0.39) | 1 | yes | chr16:80414859-80438605 + 10ex | 23.7 | 2471/4571 (0.54) | 3 | yes | - |
| NPIPB1P | - | 407 | chr18:11781592-11796212 - 7ex | 810/981 (0.83) | 1 | yes | chr18:11730340-11832558 - 22ex | 102.2 | 981/981 (1.00) | 1 | yes | - |

G3 vs G0 per copy: identical M1 nodes except LOC128966608, NPIPB6, NPIPB12 (same single exon, strand + → −) and NPIPB7 / LOC124907834 / NPIPB13 (same 8/7/7-exon node, extended by merged re-stranded single-exon pieces: two for NPIPB7, one each for the others). NPIPB2 stays opposite (its M1 node is a read-locus node, strand from reads). Full-length stays 5/27.

## What the numbers show

- **G3 (strand fix) is a small, clean-looking gain that fails the family-loss gate.** 215/371 single-exon reps flip to '−' (144 stay '+' by reads, 12 have no reads). Nodes 560 → 483, NPIP family 93 → 84, F strict 0.450 → 0.486 (+0.036), opposite-strand copies 4 → 1, full-length unchanged at 5/27. 11 small G0 families (2-12 members) fall below J 0.5; none loses its bases.
- **The read-group arms' F strict and full-length gains come from nodes that merge neighbouring copies.** In all four, the 27 copies sit on 17 nodes; 19 copies share their node (9 groups, e.g. in G1b NPIPB8+B9+B10P on one 1.03 Mb, 210-exon node (897 kb in G1a, 430 kb in G2a/G2b); NPIPB4+B5 on 743 kb in all four). Median M1 node span is 233-304 kb vs 14 kb in G0. The frozen scorer counts each copy separately in size_strict, so copy-merging nodes are not penalised at FAMILY level, and one long node "covers" several copies full-length. G1b's 0.818 (≥ A1's 0.806) and 25/27 full-length are therefore not node quality.
- **Why they merge:** `overlap_groups` on reads keyed (chrom, strand) chains adjacent copies through reads whose blocks overlap, including spliced reads whose blocks sit on both sides of 100-350 kb N gaps, so read-free stretches of up to 237 kb are bridged (inside the NPIPB4+B5 group by a single read); the AF-3 split (G2) does not cut these links (it only cuts at < 2 linking reads), and MAX_INTRON cutting fires in 0-2 groups per arm.
- **Spliced-read strand (G1b/G2b vs G1a/G2a)** moves 1,486 of 10,429 unspliced reads to the spliced majority strand (501 have no spliced overlap). It removes nodes (259 → 191; 379 → 297) and non-NPIP members of the NPIP family (27 → 10 other), raising F strict by +0.24 / +0.18. It is the one factor here that helps without copy merging being the whole story, but it was only tested inside the merged read-group arms, not on G0.
- **Footprint gates:** every read-group arm loses 59-66 of 68 G0 families by J < 0.5 and 5 by zero bases; A1 (target) also fails these gates (65 lost, 25 zero-base). The gates as declared cannot separate "different node footprints" from "families destroyed".
- SUBFAMILY levels: R 0.704 / 0.148 in every arm (no arm splits NPIP); only strict columns move, in the same order as FAMILY.

## Readings (full text in DECLARATIONS.txt)

- Reads: the ladder's 52,917 primary MAPQ >= 1 reads from `reads.bam` (the `-M -L regions.bed` extraction), strand from ts:A as `read_blocks()`; not clipped further.
- G3: single-exon dump rep strand = s if >= 2/3 of any-strand reads with a block overlapping its exon have strand s, else kept; then shipped consolidate + read-locus nodes (38 added, as G0).
- G1: one node per MAX_INTRON piece (>= 100 bp) of each >= 3-read group; n_reads = group size; rep = most-supported exact intron chain among the piece's spliced reads, ends at median_low of first-block starts / last-block ends, ties more exonic bp then leftmost then smallest chain; no spliced read → rep = piece exons (106/259 G1a nodes, 48/191 G1b).
- G1b: unspliced read (1 block) takes strand s if >= 2/3 of spliced reads with a block overlapping its block have s; spliced reads unchanged; simultaneous.
- G2: `split_linked` on the group's depth-2 segments; sub-locus kept if sup >= 3 and >= 100 bp (shipped `consider()`), then G1 steps with n_reads = sup (48 G1a groups yield >= 2 kept sub-loci (89 have >= 2 components before the sup/bp filter), 276 sub-loci dropped; G2b 44 (84) / 245).
- Queries named by sequence md5; PAF records reused from `npip_ladder/union/` and `sd_fold/map/`; 450 new tx + 441 new body queries mapped. Edges, triangle leaders, components, scorer: copies of the verified mirror.
- Footprint Jaccard: strand-blind merged member exons per chrom; "new" = arm family with best J < 0.5 against all G0 families.
- Best arm for per-copy table: no eligible arm → highest-F-strict non-G0 arm (G1b), labelled INELIGIBLE.

## Provenance

- Repo `/mnt/c/Users/jfris/Desktop/Rustle`, HEAD `63a0843c`, `src/` clean; `shared_definition.rs` unchanged since ladder HEAD `111f6727`. minimap2 2.30-r1287. Python 3.14.4.
- Work dir `/mnt/linuxdisk/home/juanfraitu/sd_readgroup/`: `DECLARATIONS.txt` (+ `.md5`), `scripts/{build,score,addendum}.py`, `arms.pkl`, `map/{tx,body}_chunk*.{fa,paf,err}`, `map/new_{tx,body}.paf`, `results.json`, `families.pkl`, `addendum1.json`, `logs/`.
- Reused: `npip_ladder/{reads.bam,regions.bed,results.json,rust/dump/x.nodes.tsv,verify/nodes.pkl,verify/edges.pkl,union/,idx/}`, `sd_fold/charz/lib.py`, `sd_fold/map/`.
- Asserts passed: G0 node list = ladder A4 (560, row for row); G0 triangle families = results.json A4 (membership, order); G0 triangle and components core rows (3 levels) = results.json A4; A1 recompute core rows = results.json A1.

| step | command | wall | peak RSS |
|---|---|---|---|
| build arms | `python3 scripts/build.py` | 63 s | 0.28 GB |
| tx mapping (450 q, 1.09 Mb, 1 chunk) | `minimap2 -c -N 50 -p 0.1 -x splice -uf -t 4 npip_ladder/idx/target.splice.mmi map/tx_chunk0.fa` | 46 s | 15.8 GB |
| body mapping (441 q, 22.8 Mb, 6 chunks) | `minimap2 -c -N 50 -p 0.1 -x asm20 -t 4 npip_ladder/idx/target.asm20.mmi map/body_chunk{0..5}.fa` | 42+35+29+29+32+29 s | 13.4 GB |
| score + blast radius | `python3 scripts/score.py` | 3 s | 0.24 GB |
| addendum 1 | `python3 scripts/addendum.py` | < 2 s | – |

## Caveats

- Development family, copy-centred windows, one CHM13 haplotype, one testis Iso-Seq library: nothing here is validation.
- The frozen FAMILY scorer does not penalise nodes that merge copies; any node rule that lengthens nodes can raise F strict and full-length without better copies (Addendum 1 shows this is what happens in G1a-G2b).
- Body queries of read-group nodes reach 1.1 Mb; asm20 chaining on such queries was not inspected.
- The footprint-Jaccard eligibility gates also reject the annotated target (A1); a NO CANDIDATE outcome from them is weak evidence against an arm by itself. For G3 the outcome would be NO CANDIDATE even without the gates (+0.036 < 0.05).
- 2/441 new body queries had no PAF record.
- Independently verified afterwards (section below); verification is a recompute, not validation.

## Verification (independent recompute)

Verifier code: `/mnt/linuxdisk/home/juanfraitu/sd_readgroup/verify/` (`v_build.py`, `v_edges2.py`, `v_score2.py`, `v_percopy.py`, `v_samtools.py`, `v_bridge.py`, `v_splitcount.py`). It does not read or import the builder's scripts. The scorer, edges, triangle-leader and component functions are reused from the independently verified ladder mirror (`npip_ladder/verify/`). `split_linked` and the G1/G2/G3 rules were written again from `shared_definition.rs` (HEAD `63a0843c`) and DECLARATIONS.txt. Development family, not validation.

- **Declarations came first.** The md5 of DECLARATIONS.txt up to the `ADDENDUM 1` banner is `f160935a` (it matches DECLARATIONS.md5, written 05:57:17). The file itself is dated 05:57:14 for that part; arms.pkl is 05:59, the PAFs 06:00-06:03, and results.json 06:04:29. The addendum was saved at 06:04:56, after results.json, and is disclosed. Its banner says "06:1x", but the file mtime is 06:04:56.
- **Nodes: 0 rows differ across all 6 arms.** My rebuild from `reads.bam` reads (52,917) and the dump matches arms.pkl row for row on chrom, strand, n_reads, exons and rep_exons: G0 560, G3 483, G1a 259, G1b 191, G2a 379, G2b 297. G0 also equals the ladder-verified A4. The build counts also match:
  - G3: 215 reps flip + → −, 144 stay +, 12 have no reads.
  - Unspliced-read restrand: 1,486 changed, 8,442 kept, 501 with no spliced overlap (10,429 unspliced reads in all).
  - Reps with no spliced read: 106/259 (G1a) and 48/191 (G1b).
  - MAX_INTRON cuts: 2 / 1 / 0 / 0.
  - Sub-loci dropped: 276 (G2a) and 245 (G2b).
  - "Groups split" is 48 / 44 only when it means groups with ≥ 2 *kept* sub-loci. Counted before the sup/bp filter, 89 / 84 groups have ≥ 2 components. The Readings line now says so.
- **PAFs.** All query FASTA names equal the sequence md5. No new query duplicates a union or sd_fold query, and no query has records in more than one source. Every arm's tx/body query is present. The concatenated chunk PAFs equal `new_{tx,body}.paf`, and 2/441 new body queries have no record. I reran minimap2 on `tx_chunk0.fa` and `body_chunk5.fa` with the declared indexes and flags, and both outputs are identical to the builder's (sorted records).
- **Edges and families.** Exon, body and pair counts match for all 6 arms: 1349/1885/1919, 1085/1583/1622, 607/730/811, 261/313/377, 1076/1236/1441, 466/575/708. Triangle families are identical to families.pkl in ids and order. G0's triangle and component families equal ladder A4 in results.json. The A1 families from the ladder verifier equal families.pkl A1 (762 pairs, 461 exon / 711 body).
- **Scores.** Every core row (triangle and components, 3 levels, 7 arms including A1) is string-identical to results.json. The following all match, apart from the one table fix below:
  - Read-first table
  - SUBFAMILY table
  - Composition, full-length, opp-strand, shared counts
  - M1 span medians: 303,547 / 286,408 / 234,057 / 232,849 bp vs 14,280 bp in G0
  - The 27-row per-copy table: all cells recomputed, 0 mismatches
- **Blast radius and decision.** The lost / zero-base / largest-J / changed / new (0-bp) counts and the lost and zero-base id lists match for every arm, A1 included. G3's lost family sizes are 12, 6, 5, 4, 3, 3, 3, 3, 2, 2, 2. No arm is eligible. F strict minus G0: G3 +0.0365, G1a +0.131, G1b +0.368, G2a +0.032, G2b +0.209. **NO CANDIDATE is confirmed.**
- **samtools spot-check (25/25 OK), done with my own CIGAR parsing of `samtools view -F 2308 -q 1`, not pysam:**
  - G3: I checked the node strand against the ≥ 2/3 read-strand rule on the node exon for LOC128966608, NPIPB6, NPIPB12 (all flipped to −) and NPIPB5, NPIPB9 (stay +). All agree with the copy strand.
  - G1a/G1b (NPIPB2, NPIPA5, NPIPB4, NPIPB14P, LOC124907807) and G2a/G2b (NPIPA8, NPIPB3, NPIPB9, NPIPB15, NPIPB1P): depth-2 exons inside the copy span and the rep chain (count, median_low ends) are equal to the node. Node strand equals the majority strand of reads on the copy exons.
  - Note: restranding needs a padded fetch. An unspliced read can reach past the node span, where a spliced read overlaps it. With a ±100 kb pad all 25 checks agree.
- **Mechanism check (report-only).** In every copy-sharing group, the read-block coverage has read-free gaps of 30-237 kb. Spliced reads with N gaps of 100-351 kb bridge them, sometimes a single read (e.g. the 237 kb gap in the NPIPB4+B5 group). The "Why they merge" sentence now says so.

**Corrections made to this report**

1. The G2b size_strict in the read-first table was [52]; it is [55] (27 + 3 + 25; P strict 27/55 = 0.491).
2. "NPIPB8+B9+B10P on one 1.03 Mb, 210-exon node" is true for G1b only; the node spans 897 kb in G1a and 430 kb in G2a/G2b. The sentence is now qualified.
3. In G3, NPIPB7's node is extended by two re-stranded single-exon pieces (dump reps 343 and 346), not by one.
4. The caveat "Not independently verified" was replaced.
5. Clarifications, where the numbers were unchanged: the definition of "groups split" (48/44 vs 89/84), and the long-N read bridging in "Why they merge".
