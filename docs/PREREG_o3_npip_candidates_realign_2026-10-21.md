# Pre-registration: do the maternal-only NPIP locus's reads find their copy once the O3 candidates join the reference? (2026-10-21)

**User (2026-10-21):** "ok lets do that". The proposal: add the two candidates from register 1363 to the reference as extra loci,
realign, and see whether C's 37 reads move onto its candidate. This is the read-recovery step of the 10-01 chain (Amendments 5-9
of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`, arms R / R+I+L / C), applied to register 1363's output.

**Question.** IsoCon left most of the locus's reads unclustered (register 1363: C 7 of 37). Does adding the candidates' sequences to
the reference give the reads a home through the aligner, without pulling in other genes' reads?

## Fixed inputs (from register 1363, frozen in `/mnt/linuxdisk/tmp/o3_npip/audit_frozen_4c64ec81/`)

- **Candidates added:** the two counted candidates (>= 2 contigs).
  - **B** = `iso_NPIP_14` + `iso_NPIP_15` (his maternal chr18 6.744-6.756 Mb);
  - **C** = `iso_NPIP_16` + `iso_NPIP_17` (6.767-6.804 Mb).
  - The single-contig `iso_NPIP_19` is under the floor and is not added.
- **Reads:** the 2,955 net reads.
- **Units** (the audit's definition, fixed here before any realignment). A target read (register 1363's R_T, 94 reads) is labelled
  by the aligned bases (>= 50) of its best alignment in T on his maternal chr18:
  - **A** = 6,720,000-6,744,300; **B** = 6,744,300-6,756,100; **C** = 6,766,000-6,822,000;
  - so a read is C, B, A, or "B into A".
  - The audit counted C 37, B 21, B into A 21, A 11.

## Arms

Both arms use the library's command, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (minimap2 2.30).
- **R0:** the net reads realigned to GGO alone (`GGO.splice.mmi`).
- **C:** the same reads aligned to GGO plus the four contigs (`minimap2 -x splice -d`, as Amendment 9 built its index).

A read's call is the locus of its primary: a candidate (B or C, either of its contigs) or a reference position.

## Bars

- **H1 (C's reads find their copy):** in arm C, at least half of unit C's reads have their primary on C's candidate.
- **H2 (B's reads):** at least half of unit B's reads have their primary on B's candidate.
- **H3 (no false moves):** at most 1% of the non-target net reads (2,955 − 94) have their primary on either candidate in arm C.

**PASS** = H1, H2 and H3; **FAIL** otherwise, naming the bars that fail.

**Reported, not bars:**
- the "B into A" and A reads;
- register 1360's 29 reads;
- the MAPQ of reads that moved;
- the target reads AS-tied (within 10) among NPIP copies in R0 vs C, i.e. the O2 domain this O3 step would shrink.

## What a PASS would claim, and what it would not

**Claim.** In his own data, a copy the reference lacks, once named by the reference-free chain, takes back its reads through ordinary
alignment, at a cost of at most 1% of other reads.

**Not claimed:**
- that the chain finds such copies in general (register 1363 is one locus);
- unit A (no candidate);
- other samples.

## Rules

- Run once; bugs → amendment + rerun.
- The index build and the alignments run in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict before it is reported.
- Gorilla only. Push only when asked.
- Driver: new `realign` and `score-realign` steps in `bench/rna_allele/npip_o3.py`, with tests.

## Results (driver 5d04d56e; run once; `/mnt/linuxdisk/tmp/o3_npip/realign/`): **PASS**

- **Index:** GGO + the 4 contigs, built in 3 min 26 s. Both arms aligned all 2,955 reads (R0 2 min 24 s, C 2 min 18 s).
- **Units as registered** (no strand): C 39, B 23, A 11, B into A 21. The audit's 37 / 21 were the main-strand reads; 2 C and 2 B
  reads lie on the opposite strand.

| bar | arm C | bar value | holds |
|---|---|---|---|
| H1: unit C's reads on C's candidate | **33 / 39** (85%; main strand 33 / 37) | >= 20 | yes |
| H2: unit B's reads on B's candidate | **21 / 23** (91%; main strand 21 / 21) | >= 12 | yes |
| H3: other reads moved onto a candidate | **0 / 2,861** | <= 28 | yes |

**Reported:**
- **Where each unit's reads go in arm C:**
  - C: 33 on C's candidate, 6 on the reference;
  - B: 21 on B's candidate, 2 on the reference;
  - A: 6 on C's candidate, 5 on the reference;
  - B into A: 2 on B's candidate, 19 on the reference.
- **The 4 opposite-strand reads** all stay on the reference.
- **Register 1360's 29 reads** (24 in unit C, 5 in A): 24 on C's candidate (21 C reads + 3 A reads), 5 on the reference (3 C
  reads + 2 A reads).
- **MAPQ of target reads on a candidate:**
  - C: 26 at >= 30, 10 at 1-29, 3 at 0;
  - B: 23 at 1-29 (its two contigs are near-identical);
  - homed reads split between C's contigs (iso_NPIP_16 18, iso_NPIP_17 15).
- **The O2 domain:** target reads AS-tied within 10 among loci fall from **34 (R0) to 10 (C)**.

**Reading.** Once the chain names the missing copies and they join the reference, ordinary alignment gives their reads a home:
89-100% of the main-strand reads of the two named units, with no other read pulled in. Reads of unit A, which got no candidate, stay
on the reference, though 6 of its 11 land on C's candidate. Most of the tied reads that O2 would have faced at this locus stop
being ties.

Audit pending.
