# PREREG SEED UNSPLICED FROM PRIMARIES (register 1518, 2026-11-03)

**Written before any number of the arm.** User (2026-11-03), after the tandem rows closed: "are there any other analysis that
could help" -> the single-exon blob finding -> "yes lets do that".

## Finding (arm A of r1495, A119b genome-wide, `/anvme/workspace/v138da13-rustle/fusedgw/A/x.gtf`)
`DN_chr4_9126141_1` (USP17L array): one single-exon model, 265 kb, `reads "3567"`, yet 130 primary records overlap the span and 227
of its 265 1-kb bins have no primary depth. Genome-wide, 707 single-exon models >= 20 kb: 152 with a zero-primary 1-kb bin, 61 with
a zero-primary run >= 3 kb, 42 >= 10 kb (`bench/mono_gaps.sh`).

## Mechanism (established before registration, not a prediction)
The default seeding pool is `good` (`tools/rustle_pipeline.sh --seed-pool good`: `RUSTLE_GTF_SECONDARY=1`, ratio 0.98, the as-table;
r1060/r1100, `docs/archive/2026-09/PREREG_locus_read_pool_2026-09-22.md`). `stream_pass1_region_with` pushes every secondary record
with `AS >= 0.98 x the molecule's genome-wide best`, spliced or not, into the same pool that `cluster_unspliced` (`src/family.rs`)
single-links by span overlap. A near-tied UNSPLICED secondary of a tandem-array molecule lands on every copy of the array, so the
placements chain into one array-wide skeleton: in the USP17L span 8,356 unspliced secondary records from 202 molecules, 3,971 at
`AS >= 0.98 x best` (+130 primaries = 4,101; the model reports 3,567 after coordinate dedupe and polish). The r1060 row scored
spliced chains and loci overlapping annotated genes/pseudogenes; it did not score unspliced seeding separately, so this row
withdraws only the part of that admission the earlier row did not test.

## Change (`RUSTLE_SEED_UNSPLICED_PRIMARY=1`, opt-in, byte-identical unset)
In the streaming pass-1 reader a SECONDARY record with an empty intron chain is dropped before the pool (`stream_pass1_region_with`,
`unspliced_primary_only`); spliced secondaries keep the r1060 rule. Out of scope (unchanged): the buffered reader
(`reads_in_region_indexed_with`, not taken by `assemble` with an as-table), `--tied-seed` (opt-in), the "split at zero-primary runs"
variant (a patch over the same cause; not registered).

## Arms (A119b, r1495 protocol; `~/fused/sgw.sbatch`)
A: `fusedgw/A` as run (RUSTLE_ASSEMBLY_LAYERS=off; the binary change is opt-in so A is not rerun).
B: same with `RUSTLE_SEED_UNSPLICED_PRIMARY=1`. Fast feedback first: `--contig chr4,chr20` assemble only (wall 30 min) for bars 1 and
2, then genome-wide assemble + families for bars 3-5.

## Bars (fixed now)
1. Mechanism check (cannot fail if the mechanism is as stated): B single-exon models >= 20 kb with a zero-primary run >= 3 kb = 0
   (A: 61), and `DN_chr4_9126141_1`'s span holds no single-exon model > 20 kb in B.
2. chr20 SQANTI (sqanti3_qc of each arm's chr20 transcripts vs `~/sqanti_sc_A119b/chr20_ref_cat.gtf`): single-exon isoforms in the
   artifact classes (intergenic, genic_intron, antisense, fusion) B <= A; single-exon FSM/ISM B >= 0.98 x A.
3. Pair precision >= 0.7158 and recall >= 0.7307 (r1495 guards); merge_audit fused loci B <= A (547).
4. Recall of single-exon genes: CAT single-exon genes (`chm13v2.0_CAT_Liftoff.slim.gff3.gz`, one exon in every transcript) overlapped
   by >= 1 same-strand single-exon model, B >= 0.98 x A (the r1060 fill-in was partly pseudogenes; this is the cost side).
5. Spliced models unchanged: transcripts with >= 2 exons identical in (chrom, start, end, introns) between A and B >= 0.99 of A.
PASS = all five -> default on. Any bar failed = FAIL, rule stays opt-in.

## Panel result (chr4 + chr20, assemble only, 2026-11-03 18:55, before the genome-wide arm)
Bar 1 PASS: single-exon models >= 20 kb with a zero-primary run >= 3 kb 4 -> 0; the USP17L span holds no model at all in B (the
85-primary bin at chr4:9.17 Mb seeded nothing; not investigated). Bar 5 PASS: 14,338 of 14,338 spliced models identical. Bar 2: the
sqanti3 call failed (`Bio` missing: the env's bin must be on PATH before python resolves, as r1515 did); rerun separately.
Bar 4 FAIL as registered: 90 -> 46 genes hit, ratio 0.51. 41 of the 44 lost genes are the USP17L array: they were "hit" by the
265-kb blob itself, so the registered metric counted the artifact as recall. The 3 real losses: FTLP, SDAD1P, AL078621.

## Amendment 1 (written before any genome-wide number; the registered bar 4 is still reported)
Bar 4 counts a gene as hit only by a same-strand single-exon model whose span is <= 2x the gene's span
(`bench/seed_unspliced_score.py`, "Amendment 1" line). Bar 4 is judged on the amended line; the as-registered line is reported beside it.

## Bar 2 (chr20 SQANTI, rerun with the env on PATH; single-exon isoforms only)
A: intergenic 117, genic_intron 95, antisense 32, genic 9 (artifact 253), FSM 16, ISM 0. B: 91 / 95 / 32 / 9 (227), FSM 14, ISM 0.
Artifact half PASS (227 <= 253, all of it intergenic). FSM half FAIL as registered: 14 / 16 = 0.875 < 0.98 (2 single-exon FSM lost).
