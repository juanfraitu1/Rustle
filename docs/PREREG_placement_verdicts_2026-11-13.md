# PREREG PLACEMENT VERDICTS AND A REDUCED SEEDING POOL (register 1701, 2026-11-13)

Why (user, 2026-11-13, after the audit of row 1699): the goal is pool reduction. Where minimap2's primary is sensible
and neither the BAM scores nor the PSVs put it in doubt, it stays. The gate widens to AS >= 0.8 x the molecule's best
and every placement gets a verdict. Then a panel assembly without the excluded placements.
Written before any verdict exists. Prototype only: `copy_assign` and every default stay unchanged until a panel PASS.

## Verdict rule (`bench/o2_contested/placement_verdicts.py`, reuses `certify_pair` of `bench/tool_difference/s1_o2.py`)
Molecules with >= 2 records and second AS >= 0.8 x best (`as_table`). Records judged: the primary and every
non-supplementary secondary with AS >= 0.8 x best. Each secondary against the primary, on the read positions both
records align as = or X and where the two reference bases differ (m columns; ka of them match the primary's
reference, kb the secondary's); p = binomial upper tail of the winner's matches at eps 0.001; alpha 1e-3.
- `excluded`: the primary is certified over the secondary (p < alpha) and the read is not mosaic.
- `over_primary`: the secondary is certified over the primary, not mosaic. Reported; NOT acted on in the judged arm
  (moves off the primary are the error-prone class of rows 1652 and 1653).
- `identical`: m = 0. `undecided`: everything else, including mosaic reads (ka >= 2 and kb >= 2).
One row per record: read, chrom, pos, strand, primary flag, AS and ms and their ratios to the best, spliced, ka, kb,
m, p, verdict. Known limits of the prototype certificate: no junction term, no per-base quality, no origin test.

## Panel (fixed before any verdict)
KB3781 (gorilla, the reference individual). The three chromosomes holding the most decidable molecules of the r1642
truth dumps: NC_086018.1, NC_073242.2, NC_073235.2 (803 of 1,570). Sub-BAM = every record of every molecule with a
record on those chromosomes, so each molecule keeps its genome-wide best score. The base run uses the same sub-BAM.

## Step 0: verdicts against haplotype truth (no assembly)
Truth: r1642 dumps (`pof_panel`, `pof_heldout`), decidable molecules, a placement is on the truth copy when its
aligned span overlaps a truth copy's span. That truth covers catalog near-tie reads only: verdicts in the
0.8-0.98 band are counted but not judged.
STOP unless both: (S1) >= 100 excluded placements of decidable molecules; (S2) excluded placements on the truth copy
<= 2% of them. Reported, not barred: `over_primary` on / off the truth copy; verdict shares per AS-ratio band
(>= 0.98, 0.9-0.98, 0.8-0.9) and per spliced / unspliced secondary.

## Assembly arms (only after step 0 passes; `rustle_pipeline.sh assemble`, same binary, panel chromosomes, wall 1 h)
- BASE: sub-BAM, default pool (0.98).
- ARM (judged): sub-BAM minus `excluded` secondary records, default pool, with BASE's score table (same
  denominators: a pure subtraction from today's pool).
- CTRL-U (control): BASE with `RUSTLE_SEED_UNSPLICED=primary` (rows 1518 / 1520), to tell the arm from that switch.
- W80 and W80-ARM (reported, not judged): pool ratio 0.8 unpruned, and 0.8 minus excluded.
Scoring against the RefSeq gorilla annotation on the panel chromosomes, spliced chains (chromosome, strand, introns):
chains BASE has and ARM lacks (lost), chains ARM has and BASE lacks (new), each split into exact annotation matches
and not; single-exon models counted separately; r1671 secondary-only loci on the panel (shadow / real) that disappear.
PASS iff all: (A1) lost exact-match chains <= 1% of BASE's exact-match chains; (A2) among lost chains the
exact-match share is below BASE's overall exact-match share (what goes is worse than average); (A3) lost + new
chains >= 50 (else the arm is inert on this panel: STOP, no claim); (A4) ARM differs from CTRL-U: fewer than 80% of
ARM's lost chains are also lost in CTRL-U. A PASS licenses a held-out sample row and the port into `copy_assign`; it
changes no default.
Expected: step 0 passes; A3 uncertain (most excluded secondaries may support models that primaries also support).
Not in this row: removing origin-rejected reads (needs `copy_assign` flags at their own gate; separate row).

## Amendment 1 (2026-11-14, before any assembly-arm number)
Job 4470937 hit its 90-minute limit while writing the second filtered BAM: my filter compressed on one thread (50 min
for 8.6 GB); the verdict pass took 30 min. Done before the limit: BASE (53 s, 7,363 transcripts), verdicts, step 0,
the >= 0.98 filtered BAM. No arm was assembled or scored. The filter now uses samtools threads (output unchanged) and
the remaining steps rerun as job arms2. The 1 h wall of the arms is read as the assemblies' wall, not the BAM
rewriting; a port would drop records inside the assembler and need no rewrite.
Filter count: 195,183 excluded records listed, 195,192 removed (9 extra records share read, position and strand with
a listed one).

## RESULT step 0 (job 4470937, `r1701_verdicts/step0.txt`)
S1 4,881 excluded placements of decidable molecules (bar >= 100). S2 33 on the truth copy = 0.0068 (bar <= 0.02).
**PASS as registered (pooled over the bands).** Per band, read after the run: >= 0.98, the only band the judged arm
removes: 30 of 1,014 excluded placements on the truth copy = 0.030, ABOVE the 0.02 bar; 0.9-0.98: 3 of 2,804;
0.8-0.9: 0 of 1,063. The pooled 0.0068 is diluted by the lower bands. Secondaries on the truth copy in the >= 0.98
band (the primary is wrong there): 140 = excluded 30, identical 37, undecided 73, over_primary 0: the reduction
destroys 21% of the true alternatives and the certificate rescues none. `over_primary`, all bands: 1 placement of a
decidable molecule, off the truth copy.

## RESULT arms (job 4471021, `r1701_verdicts/arms.txt`; each assembly 15-22 s on 64 threads)
Verdicts on the sub-BAM (15.2 M records): 1,627,466 rows. Secondaries inside today's pool (>= 0.98): excluded
195,183, identical 101,521, undecided 9,727, over_primary 54. In 0.8-0.98: excluded 798,748 of 864,116.
BASE: 7,324 spliced chains, 2,673 exact annotation matches (0.365), 39 single-exon models.
ARM (pool minus the 195,183 excluded secondaries, BASE's score table): lost 51 chains, 7 of them exact matches
(0.137); new 1; single-exon models 39 -> 31.
A1 lost exact 7 / 2,673 = 0.0026 (bar <= 0.01) PASS. A2 0.137 < 0.365 PASS. A3 lost + new 52 (bar >= 50) PASS, by 2.
A4 lost chains also lost in CTRL-U 0 of 51 PASS; CTRL-U changes no spliced chain at all (single-exon 39 -> 25), so the
separation holds for spliced chains only.
**PASS as registered.** Size of the effect: 0.7% of BASE's chains; 195,183 removed records change 52 chains, so most
excluded secondaries supported models that primaries also support (the expected reading).
Secondary-only loci of row 1671 on the panel: 13 shadow loci (`P`: the reads' haplotype origin is the primary's
locus), ARM removes 10, CTRL-U 9: on these loci the arm is no better than the existing switch; the 1 `both` and 2 `other` loci stay; no real (`L`) locus lies on the panel, so
"real loci are kept" is not tested here.
Reported, not judged. W80 (pool at 0.8, unpruned) against BASE: new 492 chains, 38 exact (0.077); lost 31, 0 exact;
single-exon 39 -> 83. W80-ARM (0.8 minus excluded) against BASE: new 117, 12 exact (0.103); lost 56, 7 exact;
single-exon 44. Against W80, W80-ARM loses 431 chains at 0.079 exact (385 of them are chains the 0.98 pool never had).
Reading: the verdicts remove most of what a 0.8 pool adds, and what is left over BASE is +12 exact chains at 10%
precision against BASE's 36.5%: widening the assembler's pool to 0.8 is not supported by this panel.
Not checked: whether the 7 lost exact-match chains are shadows that happen to match the annotation at an unexpressed
copy (the 2026-10-12 arm found that class) or real models; verdicts in 0.8-0.98 against truth.
Same binary for every arm (`~/rustle_target/release`, built 2026-11-13 19:51, before both jobs).
Next (needs the user's go-ahead): a held-out step 0 judged per band; held-out row on a second sample, then dropping excluded records inside the assembler.

## Row 1704: held-out step 0, judged per band (2026-11-14, written before any HG002 verdict)
Why: row 1701's pooled step 0 hid 0.030 wrong in the >= 0.98 band. Sample: HG002 Kinnex on its own paternal
haplotype (`hg002.pat.bam`, `hg.fam.copies.tsv`); the reads come from both haplotypes, so maternal alleles test the
certificate. Rule, script and constants frozen as in row 1701. Truth: `r1653_o2_guard/hgpat/dump.tsv` and
`r1644_hg002_unspliced/dump.tsv`, decidable molecules; only their records are fetched (no assembly in this row).
Bars, all on the >= 0.98 band (the one a pool reduction removes):
(H1) >= 100 excluded placements of decidable molecules;
(H2) excluded placements on the truth copy <= 2% of them;
(H3) of the secondaries on the truth copy, excluded <= 10% (what the reduction destroys where the primary is wrong).
PASS = all three: licenses a HG002 panel assembly arm. FAIL on H2 or H3: the certificate is tightened (a new row)
before any pool reduction; no port. Reported: the 0.9-0.98 and 0.8-0.9 bands, verdict shares, `over_primary` on / off.
Expected: H2 fails (KB3781 gave 0.030 on the reference individual; maternal alleles should not help).

## Row 1705: planted cases with known read origin (2026-11-14, written before any simulated verdict)
Why (user): a synthetic sample where every odd case is planted and every read's origin is known. Use: find and
explain the verdict rule's failure classes (row 1701 gave 0.030 wrong in the >= 0.98 band). Not accuracy evidence:
earlier simulation passes did not carry to real reads.
Tool: `bench/famsim` (new op `retro` = all introns removed; verifier updated; its tests pass) driven by
`bench/o2_contested/verdict_sim.py`. Synthetic 6-exon template (chain 1,490 bp), random background, the shipped
minimap2 command, `as_table`, then the frozen `placement_verdicts.py verdicts` at ratio 0.8. 200 reads per expressed
copy. Each case at perfect reads (err 0, indel 0) and at err 0.001 / indel 0.0003.
Home of a read = its planted copy; for a copy absent from the reference, the home named below.

| case | planted | expected verdict |
|---|---|---|
| identical | A, A2 = A | identical; nothing excluded |
| psv_0.001 / 0.005 / 0.02 | A, A2 = A + SNPs at that rate | rival placement excluded when the sites suffice, else undecided |
| retro_0.002 / 0.01 | A spliced, AR = A without introns + SNPs | rival excluded; the unspliced higher-AS record is the trap |
| unexpressed | A expressed, A2 = A + 0.005 SNPs with 0 reads | every placement on A2 excluded |
| absent_copy | A, A2 in the reference; A3 = A + 0.01 SNPs expressed, absent | no home; counted only |
| allele_private | A, A2 (0.005); A* = A + 0.003 SNPs absent, home A | placement on A kept |
| allele_minor | A, A2 (0.01); A* = A with exon 3 converted from A2, absent, home A | mosaic -> undecided, home kept |
| allele_major | as above with exons 2-6 converted | predicted failure: reads fit A2; home placement excluded |

Bars. (B1) perfect reads: 0 excluded placements on the home copy in every case except `allele_major`.
(B2) with errors: excluded on home <= 1% of home secondaries, same cases.
(B3) power, perfect reads: >= 90% of off-home secondaries excluded in psv_0.005, psv_0.02, unexpressed, retro_0.01.
Every case must pass famsim's own `verify` (a failed claim voids the case, reported). `allele_major`, `absent_copy`
and psv_0.001 are measurements. Reading registered now: if B1 and B2 hold and `allele_major` shows exclusions on
home, the real-data error is the allele class, which reads alone cannot separate; if B1 or B2 fails elsewhere, the
rule has a defect to fix before any pool reduction.

## RESULT row 1705 (2026-11-14, login node, seconds per case; `/anvme/workspace/v138da13-rustle/r1705_verdict_sim/`)
`summary.registered.tsv` = the registered run. 10 of 11 cases pass famsim's verify; `allele_minor` is VOID (A2 drew
no SNP in exon 3, so the converted allele equals A).
B1 perfect reads, excluded on home: 0 in all nine judged cases: PASS. B2 with errors: 0: PASS.
B3 off-home secondaries excluded: psv_0.005 400/400, psv_0.02 400/400, retro_0.01 400/400, unexpressed 200/200: PASS.
What the bars do and do not show: in every case except `identical` and the converted alleles the primary is on the
home copy for every read, so there is no home secondary to exclude wrongly; B1 and B2 bite only on `identical`
(206-210 home secondaries, all `identical`). The retrocopy trap did not fire with full-length reads (primary on home
400/400 at both divergences). Errors at HiFi rate change nothing (psv_0.02 loses 28 reads below the 0.8 gate).
Measurements. `absent_copy`: all 200 reads of the absent copy keep the primary on A and have A2 excluded.
`allele_private`: home kept, 600/600. `allele_major` (exons 2-6 of A2 at A's locus): 163 of 200 primaries on A2, and
of the 200 home secondaries 124 excluded, 76 undecided: the predicted failure, the same with and without errors.
Amendment after the void (labelled post hoc): `allele_minor6` (exon 6 converted, `amend/`): 163 of 200 primaries on
A2 by 6 sites against 5; all 200 home secondaries `undecided` through the mosaic rule, 0 excluded.
Reading as registered: B1 and B2 hold and `allele_major` shows exclusions on home, so the rule has no defect on the
planted clean cases and the only planted class that failed is the converted allele, which reads alone cannot separate. Not
shown: that the 30 wrong exclusions on KB3781 are that class; partial or truncated reads; the retrocopy trap.

## RESULT row 1704 (2026-11-14, job 4471032, 636 s; `/anvme/workspace/v138da13-rustle/r1704_verdicts_hg/step0.txt`)
HG002 on its own paternal haplotype, 5,873 decidable molecules with a verdict row. Band AS >= 0.98 x best:
H1 excluded placements 6,563 (bar >= 100): PASS.
H2 on the truth copy 604 of 6,563 = 9.2% (bar <= 2%): FAIL.
H3 true-copy secondaries excluded 604 of 2,286 = 26.4% (bar <= 10%): FAIL. `over_primary` 14 on / 17 off the truth copy.
Lower bands: 0.9-0.98 141 of 11,505 = 1.2%; 0.8-0.9 32 of 11,074 = 0.3%. Pooled 777 of 29,142 = 2.7%.
Verdict as registered: FAIL (the expected outcome); no port.
Audit of the 604 (post hoc, does not change the verdict): 337 have the primary on the same truth copy (a second
alignment inside the copy's span, harmless), 267 have the primary off the truth copy = true placements removed.
Corrected band >= 0.98: 267 of 6,563 = 4.1% (bar 2%), and 267 of the 1,391 true-copy secondaries whose primary is off
truth = 19.2% (bar 10%); still FAIL on both. `over_primary` there: 13 on the truth copy. The same split on KB3781
(row 1701): 29 of 1,014 = 2.9%, 29 of 138 true-copy secondaries, 0 `over_primary`; row 1701's figures stand.
Lower bands are not safe: of true-copy secondaries whose primary is off truth, 29 of 58 (0.9-0.98) and 12 of 34
(0.8-0.9) are excluded; the low shares of exclusions (1.2%, 0.3%) come from the large off-truth base.
Post hoc, band >= 0.98 as scored (ON share of exclusions): AS below the best with 1-4 deciding sites 0.10 (kb 0) and
0.22 (kb 1), 5-9 sites 0.04, >= 10 sites 0 of 75; exact AS tie with 1-4 sites and kb 0: 6 of 763. Wrong exclusions sit
where the certificate restates a small AS gap on few sites. Any threshold read off this is post hoc on two spent
samples. The truth is the paternal haplotype only. Row 1705 shows one planted class that produces such exclusions
(converted allele); nothing ties these 267 to it.
