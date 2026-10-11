# Compatible-read support in the assembler's chain admission, 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** User (2026-10-11): the advisor wants the assembler benchmarked
against StringTie and FLAIR; make the primary-only assembler as close to them as possible with their tricks, extendable later to the
seeded loci (secondaries, AS-tied primaries resolved by O2).

**Diagnosis that motivates it (dev only, gorilla OR6737).** Figure 1 rebuilt on the current assemblies (identical to 2026-09-25): we lead
on intron-chain sensitivity / precision (0.269 / 0.356 vs StringTie 0.244 / 0.346) and trail on intron sensitivity (0.593 vs 0.639).
Genome-wide, 18,961 reference introns are StringTie's and not ours (6,329 the reverse); their primary junction reads have a median of 4.
On NC_073244.2 (chr20, the dev contig): of 941 such introns only 20 exist before polish; **866 (92%) are carried by >= 2 primary reads
of which no two share an exact intron chain** — the assembler groups reads by EXACT chain and admits a chain only at `>= 2` reads
(`PASS1_MIN_READS`), so truncated or end-variant reads of one isoform never pool. StringTie's splice graph pools them.

**Not a re-proposal.** Register 848 / 855 (read-isoform widening: a chain admitted when each of its junctions has >= k reads counted
over ALL spliced reads of the region) failed (+50% transcripts for +2 NPIP junctions; refuted held-out on human chromosomes): it let a
chain borrow support from reads of OTHER isoforms. Here a chain borrows only from reads COMPATIBLE with it, and only chains that some
read actually carries are candidates — no new combination of junctions is ever built.

## The rule (one flag, `copy_assign --compat-support`, default off = byte-identical)

In pass 1 (`Pass1Acc`), after reads are grouped by exact `(chrom, intron chain)`: a group D is COMPATIBLE with a group C when D's chain
is a contiguous sub-chain of C's, `D = C[i .. i+m)`, D is not C, and D's reads stay inside C's exons around it — every read of D starts
at or after the acceptor of C's intron `i-1` (when i > 0) and ends at or before the donor of C's intron `i+m` (when i+m < len C), tested
on the group's outermost start and end. For each junction t of C, its compatible support `s_C(t)` = reads of C + reads of every
compatible D that contains t. A group with fewer than `min_reads` exact reads is ALSO admitted when `min_t s_C(t) >= min_reads`; its
skeleton then carries `n_reads = min_t s_C(t)`. Groups already admitted are untouched. Unspliced reads, the gate, polish and every
other step are unchanged; `min_reads` is the existing floor (2), so the rule has no new parameter.

## Test

**Arms.** OFF (current) and ON, each in the two seeding pools: `primary` (the benchmark arm) and `good` (the default: primaries +
secondaries with AS >= 0.98 of the molecule's best). Scored with gffcompare against the same reference as figure 1
(`assembly/<species>/eval_all/ref.gtf`).
**Dev.** Gorilla OR6737 NC_073244.2 only (used for the diagnosis): correctness and sanity, no tuning (there is nothing to tune); at
most ONE amendment, written before any held-out number.
**Held-out (decides).** (H-g) gorilla OR6737 genome-wide minus NC_073244.2; (H-h) human A119b genome-wide. Per sample and pool:
- **B1** intron sensitivity up by >= 1.0 point;
- **B2** intron-chain precision down by <= 0.5 point;
- **B3** intron-chain sensitivity not lower.
Reported: transcripts emitted, transcript-level sensitivity / precision, locus level, the share of the StringTie intron gap closed.
**Verdict.** B1-B3 hold in both held-out samples for BOTH pools => ADOPT as the default (and as the base for the seeded-loci extension);
for the primary pool only => default for the benchmark arm, opt-in otherwise; else NOT ADOPTED, kept opt-in with the numbers.

## Result (2026-10-11): NOT ADOPTED, kept opt-in

**Dev (NC_073244.2, primary):** OFF byte-identical to the previous binary; ON +17 transcripts (3,862 -> 3,879), intron sensitivity 54.8 -> 54.9,
intron-chain 28.3 / 40.8 -> 28.3 / 40.7. Why it barely moves (dev diagnosis, before any held-out number): for 864 of the 927 missed
junctions, no read chain carrying the junction reaches a compatible support of 2 at ALL its junctions — the supporting reads are
STAGGERED (each covers a different overlapping stretch), so no OBSERVED chain is long enough. StringTie stitches overlapping reads into a
longer path; this rule, by design, never builds a chain no read carries. No amendment was made.
**Held-out (gffcompare, sensitivity | precision, points):**

| arm | transcripts | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla primary OFF | 68,828 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 |
| gorilla primary ON | 70,375 | 59.9 / 82.6 | 27.2 / 35.0 | 24.4 / 34.8 | 30.3 / 71.0 |
| gorilla good OFF | 70,134 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 |
| gorilla good ON | 71,700 | 60.3 / 82.3 | 27.3 / 34.6 | 24.5 / 34.4 | 30.7 / 69.4 |
| human primary OFF | 238,406 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human primary ON | 250,598 | 61.2 / 63.9 | 22.3 / 18.1 | 20.7 / 15.6 | 34.8 / 22.2 |
| human good OFF | 252,217 | 61.2 / 62.7 | 22.0 / 17.7 | 20.5 / 15.4 | 34.8 / 22.0 |
| human good ON | 263,799 | 61.6 / 62.2 | 22.4 / 17.3 | 20.9 / 14.9 | 35.0 / 21.6 |

(gorilla scored without NC_073244.2; human on the figure 1 reference `eval_annotated/ref.gtf`.) **B1 fails in all four** (intron
sensitivity +0.3 / +0.3 / +0.4 / +0.4, bar +1.0); B2 holds (intron-chain precision -0.3 / -0.4 / -0.5 / -0.4, bar <= 0.5); B3 holds
(+0.4 / +0.3 / +0.4 / +0.4). The flag stays opt-in (`copy_assign --compat-support`, `RUSTLE_COMPAT_SUPPORT=1`): about +0.4 point of
intron-chain sensitivity for -0.3 to -0.5 of precision and +2-5% transcripts. Next mechanism to pre-register separately: stitching
overlapping, consistent reads into longer paths (StringTie's splice graph), after measuring how many missed junctions it could recover.
