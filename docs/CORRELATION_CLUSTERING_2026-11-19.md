# Correlation clustering (Strainy / SDA lineage) for assignment and for copies the reference lacks (2026-11-19, PSU)

User's request: look for ways to use correlation clustering in the assignment (O2) or in finding a copy the reference lacks
(O3); the advisor admires Strainy; use the ThesisVault for implementations. Development only (no bar): all four HG002
quarters are spent and were reused under the user's ruling. Workspace `goal_multicopy_dev/cc`, `cc_all` (PSU; linked from `r1823_correlation_clustering`). Register row 1823 (descriptive). Open items: `PARKED.md` item 4.

## The lineage in the vault
| paper | unit clustered | edges | count of copies | what the thesis took |
| Chaisson 2017 (polyploid phasing of collapsed duplications) | PSV sites | + co-read in c_min..c_max reads, - nearby and never co-read | not given (CC) | the formulation |
| Vollger 2019, SDA | PSV sites | attraction by a binomial ratio test of co-read reads, repulsion for mutually exclusive PSVs | not given | PSVs as the signal |
| Kazantseva 2024, Strainy (advisor: "extremely important") | reads | overlap with zero disagreement at informative SNPs; label propagation, re-clustered per cluster until the SNP set is stable | not given | the read graph, per-cluster refinement |
| Chen 2025, Paraphase | reads realigned to ONE archetype copy | haplotypes assembled from the reads' bases at variant sites | haplotype count = copy number | not ported |
| Shaw 2025, devider | reads, positional de Bruijn graph | tied paths abstain | not given | tried, refuted (all six clauses) |

## What is already built: `copy_phasing` (rows 1375-1401, October 2026)
`src/copy_phasing/` (Rust port of `bench/o3_graph/o3port.py`): per copy, the pile of reads aligned there; informative columns
by IsoCon's test; reads clustered by **correlation clustering** (default since row 1378, P5 0.83 vs label propagation 0.66):
reads with the same allele pattern are one node, + when they differ at < 2 shared columns, - at >= 2 (the two-observation
rule), pivot start then local moves until no move lowers the cost. Haplotypes take allele slots (2 per copy); more
haplotypes than slots = a copy the reference lacks (certified at >= 6 reads); `--semi` places reads on sample-supported
distinguishing sites, `--guard` abstains where a copy is flagged. Tested on fresh held-out panels (simulated deletions),
GRCh38 collapses, cross-individual piles; never against HG002's own-haplotype truth until today.

## New today: `copy_phasing` against HG002's own-haplotype truth (`bench/famwise_eval/cp_truth.py`)
Families = the 365 truth families (Compara recent paralogues, 1,243 CAT genes on CHM13); net = every read with a record on a
copy; `--clusterer cc --semi --guard --tie-ms 10`. Truth = `neutral_truth.py` with the repaired rule, all four quarters pooled
(170,727 molecules). Scorer pitfall found and fixed: `assigned` also holds ABSTAIN, HOPELESS (exact ties it refuses) and
`new:N`; an early table counting them as wrong showed 766 false errors.

**O2 (assignment): no gain.** Near-tied molecules the truth decides: 1,673. copy_phasing places 783 (757 right, 26 wrong =
3.3%), abstains on 867, puts 23 on a new copy; minimap2's primary on the same 783 is right 764 (2.4% wrong). Fixes 0,
breaks 7. On clear-lead molecules it places 13 (all right) and puts **807 molecules from reference copies on "new" copies**.
The seen quarter alone gave the same (fixes 0, breaks 1 of 185). Reading: where it decides it repeats the primary; on the
primary's errors it abstains. As rows 1777-1785 found with per-read rules, HG002's alleles differ from CHM13's, so the
reference's distinguishing sites do not separate HG002's copies (row 1783's lesson).

**O3 (copies the reference lacks): a signal, confounded with depth.** Positives = families with >= 6 `outside` molecules
(origin at an HG002 place that no CHM13 locus of the read represents): 13. Certified flags: in all 13; in 15 of the 276
families with no outside molecule. By the family's read count:
| net reads | >= 6 outside: flagged / n | 1-5 outside | 0 outside |
| 200-1,000 | 2 / 2 | 1 / 1 | 0 / 32 |
| 1,000-5,000 | 2 / 2 | 0 / 0 | 7 / 58 |
| >= 5,000 | 9 / 9 | 3 / 3 | 8 / 32 |
Recall 13 / 13 at every depth; false flags grow with depth (0%, 12%, 25%). Caveats: `outside` means the read's best HG002
place is not among its CHM13 loci, which is a copy CHM13 lacks OR an HG002 allele too far from CHM13's copy (not separated
here); 13 positives; development data.

## Ways forward (none run)
1. **O3 flag, judged.** The only promising use. Needs a held-out with own-haplotype truth: KB3781 on its own reference
   (outside = a copy on one haplotype the primary assembly lacks) or HG002 against its own paternal assembly
   (`hg002.pat.fa`: maternal-only copies are outside). Bars to set before looking: recall on >= 6-outside families, and false
   flags WITHIN depth bins (the depth null above). Split `outside` into copy-absent vs diverged-allele first (the truth
   records the HG002 place; check whether it is a separate copy on the haplotype).
2. **Depth-aware certificate.** False flags rise with depth because more reads give more chances for an extra cluster of
   >= 6 reads; a certificate scaled to the pile (e.g. a share of the pile, or a Poisson tail on the cluster size given depth
   and error) is the obvious repair. Development only until 1.
3. **Paraphase-style collapse (not built).** Realign every read of a family to one archetype copy and cluster there, so a
   read's cluster does not depend on which pile minimap2 put it in, and haplotypes (not reads) are matched to reference
   copies. The O2 ceiling is small (the primary is wrong on 3-6% of decidable near ties) and the cross-individual allele
   problem stays; worth it for O3 (an extra haplotype in the collapsed pile is the direct signal), not for O2.
4. **SDA's PSV-level graph (not built).** Nodes = variant sites, attraction by co-read tests, repulsion for exclusive sites.
   Pools weak evidence across reads that individually cover few sites (short exons, RNA). A variant of 3; same verdict.
5. **Not worth it for O2:** every placement rule tried on this truth (seven now, with copy_phasing) repeats the primary or
   abstains; correlation clustering does not change that.

## Files
`bench/famwise_eval/cp_truth.py` (inputs from the truth families; scorer, O2 and O3). Runs: PSU workspace
`goal_multicopy_dev/cc` (seen quarter; settings default, `--semi --guard`, `--tie-ms 10`, `--anchor`) and `cc_all` (pooled).
