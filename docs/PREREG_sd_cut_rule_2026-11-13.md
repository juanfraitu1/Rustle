# Row 1687 — SD cut rule for models fused across tandem duplicates (2026-11-13)

Goal (user, 2026-11-13): explore other ways to include segmental duplications. User's ruling after row 1686
("Register the SD cut rule"): go on although the ceiling (28 of 69) missed my bar of 35. **Opt-in only.** Assembly
layer; under the switch the families input changes, so `fam.copies.*` and `assign` output would differ from the
default; defaults unchanged; not for `advisor_share`. No read is assigned to a copy by this row.

**Rule S (blind, frozen; `bench/sd_cut.py ASSEMBLY.gtf SEDEF.bed OUT.units.tsv`).** No annotation is read.
1. SD pairs used: both mates on one contig, disjoint, in direct orientation (the two strand columns equal).
2. For a transcript with >= 2 exons and such a pair (left mate, right mate): p = the midpoint of the gap between the
   mates. The transcript is cut when (a) one of its introns contains p; (b) its exons left of p have >= 100 bp inside
   the left mate and those right of p >= 100 bp inside the right mate; (c) on each side, at least half of the
   transcript's exonic bases on that side lie inside that side's mate (the span idea of `--sd-damp`: the side IS the
   duplicated unit, which an exon duplication inside one long gene does not satisfy).
3. No cut when an exon spans p or no intron contains it. One row per cut transcript in the `--bridge-units-list`
   format (label `sdcut`); a transcript cut by several pairs gets every cut intron.

**Step 0 (list only, A119b default assembly `r1662_edgemiss/e2e/a1.gtf`, CHM13 SEDEF; stop rule).**
`bench/annot_compare/cut_audit.py` against the benchmark RefSeq GTF. GO needs both: cuts inside one RefSeq coding gene
<= 10% of audited cuts, and >= 15 assembly loci (gene_id) with a cut. The 15 is fixed now, in absolute loci, about
half of row 1686's ceiling and three times rows 1516/1517. Otherwise STOP: no arm is run.

**Full arm (only after GO; control re-run today with one pinned binary and a copied driver).** Control = `f1units`
with F1's own bridge list; arm = that list plus Rule S rows built from the control's own assembly; f1v2 default
reported beside. Bars against the control: pair precision and recall each within 0.01 (pairscore, simiiformes);
fused loci (merge_audit) not above control; spliced tandem-class fused loci (recounted on the control) down by >= 15;
cuts inside one gene <= 10%. PKD1P6-NPIPP1 reported as in row 1516.
**Samples.** A119b is development (the 28 are seen). OR6737 with the gorilla SEDEF decides: same safety bars, and
>= 5 loci cut between two different coding genes.

## Step 0 result (2026-11-13, `r1687_sdcut/step0.out`) — GO

Rule S on the A119b default assembly: 128 transcripts cut in 79 assembly loci. Cut audit: 9 of 128 cuts inside one
RefSeq coding gene (7.0%, bar <= 10%); 13 between two coding genes; 7 overlapping gene sets; 99 have a side with no
coding gene (non-coding or pseudogene side, or unannotated). Stated plainly: the 10% bar is row 1516's, over all cuts;
among the 29 cuts with coding genes on both sides, 9 (31%) are inside one gene. The full arm decides.

## Amendment 1 (2026-11-13, before any arm number)
- Three arms in one job, one pinned copy of the driver and binaries: D0 = default (f1v2); E = `f1units` with F1's own
  list (rows of D0's `x.bridge_junctions.tsv` with `bridge` True, one row per TJ transcript); S = E's list plus Rule S
  rows built from D0's `x.gtf` (a transcript in both gets the union of its introns). All with today's defaults
  (CHM13 SEDEF damping, exon-sharing edge, extra representatives).
- Tandem-class recount: the 69 spliced gene pairs of `~/fused/tandem_fused.tsv` are the class list; a pair counts as
  fused in an arm when one `fused` locus of that arm's merge audit holds both genes. The bar is S at least 15 below E.

## Result, A119b development (2026-11-13, job 4470775, workspace `r1687_sdcut`) — FAIL on reach

| arm | loci | fused | tandem-class pairs still fused | pairs predicted | hits | precision | recall |
|---|---|---|---|---|---|---|---|
| D0 default (f1v2) | 74,701 | 544 | 69 of 69 | 1,960 | 1,768 | 0.902 | 0.808 |
| E f1units, F1's list | 75,510 | 539 | 68 | 1,944 | 1,769 | 0.910 | 0.809 |
| S = E + Rule S (126 new rows) | 75,555 | 535 | 65 | 1,954 | 1,775 | 0.908 | 0.811 |

Rule S on D0's assembly: 129 transcripts cut in 81 loci; cut audit 10 of 129 inside one coding gene (7.8%).
PKD1P6-NPIPP1: no fused row in S's audit.

Bars against E: precision -0.002 and recall +0.002 (within 0.01, pass); fused 535 <= 539 (pass); inside one gene 7.8%
(pass); tandem-class fused pairs down by 3, bar 15 (**fail**). The rule is safe and nearly inert: it separates 3 of the
69 class pairs and adds 6 true and 4 false family pairs. Row 1686's ceiling of 28 counted SD pairs over RefSeq gene
spans; the blind rule needs a transcript with an intron over the gap between disjoint direct mates and most of each
side inside its mate, and few of the 28 loci offer that. OR6737 was not run (A119b did not pass). Row closed; nothing
wired; the detector `bench/sd_cut.py` stays as a bench script.
