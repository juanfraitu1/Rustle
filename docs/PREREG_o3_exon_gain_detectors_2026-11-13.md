# PREREG row 1694 (2026-11-13): which existing detector fires on the exon gains reads can see? Measurement only

User go 2026-11-13 ("O3 open list", first item: detector sensitivity on exon gains). Rows 1417 (BPCLUST) and 1418 (EXHAP)
failed at family level; row 1420 (EXCEIL) showed 41.9% of gaining exon units leave >= 2 reads only the haplotype-only
sequence explains (control 14.0%). Nobody joined the three tables per unit. This row does that and builds no detector.

## Measure
Units = the 272 gaining and 300 control exon units of row 1420. VISIBLE = EXCEIL level B (>= 2 reads closer to a
haplotype-only exon sequence than to every reference copy). For VISIBLE gaining units and VISIBLE control units, the
share on which each existing detector fires: EXHAP certified tier, EXHAP putative tier, BPCLUST (its family flag).

## Reading rule (written before the join)
- "The detector is the limit" (row 1420's reading) stands when the best existing tier fires on < 10% of VISIBLE gaining units.
- A tier is worth a follow-up row only when it fires on >= 25% of VISIBLE gaining units and at most a third as often on
  VISIBLE control units. Otherwise no follow-up is proposed.
Either way the O3 status sentence (`docs/LEDGER_O3.md`, "Defensible sentence") is unchanged by this row.

## Inputs
Per-unit outputs of rows 1416-1420: `units.tsv`, `copies.tsv`, `truth.tsv` (`bench/exhap_units.py`), the EXHAP
`*.families.tsv` / `*.putative.tsv` per unit, EXCEIL's per-unit level B table, `bench/bpclust_truth.tsv`.

## Status: NOT RUN, inputs are not on Fritz
The scripts read `/mnt/linuxdisk/...` (the workstation the rows ran on): `GGO_gwcat.copies.tsv`, `/mnt/linuxdisk/tmp/hapcn`,
the sdcn and exhap directories. On Fritz only the scripts and the two family-level truth tables exist. Rebuilding the
chain here would produce new units against a catalog that may differ from row 1416's 463 families, so its counts would
not be comparable with 41.9% / 14.0%. To run: copy those four inputs to the workspace; the join itself is minutes.

## Amendment 1 (2026-11-13, written before any rebuilt number): rebuild on Fritz (user: "Rebuild on Fritz")
The per-unit inputs are rebuilt here from the data on Fritz. What changes, stated up front:
- Catalog: the workstation's `GGO_gwcat` (494 families, 1,415 copies) is not on Fritz. The rebuild uses the driver's
  legacy `catalog` stage (today's `gw_family_catalog`) on the same fibroblast BAM (`GCA_029281585.2_flnc_mm.bam`) and
  reference (`GGO.fa`). Its family and copy counts are reported next to 494 / 1,415; family ids do not carry over, so
  the dev / held-out split of rows 1417-1418 does not either: all units are used as one set.
- Haplotypes: `KB3781.mat.fa.gz` / `KB3781.pat.fa.gz` (GenBank contig names). Sex chromosomes are excluded by the
  assembly's own chromosome assignment (X on mat, Y on pat, and the reference's X / Y / MT), checked by name before
  counting; unplaced contigs are kept, as in row 1416.
- Method otherwise as rows 1408 / 1416 / 1418 / 1420: copies mapped `splice:hq -uf -N 50 -p 0.5`; exons >= 100 bp
  mapped `-x sr -N 500 -p 0.5`; copies of an exon = merged loci at coverage >= 0.8, identity >= 0.90; unit = distinct
  set of reference loci; GAIN = mat or pat holds more copies than the reference; EXCEIL level B; EXHAP = `copy_phasing`
  per unit on reads whose primary covers >= 0.8 of a unit locus (cap 200 per unit).
- BPCLUST is dropped from the join (its per-family flags are not on Fritz and its pipeline is not rebuilt).

Validity gates (the row is INVALID, and reported as such, when one fails):
- G1: every exon's own reference locus is recovered in the reference mapping for >= 99% of exons (row 1416: 100%).
- G2: the share of families with >= 1 exon differing between mat and pat is within 24.0% +- 8 points.
- G3: EXCEIL level B on the rebuilt units is within 41.9% +- 10 points for gaining units and 14.0% +- 10 for controls.
The reading rule above is unchanged and is applied only when G1-G3 hold.

## Outcome (2026-11-13, jobs 4470935, 4470947, 4470960): the detector is the limit; no follow-up tier
Rebuilt on Fritz (workspace `r1694_exgain`, scripts `bench/exgain/`): legacy catalog 1,210 multi-copy families / 6,076 copies (workstation 494 / 1,415), 16,223 autosomal exons in 1,176 families, 15,469 units, 469 gaining. Gates: G1 99.99%, G2 21.7%, G3 45.4% gaining / 21.0% control, all inside the registered bars. Join: of 213 VISIBLE gaining units the certified tier fires on 0 and the putative tier on 5 (2.3%); of 63 VISIBLE control units 0 and 0. Over all 14,815 units with reads: certified 0 gaining / 14 non-gaining, putative 7 / 22.
Reading rule: the best existing tier fires on 2.3% of VISIBLE gaining units (< 10%), so row 1420's reading stands; no
tier reaches 25%, so no follow-up row is proposed. The O3 status sentence is unchanged.

Disclosed deviations:
- **A fix made after gate numbers were seen.** The first chain run (4470947, kept in `run1_nosecondary/`) failed G1
  (96.60%) and G2 (5.1%): minimap2 2.30's `sr` preset prints no secondary records, so `-N 500` had no effect and
  multi-copy exons had one hit each. `--secondary=yes` was added to the three exon mappings and the chain rerun. The
  cause was confirmed on a small test before the rerun, and the amended method ("-x sr -N 500") presupposes those
  records; but the fix did follow a failed gate. Read strictly, the row is INVALID from run 1; the user decides.
- G1's test was written before any count as: a reference hit at coverage >= 0.8 and identity >= 0.90 overlaps the
  exon's own coordinates (the original QC code is not in the repo).
- Sex contigs: reference NC_073247.2 / NC_073248.2 / NC_011120.1; maternal X CM054605.2, paternal Y CM054581.2 (from
  the GenBank assembly headers). No dev / held-out split; BPCLUST dropped (Amendment 1). Today's `copy_phasing` build.
- Not checked: whether the higher level B shares (45.4% / 21.0% vs 41.9% / 14.0%) come from the larger catalog.
The rebuild itself was run by a delegated agent; the join table, the gate files and the mapping command were re-read
before this record.

User ruling 2026-11-13: "Accept the fix". The run with `--secondary=yes` is the result of this row; the deviation stays disclosed above.
