# PREREG — L1 dropped members as O2 candidates + L2 the locus extent owned by O1 (2026-09-05)

Context: `docs/O1_O2_LOOSE_ENDS.md`. Both change what O1 hands O2; neither changes O2's certificates.

## L1 — dropped members become candidates
Now: a member with core = 0 under the family's depth threshold (`CoreStatus::Dropped`) is not emitted as a unit;
its reads have no candidate in O2 (§6ee/§6fg: they end as orphans or as origin-rejected reads on a paralogue).
Change: `mcl_families` emits every cluster member that has a chain (or GFF exons) and ≥ 1 read inside it, with a
new column `member_status ∈ {kept_full, kept_trimmed, dropped, ungated}`. Family membership is the flag; O2's
candidate set is the cluster's loci. `--no-units-include-dropped` reproduces the previous row set (previous
columns byte-identical; three columns appended: `member_status`, `locus_start`, `locus_end`).
Unit merging treats a dropped unit like any other (a base cannot belong to two copies); the representative of a
merged group is kept-before-dropped, then longest exon union.

## L2 — the locus extent is O1's
Now: O2 pads each unit's extent by the family's longest molecule (§6fd), a rule O2 invented and the catalog does
not know. Change: `mcl_families` writes `locus_start`/`locus_end` (0-based half-open) = the union of the chain's
extent and the reference span of EVERY BAM record overlapping the member's locus region (primary, secondary and
supplementary alike: if a record of a read lies at this locus, the locus extends to it). No constant.
`copy_assign` uses the catalog extent when the columns are present (`--read-star-pad-locus` restores the §6fd
padding; catalogs without the columns are byte-identical to before).

## Predictions (3-contig catalog `rna_units_v10` = v9 rebuilt; same BAM `npip3.bam`)
- P1 (byte identity): `rna_units_v10 --no-units-include-dropped`, first 16 columns == `rna_units_v9.units.tsv`;
  `units.fa` identical. `copy_assign` on a v10-derived sweep dir with the new columns REMOVED == the sweep_v10
  output of the same family (the pad rule unchanged).
- P2 (L1 on NPIP `MCL1_073242` and LCR16u `MCL7_073242`): `origin_rejected` falls (NPIP: by ≥ 100 reads, the
  admission run moved 558 → 333) and reads whose primary lies in a dropped unit are assigned or tied to it with
  MAPQ-60 placement agreement ≥ 0.99; the audited 62 anchors: 0 wrong.
- P3 (L2 on the paired 35 under the catalog extent vs `--read-star-pad-locus`): origin-rejected reads do not
  increase (the extent contains every record of O1's BAM by construction); assigned within ± 1 point of the pad
  rule; MAPQ-60 placement agreement ≥ 0.998.
- P4 (extent size): median locus extent ≤ median padded extent (the pad adds 2 × the longest molecule to both
  sides; reads rarely reach that far) — reported, not a pass/fail.
Fail rules: P1 fails ⟹ a defect, fix before anything else. P2 wrong anchors > 0 ⟹ dropped units are not
candidates by default (escape becomes the default) and the reason goes to the register. P3 rejections increase by
more than 1 % of unit reads ⟹ the extent definition misses records (report which) before shipping the default.
Report every number for both arms; nothing is chosen after looking.
