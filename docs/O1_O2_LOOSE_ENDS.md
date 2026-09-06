# O1 → O2 loose ends — the ordered queue (2026-09-05)

Assessment after §6fg (orphan reporting closed). Six loose ends between O1's catalog and O2's assignment;
two of them are the composition principle itself ("O1's vertex set IS O2's path set"). Tackle in this order;
each item gets a PREREG before the first run and a ledger section when closed. Status column is the live record.

| # | Loose end | Why it matters | Fix | Status |
|---|-----------|----------------|-----|--------|
| L1 | **Dropped members leave O2 blind.** The core rule's drop arm (core = 0) removes a locus from the unit table, so its reads have no candidate and surface as "missing copy" signals (NPIP's ABCC1-region records: 332 reads; the roster-admission prototype re-admitted exactly those loci). | A false O3 signal, manufactured by O1's own membership rule. MCL clustered the locus by homology, so it IS a competitor for the family's reads whatever its core status. | `mcl_families` emits dropped members as units with `member_status = dropped`; O2's candidate set = every locus of the cluster, family membership stays a flag. Escape `--no-units-include-dropped` (previous row set). | **CLOSED §6fh** (17 dropped members emitted; their reads assigned only to themselves, 0 stolen) |
| L2 | **O1's node ≠ O2's candidate.** O1 emits the unit (read-supported chain inside the annotated locus); O2 scores the locus padded by the family's longest molecule (§6fd). The anchor read whose exons lie 4–16 kb outside its unit proves the unit is short and the padded locus right. | The chapter's composition sentence is false by that margin, and the padding rule lives in O2 instead of the catalog. | O1 emits `locus_start`/`locus_end` = the read-supported extent (every record overlapping the locus, unioned with the chain); O2 consumes exactly those (`--read-star-pad-locus` = the §6fd rule, byte-identical on catalogs without the columns). One definition, owned by O1. | **CLOSED §6fh** (two extent definitions refuted, rows 713–714; paired 35: 46.5 % vs 43.4 % assigned, rejections 4,327 vs 5,508) |
| L3 | **Single-candidate ties never retested under the pairwise rule.** Row 704's `-p 0` result was measured under the intersection rule, which was itself the defect. | Under pairwise certificates a weak competitor is a pair with few shared columns; reporting every chain may convert part of NPIP's single-candidate ties (1,587) into certified assignments or honest K = 0 ties. | One pre-registered NPIP run with `RUSTLE_STAR_P=0` under the current default; compare assigned / wrong anchors / agreement. | **CLOSED §6fi/§6fj**: `-p 0` not adopted (row 716); user decision 2026-09-06: a certified sole candidate is `assigned` with `sole_candidate = 1` (column 18), a rejected single keeps `ambiguous`; escape `--no-sole-candidate`. |
| L4 | **Per-copy abundance undefined under read-star.** `quant.tsv` gets hard argmax counts; the EM sees no observations. Found while gating L1/L2: an ORPHAN row (`n_candidates 0`) is counted in copy 0's `n_reads_hard` (MCL106: 133 → 137 with 4 orphans), so the hard count is wrong too. | The quantification deliverable of O2 is empty in the default mode. | Soft counts = Σ per-read posterior per copy from the pairwise softmax; write beside `n_reads_hard`. | **CLOSED §6fj**: `abundance` from posteriors under read-star, `n_reads_soft` column, orphans out of `n_reads_hard`. |
| L5 | **Genome-wide O2 has never run.** `gw_units_v1` exists; family regions for dispersed families span whole contigs (NPIP: 90 Mb), so read loading scales badly. | The thesis number for O2 is 3-contig only. | Regions were never the cost: §6dh already gathers reads from the copies' 50-kb neighbourhoods, so a contig-wide region binds without loading it. Done: `gw_units_v2` (L1/L2 columns) + O2 over its 754 families. | **CLOSED §6fk**: 71,047 unit reads, 66.9 % assigned, agreement 44,364/44,365, 27.8 % origin-rejected (O3's material), 30 min wall. |
| L6 | **No read-level proof dump for read-star.** `--dump-psv` writes placeholders in that mode. | The chapter's assignment-proof figure (a read's own columns against its candidates) has no source. | Dump per read: candidate, column position, read base, each candidate's base, from the read-star observations. | **CLOSED §6fj**: `--dump-star` → `<out>.star_reads.tsv`. |

Not in the queue (measured, settled): the opt-in forms of §6fd (`--origin-substitutions-only`, `--read-star-junctions`,
`--read-star-unit`) — rows 707–710; MCL15/MCL28 residual artefacts (cell B, never scored); the human PMS2P14-type
exon-less block spans (Soto slice only); MCL scattering of Soto families (row 694, 9 members).

Cross-references: `docs/ROADMAP_O1_O2.md`, `docs/O1_O2_COMPOSITION.md`, `docs/CHAPTER_O1_O2_COMPOSITION_DRAFT.md`,
`docs/NEGATIVE_RESULTS_REGISTER.md` rows 689–712, ledger `docs/o1_ledger.md` §6et–§6fg.

## After the queue: the O3 flag pass (§6fl, 2026-09-06)
Done (`bench/o3_flag_pass.py`, PREREG `docs/PREREG_o3_flag_pass_2026-09-06.md`; P1 held, P2 held on the 3
contigs and failed genome-wide, P3 missed by 4 points, P4 done). The four follow-ups (2026-09-06, §6fm–§6fn):
(1) the 15 suspect controls → 13 are the assembly's haplotype at 0.5 %, 2 exposed an L2 defect (extents
containing other families' units; fixed: clip at every unit, hit must touch the unit; reruns `sweep_v14`,
`sweep_gw_v3`); (2) reconstruction over every flag → `reference_absent` by construction, a flag means "≥ 0.7 %
from every reference locus", copy vs haplotype undecidable from RNA (row 720); (3) annotated-no-unit loci →
44 short-chain members / 9 other families / 67 unclustered (partial paralogues below the coverage threshold);
folding them in is a definition decision (§6fn); (4) cut certificate → 21 genome-wide cuts at ≥ 0.90 identity,
44 nested, 27 without homology (row 721).
