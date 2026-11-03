# PREREG r1523: family-wide assembly mode (tied-seed members + fused-model cuts), judged on the r1521 panel

User's goal text (2026-11-04), verbatim:
> Mechanism: for each truth family, members with T >= 2 but no base model get assembled from their tied secondaries as seeds
> (family mode, reads stay per member). Base models spanning two members are cut at the member boundary when each side has
> >= 2 own or tied placements.
> Panel: r1521, 420 primates-level families, A119b, judged under T.
> Bars: complete families >= 204 + 15; fused families <= 49; ordinary-member hit rate not below 0.82; wall <= 30 min on the panel.
> PASS -> default on in final.gtf layers. FAIL -> tied-seed closed.

## Prior diagnostic (before any arm number; ~/famcmp/assemblable.tsv, bench on the r1521 genes table and fusedgw/A/x.gtf)
Per expressed member, the largest intron-chain support from reads whose introns lie within the gene +- 5 kb:
| | chainP >= 2 | chainT >= 2 only | one local spliced read | no local spliced read |
|---|---|---|---|---|
| hit by base (914) | 788 | 65 | 18 | 43 |
| missed by base (222) | 34 | 3 | 107 | 78 |
- The good seeding pool already feeds AS-tied spliced secondaries into pass 1 (65 hits rest on them). Step A can only add the
  3 members reachable by a tied chain alone (PRAMEF7, SPANXC, RAMACL). Step B cannot add a hit (a fused model already hits both
  members under the panel's rule); it can only lower the fused-families count.
- So bar 1 (+15 complete families) cannot be met by this mechanism on this panel; it is kept as the user wrote it.
- The real losses are base-mode: 34 members with a >= 2-read primary chain and no model. Two causes found:
  (a) single-exon models carry the '+' placeholder strand unless RUSTLE_READ_STRAND is set (GMCL2: 5,047 reverse reads, model
  emitted '+'); 30 of the 73 opposite-strand-only misses are single-exon truth under single-exon models;
  (b) PRAMEF8: two primaries, identical canonical 3-intron chain, no clipping, MAPQ 22, no model even with polish none and
  primaries only (region run chr1:12700000-12740000); the pass-1 stage that drops it is not yet identified.
  These are follow-up rows, not part of this one.

## Implementation (smallest)
- Step A: `copy_assign --assemble-only --tied-seed` (existing flag; tied skeletons where no primary skeleton overlaps) on
  region runs around the 3 members; a member counts as gained when a same-strand spliced model overlaps its exons.
- Step B: `bench/famcmp/cut_fused.py x.gtf OUT.gtf` cuts every model overlapping exons of >= 2 truth members of one family
  (each with T >= 2, i.e. every expressed member) at the midpoint between the members' exon spans; exons on each side form
  one model each (single-exon remainders kept). Scored with `family_panel.py` as a fourth tool against base.
Wall: region runs seconds; cutter + panel < 5 min.

## RESULT (2026-11-04; region runs + ~/famcmp/r1523.out; cut_fused.py fixed once before the final score: an exon spanning the cut was truncated, not split)
- Step A, `--tied-seed` region runs at the 3 tied-only members (PRAMEF7, SPANXC, RAMACL): same-strand spliced model on the
  gene 0/3 in both arms (RAMACL gets 5 models in the region run with or without the flag, none in the genome-wide run: the
  r980 sub-region effect, not this row). Gain 0.
- Step B, 55 models cut (truth-family members, each T >= 2). Panel under T, 324 judged families, 1,092 expressed members:
| arm | complete | partial | none | fused families | member hit rate | ordinary (mapq0 < 0.5) | T-only |
|---|---|---|---|---|---|---|---|
| base | 204 | 97 | 23 | 49 | 0.823 | 0.819 (n=872) | 0.690 (n=87) |
| cut | 204 | 97 | 23 | 36 | 0.823 | 0.819 | 0.690 |
Bars: complete >= 219 FAIL (204, unchanged); fused <= 49 PASS (36); ordinary hit rate >= 0.82: 0.819, identical to base
(the bar was set from the rounded base value; unchanged, not a loss); wall PASS (seconds + 90 s panel).
Verdict: FAIL as registered (bar 1). Tied-seed closed, as the row says.
Kept: the cut lowers fused families 49 -> 36 at zero cost on every hit rate, but it is annotation-guided (truth gene
spans), so it belongs with the r1495 guided f1units option, not the de novo default; not shipped here.
Follow-up rows (base mode, from the prior diagnostic): single-exon strand placeholder (RUSTLE_READ_STRAND opt-in; 30 misses)
and the PRAMEF8-class drop of a 2-read canonical chain (cause unidentified; 34 members with a >= 2-read chain and no model).

### Post hoc note (after the numbers): step A was a plumbing null, not a method null
`--tied-seed` is consulted only on the family-detection path (`tied_seed_skeletons` at family.rs ~30046, TSFAM families
appended at ~30900); `--assemble-only` skips that stage, so the region runs with and without the flag are the same binary
path (RUSTLE_TIED_SEED_DEBUG writes nothing). This is the r980 trap again. Wiring it into the assembly product was not
pursued: the measured ceiling is 3 members. "0/3 gained" above reads as "not executable from the assembly product".
