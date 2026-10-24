# Pre-registration: dropping diagonal self-hits from O1's family graph (2026-10-21)

**User (2026-10-21):** "lets try to fix the O1 defect". The defect is the one register 1368's audit found: copies grouped with an
overlapping or antisense copy of the same locus. Written before the fix is built or run.

## Root cause (systematic debugging, phase 1; evidence here, before any change)

- **The loci are genomic spans.** The de novo families stage (`mcl_families --from-gtf`) writes each locus as its GENOMIC SPAN and
  aligns the spans all-vs-all (`minimap2 -x asm20 -c -X`). `-X` only skips a record whose query and target names are identical.
- **Overlapping loci align to themselves.** Two different loci whose spans overlap carry the same genomic bases, so minimap2 aligns
  those bases onto themselves. In `human_A119b.fam.loci.paf` there are:
  - 253,893 records between different loci that lie on the main diagonal: same contig, + strand, query and target genomic
    intervals equal within 5 bp. 253,710 of them are at identity >= 0.999.
  - 37,451 records between overlapping intervals that lie off the diagonal (tandem period shifts). These are not self-hits.
- **Worked example:** the only record joining MCL1092's two loci maps chr2:192,437,263-192,525,700 onto itself at identity
  1.0000 over 88 kb. That record made the "family".
- **The existing guards do not catch it:**
  - `merge_overlapping_loci` (a node rule) is skipped under `fold_within_clusters` (default on), and merges exon-union overlaps
    only;
  - `--reject-overlapping` drops every record between span-overlapping loci. It is off by default, non-monotone (o1_ledger), and
    too broad: register 1010 refuted the span-level guard because overlapping de novo loci in duplicated regions are often real
    tandem copies.

## The fix (opt-in; default byte-identical)

`mcl_families --reject-self-hits` drops a PAF record when:
- query and target lie on the same contig;
- the record is on the + strand;
- the aligned genomic intervals coincide (start and end within 5 bp).

It is counted (`rejected_self_hits`) and reported in `params.tsv` and the log. An alignment of DNA onto itself is never evidence of
homology between two copies. Tandem copies keep their off-diagonal records, and inverted copies keep their − strand records.

## Evaluation

- **Catalogs:**
  - human_A119b: development; the defect was measured there.
  - **human_testis: held out.**
  - gorilla_KB3781: the second species.
- **Runs:** each is re-run twice with the current binary, from the same assembled GTF and the cached all-vs-all: a **control**
  without the flag and a **fix** with it. Comparisons are fix vs control (same binary), not vs the 10-10 catalogs.
- **Self-only copy:** a copy of the control catalog whose best same-family record at >= 0.90 identity is a diagonal self-hit,
  with no non-self same-family record at >= 0.90. This is the audit's definition, computed from the control's PAF.

**Bars, each in all three catalogs:**
- **B1 (the defect is removed):** >= 90% of the control's self-only copies, in the fix catalog, are either not in a family, or have
  a non-self same-family record at >= 0.90 identity.
- **B2 (no collateral):** of the control's copies SD-paired in their family (register 1368's test and tool, `bench/o1_sd`), >= 98%
  are still SD-paired in their family in the fix catalog.

**Reported, not bars:**
- register 1368's fraction of >= 0.90 copies SD-paired, control vs fix;
- families and copies;
- largest family size (the span-level guard was non-monotone);
- mainstay copies on CAT (human) or RefSeq (gorilla) in families.

**PASS** = B1 and B2 in all three catalogs. A PASS supports proposing `--reject-self-hits` as the families-stage default (the
user's decision); nothing is flipped by this pre-registration.

## Rules

- TDD for the Rust change: the unit test fails first.
- Build with `CARGO_TARGET_DIR` on `/mnt/linuxdisk`, `--release`, with the log captured to a file.
- Heavy runs in the foreground, one at a time, each under 10 minutes; a run that cannot fit gets an amendment.
- The metric-auditor checks the result. Never pool human and gorilla. Push only when asked.

## Results, v1 (code 72fdfbb6; evaluation e0e43642; `/mnt/linuxdisk/tmp/o1_selfhits/`): **FAIL** (gorilla B2)

Each control rerun reproduces its shipped catalog byte for byte (clusters, copies; the PAF replayed from cache).

| catalog | self-hit records dropped | B1 self-only copies fixed | B2 SD-paired copies kept | copies / families | >= 0.90 copies SD-paired |
|---|---|---|---|---|---|
| human_A119b | 253,893 | 549 / 561 (97.9%) ✓ | 5,857 / 5,974 (98.0%) ✓ | 9,907 / 2,150 → 9,495 / 1,913 | 72.1% → 76.8% |
| **human_testis (held out)** | 2,159 | 46 / 46 (100%) ✓ | 834 / 850 (98.1%) ✓ | 1,115 / 342 → 1,061 / 316 | 81.1% → 84.8% |
| gorilla_KB3781 | 6,645 | 69 / 70 (98.6%) ✓ | **687 / 705 (97.4%) ✗** | 918 / 326 → 834 / 283 | 82.4% → 90.4% |

- **Reported:**
  - largest family 260 → 262 (A119b) and unchanged elsewhere;
  - mainstay copies in families: 155 → 156 (A119b), 56 → 56 (testis), 2 → 2 (gorilla).

**Why B2 fails (post-hoc diagnosis).**
- The lost SD-paired copies are real pairs: in gorilla, all 18 come from 2-copy families at 0.977-1.000 identity, with a direct
  non-self record (e.g. MCL71: an inverted pair 49 kb apart, its − strand record at identity 1.000, 11.6 kb).
- Most of them carry a self-hit partner on the **same strand that shares exon bases**: a second model of the same transcription
  unit nested in it (human_testis 14 / 16, gorilla 13 / 18, A119b 35 pure and 59 mixed of 117).
- In the control, that self-hit kept the two models co-clustered, and `fold_within_clusters` merged them into one copy. Without
  it, host and nested model compete as separate nodes (gorilla's MCL71 becomes a 4-node bipartite graph), and MCL splits the
  pair apart.
- The self-only artefacts are the other kind: their self-hit partner is on the **opposite strand**, or on the same strand with
  **no shared exon bases** (a gene nested in another's intron). A pure "same strand, shared exons" partner is rare among them.
- `--pair-rule weight-blind` (post hoc) joins 19 two-node components but leaves B2 at 97.4%: the two-node rule is not the cause.

## Amendment 1 (2026-10-21, after v1's results; post hoc): v2 rule

- **v2** drops a diagonal self-hit only when the two loci are on **opposite strands**, or share **no exon bases**.
- A self-hit between same-strand loci that share exon bases is kept: they are two models of one transcription unit, and the fold
  step merges them into one copy.
- `mcl_families` passes the loci strands (`gene_strands` of `--gff`) to the graph builder. Without strands the rule is v1's.
- **Same bars (B1, B2) on the same three catalogs.** They are no longer held out: v2 was designed after seeing v1's results on all
  three, so a v2 PASS is a development result and needs a fresh catalog before any default flip.

## Results, v2 (Amendment 1; code e90d63f2; same evaluation): **PASS on all three catalogs (development, not held out)**

| catalog | self-hit records dropped | B1 self-only copies fixed | B2 SD-paired copies kept | copies / families | >= 0.90 copies SD-paired | mainstay copies |
|---|---|---|---|---|---|---|
| human_A119b | 170,084 | 548 / 561 (97.7%) ✓ | 5,931 / 5,974 (99.3%) ✓ | 9,907 / 2,150 → 9,355 / 1,872 | 72.1% → 77.4% | 155 → 152 |
| human_testis | 812 | 46 / 46 (100%) ✓ | 849 / 850 (99.9%) ✓ | 1,115 / 342 → 1,069 / 319 | 81.1% → 84.9% | 56 → 56 |
| gorilla_KB3781 | 3,357 | 69 / 70 (98.6%) ✓ | 704 / 705 (99.9%) ✓ | 918 / 326 → 843 / 290 | 82.4% → 91.5% | 2 → 2 |

- **Largest family:** 260 → 261 (A119b); unchanged elsewhere.
- **Collateral in A119b:** 43 SD-paired copies lose their family. 6 of them overlap mainstays.
  - Five are 2-exon models at NPIPB3 / B4 / B5 (MCL164, MCL112, MCL50) and one is at GOLGA8Q.
  - All are SD-paired at 0.996-0.999 non-self identity, and most have no self-hit of their own.
  - MCL re-routed once self-hit edges elsewhere in their component were gone: the non-monotone effect the o1_ledger recorded for
    edge-dropping rules.

**Reading.**
- v1, the pre-registered rule, FAILS: gorilla B2 is 97.4%.
- v2, designed after v1's results, removes 97.7-100% of the self-overlap artefacts and keeps 99.3-99.9% of the SD-paired copies
  in all three catalogs.
- Because v2 was designed on these catalogs, its PASS is a development result. A default flip needs a fresh catalog, e.g. a
  families run of another sample, pre-registered with v2 as is.

Audit: below, after Amendment 2. It holds with caveats; the A119b collateral reading above is corrected there.

## Amendment 2 (2026-10-21, before any run on it): confirmation on a fresh catalog, gorilla_OR6737

**User:** "ok lets ensure this works, what other dataset can we use?"

- **Fresh catalog: gorilla_OR6737.** A different gorilla individual and tissue (testis Iso-Seq on the same GGO reference). Its
  assembled GTF (`runs/gorilla_OR6737/gorilla_OR6737.families.gtf`, 73,965 transcripts) has never been through the families
  stage. Neither v1 nor v2 has seen it.
- **No fresh human catalog exists locally.** human_A119b and human_testis were used to design v2. A new human library (HG002
  Kinnex, the SCRI cerebellum) would first need a whole-genome alignment.
- **Rule: v2 exactly as committed** (code e90d63f2, `--reject-self-hits` with the loci strands); no change after this amendment.
- **Runs:** the families stage command, a control without the flag and a fix with it, the same binary, GTF and FASTA (GGO).
  The all-vs-all is computed once through `tools/mm2_shard.sh`, in resumable calls of at most 540 s, and cached; the fix
  replays it.
- **Bars:** unchanged.
  - **B1:** >= 90% of the control's self-only copies fixed.
  - **B2:** >= 98% of the control's SD-paired copies (SEDEF >= 0.90) still SD-paired.
  - **PASS** = both.
- **Reported:** the same descriptive figures as for the other catalogs.
- **Reading, fixed now:** a PASS here is the first held-out confirmation of v2, and supports proposing it as the families-stage
  default (the user's decision). A FAIL means v2 is fitted to the development catalogs.

## Audit of v1 and v2 (metric-auditor, 2026-10-21, on `audit_frozen_99782d47`): holds with caveats; readings corrected

Every number in both results tables reproduces to the digit in an independent recompute keyed on coordinates:
- diagonal records 253,893 (253,710 at >= 0.999), off-diagonal 37,451;
- records dropped, B1, B2, mainstay copies and largest family.

The code is correct:
- the loci FASTA is the forward strand, so the 5 bp slack is never used;
- exon blocks are 0-based half-open;
- the strand keys match;
- no record between different positions is dropped.

Corrections:
- **The A119b collateral is described wrongly above.**
  - 28 of the 43 "lost" SD-paired copies are representative switches: a model of the same locus that shares exon bases is
    still SD-paired in a fix family.
  - All 6 mainstay copies are among the 28 (e.g. DN_chr16_21697777_2 → DN_chr16_21699518_2, MCL136). In the control, those
    NPIPB loci sat in two families at once (MCL164 and MCL50).
  - Counted by locus, B2 is 5,959 / 5,974 (99.75%).
  - "Most have no self-hit of their own" is false: 37 of the 43, the 6 mainstays included, lost their own diagonal records.
  - The real collateral that no bar sees: 8 distant 2-copy families dissolved (identity 0.87-0.999). 4 of them are outside SD90
    (MCL153, 207, 963, 940).
- **v2 is not a pure deletion in A119b** (trap γ).
  - 43 fix families combine members of two or more control families (621 new co-member pairs).
  - Copy pairs in different families that share exon bases rise from 1,368 to 1,545.
  - 16 new same-family overlapping pairs appear (62 remain), and 4 new self-only copies (chr13).
  - Negligible in human_testis and gorilla_KB3781.
- **B1 is close to a tautology.** "Fixed" almost always means "left the families": 545 / 548, 45 / 46, 69 / 69. Only 3 / 1 / 0
  were repaired. By register 1368's definition the self-only copies go 497 → 0 (human), 70 → 0 (gorilla), 42 → 0 (testis), but
  485 and 69 of them left the families.
- **v1's FAIL is fragile.** It is 4 copies short under a gene_id join that the scorer, committed at 16:05, adopted after the v1
  runs. Counted by locus, v1 is 692-693 / 705 and would pass.
- **The ">= 0.90 SD-paired" gain is a shrinking denominator.** A119b: paired copies +9, >= 0.90 copies −523. Testis +2 / −42,
  gorilla 0 / −77. Quote counts, not percentages.
- **Default byte-identity was shown only for v1.** The controls ran on the v1 binary. The v2 build without the flag matched
  them by code reading only.
- **The fold does not depend on strand.** v2's rationale, "the fold merges them", would hold for opposite-strand pairs that
  share exon bases too. v2 drops those self-hits on purpose: a sense model and an antisense model of one locus are not two
  copies.

**What stands:**
- v1, as pre-registered, failed gorilla B2 by 4 copies under a gene_id join.
- v2, designed after seeing all three catalogs (including the testis held out for v1), passes them. That is a development
  result.

**What would settle it:**
- an evaluation fixed in advance: B2 counted by locus, plus guards on cross-family exon sharing, dissolved distant pairs and
  leftover self-overlap;
- run on gorilla_OR6737, plus a CAT-annotation arm.

A failure keeps the flag opt-in.

## Amendment 3 (2026-10-21, after the audit, before any OR6737 output exists): the evaluation, fixed in advance

Written while the OR6737 all-vs-all is still running. Its PAF and both catalogs do not exist yet. v2 stays as committed
(e90d63f2). Scorer: `bench/o1_sd/self_hits_eval.py` at the commit of this amendment.

**Bars:**
- **B1 (unchanged):** >= 90% of the control's self-only copies fixed.
  - It needs >= 10 control self-only copies; with fewer, the result is UNDERPOWERED, not PASS.
  - Reported: how many fixed copies left the families and how many were repaired.
- **B2, counted by locus** (replaces the gene_id join; audit findings 1 and 4): a control SD-paired copy is kept when it, or a
  fix copy sharing exon bases with it (either strand), is SD-paired in its fix family. >= 98%. The gene_id count is reported.

**Guards**, for what B1 and B2 cannot see; all must hold:
- **G1, one transcription unit split across families:** same-strand copy pairs in different families that share exon bases.
  - Fix <= control + max(5, 5% of the control's).
  - Opposite-strand pairs are reported only: separating sense from antisense is what the fix is for.
- **G2, distant pairs kept** (independent of SD tracks):
  - A control copy is distant-paired when its family holds a copy at a span-disjoint locus, joined to it by a non-self record
    at >= 0.90 identity.
  - The pair is kept when fix copies at its two loci (the copies themselves, or copies sharing exon bases with them) share a
    family.
  - >= 98% of distant-paired copies kept.
- **G3, the defect is gone from the output:** self-only copies under the fix's own families <= max(2, 10% of the control's
  self-only count).

**Reported:**
- SD-paired counts (not percentages);
- new and lost co-member pairs among copies present in both catalogs;
- fix families that combine members of two or more control families;
- dissolved 2-copy control families;
- the other figures as before.

**Checks before the verdict:**
- **Default unchanged:** the v2 binary without the flag must reproduce the development controls byte for byte (clusters and
  copies). Checked on the replays of gorilla_KB3781 and human_testis.
- **Development figures:** G1-G3 are computed on the three development catalogs and reported. They do not set the thresholds.
  - Known when this was written (from the audit), with no strand split seen: A119b's cross-family exon-sharing pairs
    1,368 → 1,545, 8 distant 2-copy families dissolved, and 4 new self-only copies.

**Verdict on OR6737:**
- PASS = B1, B2, G1, G2 and G3 all hold. A PASS supports proposing v2 as the families-stage default (the user's decision).
- Any FAIL keeps it opt-in.

## Results, Amendment 3 (scorer 7cda7ffd; frozen copy `audit_frozen_or6737_7cda7ffd`): **OR6737 (held out) PASS**; human_A119b (development) fails G1

- **The default is unchanged.** The v2 binary without the flag reproduces the KB3781 and testis controls byte for byte (clusters,
  copies, merged copies, PAF).
- **OR6737 control:** 1,607 copies in 498 families. The fix replays the same PAF from cache and drops 2,963 self-hit records.

| catalog | B1 self-only fixed (left / repaired) | B2 by locus (gene_id join) | G1 same-strand cross-family pairs (opposite strand) | G2 distant-paired copies kept | G3 self-only left | verdict |
|---|---|---|---|---|---|---|
| **gorilla_OR6737 (held out)** | 71 / 72 (71 / 0) ✓ | 1,272 / 1,272 (1,270) ✓ | 25 → 24 (36 → 36) ✓ | 1,364 / 1,365 ✓ | 0 ✓ | **PASS** |
| human_A119b (development) | 548 / 561 (545 / 3) ✓ | 5,959 / 5,974 (5,931) ✓ | **1,036 → 1,140 > 1,088 ✗** (332 → 405) | 7,560 / 7,598 ✓ | 10 ✓ | FAIL (G1) |
| human_testis (development) | 46 / 46 (45 / 1) ✓ | 849 / 850 (849) ✓ | 8 → 8 (4 → 4) ✓ | 961 / 962 ✓ | 0 ✓ | PASS |
| gorilla_KB3781 (development) | 69 / 70 (69 / 0) ✓ | 705 / 705 (704) ✓ | 9 → 8 (12 → 13) ✓ | 746 / 749 ✓ | 0 ✓ | PASS |

**Reported for OR6737:**
- copies 1,607 → 1,521; families 498 → 456; largest family 86 → 86;
- co-member pairs +2 / −6; 1 fix family merges control families;
- 43 two-copy families dissolved, 8 of them with members that are not self-only;
- >= 0.90 copies SD-paired: 1,159 / 1,383 → 1,160 / 1,293;
- mainstay copies 34 → 34, with the same gene_ids.

### Audit (metric-auditor, 2026-10-21, on `audit_frozen_or6737_7cda7ffd`): holds with caveats

Every OR6737 number reproduces to the digit without the scorer. The PASS survives every stricter count:
- B2 with a same-strand-only join: 1,270 / 1,272;
- G2 requiring every pair: 1,358 / 1,365 (by pair, 6,595 / 6,599);
- G1 counting folded models: 25 → 26 (limit 30).

**The timeline is clean.**
- Amendment 2 was committed at 16:34:29 and Amendment 3 at 17:39:47. At that point 49 of the 100 all-vs-all shards were done
  and no family output existed.
- The family outputs date from 18:25.
- The scorer is identical to 7cda7ffd.
- The binary was built at 16:19:44, and `src/` is unchanged since e90d63f2.

**Caveats:**
- **B1 and G3 hold by construction.**
  - All 72 self-only copies form 36 two-copy families. Each pairs opposite-strand models that overlap in span, share no exon
    bases, and are joined only by a diagonal record (4.7-87.7 kb). v2 drops every such record.
  - The deletion is consistent with the catalog's own edge rules: the pairs' other records are internal repeats at 0.83-0.894.
  - The fix removes more than the self-only definition counts. Same-family opposite-strand pairs with overlapping spans go
    47 → 0. That includes the 8 "not self-only" dissolved families (the same artefact, hidden from the definition by
    off-diagonal internal-repeat records) and MCL233.
  - All 90 lost copies are accounted for: 71 self-only, 16 from those 8 families, 1 from MCL233, and 2 antisense
    representative switches.
- **One real collateral loss, which G1 as written does not see.**
  - MCL70's locus (NC_073240.2:21.11-21.23 Mb) is split between MCL175 and MCL368.
  - Cause: v2's same-unit test is pairwise, but the fold is transitive. 17 of the 334 dropped same-strand records join loci
    that are one unit through a third model.
  - This is G2's one lost copy. Its inverted 27 kb pairing survives through the locus's other model.
- **OR6737 is a held-out sample, not a held-out genome.**
  - It is on the same GGO reference as KB3781, where v1 failed and v2 was shaped.
  - 550 / 1,607 of its copies share exon bases with KB3781 copies.
  - 14 / 72 of its self-only copies and 8 / 43 of its dissolved families sit at loci also seen in KB3781.
  - Its 10-10 GTF had never been through the families stage. Earlier assemblies of the same library had been (`rt_arms`
    09-27, `rg3_run` 09-28), but none was used for v1 or v2.
- **A119b's G1 failure is mostly relabelling, but not only that.** My first label, "all switches", was not a test. Checked by
  locus:
  - 289 of the 313 new pairs lie at loci already split across control families. 291 of the 313 are on chr13, 13.9-15.5 Mb.
  - 17 are new splits of loci the control held in one family.
  - 7 involve a locus with no same-strand control copy.
  - Counted by unit, families beyond the first go 711 → 742 (+4.4%). That would pass a 5% locus bar, but only post hoc.

**Reading** (the auditor's wording):
- v2 (e90d63f2) passes every Amendment 3 bar and guard on one held-out gorilla catalog, OR6737, which is on the same GGO
  reference as KB3781. The result was reproduced independently.
- B1 and G3 hold by construction. The collateral is one locus split across two families.
- There is no held-out human catalog. human_A119b fails G1 as pre-registered (1,036 → 1,140 > 1,088). The guided (annotation)
  arm is untested.
- As pre-registered, this supports proposing v2 as the default of the de novo families stage **for gorilla**. For human, and
  for catalogs built from an annotation, it is not shown.
- Nothing is flipped; the user decides.

**What would settle the rest:**
- a held-out human catalog, e.g. HG002 Kinnex after a whole-genome alignment;
- control and fix runs on a guided CAT arm, both under Amendment 3;
- for MCL70-type splits, a same-unit test that follows the fold's transitive components rather than pairwise exon sharing.
  That would be a v3, with its own pre-registration.

## Decision (user, 2026-10-21): `--reject-self-hits` is the default

**User:** "lets make reject self hits the default". The user decided this knowing the scope of the evidence above:
- confirmed on one held-out gorilla catalog;
- human_A119b fails G1 in development;
- annotation-built catalogs are untested.

**What changed:**
- `mcl_families` now ships `--reject-self-hits` ON, with v2's rule (strands from the loci GFF).
- `--no-reject-self-hits` restores the graph of every catalog built before the flip. The control arms above are now run with
  that flag.
- `--reject-self-hits` still parses and is a no-op.
- The library default (`GraphParams::default()`) stays OFF, as with `min_shared_exon_frac`. `params.tsv` records the
  effective value.

**Checked after the flip:**
- **gorilla_KB3781 and gorilla_OR6737, new default vs the fix arm:** clusters, copies and merged copies are byte-identical,
  with the same records dropped (3,357; 2,963).
- **The same catalogs, `--no-reject-self-hits` vs the controls:** byte-identical.
- **Tests:** 19 / 19 `mcl_families`, 1,106 library, 3 driver.

**Run without `--gff`:** when `--paf` is given without `--gff`, no strands or exon blocks are known. The rule then falls back to
v1, which drops every diagonal self-hit. Every families-stage run in this project uses `--from-gtf`, which writes and passes its
loci GFF.

**What changes downstream:**
- Every families catalog built from now on drops self-overlap families. On the catalogs scored here:
  - human_A119b: 9,907 → 9,355 copies;
  - human_testis: 1,115 → 1,069;
  - gorilla_KB3781: 918 → 843;
  - gorilla_OR6737: 1,607 → 1,521.
- The shipped catalogs in `rustle_figures/runs` are unchanged until they are rebuilt.
