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
