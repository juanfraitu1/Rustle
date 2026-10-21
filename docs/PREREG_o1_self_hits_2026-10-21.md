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
