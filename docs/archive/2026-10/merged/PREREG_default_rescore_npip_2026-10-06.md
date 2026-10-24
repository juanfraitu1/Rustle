# PREREG: strict FOUND, U2 and Compara CF153 for the current default pipeline at NPIP (2026-10-06)

Written and committed before any product of the HEAD binaries on this protocol exists (`DEF` and `PRE` below are generated after this commit).
Everything here is DEV: human A119b chr16 is the development block of NPIP and of every fusion rule. A PASS is "no regression against the registered
development values", not a validation.

## Question

The registered strict-FOUND numbers at NPIP (`docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, Amendment E: 10 / 6 / 2 found within own nodes 23 / 21 / 24, locus
level 19 / 18 / 16, for P / GOOD / ALL) were measured on pre-f1v2 arms from a frozen 2026-10-01 binary. The family scores of the current default
(f1v2 + `--min-cov-shorter 0.70`, U2 F .645, Compara chr16 F .667, CF153 F .812) were measured on 2026-09-30 at e163d955. No product has the current default's
loci scored for strict FOUND, and none of these numbers was produced by the HEAD binaries (the 2026-10-06 consolidation changed `src/`, and its real-data
slice is not recorded as run). Does the default at HEAD do no worse than the registered development values?

## Arms (human A119b, chr16, the 25 CAT/Liftoff NPIP copies of `copy_recovery_tools_cat/ann`)

- **DEF**, the default: HEAD `copy_assign --assemble-only` with the driver's `assemble` flags (strict junctions, shipped polish, `--gtf-tpm`,
  `--bridge-regroup f1v2`), seeding from the stored genome-wide best-AS table (secondaries >= 0.98 of the best AS, the default `GOOD` pool), on
  `--region chr16` in place of `--genome-wide` (the two are equal on the contig: the gates of `docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md`); then the driver's
  `families` stage with no override (`--min-cov-shorter 0.70`, most-reads representative, `--min-shared-exon-frac 0.60`).
- **PRE**, the control: the same HEAD binaries with `--bridge-regroup off` and `RUSTLE_MIN_COV_SHORTER=0` (the 2026-09-25 pipeline). It is what the
  registered `GOOD` arm was meant to be, rebuilt by the binaries under test.
- The frozen read-pool arms `P`, `GOOD`, `ALL` (`/mnt/linuxdisk/tmp/readpool_npip`) are used only by the gates.

## Instruments

- `bench/copy_support.py` (Amendment E, unchanged since 2026-10-04; sha1 in every log), arms given as `--loci NAME=loci.gff3,transcripts.gtf`
  (DEF: `PREFIX.families.gtf`, the GTF the families were built from; PRE: `PREFIX.gtf`) and `--nodes`.
- `bench/default_rescore/nodes.py`: the own-node rule of `bench/npip_read_pool/pagedata.py`, with the copy exons of the registered scoring (`npip_read_pool.json`) (an NPIP cluster = a cluster that holds a locus overlapping a
  same-strand NPIP copy; a copy has an own node iff such a locus overlaps it).
- `family_score` of the HEAD build through `bench/default_rescore/score_families.py` (clusters and truth restricted to chr16, `--chrom ALL --pairwise
  --per-family`, truths compara / u2 / soto of `families_gw/species/human`).
- `bench/default_rescore/verdict.py` computes the gates, the rules and the verdict from the products. Nothing is judged by hand.
- Denominator: 25 copies, fixed. The own-node count is conditioned on the prediction and is reported beside, never used as a denominator.

## Gates (instrument checks; a failed gate makes the verdict INVALID, the arms are then reported without a verdict)

- **G0** the HEAD scorer on the frozen arms reproduces `support_hsa_E.json` exactly for each arm: `old_overlap_in_npip_nodes` 23 / 21 / 24,
  `tc_found_in_npip_nodes` 10 / 6 / 2, `locus_tc_found` 19 / 18 / 16, `tc_found` 10 / 6 / 2.
- **G1** `nodes.py` on the frozen arms P and GOOD reproduces every stored `pagedata.json` node flag (25 copies x 2 arms; copy exons from `npip_read_pool.json`, see Amendment 1).
- **G3** the HEAD `family_score` on the stored e163d955 default clusters (`container_headroom/data/human_A119b/fam/chr16/D.fam.clusters.tsv`) reproduces the
  stored U2 (sens .588, prec .714, F .645) and Compara chr16 (sens .500, prec 1.000, F .667).
  (Before this registration I ran the HEAD `family_score` once on those stored clusters, U2 only, to learn its output format: it printed .588 / .714 / .645.
  That is a reading of stored data, not of DEF or PRE.)
- **G2** (a diagnostic, not a gate): PRE reproduces the registered `GOOD` values (own node 21, E-found in own nodes 6, locus level 18, E-found chr16-wide 6).
  If it does not, the pipeline has drifted since 2026-10-01 and the rules below still compare DEF with the registered values, with the drift reported.
- The scorer is run twice under different `PYTHONHASHSEED` values; the two summaries must be equal (the Python `set` trap of register 1044).

## Decision rules on DEF (fixed before the run)

- **R1** strict E-found within own nodes (`tc_found_in_npip_nodes`) >= 6 of 25 (the registered `GOOD` value).
- **R2** locus-level E-found (`locus_tc_found`) >= 18 of 25.
- **R3** U2 bipartite F (`family_score` pooled line) >= 0.645, at three decimals.
- **R4** NPIPB2 and NPIPB6 each have an own node in DEF (the registered `GOOD` arm had neither; the 09-29 flip separated them from GSPT1 and EIF3CL).
- **Verdict:** INVALID if G0, G1 or G3 fails; otherwise PASS iff R1 to R4 all hold; otherwise FAIL, naming the rules that fail. No rescue arm and no
  redefinition follows a FAIL: it is recorded as the finding (the representative problem of `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` would then be a
  property of the default, not of an older arm).

## Reported beside, no bar

DEF and PRE: own-node count, E-found in own nodes, locus level, chr16-wide E-found; the copies that change between PRE and DEF (found and own-node flags);
U2, Compara and Soto sens / prec / F and pairwise tp / truth / predicted pairs; the Compara CF153 row (hit / truth genes, sens, prec, F; stored value 13 of 19
genes, F .812); `cmp` of the HEAD products against the e163d955 products (assembled GTF, families GTF, clusters, loci, copy table; DEF against `D`, PRE against `B0`).
Liftoff recall is not recomputed (its scorer is not in the repo's runner).

## Not claimed

Anything about held-out NPIP substrates, gorilla, O2 or O3. Nothing about f1v2 or `--min-cov-shorter` as such: DEF against PRE is a descriptive
comparison on the development block. Strict FOUND is the registered rule of `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` (Amendment E); it depends on the
cap signal of this library, which gorilla lacks.

Runner: `bench/default_rescore/run.sh gates | asm DEF | fam DEF | asm PRE | fam PRE | score | verdict`.

## Amendment 1 (2026-10-06, written after the first gate run and before any DEF or PRE product exists)

The first run of G1 failed, on the frozen arms only (no DEF or PRE product was generated or read): `nodes.py` gave own-node counts 24 / 22 for P / GOOD
against the stored 23 / 21. Two causes, both in how the stored rule was run, not in the arms: (a) the stored rule took each copy's exons from the CAT gene's
exon blocks (`npip_read_pool.json`, `copies[].exons`), while `nodes.py` took the union of the copy's transcripts in `truth.hsa.gtf`; they differ at
PKD1P6-NPIPP1, where the stored flags are False in both arms; (b) the frozen `ALL` arm has two loci on one span (`DN_chr16_33611806_2` and `_3`), which
`pagedata.py`'s span join resolved silently and `nodes.py` refuses. Changes: `nodes.py` takes the copy exons from `npip_read_pool.json` (`--exons-json`; used for
DEF and PRE too, so the own-node rule is the registered one), and G1 is checked on P and GOOD, the two arms whose loci are unique. With the change G1 compares
50 flags and finds 0 mismatches. G0 still uses all three frozen arms (it reads the stored `pagedata.json`, not `nodes.py`). The rules R1 to R4 are unchanged.
