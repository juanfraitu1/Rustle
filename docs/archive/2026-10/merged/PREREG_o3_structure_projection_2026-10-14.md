# PREREG — O3 structure projection: held-out simulation test (2026-10-14)

**Written after development on seeds 11, 12, 13 and before any held-out world is built.** Design and the frozen rule:
`docs/superpowers/specs/2026-10-14-o3-structure-projection-design.md` §3 (reference implementation `bench/o3_structure/proto.py`,
md5 511562b804907e1f3cb54b008000dc68; worlds `bench/o3_structure/world.py`, md5 fd2bc142f13a18ee4acfb63f95d05982).

## Claim under test

`candidate_copies`, given the reads of a family with one expressed copy absent from the reference, outputs a candidate that carries
that copy's **exon structure**: every exon boundary of the true copy inside the candidate's span, and no other, within ±5 bp — for
copies that differ from their nearest reference paralog by SNPs (≤ 4%) or by a structural change (exon loss, gain, duplication,
shuffle, 5' truncation, resized intron), on either strand.

## P0 — implementation parity (gate, before the held-out worlds are built)

The Rust implementation in `candidate_copies` (product `<out>.structure.tsv`) must give the **same boundary list as the frozen
prototype for every flagged candidate** of the dev worlds s11, s12, s13 × E0, E2 (148 candidates). A difference is a defect of the
implementation, fixed in the implementation (never in the rule), and P0 re-run. The held-out worlds are built only after P0 passes.
Also required: `cargo test --release` fails only the known tests (the 3 `pipeline_driver_dispatch` tests and, on this laptop, the 10
missing-fixture tests), and the products that existed before (`candidates.tsv`, `clusters.tsv`, `contigs.fa`, `families.tsv`,
`reads.tsv`, `clusters.fa`, `nets.fa`) are byte-identical to the pre-change binary on the dev worlds.

## Held-out worlds

`world.py build OUT --seed S --strands` for **S = 21, 22, 23, 24, 25**; read variants **E0** and **E2**; `world.py run`
(minimap2 2.30, `candidate_copies` at its defaults). 27 families with an absent copy E per world (DIV@0.5/1/2/4/8% × 2, LOSS × 3,
GAIN × 3, DUP × 2, SHUFFLE × 2, TRUNC5 × 2, INTRON × 2, LOWEXPR 8/12/16 reads) and 4 NULL families.

## Unit and scoring (`proto.py WORLD VARIANT --score-rust`, the prototype's scorer applied to the Rust product)

Unit = a flagged candidate of an E family whose union's best true-transcript hit (`minimap2 -c --eqx -x map-hifi`) is that family's E
on `+` ("scored"). True boundaries = E's exon junctions in transcript order inside the union's aligned span (≥ 20 bp from its ends),
mapped to union coordinates through that alignment. A candidate is **structure-correct** iff every true boundary has a projected
boundary within ±5 bp and every projected boundary is within ±5 bp of a true one. Boundary recall / precision pooled over boundaries.

## Bars (pooled over 5 seeds × 2 variants)

| id | measure | bar | dev value |
|---|---|---|---|
| **S1** | structure-correct / scored, all scenarios except DIV@0.08 | ≥ 0.90 | 0.986 (all scenarios incl. 8%) |
| **S2** | structure-correct / scored, per structural scenario (LOSS, GAIN, DUP, SHUFFLE, TRUNC5, INTRON) | ≥ 0.80 each | 1.00 each except SHUFFLE 5/6 |
| **S3** | boundary recall and boundary precision (±5), all scenarios except DIV@0.08 | ≥ 0.97 each | — (≥ 0.99) |

DIV@0.08 is reported, not barred (8% is beyond the family-edge identity the stage's families are built at; dev 6/6).

## Reported beside (no bar)

Exact (±0) structure-correct rate; every scenario and strand separately; candidates not scored and why; detection (E flagged) per
scenario and expression; flags in NULL families; **reconstruction** — identity x coverage of the union to E (≥ 0.999 = exact), and the
same for **IsoCon 0.3.3** (`pixi run -e isocon IsoCon pipeline`, defaults, `--nr_cores 4`) given all reads of the family (A + B + E,
FASTA, as Amendment 47's targeted input): its best final candidate against E. This is the baseline for the read-correction step, which
gets its own pre-registration.

## Reading, fixed now

- **S1, S2, S3 pass** → the projection is the stage's structure output: `<out>.structure.tsv` and `<out>.contigs.gtf` are written by
  default (they add files and change none), `o3_export.py` uses them, and the claim becomes "structure-correct in X of Y simulated
  reference-absent copies". The next level (real reads, bed H, annotated erased copies) is pre-registered before it runs.
- **S1 or S3 fails** → the projection is reported as a diagnostic, not a structure claim; the failing scenarios are named. No rule
  constant is changed in response; a new rule needs a new development round on new dev seeds.
- **S2 fails for some scenarios only** → structure is claimed for the passing scenarios only; the others are listed as limits.

## P0 result (2026-10-14, before any held-out world was built) — PASS

- Rust `candidates::structure` + `candidate_copies` phase 4 (binary sha1 bf1f968f, working tree on 49e62be3); prototype md5 unchanged
  (511562b8…); scorer `bench/o3_structure/score.py` md5 0a4315ef…, comparator `isocon_cmp.py` md5 6d7c141e….
- **Parity: 148 / 148 flagged dev candidates give the prototype's boundary list exactly** (s11, s12, s13 × E0, E2;
  `score.py parity`). One implementation defect was found and fixed on the way, in the implementation: the first port named the
  re-aligned segments `<candidate>|<start>_<end>` and gave 147 / 148. minimap2 breaks equally good hits by a hash of the read name
  (its `--seed`), so a terminal stretch with two equal hits (dev s12 E2 `cand_f12_0`, 2 mismatches on either copy) landed on the other
  copy. The segment names are therefore part of the reference behaviour; the binary now names them `<start>_<end>` as the prototype
  does, in rounds that never repeat a name within one minimap2 call.
- The seven earlier products are byte-identical to the pre-change binary's on all six dev runs (`cmp`).
- `cargo test --release`: 1,295 passed; the 13 failures are the known ones (3 `pipeline_driver_dispatch`, 10 missing fixtures on
  this laptop); the 9 new `candidates::structure` unit tests pass.

# OUTCOME (2026-10-14) — S1, S2, S3 all PASS on the held-out worlds

`pixi run bash bench/o3_structure/heldout.sh ~/rustle_data/o3_structure/heldout` (seeds 21-25 `--strands`, E0 and E2; binary sha1
bf1f968f): 270 absent copies, 242 flagged candidates, all 242 scored (each union's best true transcript is its own E).

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| S1 | structure-correct (±5), all but DIV@0.08 | 210 / 222 = 0.946 | ≥ 0.90 | ✅ PASS |
| S2 | LOSS / GAIN / DUP / SHUFFLE / TRUNC5 / INTRON | 30/30, 30/30, 17/20, **16/20 = 0.800**, 20/20, 19/20 | ≥ 0.80 each | ✅ PASS (SHUFFLE exactly at the bar) |
| S3 | boundary recall / precision (±5), all but DIV@0.08 | 1,553 / 1,561 = 0.995 / 1,553 / 1,561 = 0.995 | ≥ 0.97 each | ✅ PASS |

**Reported beside.**
- Exact (±0) structure: 185 / 222 = 0.833; GAIN 14/30 and SHUFFLE 10/20 exact (their novel-segment edges sit within a few bp), every
  other scenario at or near its ±5 count. DIV@0.08: 18 / 20 structure-correct. Strand: + 110/115, − 100/107.
- **The 14 misses, named:** 10 carry exactly one extra boundary and miss none (DUP s21 E2 f17, s24 f16 ×2; INTRON s22 E0 f22; DIV@0.02
  s23 E2 f04, s25 E0 f04; LOWEXPR s25 f25 ×2; DIV@0.08 s22 f09 ×2). The 4 SHUFFLE misses are all family f19 (s23, s24 × E0, E2), the
  layout with the absent copy on the strand OPPOSITE its home paralog (`strand_of`: k % 4 == 3): they miss 1-3 boundaries around the
  moved exons. No miss is a wrong placement of a recovered boundary.
- **Detection** (E flagged): DIV@0.005 0/20 (linked as an allele, by design: delta = 0.958%), DIV@0.01 12/20, every other scenario 20/20
  or 30/30, LOWEXPR (8-16 reads) 30/30. NULL families: no flag.
- **Reconstruction** (identity x coverage of the candidate against E ≥ 0.999), all 270 absent copies: **candidate_copies 242/270 = 0.896,
  IsoCon 256/270 = 0.948**. The whole difference is detection at ≤ 1%: DIV@0.005 0/20 vs IsoCon 17/20, DIV@0.01 12/20 vs 20/20.
  Wherever the copy is flagged, our candidate is exact (242/242); IsoCon is ≥ 0.999 in 239 of those 242 cases. IsoCon reports
  transcripts and has no test of reference absence; the 0.5% copies are what the allele cutoff delta excludes by design.

**Reading (as fixed).** The projection is the stage's structure output: `<out>.structure.tsv` and `<out>.contigs.gtf` are written by
default (new files; the earlier products are byte-identical). The claim: *in simulated families, `candidate_copies` gives a reference-
absent copy its exon structure within ±5 bp in 210 of 222 flagged cases (95%), across SNP divergence up to 4%, exon loss, gain,
duplication, shuffle, 5' truncation and resized introns, on both strands.* The named limit: shuffled exons on a copy oriented
opposite its home paralog. Not yet done: `o3_export.py` exports the unmapped-rescue chain's candidates and does not yet read
`contigs.gtf`; the next level (real reads, bed H, the erased copies' annotated structure) gets its own pre-registration on the main
machine. In simulation, reconstruction is not the bottleneck where a copy is detected (242/242 exact), unlike real reads (14/40 exact
on bed H), so the read-correction step must be developed against real-read error, not famsim's.
