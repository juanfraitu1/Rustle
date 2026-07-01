# SUN identifiability catalog — the single-position private-allele witness (per-copy, O2)

Results note for the **Singly Unique Nucleotide (SUN)** layer of the copy-assignment identifiability
framework. Formal treatment, proofs, and machine-check statements live in `bench/THEORY.md` **§5·SUN**;
this note reports the genome-wide per-copy catalog and named examples. This is an **O2** object
(within-family, per-copy *copy-assignment* identifiability), complementary to the four O1
family-**definition** objects (segdup / DNA / RNA / protein) already shipped.

## What a SUN is (and the one distinction that is the whole point)

A **SUN** (Sudmant et al. 2010, *Science* 330:641, doi:10.1126/science.1197005; PMID 21030649) is a
**single column** where **one** copy's base is **unique** among all copies of the family — a
private/singleton allele at that one position. In Sudmant's pipeline (→ QuicK-mer2 / fastCN, Eichler lab)
short-read depth over the genome's 4.1M SUNs yields **paralog-specific copy number** (parCN). We import it
as the concrete, checkable, **single-read** witness for per-copy identifiability.

The load-bearing distinction, versus our existing unit:

- **PSV** (`psv_graph_genomewide.py`): a column with ≥2 distinct copy-alleles. A copy is "resolvable" iff
  its **full hap-vector** across all PSV columns is unique (a singleton group).
- **SUN**: a column that additionally isolates **one** copy against **all** the rest.
- Therefore **{copies with a SUN} ⊊ {copies with a unique hap-vector}**. A copy can be hap-vector-unique
  with **no** single-position SUN — its uniqueness is a *combination* of otherwise-shared alleles, so
  resolving it needs a read spanning **≥2 PSVs**, and it is exactly the copy a **recombinant** read can
  spoof (the K≥3 obstruction, THEORY.md §5). SUN is single-read deterministic; hap-vector-unique-only is
  co-observation-dependent.

On RNA this is operational: **a single read covering one SUN column assigns deterministically to that
copy at the per-read gate** — no phasing, no linkage, `N(r)={c}`, `|N(r)|=1`, and the read is never
misassigned to another *true* copy.

**Scope of "unique" (honest).** Uniqueness here is **within-family and conditional on catalog completeness**:
unique among the family's *enumerated* validated copies, not genome-wide (Sudmant's SUNs are genome-wide-unique;
a missing paralog could share the allele and demote a called SUN). Columns come from an all-to-longest-copy
(**star**) projection, not a true MSA, recording only aligned substitutions — this only ever **under-counts** SUNs
(drops gaps/indels, themselves good paralog markers).

**One thing a SUN does NOT buy.** SUN immunity is a **per-read gate** property, *not* cover-level immunity of the
copy. Even a SUN-rich copy can be dissolved by an alternative minimum cover in the NP-hard `RECOVER`/MCC problem:
the theory's canonical K=3 witness has a second minimum cover that dissolves a copy holding SUNs at 2 of 3 columns
(THEORY.md §5, `sun_theory_check.py::S3_cover`). That is exactly why the shipped pipeline runs the per-read gate
(Thm 4), which inherits SUN determinism, and not cover-recovery.

## The 3-tier per-copy ladder (a strict refinement of psv_graph's family K)

Per copy `c`, from copy-only asm20 alleles (no reads, no SEDEF): let `group(c)` be the copies sharing
`c`'s full hap-vector, and `SUN(c)` its private columns.

| Tier | Predicate | Meaning |
|------|-----------|---------|
| **1 — SUN-identifiable** | `|group(c)|==1` **and** `SUN(c)≠∅` | single-read **gate-deterministic** (a read over the SUN column is never misassigned to another *true* copy; `|N(r)|=1`, Thm 4/7). **Per-read** immunity — **not** cover-level: even a SUN-rich copy can be dissolved by an alternative minimum cover (§5, `S3_cover`) |
| **2 — hap-vector-unique-only** | `|group(c)|==1` **and** `SUN(c)=∅` | uniqueness is combination-based; **no single read pins it** — needs a read co-observing ≥2 PSVs (linkage); at the per-read gate a single-column read is ambiguous. This is the K≥3 recombinant-spoof regime for read assembly |
| **3 — K-frontier / unresolvable** | `|group(c)|≥2` | shares full hap-vector with another copy; `NM:i:0` collapse; Strong Sep fails on that pair; gate certifies `min_p=1` (tied) |

`{psv_graph resolvable copies} = {singleton groups} = Tier 1 ⊎ Tier 2`;
`{psv_graph frontier} = {non-singleton groups} = Tier 3`. The ladder **splits** psv_graph's family-level
"resolvable" set into the single-read core (Tier 1) vs the co-observation-dependent shell (Tier 2).

## Genome-wide catalog (154 GGO validated multi-copy families / 412 copies)

Source: `bench/sun_identifiability.py` → `bench/sun_identifiability.{tsv,json}` (copy-only asm20
aligned-pairs; drops psv_graph's `MIN_READS_PSV≥3` read-support gate). Determinism: byte-identical across
re-runs (single-threaded asm20 on a handful of short copy sequences; no RNG, no timestamps).

| Tier | Copies | % |
|------|-------:|----:|
| **1 — SUN-identifiable** (single-read gate-deterministic) | **338** | **82.0%** |
| **2 — hap-vector-unique-only** (no single read pins it; needs ≥2-PSV co-observation) | **1** | 0.2% |
| **3 — frontier / collapsed** (`NM:i:0`) | **73** | 17.7% |
| **unique-hap copies** (Tier 1 + Tier 2) | 339 | — |

**SUN-identifiable (338) ⊊ unique-hap (339)** — the strict subset is real and machine-provable in the
abstract (THEORY.md §5, witness S4), but on this substrate it is witnessed by **exactly 1** Tier-2 copy.
Honest headline: on gorilla, essentially every resolvable copy earns its uniqueness through **≥1
single-position private allele**, so **82% of copies are single-read gate-deterministic** (a read over the private
column is never misassigned to another true copy) and the K≥3 Tier-2 regime — where no single read pins the copy
and ≥2-PSV co-observation is needed — is empirically **rare** (1/339).

Per-family: **135/154** carry ≥1 SUN-identifiable copy; **132/154** are fully single-read taggable (every
copy Tier 1); **19/154** are all-Tier-3 frontier (`135 + 19 = 154`; the lone Tier-2 family is inside the
135). Copy-only PSV set is a **superset** of the read-gated psv_graph catalog, so 135 SUN families ≥ the
126 read-supported unique-hap families of `psv_graph_genomewide.json` — identifiability is a property of
the **copies**, observability of the **reads**.

psv_graph class cross-reference over the same 154 families: **124** `fully_resolvable`, **3** `partial`,
**27** `no_psv` (of which 18 have zero copy-only PSVs).

## On-real-data Strong-Separation machine-check (all green) — a self-consistency check

`bench/sun_identifiability.py::strong_sep_witness` recomputes, from the raw per-column allele dicts (**not** from
the `has_sun` flag), the single-position witness set `W(c) = { p : ∀c'≠c, (c)_p ≠ (c')_p } = ⋂_{c'≠c} D_{c,c'}`
for all 412 copies:

- `n_sun_identifiable_copies = 338` **=** `n_copies_with_single_position_strongsep_witness = 338`
  (`SUN_equals_witness = true`); `SUN_witness_equivalence_violations = 0`, `Nr_equals_1_violations = 0`,
  `SUN_implies_unique_hap_violations = 0` → **`all_green = true`**.

**Honest framing (not independent).** The witness recompute uses `all((c)_p ≠ (c')_p)` while `has_sun` uses
`count((c)_p) == 1` — these are the **same predicate** (copy `c` contributes exactly one occurrence to a column),
so `SUN_equals_witness` is true **by construction**. It is a valuable **no-coding-bug / self-consistency** check,
**not** independent corroboration of SUN. Genuine independence lives in `bench/sun_theory_check.py` (abstract,
exhaustive) and in re-deriving the tier counts from raw sequence: **S1** SUN⇒unique-hap 0/1,252,380 SUN-copies;
**S2** per-read gate immunity 0/6,675,294 reads; **S3** canonical K≥3 witness (factual, neutral); **S3_cover** the
**load-bearing counterexample** — a SUN-rich copy dissolves in an alternative minimum cover, so cover-level
copy-immunity is FALSE while per-read gate immunity still holds; **S4** NOT-iff both directions (incl. a K=4
all-no-SUN Strong-Separated family). **Honest scope**: SUN is a single-position, single-read **sufficient** witness
for the **per-read** side — **not** equivalent to family-level Strong Separation (all-pairs + coverage) in either
direction, and **not** a cover-survival guarantee for the copy.

## K-refinement: the one family where K over-counts (the strict-subset witness)

Exactly **one** family is `fully_resolvable` by psv_graph's `K` yet harbors a non-single-read-taggable copy:

- **family 42** — an 8-copy tandem block (`LOC129529768`/`LOC129529765..775` + `LOC101132445`,
  `NC_073247.2:~59.7 Mb`, 49 copy-only PSVs), psv_graph `K=8=n_copies`, `fully_resolvable`. SUN lens:
  **t1/t2/t3 = 7/1/0**. The Tier-2 copy is **`copy4 = LOC129529774`** (`NC_073247.2:59,772,963-59,781,277`):
  `group_size=1` (unique full hap-vector) but **`n_sun=0`, `n_witness=0`** — at *every* PSV column its
  allele is shared with some sibling, so it has no single-position private allele. **No single read pins it**;
  resolving it needs a read co-observing ≥2 PSVs (linkage), the exact fragment a K≥3 recombinant read can spoof
  at read assembly. The other 7 copies are all Tier-1. So `K` over-counts the single-read gate-taggable copies by
  exactly **1** genome-wide (`distinct_hap_vectors=339`, `single_read_taggable=338`).

## Named examples

**SUN-rich** (every copy Tier-1, dense private alleles):

| Family | Gene(s) | Copies | PSVs | Tiers |
|--------|---------|-------:|-----:|-------|
| 12 | **SMG1** cluster (SMG1 backbone `n_sun=652`) | 5 | **1479** | 5·T1 |
| 4 | LRPAP1-block (`LOC129523503` …) | 7 | 834 | 7·T1 |
| 38 | `LOC129530205`-block | 8 | 662 | 8·T1 |
| 8 | **DAPK1** + `LOC115935778` (both `n_sun=610`, fully private) | 2 | 610 | 2·T1 |
| 18 | `LOC134757299`-block | 9 | 401 | 9·T1 |
| 39 | **RABL2** cluster: **RABL2A** (`n_sun=25`) + **RABL2B** (`n_sun=153`) + 3 LOC paralogs | 5 | 357 | 5·T1 |
| 0 | **GSTM2** cluster (3×GSTM2 + LOC paralogs) | 7 | 246 | 7·T1 |

RABL2 is the clean named positive: fully single-read taggable, `K` does **not** over-count it. (APOBEC /
DAZ / RFPL4A are *absent* as validated ≥2-copy families in this substrate — no K-over-count case exists to
show for them, because they never formed a family here, not because they were clean.)

**The sole Tier-2 copy** (the strict-subset witness): **family 42 / `copy4 = LOC129529774`** (above).

**K-frontier no-SUN** (0 copy-only PSV, all copies collapse; genuinely RNA-unresolvable, `NM:i:0`):

| Family | Gene(s) | Copies | PSVs | psv_graph |
|--------|---------|-------:|-----:|-----------|
| 34 | `LOC115930538` | 8 | 0 | `no_psv`* |
| 1 | RGPD8-block (`LOC109025447` …) | 7 | 0 | `no_psv` |
| 7 | `LOC129534585` | 7 | 0 | `no_psv` |
| 22 | **ANKRD18A/ANKRD18B** block (+`LOC115933254`/`LOC101142783`/…) | 6 | 0 | `partial` → copy-only fully collapses |
| 21 | `LOC101141440` | 5 | 0 | `no_psv` |
| 14 | `LOC115930164` | 4 | 0 | `no_psv` |

## Honest caveats

1. **Structural, not observed.** This catalog is copy-only (asm20 self-alignment of the reference copies):
   it reports **potential** / upper-bound identifiability (what the reference copies permit). The RNA-
   **achievable** subset is SUNs that are **both read-covered AND exonic** — intronic SUNs are
   unwitnessable by RNA, and Sudmant's SUNs are genomic DNA-parCN markers. The observed/exonic intersection
   is a noted extension, not yet in this table (the per-copy JSON carries `n_sun` structural only).
2. **Substitution-only.** `column_alleles` records only aligned M/=/X bases; a copy-private **indel** is an
   excellent paralog marker but is **not** counted here (conservative → SUN undercounts identifiability).
3. **Divergent-copy alignment gotcha.** A copy too divergent to asm20-align has no columns, collapsing the
   *whole* family to 0 PSVs / all-Tier-3 even though that copy is trivially identifiable. **3 families
   (34, 75, 175) have `n_aligned < n_copies`** and should route to the **O4 divergent / reference-absent**
   path, **not** be read as pure `NM:i:0` collapse — e.g. family 34's all-Tier-3 label (6/8 aligned) is
   partly an alignment artifact. `n_aligned` is emitted per family precisely to flag this.
4. **The Tier-1/Tier-2 gap is empirically near-vacuous on gorilla** (1/339). The distinction is
   theoretically real and machine-provable in the abstract (S4), but only 1 real GGO copy realizes it — a
   result worth stating plainly rather than dressing up.
5. **Per-read, not cover-level.** "SUN-identifiable / gate-deterministic" is a **per-read gate** guarantee
   (`|N(r)|=1`, no misassignment to another true copy). It is **not** immunity of the *copy* in the NP-hard
   `RECOVER`/MCC cover: even a SUN-rich copy can be dissolved by an alternative minimum cover (§5, `S3_cover`).
   The earlier "recombination-immune copy / unspoofable / Tier-2→Tier-1 immune" phrasing overclaimed this and
   has been retracted; the tier numbers are unaffected.

## Files / reproduction

- `bench/sun_identifiability.py` — genome-wide per-copy catalog + a **self-consistency** (no-coding-bug)
  on-real-data Strong-Sep witness recompute (copy-only asm20, no reads / no SEDEF; ~14 s, deterministic).
- `bench/sun_identifiability.tsv` — per-copy: family, copy, chrom, start, end, gene, n_psv, n_sun, tier,
  hap_unique, group_size.
- `bench/sun_identifiability.json` — genome-wide tier breakdown, K-refinement, examples, `machine_check`.
- `bench/sun_theory_check.py` — exhaustive abstract checks S1, S2, S3, **S3_cover** (the cover-level
  counterexample), S4 (`ALL_GREEN=True`).
- `bench/THEORY.md` **§5·SUN** — formal definition, lemma, NOT-iff relation to Strong Separation, the
  per-read-vs-cover recombination link (with counterexample), and Sudmant 2010 provenance (References).
