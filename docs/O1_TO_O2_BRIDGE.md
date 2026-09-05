# O1 → O2: feeding annotation-defined families into copy assignment

**Status 2026-09-03: WORKING, including on dispersed families (§6dh fixed the blocker).**
Ledger §6da–§6dg. This document is the runnable procedure, its structural limits, and the one bug
that had to be fixed to make the correction leg usable.

---

## 1. What connects to what

| stage | object | implementation |
|---|---|---|
| **O1** | families = MCL clusters of annotated genes, **by sequence** | `bench/mcl_annotation.py` |
| bridge | clusters → `copies.tsv` + `copies.fa` | §3 below |
| **O2** | PSV assignment of reads among copies | `copy_assign --families --copies-fa` |
| (optional) | re-align hard reads to the family's copy-paths | `--vg-realign-correct` |

⭐**The non-linear reference the advisor asks for is `--vg-realign-correct`**: it re-threads low-MAPQ /
heavily-clipped / divergent reads through the family's OWN copy sequences instead of the linear genome,
and reassigns them among the existing copies. It is **default OFF** and must stay a stated arm.

---

## 2. Two structural limits — read these before blaming a run

### 2a. `--families` REFUSES every roster-changing leg
`copy_assign.rs:987-1000` hard-refuses `--absent-copies`, `--vg-realign`, `--iterative-prune`,
`--collapse-gate`, `--tied-seed`, `--recover-copies` when `--families` is supplied: with a supplied
catalog the copy set must remain exactly the supplied one.
⭐**`--vg-realign-correct` is deliberately NOT on that list** — it re-threads reads among the GIVEN
copies and never touches the roster. ⟹ **O2 can never silently become O3.** Keep it that way: the
admission leg (`--vg-realign`, which admits novel read pools as new copies) is a SEPARATE experiment
that must run without `--families`.

### 2b. ⭐ DISPERSED FAMILIES — WAS A BLOCKER, FIXED IN §6dh
`copy_assign` is region-scoped and binds the **whole family span** to a single `--region`/`--regions`
entry. Its own rationale: a cross-chromosome family is *"structurally unassignable by a region-scoped
binary and is refused, never truncated."*

⟹ A DISPERSED FAMILY ONCE HIT THE SAME WALL. MCLFAM1 (38 NPIP copies on NC_073242.2) spans **89.5 Mb**;
covering it loaded **254,726 primaries** to assign copies occupying a few hundred kb, and the run OOMed.
⭐**FIXED (§6dh)**: reads are now gathered from the union of the COPIES' neighbourhoods
(`COPY_READ_PAD = 50 kb`, merged and deduplicated) rather than the whole hull — a read overlapping no copy
can never be assigned to one, so nothing assignable is lost and every copy's reads are still loaded in
full. **89.5 Mb -> ~3.2 Mb.** The 38-copy family now completes: **5,278 assignments, unique-mapper
agreement 292/293 = 99.7%.**
⚠**PROCESSABLE IS NOT RESOLVED**: only **8 of 38** copies receive a decisive read and **84.4%** of that
family's assignments are TIED. Report both numbers.
⚠**`--vg-realign-correct` STILL DOES NOT SCALE** to 38 copies (full NW per read against every copy,
~900M DP cells/read); the completing run had it OFF. Small families only until it gets a prefilter.

⭐**THE UNDERLYING TENSION REMAINS WORTH KNOWING**: O1 defines families by homology, which is scale-free
and produces genuinely dispersed families; O2's architecture assumes locality. NPIP is dispersed — that is
the biology, not an inconvenience.
⛔**DO NOT split a dispersed family into local sub-families to fit a tool.** That reshapes O1's answer to
suit the instrument.

---

## 3. The bridge, exactly

From an MCL cluster, emit `gw_family_catalog`-format inputs:

**`copies.tsv`** — header `family_id copy_idx tid chrom start end n_exon strand n_reads exons`,
`exons` as comma-separated 0-based half-open `start-end` blocks.
Contract (all loud, none silent): every copy must (a) be named by those columns, (b) belong to a
**same-chromosome** family, (c) fall inside exactly one region, (d) have a sequence, (e) have ≥1
overlapping read. A violation aborts.

**`copies.fa`** — header **`fam|idx|chrom:start-end|strand|nexon=N`**, verbatim. The sequence is the
spliced product: exon blocks concatenated, reverse-complemented on `-`.
⭐Always pass `--copies-fa`: it makes O2 assign against **exactly the bytes O1 defined the family with**,
instead of a reconstruction that could differ.

⚠**COORDINATE CONVENTION**: `seedmode/allgenes.fa` headers carry **GFF 1-based coordinates verbatim**
(`NC_073241.2:31346-41669` is exactly GFF `31346 41669`), NOT 0-based. A `start-1` key joins **0/4,477**
and silently falls back — ⭐**always print the join rate; a no-op result is the signature of a failed join.**

### What reached O2 on the NPIP pilot
44 cluster members → **43** on the 3 BAM contigs → **41** with ≥1 read in a same-chromosome family of
≥2 copies → 2 families (38 + 3 copies). **The loss is structural — cross-contig membership and read
coverage — not a scoring choice.**

---

## 4. Result on MCLFAM2 (3 copies, 292 kb)

1,711 mapped reads → **1,242 assignments**; **498 assigned / 726 tied / 18 ambiguous**;
median **2,884 decisive PSVs** on assigned reads; ⭐**genome-wide unique-mapper agreement 68/68 = 100%.**

⚠**The majority are `tied`** — near-identical copies are not separable, and the abstain behaviour is
working as designed rather than forcing a call. Report `tied` alongside `assigned`, never precision alone.
⚠`assigned_copy` is written for EVERY row regardless of status (`copy_assign.rs:1622`); on a `tied` row
it is only the earliest-index tie-break. **Filter on `status == assigned` before interpreting it.**

---

## 5. ⛔ The bug that had to be fixed first (§6df/§6dg)

`--vg-realign-correct` was length-confounded and reassigned **362/362 = 100%** of the reads on a 721 bp
fragment copy, unconditionally.

- `aln_id` is `1 - hw_distance(q,t)/len(q)` with free end-gaps on the TARGET only ⟹ a target shorter than
  the read caps identity at `len(t)/len(q)`. For a 721 bp copy against 2,685 bp reads that is **0.269** —
  the copy was excluded **by arithmetic, not homology**.
- `n_decisive = (id_best - id_linear) * read_len` is a whole-read EDIT DISTANCE, then
  `min_p = (error_rate/3)^n_decisive` treats each as an independent allele observation ⟹ `n_decisive ≈ 754`,
  `min_p` underflows to 0, reassignment is unconditional.

**FIXED**: `aln_id_len_safe` aligns the shorter as query (the existing idiom at `bridge_detector.rs:363`);
`psv_decisive_count` counts PSV COLUMNS where the read's allele matches `best` and contradicts `linear`.

| | before | after |
|---|---|---|
| reads on the 721 bp copy reassigned away | **362/362 = 100.0%** | **9/362 = 2.5%** |
| reassigned / rejected | 626 / 907 | 294 / 1,239 |
| unique-mapper agreement | 68/68 | **68/68** |

⚠**Any family with copies of unequal length was affected**, and §6db showed fragments are common.

### Two traps this exposed, worth carrying forward
1. ⛔**A same-typed parameter swap has no compiler safety net.** Replacing `read_len: usize` with
   `n_decisive: usize` left a failed patch **compiling and wrong** — `cargo build` returned 0 on code
   that would have reassigned everything. Only the tests caught it.
2. ⛔**`path_obs_at` returns the READ's base, not the copy's.** A counter comparing two `path_obs_at`
   vectors compares the read to itself and always returns 0. Compare the read's allele to each COPY's base.

---

## 6. Runnable

```bash
copy_assign \
  --bam   npip3.bam \
  --fasta npip3_contigs.fa \
  --families  fam.copies.tsv \
  --copies-fa fam.copies.fa \
  --regions   fam.regions \
  --vg-realign-correct \
  --out o2_fam
```
⚠**Run `copy_assign` FOREGROUND, serial, small batches** (standing WSL2 rule — it has OOM-killed this
machine before). Outputs: `.assignments.tsv` (one row per MOLECULE) and `.vg_realign.tsv`
(**one row per ALIGNMENT RECORD** — 1,533 rows over 1,092 names on the pilot).
⛔**A `read_name` join between them is MANY-TO-ONE** and silently collapses up to 5 verdicts onto one
assignment. Deduplicate deliberately, and say which grain a number is on.
