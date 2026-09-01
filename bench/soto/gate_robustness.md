# Soto gate robustness — do the readthrough / mis-chain gates distort detection? (2026-07-23)

> ⚠ **PROVENANCE (2026-09-01).** The script that produced this table,
> `bench/soto/gate_sensitivity_sweep.sh`, was found to be a FALSE-PASS instrument (`o1_ledger.md` §6am):
> its binary path was dead, it truncated this result file *before* attempting any binary, and nothing
> aborted — so it could have scored stale files and appended a fabricated row. It has since been repaired
> (§6an) to abort above the truncation and to preserve the prior file.
>
> ⭐**THESE NUMBERS NEVERTHELESS STAND, on an internal positive control**: the four configs report **four
> DISTINCT copy counts — 467 / 376 / 483 / 513, a 137-copy spread** — each moving in the direction its
> gate predicts (`gates_off` fewest, `mischain_aggressive` most). A dead binary scores the SAME stale
> files for every config, so a spread of that size requires the catalog to have been genuinely rebuilt
> per config. The recall column is the weaker half: `soto_cache_score.py`'s `hit()` is a union over four
> legs, three of them static files independent of the rebuild, so the 81.2–82.3% band would move little
> even if one leg were empty. **Quote the copy swing as measured; treat the 1.1-point recall band as
> insensitive by construction rather than as evidence of stability.**


**The question (advisor):** the pipeline drops two kinds of spurious de-novo transcripts before it builds
families — *unspliced readthroughs* and *giant-intron mis-chains* — both of which can connect otherwise-distinct
copies. If those gates are load-bearing, the Soto recall number would be an artifact of the thresholds. So:
**does member detection move when we vary the gates?**

**Answer: no.** Across the full range — from both gates *off* to *maximally aggressive* — member recall stays in
a **1.1-point band (81.2 %–82.3 %)**. The gates change how many *copies* we report (376 → 513, a ~36 % swing) but
not which Soto *members* we recover. They are a copy-**separation** knob, not a recall knob.

Reproducible on the cached Soto-region subset, scored on the real `80_fams.chr.bed` (362 members, 83 families)
exactly like the [precision/recall audit](precision_recall_audit.md).

## The two gates (definitions from source, `denovo_pipeline.rs`)

Both run on assembled **transcripts, before `collapse_loci_span_aware`** — deliberately, because a readthrough is
the *longest* object in its span, so after collapse it becomes the representative of every locus it covers and
dropping it would delete the real copies with it (the DAZ 298 kb span absorbed DAZ2 into DAZ1's group). Gating
before collapse is what lets O1 (family catalogs) and O2 (assignment) agree on what a locus *is*.

| gate | fires on a transcript iff… | compiled default | env knob |
|---|---|---|---|
| **readthrough** (`retain_non_readthrough`) | it is **single-exon** AND ≥ *N* **distinct** junctions (each ≥ 2 primary reads) lie **entirely inside** its span — an unspliced pre-mRNA engulfing whole gene structures | *N* = 5 distinct | `RUSTLE_READTHROUGH_MIN_DISTINCT` |
| **mis-chain** (`retain_non_mischain`) | any of its **own** introns is **giant** (> *G* bp) yet carried by **< *R*** primary reads with that exact junction — a spurious splice minimap2 made by mis-aligning read ends across a large gene | *G* = 50 kb, *R* = 3 | `RUSTLE_MISCHAIN_GIANT_BP`, `RUSTLE_MISCHAIN_MIN_READS` |
| **master toggle** | — | both gates ON | `RUSTLE_KEEP_READTHROUGH=1` disables **both** |

Neither threshold was tuned to Soto: *N*=5 is read off 30 regions (every single-exon de-novo transcript engulfs
≥ 14 distinct junctions; no expressed intronless annotated gene of 260 exceeds 4). *G*=50 kb sits above the
largest intron in any real recovered paralog on this substrate (POTE, 48 kb).

## The sweep

Four configs, each a full per-chrom Soto recompute on the cache; thresholds injected by env (no recompile),
so the `default` row is byte-identical to a normal run. Driver: `bench/soto/gate_sensitivity_sweep.sh`.

| config | gate setting | copies | on-member % | **member recall** |
|---|---|---|---|---|
| `gates_off` | both gates OFF | 376 | 33 % | 294/362 = **81.2 %** |
| `default` | ON — distinct = 5, giant = 50 kb / min = 3 | 467 | 35 % | 295/362 = **81.5 %** |
| `readthrough_strict` | distinct = **3** (lower bar → drops more readthroughs) | 483 | 34 % | 298/362 = **82.3 %** |
| `mischain_aggressive` | giant = **20 kb**, min = **5** (drops more mis-chains) | 513 | 37 % | 296/362 = **81.8 %** |

*(rows ordered by increasing gate aggressiveness)*

## What it means

- **Recall is robust.** 81.2 % → 82.3 %, a **1.1-point band** spanning gates-fully-off to maximally-strict. The
  readthrough / mis-chain gates are **not** what sets the Soto recall figure — so the recall result cannot be
  dismissed as a threshold artifact.
- **The gates control copy separation, monotonically.** Copies climb with aggressiveness
  (376 → 467 → 483 → 513) exactly as the mechanism predicts: keeping readthroughs/mis-chains lets giant spurious
  spans bridge and collapse distinct loci into *fewer, dirtier* copies; gating them apart yields *more, cleaner*
  copies. This lever is **orthogonal** to member recall.
- **The extra copies are finer splits, not new biology.** As copies rise 467 → 513, the distinct members
  recovered stay flat (~296) while on-member % rises (35 → 37 %) — i.e. the aggressive config re-splits
  already-detected members rather than discovering new ones. That reads as mild **over-splitting** (a precision
  cost), and is a point *for* the default gates: they are recall-neutral and avoid it.

**Scope / honesty caveat.** This sweep measures **recall** robustness, not precision. It shows the gates do not
create or destroy member detection; it does **not** certify that the copy-count swing is all-real vs. partly-FP.
The precision of off-benchmark copies is characterized separately in
[`precision_recall_audit.md`](precision_recall_audit.md) (119/125 off-benchmark families are real non-Soto SDs;
~2 spurious). A dedicated precision check on the 467 → 513 swing would be its own analysis.

## Reproduce

```bash
bash bench/soto/gate_sensitivity_sweep.sh      # 4 configs, serial, ~3 h on 5 cores
cat /home/juanfra/winloci_scratch/soto_cache/gate_sensitivity.tsv
```

Each config clears and regenerates the per-chrom `*.copies.tsv` (cached mini-BAMs are reused), scores member
recall with `soto_cache_score.py`, and appends one row to `gate_sensitivity.tsv`.
