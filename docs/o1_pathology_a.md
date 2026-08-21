# Census pathology (a) decomposes — and its largest class was misattributed

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status 2026-08-19. Offline (T8), whole catalog, nothing through the shipped binary.**
Scratch: `/mnt/linuxdisk/home/juanfraitu/o1_gmult/patha.py`.

## 1. Why this was worth running

Today's read-based guards all died, and the reason converged: in the shipped catalog the 47
node-construction failures are dominated by **pathology (a), "one locus emitted as two"**, which is a
**coordinate/identity** signature, not a read signature. The census found it but nobody had checked
what merging the flagged pairs would do.

Four signatures, applied within each shipped family, over all 3,751 same-family copy pairs:

| signature | definition | pairs | families | census |
|---|---|---:|---:|---|
| `A_OVERLAP` | genomic intervals intersect | 35 | 35 = 7.09% | — |
| `A_ADJACENT` | same chrom, gap ≤ 10 kb | 50 | **28 = 5.67%** | **matches 28/494** |
| `A_IDENT` | rep alignment ≥ 0.99 id over ≥ 500 bp | 204 | 106 = 21.46% | broader than the census's |
| `A_SAMEGENE` | both best-overlap the **same GFF gene instance** | 29 | **20 = 4.05%** | **matches 20/494 exactly** |

Two signatures reproduce the census's numbers exactly, from an independent re-derivation.

## 2. ⭐ The largest class is NOT pathology (a)

**All 35 `A_OVERLAP` pairs are on OPPOSITE STRANDS — 35/35.** They are not one locus split in two.
They are **overlapping antisense / nested gene pairs**, a real and common arrangement, and the
`distinct_locus_reps` predicate is right to treat opposite-strand loci as distinct.

The question then becomes why they are in one *family* at all. They share the same DNA read in
opposite directions, so they align — **on the minus strand**:

| | minus-only edges | rate |
|---|---:|---:|
| the 35 overlapping antisense pairs | **33** | **0.9429** |
| all shipped edges | 56 / 2,727 | 0.0205 |

**A 46× enrichment.** These are reverse-complement artifacts, which is precisely the class the
**transcript-orientation guard** rejects. Merging them would be flatly wrong — it would fuse two
genuinely different genes. They should lose their *edge*, not their node.

### Corroboration for the orientation guard, from a new direction

The guard's evidence so far was the 74-pair FP arm (29/74 rejected, 4 lost edges of 9,032). This is
independent — it comes from **coordinates and strand**, not from the FP arm — and it shows:

* the guard's genome-wide reach is bounded: minus-only edges are **56/2,727 = 2.05%** of the catalog,
  so **that is its maximum cost**, against junction-crossing's 12.80% and the genome-anchored veto's
  3.67%. **It is the cheapest of the three guards by a factor of ~2–6.**
* **33 of those 56** are overlapping antisense pairs ⟹ **59% of what the guard removes carries no
  homology evidence.** ⚠ **Stated precisely (the earlier "provably artifactual" overstated it):** a
  gene and its antisense partner share the same DNA *by construction*, so a `-` alignment between their
  transcript-oriented reps is **entailed by the overlap** and is therefore not evidence of homology.
  That is not the same as proving the two are non-homologous — a duplicated region could in principle
  contain both. The edge is uninformative, not disproven.

## 3. What real pathology (a) costs to fix

Merging every flagged pair, judged on **family-level** outcomes only (**T7** — three prior node
changes passed on node metrics and failed end to end):

| signature | families touched | copies lost | families DISSOLVED (n < 2) |
|---|---:|---:|---:|
| `A_SAMEGENE` | 20 | 25 | **6** |
| `A_ADJACENT` | 28 | 50 | 10 |
| `A_OVERLAP` | 35 | 35 | 30 ⚠ *should not be merged at all* |
| `A_IDENT` | 106 | 148 | 75 ⚠ *too broad* |
| ANY | 133 | 201 | 86 = 17.4% of the catalog |

**`A_SAMEGENE` is the defensible detector.** It reproduces the census exactly, it is targeted (6
families dissolve, 1.21% of the catalog), and merging is the *right* operation for it: two nodes whose
best-overlapping gene is the **same gene instance** are one locus by construction.

⚠⚠ **But it reads the ANNOTATION**, and O1's standing position is *annotation as SEED, not DEFINITION*.
So `A_SAMEGENE` may be a **catalog QC flag**, never a membership condition — otherwise O1 stops being
annotation-independent, which is a larger loss than 20 families.

⚠ **`A_IDENT` must not be used.** At 106/494 = 21.46% it cannot separate a split locus from a **real
tandem duplicate at 99.9% identity** — the two present identically on identity alone. This was stated
before the numbers were computed and the numbers did not resolve it.

## 4. Net

Pathology (a) is smaller and better-behaved than "47 node-construction failures" suggests:

* **35 families** are misattributed — antisense overlap, an **edge** problem the orientation guard
  already covers, not a node problem;
* **20 families** (4.05%) are genuine same-gene-instance splits, fixable at a cost of 6 dissolved
  families, but only with the annotation;
* **28 families** (5.67%) are adjacent-locus candidates, unresolved;
* the identity-only route is unusable.

**Nothing here changes the definition.**

## 5. ⭐ CONFIRMED END-TO-END ON THE REAL BINARY (2026-08-20) — this is no longer T8

The prediction above was offline. The genome-wide catalog built with the guard **on by default** now
measures it directly:

| | families with an overlapping same-family pair | opposite-strand |
|---|---:|---:|
| 494 catalog, guard OFF | **35/494 = 0.0709** | **35/35** |
| **627 catalog, guard ON** | **4/627 = 0.0064** | 3/4 |

**A 91% reduction in the antisense-overlap class, measured by the shipped binary.** Every one of the
35 pre-guard cases was opposite-strand, which is what made the mechanism identifiable in the first
place. Together with the human negative panel (spurious E_r edges **28 → 3**,
`o1_false_merge_remeasured.md`) the guard's precision benefit is now measured on **two independent
substrates and the real binary**, not on the GGO FP arm it was derived from.
