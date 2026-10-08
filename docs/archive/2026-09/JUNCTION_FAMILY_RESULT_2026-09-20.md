# A family definition from splice junctions — it exists, it works, and it cannot be the definition

Run 2026-09-20 against `docs/PREREG_junction_family_2026-09-20.md` (committed `ea853f00` before any
junction edge was built). Tool: **`bench/junction_family_edges.py`**. Held-out chr2/chr8/chr10,
Soto S1C truth, `bench/heldout_family_score.py` — identical to every other arm.

## The option, stated

A junction cannot be compared by coordinate across loci, so it is compared **through the pairwise
alignment**:

```
for a pair (A,B) with a PAF record carrying a CIGAR:
    project each junction of A into B's frame through the alignment
    a junction MATCHES if BOTH endpoints land within ±10 bp of a junction of B
    shared(A,B) = number of matched junctions      EDGE iff shared(A,B) >= k
```

This is **drift-tolerant by construction**, which is what r344's refutation of *strict* concordance
demanded: intron **length** drift does not matter because junctions are matched by projected position,
and intron **count** drift does not matter because the rule asks for *k shared*, never *all shared*.
It runs, it is deterministic, and it produces families.

## As a definition it fails, on both axes

| k | pooled F | sens | prec | nodes | coverage of spliced genes |
|---|---|---|---|---|---|
| **1** | **0.4513** | 0.441 | 0.541 | 369 | **6.5%** |
| 2 | 0.3915 | 0.378 | 0.469 | 252 | 4.5% |
| 3 | 0.3287 | 0.313 | 0.395 | 198 | 3.5% |
| 5 | 0.2187 | 0.212 | 0.246 | 93 | 1.6% |

Against the shipped rule on the same substrate — **F 0.7123 at 95.4% coverage** — every k is below the
pre-registered 0.60 floor **and** below the 80% coverage requirement. ⛔ **Verdict: NO.**

## Why, exactly — a two-stage loss

Taking the **shipped graph's own edges** (known homologous pairs) and asking where the junction rule
drops them:

| | pooled |
|---|---|
| shipped edges | 3,785 |
| …with **both** genes spliced at all | 1,841 (**48.6%**) |
| …of those, sharing ≥ 1 projectable junction | 934 (**50.7%**) |
| **net surviving a ≥1-shared-junction rule** | **934 / 3,785 = 24.7%** |

**Half the family relationships have an intronless member**, and half of what remains has junctions that
do not project onto each other even at ±10 bp. A junction rule therefore reaches at most a quarter of the
edges the shipped rule uses — which is the 6.5% node coverage, restated.

⭐ This is r344's *"paralog intron lengths and counts DRIFT after duplication"* re-derived with a rule
built specifically to tolerate that drift. Tolerance does not rescue it, because **the dominant loss is
not drift — it is members with no introns to compare at all** (19.5% of genes on these chromosomes have
zero junctions, and pseudogene fragments are exactly the population §6t7–§6u0 is about).

## Where it is worth something: as a conjunct, weakly

Splitting the shipped edges by whether they share junctions:

| shipped edges | TRUE | total | precision |
|---|---|---|---|
| share ≥ 2 junctions | 50 | 53 | **0.943** |
| share 1 junction | 3 | 4 | 0.750 |
| share no junction | 158 | 178 | 0.888 |
| all | 211 | 235 | 0.898 |

⭐ Edges sharing ≥ 2 junctions are the most precise bucket, **+0.055 over the no-junction edges**.
⚠ But it is **53 of 235 labelled edges**, and 50/53 vs 158/178 is a 5.5-point difference on small
numbers — well inside what this session has repeatedly seen evaporate (§6t6's +0.0133 was two families).
Treat it as a **flag worth reporting on an edge**, not as a conjunct to admit or reject on, and confirm
it on more chromosomes before quoting the number.

## The honest summary

**The option exists and is now implemented.** It is a clean combinatorial object — shared junctions
counted through an alignment, one tolerance parameter, no weighting — which is the shape the advisor
prefers. It is simply scoped to a minority: it can define families **only among members that are both
spliced and have non-drifted junctions**, which is a quarter of the relationships the current rule
carries. It is a *sub-definition*, not a replacement, and the sub-population it covers is precisely the
one the existing rule already handles well.

## Limits

- ±10 bp tolerance, fixed in advance; a much looser tolerance would admit more but stops being a
  junction identity.
- Annotated junctions from the max-exonic transcript per gene; read-derived junctions on these loci
  would change the census and are the obvious follow-up (§6o8's node-construction priority).
- 27 Soto families / 235 labelled edges over three chromosomes.
