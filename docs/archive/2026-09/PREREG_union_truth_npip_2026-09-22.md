# Pre-registration — does the guided arm's lead survive a UNION truth, and how much rides on SEDEF?

**Written 2026-09-22, §6x1, before any arm is re-scored.** User: *"re-score against a union truth ... also
can we ensure we are not over-relying on SEDEF results? like great if we have them we could even replicate
most of Soto, but if we dont that does not weaken anything."*

## Why — the defect §6x0 left behind

`bench/mode_family_score.py:121` intersects every predicted cluster with the truth universe
(`gs &= set().union(*truth.values())`). A predicted member that is a real paralog but carries no Soto
label is therefore **deleted from the numerator AND the denominator** — it is not a false positive, it is
invisible. So §6x0's precision is **inflated** (register 770, universe-conditioned-on-comparator), not
"kneecapped" as first suspected. Measured on the de novo chr16 NPIP cluster: of 26 members, **1** is a
cross-family FP (PDXDC1, the PKD1P6-NPIPP1 readthrough locus's name) and **9 unlabeled members align at
100% of their own length** to the family. The question is whether the arm RANKING survives counting them.

## The three truths — fixed now

- **U0 = Soto NPIP**, `bench/soto/soto_famCN_S1C.tsv`, families containing an NPIP-named chr16 gene:
  **24 members / 6 family IDs**. The §6x0 comparator, unchanged.
- **U1 = U0 ∪ chr16 genes whose RefSeq `Name` contains NPIP**: **+4** (`NPIPA6`, `NPIPA9`,
  `PKD1P3-NPIPA1`, `PKD1P4-NPIPA8`). ⚠**r902-adjacent and declared so**: that row voided symbol roots as a
  *genome-wide truth generator* (they manufacture MIR/RNA5S/ZNF/C16orf "families"). Here the family is
  already defined by Soto and RefSeq only adds curated members to it. Reported separately, never pooled
  into a headline, so a reader can discard U1 alone.
- **U2 = U1 ∪ alignment-admitted members**, prediction-independent: every chr16 annotated gene body is
  aligned with `minimap2 -c -N 50 -p 0.1 -x asm20` against the **U0 member bodies only** — no arm's
  clusters are consulted — and a gene joins if **≥95% of its OWN length** aligns to ≥1 U0 member.

⚠⚠**U2 needs a length floor and the floor is fixed BEFORE scoring, from the truth, not from the list.**
Without one, `CLN3` (253 bp), `UBL5P4` (236 bp) and `PAWRP2` (634 bp) join at 100% — short elements sitting
*inside* a duplicated block, which is register 913's "short genes become hubs" arriving through the truth
instead of the metric. Floor = **the shortest U0 member, 10,633 bp (`NPIPB8`)**: a family member is at
least as long as the shortest known member. No other value is tried.

## Arms — unchanged from §6x0, re-scored only

de novo (`dn16_fam3`) · guided (`chr16_guided`) · semi-guided `exonic` · semi-guided `whole`.
Same `bench/mode_family_score.py`, same one-to-one bipartite matching, sensitivity/precision/F/collapse.

## The bar — committed now

| outcome | verdict |
|---|---|
| the §6x0 ranking (guided > semi > de novo) holds under **both** U1 and U2 | ⭐ **ROBUST** — §6x0 stands as written |
| ranking holds under one and flips under the other | ⚠ **TRUTH-SENSITIVE** — §6x0's verdict must carry the caveat |
| guided loses its lead under both | ⛔ §6x0's verdict is an artifact of an incomplete truth |

**Predicted, before looking**: ⭐ ROBUST *on ranking*, with **every arm's precision falling** — the
guided arm proposes annotated bodies, so the unlabeled paralogs it finds were already annotated genes and
enter all three arms' numerators together. What should move most is **de novo**, which is the arm whose
extra members were being deleted.

## SEDEF — the second question, answered by audit not by score

Claim under test: **SEDEF is a nice-to-have, not a load-bearing dependency.** Falsified if any arm other
than semi-guided consumes an SD/SEDEF input. Checked by reading every SEDEF call site in `src/` and the
`params.tsv` each arm actually emitted.

I will not change the truths, the floor, the arms or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⚠ **TRUTH-SENSITIVE. Guided's lead is ROBUST; everything else in §6x0 is not.**

Human chr16, one-to-one bipartite, identical downstream in every arm. Only the truth changes.

| arm | U0 Soto (24) | U1 +NPIP-named (28) | U2 +aligned (37) | prec U0 → U2 | collapsed U0 → U2 |
|---|---|---|---|---|---|
| **guided** | **0.833** | **0.762** | **0.721** | 0.938 → 0.815 | 0 → 1 |
| de novo | 0.727 | 0.667 | 0.610 | 0.923 → 0.720 | 2 → 8 |
| semi-guided `exonic` | 0.750 | 0.684 | 0.593 | **1.000 → 0.800** | 2 → **11** |
| semi-guided `whole` | 0.750 | 0.565 | 0.522 | **1.000 → 0.514** | 1 → 7 |

⭐**The headline holds and STRENGTHENS.** Guided is best under all three truths and its margin over second
place *widens* with truth completeness: **+0.083 → +0.078 → +0.111**. Guided is also the only arm whose
precision degrades gently (−0.123 against de novo's −0.203 and semi-`whole`'s −0.486) — expected, since its
nodes *are* annotated genes, so it has the fewest unlabeled members for the universe intersection to delete.

⛔⛔**But two §6x0 claims are now refuted, both by the same mechanism.**

**1. "Semi-guided lands BETWEEN de novo and guided" — an artifact of the small truth.** It is 2nd under U0
and **last (both variants) under U2**. De novo moves the other way, last → 2nd.

**2. "Precision 1.000 on NPIP in both variants, the highest of any arm" — that number was the universe
intersection, not precision.** Under U2 it is **0.800 / 0.514**. Register 770 arriving exactly as written:
the semi-guided arm has the MOST unlabeled members, so it had the most deleted, so it scored the highest.
⚠**The arm that looks most precise under an incomplete truth is the arm with the most unscored predictions.**

**3. §6x0's own correction #1 is itself wrong.** It recorded that r817's collapse failure "did not
reproduce on this substrate" (2/1 collapsed). Under U2 semi-`exonic` collapses **11 of 34** truth genes
against guided's **1**. It did not reproduce because the truth did not yet contain the genes being
swallowed. ⭐**A collapse metric cannot see a collision between two genes when only one of them is in the
truth** — collapse counts are only interpretable at a stated truth completeness.

## What the union truth cost to build, and the trap it walked into

⚠**`CLN3` joined U2 at 100% alignment before the length floor was applied** — and the reason is a
**duplicate gene `Name` in the GFF** (2 on chr16: `CLN3`, `LOC102724181`). One `CLN3` record is 253 bp and
aligns 100%; the real gene is 16,713 bp and aligns **0%**. Keying by `Name` silently merged them and gave
the short element the long gene's length. Fixed by making every record `Name|chrom:start-end`.
⭐**Any truth built by joining an alignment to an annotation on gene NAME is exposed to this**; the
symptom is a biologically absurd member that survives a length filter.

The floor itself did the work it was pre-registered for: without it `CLN3` (253 bp), `UBL5P4` (236 bp) and
`PAWRP2` (634 bp) all join at 100% — register 913's "short genes become hubs", reaching the scorer through
the truth rather than through the metric.

## SEDEF — ⭐ **not load-bearing. Exactly one arm consumes it, and that arm is the one U2 demotes to last.**

Audited every call site in `src/` and every arm's emitted `params.tsv`:

| consumer | role | SEDEF-free path? |
|---|---|---|
| `gw_family_catalog --from-genome-sd` | the **semi-guided** node source | — it *is* the mode |
| `mcl_families --sedef` | optional post-MCL core refinement, gated on `--core-refine` | ✅ `--core-from-paf` derives the same SD pairs **from the input PAF** |
| `from_genome.rs:346,572`, `annotation_families.rs:1199` | hardcoded `GGO_sedef_final.bed` | ✅ all three are `#[ignore]`d real-data tests behind `std::fs::metadata(p).is_err() → skip` |

- **De novo and guided both emitted `core_refine=false sedef=<unset>`** — they consumed zero SEDEF.
- SEDEF appears in **no** assembly or O2 source file, and in **neither `REPRODUCE.md` nor `docs/DATA.md`** —
  it is not a declared input for reproducing anything.

⭐**So the answer to "what if we don't have SEDEF" is: the best arm is unaffected.** What is lost is the
semi-guided arm, which under the fuller truth scores **last**. SEDEF remains useful for what §6u5 used it
for — *agreeing* with Soto, who built their families on segmental duplications — but agreement with a
comparator is not the objective, and O1's definition does not route through it.

## Correction to my own prediction

I predicted ⭐ ROBUST with **de novo moving most, upward**, on the reasoning that de novo was the arm whose
extra members were being deleted. Wrong on the direction: **every arm's F falls monotonically with truth
completeness**, because a bigger truth adds genes to the denominator faster than any arm recovers them.
De novo did move most in *rank* (last → 2nd), just not in score. The prediction confused *recovering a
deleted true positive* (raises precision) with *adding an unrecovered truth member* (lowers sensitivity);
the second is larger because these arms miss more of the family than they mis-assign.

> **Scorer (2026-09-22 port):** every `bench/mode_family_score.py` number above is reproduced byte-for-byte by `target/release/family_score` (same flags; 732/732 parity, scipy's assignment tie-breaking included — r1045/r1046). The Python was retired in §6z2.
