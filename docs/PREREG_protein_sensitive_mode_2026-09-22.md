# Pre-registration — a SENSITIVE MODE: protein edges from translated de novo locus representatives

**Written 2026-09-22, §6y3, before any arm is scored.** User: *"what if we add a sensitive mode that takes
the loci representatives, translates them to protein and performs the graph based on the protein
sequences"*.

## What is already known, and what is genuinely new here

**r906 (§6t0) already measured protein edges and they WORK** — from **annotated CDS**, §6ko's rule: on
held-out chr2/chr8/chr10, within-family pair coverage **80.0% → 93.8% (+13.8 pts)**, no-edge families
**halve 14.8% → 7.4%**, **0 cross-family protein edges on all three held-out chromosomes**, pooled
precision **0.919**. It was recorded ⚠PARTIAL and **NOT adopted**, blocked on three things: *its own
prereg, a §6bt-style false-merge measurement, and a decision on whether a protein-space edge belongs in a
definition called topological.* ⭐**The user's "sensitive mode" framing dissolves the third** — an opt-in
mode is not the definition. This file supplies the first, and the cross-family count supplies the second.

**r909 (§6t2) says protein edges "buy nothing"** because LCS-on-protein recovers blastp's pairs only in
bands the nucleotide edges already reach, and blastp's unique band (<0.60 identity) is out of RNA scope.
⚠**That verdict is about the blast-free LCS core, not about whether a protein GRAPH helps**, and r906
measured a real +13.8 pts on the same substrates. The two rows are in tension and this arm settles it.

**Genuinely new: translating the DE NOVO locus representative, not an annotated CDS.** Measured today:
**0 of 457 chr16 pseudogenes have an annotated CDS**, and 6 of the 24 Soto NPIP-family genes have none —
so the annotated route is blind exactly where the thesis lives. De novo 6-frame ORF-finding over the
spliced locus yields **≥100 aa for 92.4% of 2,550 chr16 loci (median 173 aa)** and **rescues 2 of those 6**
(`NPIPB10P` 232 aa, `PDXDC2P-NPIPB14P` 328 aa). The other 4 have no de novo locus at all — a node gap,
not a translation gap.

## ⚠⚠ The truth must NOT be the protein referee

Every other arm this session used the protein-family referee as the neutral truth. **Here it is CIRCULAR
by construction** — the referee IS built from translated CDS clustered by protein homology, so scoring
protein-derived edges against it measures agreement with its own construction. **Primary truth is Soto**
(segmental-duplication-derived, independent of protein space), with the referee reported only as a
consistency check and labelled circular. ⚠Soto carries its own SD-circularity for SEDEF-derived arms
(r741/r1085), but this arm consumes no SD, so Soto is the clean choice here and the referee is not.

## ⚠ The chance-ORF hazard, stated before measuring

**100% of de novo loci yield an ORF ≥50 aa** and the median is 173 aa — a six-frame longest-ORF is
non-zero for any sequence, so **ORF existence is not evidence of coding**. The arm is therefore judged on
whether protein edges add TRUE pairs at acceptable precision, never on ORF yield. A NULL arm is included:
the same pipeline on **shuffled** ORFs (codon-preserving), which must produce near-zero edges. If the null
produces edges, the signal is chance similarity and the arm is void.

## Arms

| arm | edges |
|---|---|
| **N0** baseline | shipped nucleotide edges only (`--dump-graph`) |
| **P** sensitive | N0 ∪ protein edges between de novo locus ORFs (≥100 aa, longest per locus) |
| **NULL** | N0 ∪ protein edges on codon-shuffled ORFs — must add ~nothing |

Protein edge rule copied from §6ko unchanged (identity ≥0.70, coverage ≥0.30 of the longer protein,
≥50 aa aligned) so no new constant is introduced.

## Substrates · metrics · the bar

Development **chr16**; held out **chr2 / chr8 / chr10**, no re-tuning.
Reported: within-family pair coverage · **cross-family edges (the false-merge measurement r906 was blocked
on)** · pairwise and bipartite family F.

| outcome | verdict |
|---|---|
| held-out pair coverage **+5 pts or more** over N0, cross-family edges **≤ N0's**, NULL adds <1% of P's edges | ⭐ **SHIP AS AN OPT-IN SENSITIVE MODE** |
| +1 to +5 pts, cross-family ≤ N0's | ⚠ **PARTIAL** — document, leave off |
| **<+1 pt, or cross-family edges rise, or the NULL fires** | ⛔ **NO** |

**Predicted, before looking — ⚠ PARTIAL.** r906's +13.8 pts came from ANNOTATED CDS on genes that have
one; de novo ORFs are noisier and the loci that most need rescuing (pseudogenes, `PKD1P6-NPIPP1`) are the
ones whose ORFs are most degraded or whose loci do not exist. I expect a real but smaller gain, and the
cross-family count to be the deciding number rather than the coverage.

I will not change the arms, the truth, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO on the bar. The signal is REAL but too small; and the de novo variant was not held-out testable.**

| arm | ORF ≥100aa | protein edges | NULL edges | pair coverage | cross-family |
|---|---|---|---|---|---|
| **chr16 de novo (development)** | 802/864 = **92.8%** | 907 | **0** | 34.9% → **37.5%** (+2.6) | 33 → 38 (+5) |
| chr2 held out (annotated nodes) | 212/639 = 33.2% | 98 | **0** | 35.6% → 36.2% (+0.6) | 3 → 3 |
| chr8 held out | 176/393 = 44.8% | 1,425 | **0** | 20.3% → 20.3% (+0.0) | 5 → 6 |
| chr10 held out | 123/352 = 34.9% | 83 | **0** | 42.6% → 43.0% (+0.4) | 49 → 49 |

⛔**All three held-out gains are under +1 pt**, which is the bar's ⛔ row. ⭐**But the NULL arm is clean
everywhere — 0 edges from codon-shuffled ORFs on all four chromosomes** — so the chance-ORF hazard did not
materialise and the 907/98/1,425/83 protein edges are real homology, just largely redundant.

## ⚠⚠ The confound, stated plainly: the user's actual proposal got DEVELOPMENT-ONLY evidence

De novo assemblies exist only for **chr16**, so the held-out arms had to translate **annotated** nodes.
Those are not the same experiment, and the difference is large: **ORF yield is 92.8% on de novo loci versus
33.2–44.8% on annotated nodes** — de novo loci are read-supported and longer, so the de novo route has
2–3× more material to work with. **The verdict ⛔ is therefore properly scoped to the annotated-node
variant**; the de novo variant shows +2.6 pts at a cost of +5 cross-family edges on development, and
testing it properly needs de novo assemblies of chr2/chr8/chr10.

## ⚠ This does NOT refute r906

r906 measured **+13.8 pts** using **annotated CDS** and §6ko's rule, from a baseline of **80.0%**
within-family pair coverage. This arm uses **six-frame longest ORFs** (a much weaker protein source — often
the wrong frame or a partial region) from a baseline of **20.3–42.6%**. Different protein source, different
rule, different baseline population ⇒ **the two are not comparable and r906's ⚠PARTIAL stands unchanged**,
still unadopted and still blocked on its own false-merge measurement.

⭐**What this arm does settle** is the narrower question the user asked: translating the *locus
representative* — with no annotated CDS, which is the only route that reaches the 457 chr16 pseudogenes
that have none — **is not by itself a sensitive mode.** A longest-ORF is too lossy a protein to beat the
nucleotide edges that already exist.

## Prediction scorecard

Predicted ⚠ PARTIAL, "a real but smaller gain, with the cross-family count deciding rather than the
coverage." **Wrong on which number decided**: cross-family stayed flat on two of three held-out chromosomes
(3→3, 49→49) and it was the *coverage* that failed, at +0.0 to +0.6 pts. The prediction that de novo ORFs
would be noisier than annotated CDS was right in direction but I attributed the cost to false merges when
it is simply that the edges are redundant with nucleotide ones.

> **Generator (2026-09-22 consolidation):** `python3 bench/edge_probes.py protein-denovo ...` — the original `bench/protein_sensitive_mode.py` was folded in verbatim (§6z2); the register rows above cite this file.
