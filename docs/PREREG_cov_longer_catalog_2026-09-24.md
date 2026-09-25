# Pre-registration — the longer-side coverage floor on the shipped catalog, against external truths (§6zj)

**Written 2026-09-24, before either arm is run.** Rows 1097–1098 found that the human chr16 catalog's within-family
pair precision against Compara is 0.268 (multi-exon copies 0.864), because 85% of its copies are single-exon loci and
the rejected families are 80% repeat by base (Alu) — §6cs's unspliced pre-mRNA stubs. Every direct lever is closed:
repeat-masked edges (§6cw, kills 16/67 true NPIP edges), single-exon drops (r356/r1049), depth floors (§6cx: any stub
lever is a low-depth filter). Node-filter measurement today (families recomputed from direct edges): dropping
single-exon nodes with < 5 reads gives chr16 precision 0.54 at −2/103 true pairs, but is precision-neutral (0.983) and
recall-negative (−4%) on held-out gorilla. §6cw's verdict pointed instead at containment on the LONGER side.

## The arm

`RUSTLE_ER_COVERAGE_LONGER_FLOOR=0.30` — the value validated on 09-02 (§6bp/§6bv: 5/5 criteria, 9/9 edge outcomes
cross-species; kept opt-in because its copy-level cost could not be judged without an external truth). ADDITIVE clause:
the shorter side must still clear 0.50 in one colinear record; the longer must clear 0.30. Nothing else changes. The
mechanism: a shared Alu covers ~300 bp of each representative and fails on the longer side; a real duplicate covers
both. It is NOT tuned here — 0.30 is the only value run.

Catalogs: `gw_family_catalog` on human chr16 (A119b, development; truth Ensembl Compara, 275 expressed pairs,
`spectrum/chr16.truth_pairs.tsv`) and on gorilla NC_073244.2 (OR6737 testis, held-out; truth = protein-referee
families, universe = referee pairs with both genes in the SHIPPED catalog, 17,047 pairs, `rule/g44.universe.tsv`; shipped baseline 1,874 recovered at precision 0.983).
Scored at the pair level exactly as rows 1097–1098 (`bench/identity_spectrum.py --catalog --universe` for chr16;
the same mapping for gorilla). Baselines are the shipped catalogs already built (`c16_shipped`, `g44_shipped`).

## Bar — committed now

| outcome | verdict |
|---|---|
| chr16 all-pair precision ≥ 0.50 (from 0.268) AND chr16 recall over the 275 universe ≥ 100/275 (loses ≤ 3 of the 103 pairs) with the ≥ 90% band ≥ 26/31 AND gorilla precision ≥ 0.978 (drop ≤ 0.005) AND gorilla recall ≥ 1,781/17,047 (drop ≤ 5% of the shipped 1,874) | ⭐ **flip the floor ON by default in the catalog builder** |
| chr16 precision ≥ 0.50 but any of the recall/held-out conditions fails | ⚠ **stays opt-in; documented as the switch for deep human libraries, with its measured cost** |
| chr16 precision < 0.50, or the ≥ 90% band falls below 26/31 | ⛔ **the floor does not reach the stub families at the catalog level; report and close** |

**Predicted:** ⚠. Stub↔gene pairs die (short stub, long gene ⇒ cov_longer ≪ 0.30) but stub↔stub pairs of similar
length survive the longer-side floor, so precision rises to 0.4–0.6, not to the multi-exon 0.86; on gorilla the
09-02 arm lost 28% of copies that were 69.5% single-exon and ~1% gene-matched, so referee recall should move < 5%.

I will not change the value, the truths or the bar after seeing any number.

---

# OUTCOME (2026-09-24) — catalogs `rule/g44_cl30.*`, `rule/c16_cl30.*`

| substrate | arm | copies | families | largest (genes) | pair precision (all) | multi-exon | recall (fixed universe) |
|---|---|---|---|---|---|---|---|
| gorilla NC_073244.2, referee | shipped | 357 | 54 | 37 | 0.983 | — | 1,874 / 17,047 |
| | `COVERAGE_LONGER_FLOOR=0.30` | 303 | 47 | 34 | 0.979 | — | **1,421 (−24%)** |
| human chr16, Compara | shipped | 1,418 | 258 | 25 | 0.268 | 0.864 | 103 / 275 |
| | `COVERAGE_LONGER_FLOOR=0.30` | 881 | 176 | 22 | **0.529** | 0.859 | **95 / 275** |

chr16 by band (shipped → floor): ≥ 90 28 → 26/31 · 80–90 11/12 = · 70–80 11/14 = · 60–70 16/29 = · 50–60 23 → 20/30 ·
30–50 12 → 9/61 · < 30 2/98 =. Edges on gorilla 1,675 → 1,469; the 54 copies removed there are 39 single-exon, 32 in
referee genes, all inside PF43's large families, where losing a copy loses its whole row of pairs (537 of the 538 lost
true pairs are PF43; 57 gained).

## Verdict — ⚠ by the pre-registered table: stays opt-in, documented as the switch for deep human libraries

- chr16 precision 0.268 → 0.529 clears 0.50 and the ≥ 90% band holds at 26/31 (exactly the bar); the 60–90% bands are
  untouched. ⭐ the prediction of the mechanism was right: stub↔gene pairs die, stub↔stub pairs of similar length
  survive (GWFAM40: 43 copies / 16 genes / 0 Compara pairs is still there), so precision stops at 0.53, not 0.86.
- chr16 loses 8 of 103 true pairs (bar: ≤ 3), all below 60% identity (50–60: −3, 30–50: −3, ≥ 90: −2).
- ⛔ the held-out cost was mis-predicted: "< 5%" became −24% of referee pairs at −0.004 precision. The 09-02 arm's
  "lost copies are 69.5% single-exon and ~1% gene-matched" did not carry over to this contig, where the removed
  single-exon copies sit inside PF43 genes.

So: `RUSTLE_ER_COVERAGE_LONGER_FLOOR=0.30` is the right switch for a deep human library whose node set is
dominated by unspliced Alu stubs (it doubles chr16 all-pair precision for 8 pairs), and the wrong default for the
gorilla substrate the thesis runs on. It is NOT flipped. The remaining chr16 false families are stub↔stub and are
not reachable by any edge clause; they are a node-admission question whose every lever is closed (r356/r1049/§6cx).
