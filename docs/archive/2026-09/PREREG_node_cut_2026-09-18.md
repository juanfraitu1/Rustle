# Pre-registration — the NODE CUT rule (NC): split a known chimeric record at its parent boundary

**Written 2026-09-18, BEFORE any metric under the cut rule existed.** Follows §6m2 (`--chimera-policy`
NOT ADOPTED — a truth-side lever cannot fix this) and the user's question *"What would entail that a node
is in two families? Can we later cut them somehow?"*

The answer this pre-registers: a record is in two families because **its pieces are**, one each. Levels stay
connected components of one graph, so T1 (nesting), T1′ (laminar) and T2 (expression view) survive
unchanged. This is a **node rule**, which is what §0★★★.7e already identified as the only thing that could
move the DNA levels — not a threshold, not a scoring policy.

## 0. What is ALREADY KNOWN (declared, so nothing below is presented as fresh)

Measured earlier today, before this document:
- **95.6%** (173/181) of readthrough records whose halves are both gene records have exons falling into
  exactly **two contiguous parent runs** — a single cut point exists. 8 are interleaved or single-parent.
- The chimeras are **NOT articulation points** in NPIP's component at any level (only `PDXDC2P-NPIPB14P`,
  and `BOLA2-SMG1P6` at L3). Deleting every non-member cut vertex takes L2 from 91 → 47 nodes and still
  leaves 21 outsiders while losing 1 member. Route "biconnected components" is already refuted.
- §6m2: NPIP h_join = 1.0000 is held by CLN3, EIF3CL and the LOC lncRNAs — **none of them chimeras**.
- Only **1 of 27** NPIP members and **0 of 19** TBC1D3 members carry the readthrough flag.

**Stated expectation, before measuring: the cut will improve definitional coherence and is unlikely to move
a family metric much, and NC-3 is expected to fail.** This is recorded so a null result cannot be
re-narrated afterwards.

## 1. The rule (deterministic, no threshold, no new alignment)

Chimera set = the frozen RefSeq `description=readthrough` set of the §6m2 prereg (209 nodes).

For each such **primary** node:
1. **Parents.** Split `Name` as `A-B` at the rightmost hyphen for which at least one half is itself a gene
   record. If neither half is, the node is `uncut`.
2. **Exon labels.** For the node's longest transcript, label each exon by parent-span overlap: `A`, `B`,
   `X` (both), `.` (neither). Where only one parent has a record, exons overlapping it take that label and
   all others take the complementary label.
3. **Admissible iff** the `A`/`B` labels form exactly **two contiguous runs** (ignoring `X`/`.`). Otherwise
   `uncut`.
4. **Cut coordinate** c = midpoint of the intron between the last exon of run 1 and the first exon of run 2.
   Piece P1 = exons ending ≤ c; P2 = exons starting ≥ c. Same chrom and strand; `X`/`.` exons follow the
   side of c they lie on.
5. **Edge re-attribution, from the PAF batches `family_cert` already wrote — no realignment.**
   - chimera as TARGET: assign each record to the piece holding the **majority of its aligned target bp**
     (`tx_exon_blocks`); ties to P1.
   - chimera as QUERY: map c to transcript coordinate t_c = cumulative exon length of P1; assign by the
     majority of `[qs, qe)` on either side of t_c; ties to P1.
   - Every record lands on exactly one piece. Aggregates (identity, coverage, f_ex, w_98) are recomputed
     per piece using **that piece's own** exon and body lengths as denominators.
6. `uncut` nodes are left **whole** and counted, never dropped.

## 2. Arms — and the limitation, stated before measuring

⚠ **A strong held-out arm does not exist for this rule.** Of 137 CAT chimeras only **10** are in any Soto
family, and exactly **one** chimera-containing Soto family sits on a chromosome that is neither already used
nor reserved for the lattice hold-out (chr8–11). The held-out arm below is therefore small, and
**chr10 ID_233 is deliberately EXCLUDED to protect the chr8–11 reservation** (`PENDING_2026-09-17` item 3).

| arm | content | status |
|---|---|---|
| **DEV** | NPIP (27 members), TBC1D3 (19) — levels L1a/L1b/L2/L3, `lattice_rules/records.tsv` | development families; the rule was designed with knowledge of these |
| **FRESH** | 5 Soto families named here before any score: **ID_453**, **ID_454** (SLX1B-SULT1A4, chr16 ~29.7 Mb — a different locus from NPIP at 15–18 Mb), **ID_380** (LINC02210-CRHR1, chr17), **ID_480** (TVP23C-CDRT4, chr17), **ID_121** (RNF103-CHMP3, chr2), **ID_305** (C1QTNF3-AMACR, chr5) | never looked at; small (15–17 members total). Weak, and reported as weak |
| **APPLIC** | all 181 records genome-wide: cut admissible / `uncut`, and the piece-size split | the 95.6% is already known (§0) and is reported as such |

## 3. Pre-registered decision rules

**NC-1 (safety, DECISIVE).** No family may lose a member. Bar: on BOTH the DEV and FRESH arms, pairwise
sensitivity does not decrease at any level, and no truth member ends up outside the group holding its
family's plurality. **A single lost member fails the rule outright**, whatever else improves.

**NC-2 (benefit).** On the **FRESH** arm, pairwise precision improves by **≥ 0.05** with no recall loss —
the same bar §6ks used for AP-1. DEV-arm precision is reported but cannot carry the decision.

**NC-3 (certificate).** Does any (h_join, h_split] become non-empty at any level after the cut?
**PRE-DECLARED EXPECTATION: NO** (§0: the boundary is CLN3/EIF3CL/LOC lncRNAs, not chimeras). A
non-empty interval requires independent confirmation before being claimed anywhere.

**NC-4 (denominator guard, inherited from §6m2 CP-5).** The cut **adds** nodes. A truth record maps to the
piece carrying the **majority of its exon bp**; the other piece counts as an outsider unless it matches a
different truth record — so one truth copy can never be credited twice. Any precision gain accompanied by
universe growth > 10% is reported as **universe growth, not a gain**.

**NC-5 (no silent drops).** The `uncut` records are counted and reported at every level. A rule that
improves a number by quietly removing records from scope fails.

## 4. Out of scope

- No edge test, grouping operator or threshold changes.
- No claim that the cut fixes NPIP's certificate; §0 says it should not.
- No new alignment is run. If the re-attribution in §1.5 turns out to be impossible from the stored PAFs
  for some record class, that class is reported as `uncut`, not approximated.
