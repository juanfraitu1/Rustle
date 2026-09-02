# O1 precision levers — what was tried, what survived, what shipped

> ## ⛔ NOTHING IN THIS DOCUMENT IS A DEFAULT.
> **UPDATE 09-02 (§6bv): two-sided coverage has now been formally ADJUDICATED against criteria fixed
> in advance (§6bu) and the verdict is KEEP OPT-IN.** It passed a cross-species pre-registered test
> 5/5 and was still not shipped, because its copy-level cost is annotation-neutral, non-monotone,
> and unadjudicable without a positive stratum. **That combination — passing the test and declining
> the default — is the point, not an inconsistency.**
> **UPDATE 09-02 (§6bp):** the two-sided coverage clause has now been run **END TO END through the real
> binary** and passes all three criteria — the OFF arm is byte-identical to `arm_f2`, the params
> certificate distinguishes the arms, and NPIP recall **improves 14/31 → 15/31**. It ships as
> `RUSTLE_ER_COVERAGE_LONGER_FLOOR`, **default OFF**. Everything else below remains offline-only.
> Every result here is an **offline reconstruction** scored against a **proxy** (CDS concordance).
> By the T8 rule an offline re-derivation is a **hypothesis generator, never a test**. No default has
> been flipped, no edge rule changed, and the only `src/` addition from this line of work is
> `src/bin/gamma_refine.rs`, a bench binary that calls the shipped γ so offline arms use the real
> partitioner instead of a re-implementation.

Session of 2026-09-01. Full detail in `o1_ledger.md` §6bg–§6bn; negative results in
`NEGATIVE_RESULTS_REGISTER.md`.

## The scoreboard

| lever | verdict | why | §
|---|---|---|---|
| **Two-sided coverage** (`cov_longer` floor) | ⭐⭐ **ADJUDICATED 09-02: KEEP OPT-IN** | edge filter validated cross-species (**5/5 pre-registered, 9/9 edges**, §6bt.2); ⛔ copy-level loss is **annotation-neutral and non-monotone** — 18 families deleted, 29 loci created; the **14→15 gain is one 274 bp 2-read node** | §6bp/§6bt.2/**§6bv** |
| **Dense clusters inside families** (γ′=0.70, families ≥3) | ⭐ **second candidate** | within-biotype separation 31× / 15×, reproduces | §6bn |
| E_c within-family refinement | ⚖️ scoped only | safe and real, but reaches only near-identical arrays | §6bi/§6bj |
| Direct-edge requirement | ⛔ no default form | as a pair rule it cannot apply to a set; as k-core it loses 2/3 of families | §6bg/§6bh |
| Connectivity (triangle, k-core, articulation, min-degree) | ⛔ subsumed / artefact | adds nothing over `direct`; global form deletes every 2-member family | §6bh/§6bm |
| Coverage **+** connectivity composite | ⛔ does not reproduce | +10% F1 genome-wide, **+0.004** on arm_f2, cost on both | §6bm |
| **Cliques** (γ′=1.00) | ⛔ too strict | worse F1 on both; its split set starts cutting GOOD pairs | §6bn |
| Raising γ globally | ⛔ wrong lever | F1 **peaks at the shipped 0.20** for the seeded catalog | §6bg |
| Same-chromosome prior | ⛔ worst tried | F1 0.7536 → 0.6880; dispersed paralogy is real | §6bg |

## The two candidates, and what each still needs

**1. Two-sided coverage.** Charge coverage on the longer sequence as well as the shorter.
Non-ZNF F1 `arm_f2` 0.3182 → **0.4043**, genome-wide 0.4306 → **0.4595**, NPIP recall unchanged at
14/31. Costs 23–31% of copies on `arm_f2`, 13–15% genome-wide.

⚠**The floor must be ASYMMETRIC.** Shorter stays 0.50; the longer floor should be **0.30**. The
symmetric 0.50 costs NPIP (14→12) and 44% of families genome-wide.
⚠**CORRECTION — the PRE-EXISTING flag does not implement this.** `RUSTLE_ER_COVERAGE_LONGER`
**replaces** the denominator, applying the *same* floor to the longer side — the symmetric 0.50 variant
the measurements reject. ✅**RESOLVED 09-02**: the asymmetric form is now implemented as
`RUSTLE_ER_COVERAGE_LONGER_FLOOR` (§6bp), default OFF, OFF byte-identical.

⭐**The code anticipated all of this.** `denovo_pipeline.rs:4900-4926` already records: coverage-of-
shorter is *"STRUCTURALLY BLIND to truncation — a 10% fragment that aligns fully into a complete
sibling scores 1.00"*; a certificate on RNA pairs giving **129 TRUE vs 2 FALSE (precision 0.985)**; a
named failure — *a 2,037 bp NPIPB6 fragment reaches coverage 0.948 against a 38,653 bp chimeric
read-through node while touching 5% of it, dragging EIF3CL into the NPIP family*; and the argument the
advisor will care about most — it is what would make component merging **MONOTONE under improving
evidence**, so that *"better reads can only merge components"* becomes true.
It also predicted the cost mechanism: *"dividing by the longer span penalises pairs whose ANNOTATION
differs rather than whose SEQUENCE does"*, with a measured NPIP ceiling of **134/171 true pairs able to
reach 0.50 at all** and NPIPB8–NPIPB2 capping at **0.215**. ⟹**that is precisely why 0.50 is too
strict and 0.30 is not**, and it is the question the knob was created to answer:
*"whether the precision gain pays for it end to end."* This session answers it **offline**; the
end-to-end half is still open.

**2. Dense clusters inside families.** After building at the shipped γ=0.20, re-refine each family of
≥3 members at γ′=0.70. Adds no new concept — the shipped γ at a second value. NPIP unchanged, costs
4–6.5% of copies, and *increases* family count rather than deleting families.
⚠ its headline F1 movement is ~1 pair; the evidence is the **within-biotype** contrast (ZNF vs ZNF:
splits 0.0126/0.0158 against keeps 0.3866/0.2429).

## Before either can ship — updated 2026-09-02

**Two-sided coverage is now ADJUDICATED (§6bu rule, §6bv verdict): KEEP OPT-IN.**

1. ~~Run through the **real binary**~~ ✅ **DONE (§6bp)**; the offline model proved faithful
   (predicted 1,643 edges, the binary gave 1,652). ⚠ Still needed for γ′=0.70.
2. ~~The **HUMAN 150-window false-merge panel**~~ ✅ **DONE (§6bq → §6bt.2): 5/5 pre-registered
   criteria, 9/9 edge outcomes, prediction committed to git before the run.**
3. NPIP recall on the **3-contig** catalog — done, **14/31 → 15/31** ⚠**but see §6bv: the gain is one
   274 bp, 2-read, single-exon copy**, in a family that exists only because γ re-partitioned.
4. A substrate whose large families are **not** zinc-finger clusters — **still open**.
5. ⭐**NEW, and it is now the binding one: a POSITIVE STRATUM** labelling a real share of the
   **220 lost copies**. NPIP labels **21 of 678**. Without it, §6bv's D4 null is the ceiling of what
   can be said about the copy cost, and the clause cannot earn a default.

### ⚖️ Why it did not ship, in one paragraph

It is a **validated EDGE filter and an INDISCRIMINATE COPY filter**, simultaneously. Edge level:
9/9 across a species boundary. Copy level: **annotation-neutral** — within exon-structure strata the
lost and kept copies carry a reciprocal RefSeq match at the *same* rate (single-exon 0.013 vs 0.009;
multi-exon 0.791 vs 0.784), so it deletes real and spurious copies alike. It is also **not monotone
at the catalog level**: the node set is identical (3,598, 0 diff) yet **220 copies are lost and 32
gained, 29 at loci that had none**, because removing edges lets γ split a sparse component into
denser passing pieces. And **18 families are deleted outright** (sizes 2–6). ⭐**One benefit is not in
dispute and is quotable alone: largest family 54 → 32 — the ZNF blob splits.**

⚠⚠**METHOD LESSON, and it cost the D4 criterion.** §6bu pre-declared **length** as D4's confound.
Within length quartiles the clause looked like it preferentially deleted uncorroborated copies, in
all four strata. The operative confound was **exon structure**, and controlling for that instead
erased the effect entirely. **Naming a confound in advance does not protect against naming the wrong
one.**

## Instruments (all in `bench/`, all runnable)

`er_both_coverage.py` · `er_both_coverage_gw.py` · `er_scoped_density.py` ·
`er_coverage_plus_connectivity.py` · `fp_sources_seeded.py` · `fp_sources_read.py` ·
`fp_connectivity_levers.py` · `ec_blob_split.py` · `catalog_partition_vs_cds*.py`

**Added 2026-09-02** — `adjudicate_covlonger.py` (prices the copy loss against §6bu's five criteria) ·
`o1neg_covlonger_predict.py` (the pre-registration predictor) · `o1neg_score_arms.py` (two-arm panel
scorer) · `catalog_self_overlap_audit.py` (§6bs) · `er_identity_band_external.py` (§6br).

⚠**The proxy's limits are part of the result.** CDS concordance gives read-catalog precision **0.8411**
on the 118 loci both catalogs see and **0.1234** genome-wide; **74.5%** of genome-wide pairs are
ZNF–ZNF scoring 0.0778. It is a **relative** instrument — valid for contrasts inside one run and
within one biotype, never as an absolute false-positive rate.
