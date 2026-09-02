# O1 precision levers — what was tried, what survived, what shipped

> ## ⛔ NOTHING IN THIS DOCUMENT HAS SHIPPED.
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
| **Two-sided coverage** (`cov_longer` floor) | ⭐ **best candidate** | gains on BOTH substrates, NPIP-neutral at 0.30 | §6bk/§6bl |
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
⚠**CORRECTION — the existing flag does NOT implement this.** `RUSTLE_ER_COVERAGE_LONGER` **replaces**
the denominator, applying the *same* floor to the longer side, i.e. it is exactly the symmetric 0.50
variant the measurements reject. The asymmetric form needs a **second constant that does not exist
yet**. Earlier statements in §6bk/§6bl that "the flag already exists" are corrected here.

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

## Before either can ship

1. Run through the **real binary**, not an offline reconstruction.
2. The **HUMAN 150-window false-merge panel** (E_r-only levers can take it; E_c cannot — no human RNA).
3. NPIP recall on the **3-contig** catalog — `fibro_gwcat` covers the truth set at only **1/31**.
4. A substrate whose large families are **not** zinc-finger clusters.

## Instruments (all in `bench/`, all runnable)

`er_both_coverage.py` · `er_both_coverage_gw.py` · `er_scoped_density.py` ·
`er_coverage_plus_connectivity.py` · `fp_sources_seeded.py` · `fp_sources_read.py` ·
`fp_connectivity_levers.py` · `ec_blob_split.py` · `catalog_partition_vs_cds*.py`

⚠**The proxy's limits are part of the result.** CDS concordance gives read-catalog precision **0.8411**
on the 118 loci both catalogs see and **0.1234** genome-wide; **74.5%** of genome-wide pairs are
ZNF–ZNF scoring 0.0778. It is a **relative** instrument — valid for contrasts inside one run and
within one biotype, never as an absolute false-positive rate.
