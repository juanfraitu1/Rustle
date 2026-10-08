# Pre-registration — `--chimera-policy`: a truth-side lever for known readthrough/chimeric records

**Written 2026-09-18, BEFORE any metric under any policy was computed.** Prompted by the user after
§6m1 + the 09-18 addendum established that PKD1P6-NPIPP1 is a REAL fused transcript (110 MAPQ-60 Iso-Seq
reads over a canonical GT-AG junction), so the NPIP DNA-certificate blocker cannot be removed by deleting a
mis-annotation.

User's request, verbatim: *"maybe could we have some or more chimeras accepted if they are known and a lever
to determine if we count them or not for precision, sensitivity and bipartite matching"*.

This is an **evaluation lever**, not a change to the family definition. No edge test, no grouping operator
and no threshold changes. Only how a *known* chimeric record is scored.

---

## 1. The chimera set — frozen before scoring, from the annotation only

The set is built from the ANNOTATION and is independent of any prediction this project makes. This is the
condition that keeps it clear of the standing metric trap "a denominator conditioned on the prediction"
(`feedback_metric_traps`, 7 retractions).

| substrate | annotation | rule | count (declared now) |
|---|---|---|---|
| human RefSeq (NPIP / TBC1D3, chr16 / chr17) | `winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz` | `gene`/`pseudogene` record whose **`description=` contains "readthrough"** | **209 genome-wide; 12 on chr16; 5 touching NPIP** |
| CAT/GENCODE (held-out chr5 / chr7 / chr21) | `lit/soto_fams/cat_genes.tsv` | record **named `A-B` where both `A` and `B` are themselves gene records** in the same annotation (no `description` field exists here) | **137 genome-wide; 15 on chr5/7/21** |

Cross-check already on disk: `/mnt/linuxdisk/home/juanfraitu/lattice_rules/S5_readthroughs.tsv` is the same
RefSeq set (209 rows) with both `by_description` and `by_name` columns; the name rule is a strict subset
(181 of 209) and every one of the 181 also carries the description flag.

The 5 NPIP-touching chimeras, all flagged, declared now: `PKD1P6-NPIPP1`, `PKD1P3-NPIPA1`, `PKD1P4-NPIPA8`,
`PKD1P5-LOC105376752`, `PDXDC2P-NPIPB14P`. The last is the record §0★★★.7e noted a containment rule *cannot*
see, because `NPIPB14P` never became a node — the curated flag sees it anyway.

⚠ The two substrates use DIFFERENT detection rules because the annotations differ. Numbers from the two are
never pooled. The RefSeq rule is curated (stronger); the CAT rule is structural (weaker) and its recall
against a curated truth is unknown on that substrate.

## 2. The three policies

Exactly one policy applies to an entire run — every family, every catalog, every arm of a comparison. A
policy is never chosen per family, per catalog, or after seeing a number. All three are always reported
side by side; reporting one alone is a protocol violation.

| policy | truth side | prediction side | question it answers |
|---|---|---|---|
| **`strict`** | chimera is one ordinary truth record in one family (whatever the truth source says) | unchanged | the current baseline — every published number to date is `strict` |
| **`exclude`** | chimera removed from truth | chimera removed from the predicted clustering **and** from the "extra genes pulled in by a touching cluster" set | what the family looks like with known fusions set aside; symmetric removal, so no denominator is conditioned on the prediction |
| **`multi`** | chimera carries the truth label of **both** halves (where a half's family is known); implemented as two pseudo-records sharing the chimera's predicted cluster | unchanged | the overlapping-membership route: a copy that genuinely belongs to two parent families |

**`multi` metric definitions, fixed now:**
- *Bipartite matching*: the standard reduction — each chimera contributes one row per parent label, each
  row carrying the chimera's predicted cluster. R and P are then computed unchanged on the expanded table.
- *Pairwise precision*: a pair (chimera, x) is a TP if x is in **either** parent family; it is never an FP
  on the grounds of the other parent. A pair (chimera_A, chimera_B) is a TP if the two share any parent.
- *Pairwise sensitivity*: the chimera's true pairs are the union over both parent families.
- `multi` ADDS truth rows; it never removes them. The added rows do not depend on the prediction.

## 3. Arms

| arm | substrate | catalogs | status |
|---|---|---|---|
| A | held-out chr5/7/21 (CAT) | `lit/ap_ho/e1` vs `lit/ap_ho/e1s`, 11 pre-chosen Soto families (§6ks `FAMS`) | substrate already USED by §6ks; this re-scores existing catalogs, it does not re-tune anything |
| B | human RefSeq chr16 / chr17 | NPIP and TBC1D3 lattice levels L1/L2/L3 (`lattice_rules/records.tsv`) | development families; descriptive only |

## 4. Pre-registered decision rules

**CP-1 (headline, decisive).** Does chimera handling change the §6ks/§6kt conclusion?
On arm A, E1S beat E1 under `strict` (bipartite F universe 0.831 → 0.881, pairwise precision 0.815 → 1.000).
- **ROBUST** if E1S ≥ E1 on BOTH bipartite F and pairwise precision under ALL THREE policies.
- **CONDITIONAL** if any policy reverses either comparison. Then §6ks's AP-1 result must be restated as
  holding only under `strict`, and the shipped `--min-shared-exon-frac 0.30` default is re-opened.

**CP-2 (the lever's value).** On arm A, does any policy move bipartite F by **≥ 0.02** from `strict` in
either catalog? If no policy moves any arm-A number by ≥ 0.02, the lever is declared a **NO-OP on this
substrate** and its scope is restricted, in writing, to families that actually contain a chimera.

**CP-3 (arm B, the overlapping-membership route).** Under `multi`, does NPIP's or TBC1D3's bipartite F rise
by ≥ 0.02 over `strict`? Reported, not decisive — these are development families.

**CP-4 (certificate, reported with a pre-declared expectation).** Under `multi`, does NPIP's DNA-level
certificate interval (h_join, h_split] become non-empty at L1, L2 or L3?
**PRE-DECLARED EXPECTATION: NO.** §6m1 measured that removing readthroughs moves the boundary to CLN3,
EIF3CL and the LOC lncRNAs, none of which are chimeras and none of which carry the flag. If an interval DOES
open, that is a surprise and requires an independent confirmation before it is claimed anywhere.

**CP-5 (guard against a shrinking denominator).** Every reported cell carries the number of genes in its
universe. If `exclude` improves a number while shrinking that universe by more than 10%, the improvement is
reported as **denominator shrinkage, not a gain** — the §6m0 addendum-3 trap, which already cost one
retraction this month.

## 5. What is NOT being tested here

- No change to any edge test, grouping operator, or threshold.
- No claim that a chimera is or is not "a copy". That definitional question stays open in
  `docs/PENDING_2026-09-17.md` item 1; this lever only makes the two answers comparable.
- No new substrate is consumed. Arm A re-scores catalogs built on 2026-09-14; arm B is development data.
