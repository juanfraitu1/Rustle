# Pre-registration — admitting SINGLE-EXON loci whose reads fall entirely inside one exon of a family member

**Written 2026-09-22, §6z1, before any candidate is admitted.** User: *"lets only select single exon
transcripts if the entirety of the long read only maps to a single exon of a bigger transcript."*

## Why this population is worth admitting — three measured priors

- **r356**: dropping single-exon candidate loci costs **precision as well as recall** —
  `P 0.709→0.664, R 0.698→0.560, F1 0.704→0.608`, copies **617→166**. *"Deleting them costs PRECISION
  too; the ≥2-loci gate dissolves the small pure families."* They are load-bearing.
- **r345**: dropping both-intronless edges (3,267, ~63% of the graph) *"would destroy real single-exon
  families at ~10:1 false:true"*, but the **one-side form flags 2/2 retros and 0/5 reals**. ⭐The
  one-sided, family-anchored shape is the one already shown to work.
- **r508 / r84**: **303 of 328** discovered loci are single-exon, *"the class our own isoform guard treats
  as suspect"*, and the shipped catalog contains **0 of 1,415** single-exon nodes. The population is
  currently excluded wholesale.

## ⚠⚠ SCOPE — what this criterion catches, and what it does NOT

A **processed pseudogene** is a retrocopy of the spliced **mRNA**, so it corresponds to *several fused
parent exons*; a read from it spans multiple exons of the parent and **this criterion rejects it**.
So the rule admits **DUPLICATED EXONS**, not processed pseudogenes — a different population from the one
§6u1/§6u2/r927 identified (*"half of family pairs have a genuinely intronless member"*, biological and
unliftable). ⭐**Stated before running so the result is not later read as closing the §6u1 gap.** The count
of candidates rejected *because* their reads span >1 parent exon is reported, so we learn the size of the
population this criterion deliberately leaves out.

## The rule

- **Candidate C**: an intronless (single-exon) de novo locus, currently excluded from the catalog.
- **Anchor F**: a family **already accepted from spliced members and FROZEN before this pass**.
  ⚠Anti-circularity: F must not be re-derived using C, or the family would define its own members.
- **Admission test (the user's criterion, verbatim in operational form)**: for **every** read supporting C,
  the read's alignment must lie **entirely within ONE exon** of a member of F — it may not span an
  exon-exon junction and may not extend past that exon's boundaries.
- C joins F only if it also clears the shipped identity floor against F's member.

⚠ The anchor is F's **SPLICED representative**, not its gene body. r394 measured that the genomic-span
substrate lifts both-spliced pairs 0.259→0.914 but a one-single-exon pair only **0.000→0.022** — the
genomic substrate is the wrong one for this class, and §6y8 showed the spliced substrate is a strict subset
*for spliced↔spliced*, which does not apply here because the candidate has no introns to align.

## Metrics · substrates · bar

Per arm: candidates admitted · **cross-family admissions** (a candidate matching ≥2 frozen families is the
false-merge signature) · families whose membership changed · family pairwise precision / recall / F against
the protein referee and Soto. Plus the **scope counter**: candidates rejected for spanning >1 parent exon.

⚠**Admission can only ADD members**, so recall can only rise and precision can only fall. Precision is
therefore the decisive number, exactly as in r345's retro/real split.

Development **chr20** (the substrate with a de novo GTF on disk). Held out **chr16**.

| outcome | verdict |
|---|---|
| held-out recall up **≥0.02** with precision down **≤0.01** and **0 cross-family admissions** | ⭐ **ADOPT** |
| recall up with precision down 0.01–0.03, cross-family ≤1 per chromosome | ⚠ **PARTIAL** — flag, leave off |
| precision down >0.03, **or** any systematic cross-family admission | ⛔ **NO** |

**Predicted, before looking — ⚠ PARTIAL.** r345's one-sided form scored 2/2 retros and 0/5 reals, which is
the right shape, but on 7 cases. r356 says the population is load-bearing, so admitting some of it should
raise recall. My concern is the scope note above: the strict "entirely within one exon" test is a *narrow*
filter, and I expect it to admit few candidates — most single-exon loci that belong to a family are
retrocopies spanning several parent exons, which this rejects by design. **If admissions are near zero,
that is a finding about the criterion's reach, not about whether single-exon loci belong in families.**

I will not change the rule, the anchor, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **0 admissions. The criterion is arithmetically unsatisfiable for this population.**

chr20, 1,779 de novo loci → **1,289 multi-exon anchors** (clustered first and FROZEN: 17 families,
40 members, anti-circularity honoured) and **490 single-exon candidates**.

| | |
|---|---|
| candidate→frozen-member records passing identity ≥0.70, alen ≥300 | 231 (53 distinct candidates of 490) |
| ⚠ rejected: candidate not fully mapped (qcov < 0.95) | **227** |
| ⚠ rejected: aligned region spans >1 anchor exon (the pre-registered SCOPE counter) | **4** |
| ⭐**ADMITTED** | **0** |
| cross-family admissions | 0 |

## Why — and it is not the reason I pre-registered

I predicted the binding filter would be the scope one: *"most single-exon loci that belong to a family are
retrocopies spanning several parent exons, which this rejects by design."* **That filter rejected 4 of
231.** The binding filter is the *other* half of the criterion — **"the entirety of the read maps"** —
and the arithmetic behind it is decisive:

- median candidate length **7,974 bp**
- median largest exon of the anchor it matches **2,141 bp**
- ⭐**230 of 231 records (99.6%) have a candidate LONGER than any exon the anchor possesses**
- median query coverage **0.169**; only **4 of 231** records reach qcov ≥ 0.95

**An 8 kb candidate cannot lie inside a 2 kb exon.** The criterion is not failing on biology, it is failing
on size.

## ⭐⭐ The real finding: "single-exon locus" does not mean "one-exon transcript"

A single-exon de novo locus has a median span of **7,974 bp**, where a real human exon is ~150 bp. These
are not duplicated exons — they are **regions with no splicing evidence**, whose extent is genomic
(an unspliced read pileup over a span that *contains* introns rather than lacking them). That is exactly
why r508 called them *"the class our own isoform guard treats as suspect"*, and it means the population
this rule was aimed at is not the population the class actually contains.

⚠**Consequence for the idea, not a refutation of it.** The user's criterion is a sound precision guard and
would do what it promises — against candidates that are genuinely exon-sized. It cannot reach the chr20
single-exon class because that class is built from genomic spans. Any future attempt must first ask what a
single-exon locus's span *is*, before testing whether it fits inside an exon.

⚠ Held-out chr16 **not spent** — a criterion that admits 0 of 490 on development does not get a held-out
substrate (the §6w6 precedent).

⚠ No threshold was relaxed after seeing the numbers. The qcov distribution is reported as a diagnostic
(median 0.169, 23 of 231 records at ≥0.50) and is **not** an invitation to lower the bar: at qcov 0.50 the
"entirety of the read" clause is no longer the user's criterion.
