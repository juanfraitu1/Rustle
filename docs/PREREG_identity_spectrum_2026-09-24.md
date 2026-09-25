# Pre-registration — the identity spectrum: which alignment tier finds which paralogue pairs, against Ensembl Compara

**Written 2026-09-24 (§6zh), before any tier is run against the truth.** User's concern: the advisor may say the
family definition keeps only the most similar copies and loses real biology. The edge builder already runs three
tiers (asm20 k = 19 · sensitive `-k11 -w5` · protein tail via mmseqs). This arm measures, band by band of
sequence identity, what each tier recovers of a truth the pipeline never saw.

## Truth — Ensembl Compara paralogues (BioMart, `hsapiens_paralog_*`)

Every human within-species paralogue pair with both genes on the test chromosome (chr16 = the family-rich
development chromosome, chr20 = the second substrate), with Compara's protein percent identity
(`paralog_perc_id`, the query-side value; pairs are symmetrised with the max of the two sides) and paralogy
subtype. Genes are matched to our loci by HGNC symbol through the RefSeq CHM13 annotation. The universe is the
set of Compara genes that have an expressed locus in the assembled GTF of that chromosome; a truth pair counts
only when BOTH genes are in the universe (a pair with an unexpressed member is not recoverable by an RNA
pipeline and is reported as such, not as a miss).

## Objects

- **Nodes**: one representative spliced sequence per expressed locus (`gene_id` group of the assembled GTF,
  representative = the transcript with most reads, tie → longer), exactly the sequence the catalog's edge
  builder aligns.
- **Tiers**, each run all-vs-all on the same node set with the builder's own flags:
  - T1 `minimap2 -x asm20 -c -X --no-long-join` (k = 19): edge iff identity ≥ 0.80, coverage of the shorter ≥ 0.50;
  - T2 `minimap2 -x asm20 -k11 -w5 -c -X --no-long-join`: edge iff identity ≥ 0.60, coverage ≥ 0.50;
  - T3 `mmseqs easy-search --search-type 2` (translated, TBLASTX-like): edge iff protein identity ≥ 0.30,
    e-value ≤ 1e-5, query coverage ≥ 0.50.
  A pair's identity is the best record's; T1/T2 identities are nucleotide, T3's is protein.

## Metrics — committed now

Per chromosome, per tier (T1; T1 ∪ T2; T1 ∪ T2 ∪ T3), per Compara protein-identity band (≥ 90, 80–90, 70–80,
60–70, 50–60, 30–50, < 30%):
- **recall** = truth pairs (both in the universe) that the tier aligns above its floor / all truth pairs in the band;
- **precision by our identity band** = among pairs the tier aligns above its floor, the fraction that are Compara
  paralogues, in bands of OUR identity (nucleotide for T1/T2, protein for T3), restricted to pairs whose two genes
  both have Compara data (a pair with an unlisted gene cannot be judged);
- the count of truth pairs recovered by NO tier, with their Compara identity.

| outcome (chr16 and chr20 agree) | verdict |
|---|---|
| T1 ∪ T2 recall ≥ 0.80 on truth pairs at ≥ 70% protein identity, and T2's added pairs at 60–80% nucleotide identity have precision ≥ 0.50 | ⭐ **the default tiers keep the biology down to the point where reads stop cross-mapping; the loss is below 70%, where only the protein tier looks** |
| recall ≥ 0.80 only at ≥ 80% | ⚠ **the sensitive tier is not doing its job in the 70–80% band — say so and quantify** |
| recall < 0.80 at ≥ 80% | ⛔ **real copies are being lost by seeding; the advisor is right** |

T3 is reported for the < 70% bands with its precision: the prediction (r1027–r1031) is that it recovers most
Compara pairs there but at low precision against Compara, because Compara itself separates old paralogues from
fold-sharers by gene trees, not by identity.

**Predicted, before looking:** ⭐ — T1 alone ≥ 0.9 at ≥ 90% and falling through 80–90%; T2 lifts 70–80% above
0.7; both near zero below 60% where T3 takes over at precision 0.3–0.6. Pairs recovered by no tier: dominated by
< 50% protein identity and by pairs with one unexpressed member.

I will not change tiers, floors, bands or the bar after seeing any number.

---

# OUTCOME (2026-09-24) — `bench/identity_spectrum.py`, runs in `/mnt/linuxdisk/tmp/gw22/spectrum/`

Truth from BioMart (`useast.ensembl.org`; the main site and REST were down): chr16 1,324 same-chromosome Compara
pairs, **275 with both genes expressed** (1,049 have an unexpressed member); chr20 965 pairs, 94 with both
expressed, none above 80% protein identity — chr20 is uninformative for the bar and is reported only.

## chr16 — recall of Compara pairs (both genes expressed), by Compara protein identity

| band | pairs | T1 asm20 | T1 ∪ T2 sensitive | T1 ∪ T2 ∪ T3 protein |
|---|---|---|---|---|
| ≥ 90% | 31 | 0.774 | **0.806** | 0.806 |
| 80–90% | 12 | 0.333 | **0.333** | 0.333 |
| 70–80% | 14 | 0.286 | 0.357 | 0.357 |
| 60–70% | 29 | 0.207 | 0.241 | 0.241 |
| 50–60% | 30 | 0.167 | 0.267 | 0.267 |
| 30–50% | 61 | 0.033 | 0.033 | 0.131 |
| < 30% | 98 | 0.031 | 0.031 | 0.071 |

Precision of the aligned pairs that Compara can judge: T1 0.957 (n 50), T2 0.957–1.000 across its bands (n 55;
its added 70–80% nucleotide-identity pairs 5/5 correct), T3 0.29–0.67 (n 38).

**By the pre-registered bar: ⛔ on chr16** — ≥ 90% clears 0.80 (0.806) but 80–90% is 0.333. **The prediction that
the sensitive tier lifts 70–80% above 0.7 was wrong** (0.357), and the reason is not the one the bar assumed.

## Why the pairs are missed (T1 ∪ T2, pairs at ≥ 60% Compara identity)

| band | no nucleotide record | record, coverage < 0.50 |
|---|---|---|
| ≥ 90% | 4 | 2 |
| 80–90% | 1 | **7** |
| 70–80% | 6 | 3 |
| 60–70% | 13 | 9 |

The 80–90% band is lost to the **coverage clause, not to seeding**: the representatives align at 0.86–1.00
nucleotide identity but over 8–44% of the shorter one (NPIPB6–B8 0.858/0.39, NPIPB7–B8 0.997/0.31, NPIPA2–A9
0.882/0.28, CHST5–CHST6 0.863/0.42, ZNF747–ZNF764 0.886/0.15). Two multi-exon paralogues whose expressed
representatives are different isoforms, or which share only part of their exons, never reach 50% of the shorter
sequence in one colinear record. The "no record" misses at ≥ 90% are the same phenomenon at its limit (NPIPB3–B5,
100% protein identity, representatives overlapping 3%); genuine seeding losses begin at 60–70% (13 pairs, e.g.
MT1X–MT2A at ~70% nucleotide identity over 400 bp, where k = 11 finds nothing).

**Post hoc, not part of the bar:** with the coverage clause computed as the union of a pair's records instead of
the best record (the change discussed for rearranged copies), and its floor at 0.30, the 80–90% band goes
0.25 → 0.75 and 70–80% 0.36 → 0.50 at precision 0.869 (vs 0.968 for the shipped clause, 84 vs 63 judgeable
pairs). That is a real lever and a real trade; it needs its own pre-registration on the shipped catalog (both
substrates, the register's family-level metrics), not a switch here.

## What this settles for the advisor

The loss is not "we only keep the most similar copies": nucleotide precision is 0.96–1.00 down to the sensitive
tier's floor, and the sensitive tier's additions are correct. What the pipeline loses at 60–90% protein identity is
(i) paralogues whose expressed transcripts share too little of each other — a coverage question, addressable and
now measured — and (ii) below ~70% nucleotide identity, seeds. Below 50% protein identity the protein tier is the
only tier that sees anything and Compara's own precision there is 0.3, because at that depth membership is a
gene-tree question. The three-regime framing stands, with one correction: the boundary of the first regime is
set by exon sharing as much as by identity.
