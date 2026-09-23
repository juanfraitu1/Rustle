# Pre-registration — the FALSE-MERGE measurement that blocks r906's adoption

**Written 2026-09-22, §6y4, before any number is produced.** User: *"ok lets finish that."*

## The one thing standing between r906 and a shippable sensitive mode

r906 (§6t0) measured protein edges from annotated CDS, §6ko's rule unchanged, and found on **held-out**
chr2/chr8/chr10: within-family pair coverage **80.0% → 93.8% (+13.8 pts)**, no-edge families **halved
14.8% → 7.4%**, **0 cross-family protein edges on all three**, pooled precision **0.919**. It was recorded
⚠PARTIAL and NOT adopted, blocked on exactly three things:

1. **its own pre-registration** — this file, plus r906's own `PREREG_protein_edges_2026-09-20.md`;
2. **a §6bt-style false-merge measurement** — ⚠**the actual blocker, and the subject of this file**;
3. **whether a protein-space edge belongs in a definition called topological** — ⭐**dissolved by the
   user's sensitive-mode framing (§6y3): an opt-in mode is not the definition.**

⚠**Why (2) is not already answered.** r906's own text: *"Precision measured only over Soto-labelled genes,
so this is NOT a genome-wide false-merge rate."* Soto labels a minority of genes; every protein edge with
an unlabelled endpoint is currently **unmeasured**, and a sensitive mode's whole risk lives there.
§6bt.1's analogue was **2/150 = 1.33%** [0.37, 4.73] on single-locus windows at `RUSTLE_GATE_MIN_READS=3`
(**3/150 = 2.00% at the shipped node floor 2** — ⚠quote the floor with the rate). That panel is not on
disk, so an equivalent is constructed below rather than re-run.

## The measurement

Protein edges by §6ko's rule, unchanged, from `bench/protein_edge_gap.py` (blastp `-evalue 1e-5`, edge iff
non-overlapping HSPs cover **≥0.30 of the longer protein**; CDS translated in-house). Then every
**PROTEIN-ONLY** edge — one with no nucleotide edge in the shipped graph, i.e. exactly what the sensitive
mode would ADD — is classified:

| class | definition |
|---|---|
| **TRUE** | both endpoints Soto-labelled, same family |
| **FALSE** | both endpoints Soto-labelled, different families |
| **UNKNOWN** | at least one endpoint unlabelled — **the population r906 could not see** |

and each UNKNOWN is further split by **structural corroboration**, which uses no label: does the pair have
**any** nucleotide PAF record at all, even one that failed the gate?
- **corroborated** — a sub-threshold alignment exists, so the protein edge agrees with weak nucleotide evidence;
- **protein-only, uncorroborated** — literally no alignment. **This is the worst case and the number that decides.**

Reported: the false-merge rate bounded two ways — **optimistic** `FALSE/(TRUE+FALSE)` (r906's view) and
**pessimistic** `(FALSE + uncorroborated UNKNOWN)/all protein-only edges`.

## Substrates and the bar

Development **chr16**; held out **chr2 / chr8 / chr10**, no re-tuning, verdict taken there.

| outcome | verdict |
|---|---|
| held-out **pessimistic** false-merge ≤ **2.00%** (§6bt.1's rate at the shipped node floor) | ⭐ **ADOPT AS AN OPT-IN SENSITIVE MODE** — r906's blocker is cleared |
| pessimistic ≤ 10% and **optimistic** ≤ 2.00% | ⚠ **PARTIAL** — ship behind the flag, document the unmeasured tail |
| pessimistic > 10%, **or** any held-out cross-family (FALSE) edge appears | ⛔ **NO** — r906 stays unadopted |

**Predicted, before looking — ⚠ PARTIAL.** r906 already measured 0 cross-family edges on all three
held-out chromosomes, so I expect the optimistic rate to be 0.00% and the decision to rest entirely on the
uncorroborated-UNKNOWN fraction. Proteins are more sensitive than nucleotide by design, so a substantial
share of protein-only edges having NO nucleotide alignment is the expected outcome, not a surprise — the
question is whether it is a tail or the bulk. ⚠If the pessimistic rate is large, that is not automatically
a refutation: it may be the sensitive mode working as intended on genuinely divergent paralogues. I commit
in advance to reporting it as **unmeasured**, not as **false**, and to letting the ⚠PARTIAL row stand
rather than reinterpreting the bar.

I will not change the rule, the classes, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. r906's 0.919 precision was a restricted-universe artefact; the mode adds DOMAIN SHARERS.**

§6ko's rule unchanged, blastp `-evalue 1e-5`, edge iff non-overlapping HSPs ≥0.30 of the longer protein.

| | chr16 (dev) | chr2 | chr8 | chr10 |
|---|---|---|---|---|
| protein edges | 1,171 | 2,927 | 1,741 | 394 |
| already nucleotide | 162 | 275 | 1,323 | 26 |
| ⭐**PROTEIN-ONLY (what the mode ADDS)** | **1,009** | **2,652** | **418** | **368** |
| TRUE (both labelled, same family) | 26 | 1 | 0 | 2 |
| FALSE (both labelled, diff family) | 31 | 1 | 3 | 0 |
| UNKNOWN, sub-threshold alignment | 62 | 65 | 34 | 35 |
| **UNKNOWN, NO alignment at all** | **890** | **2,585** | **381** | **331** |
| optimistic false-merge `F/(T+F)` | 0.544 | 0.500 | **1.000** | 0.000 |
| ⭐**pessimistic `(F+bare)/all`** | **0.913** | **0.975** | **0.919** | **0.900** |

⛔**Both ⛔ conditions fire on held-out**: pessimistic false-merge is **90.0–97.5%** against a >10% ⛔
threshold, and **FALSE edges appear** on chr2 and chr8. Against §6bt.1's 2.00% this is not close.

## ⭐⭐ What the mode actually adds: domain sharing, not paralogy

Of the protein-only, uncorroborated edges, the two genes carry **different symbol roots** in
**54.3% (chr16) / 85.3% (chr2) / 64.6% (chr8) / 71.9% (chr10)** of cases. ⚠Symbol root is used here only
as a DIAGNOSTIC — r902 voided it as truth — but the examples are unambiguous: `ADRA2A~HTR7`,
`ADRA2A~NPY4R`, `ADRA2A~NPFFR1` (GPCRs sharing a 7-transmembrane domain), `AAMP~CIAO1`, `AAMP~WDSUB1`
(WD40), `ACTA2~ACTR1A` (actin fold), `ABCC1~ABCC6`, `ADAM18~ADAMDEC1`.

⭐**This is §6ko's rule behaving exactly as written**: a shared domain easily covers ≥0.30 of the longer
protein, so the rule admits domain homology. That is real homology — it is simply not *gene family* in the
sense O1 defines, and it is what 90%+ of the mode's additions are.

## ⚠⚠ The correction to r906, and it is register 770 again

r906 reported pooled precision **0.919** and **0 cross-family edges on all three held-out chromosomes**.
Both are reproduced here **within the Soto-labelled universe** — chr10's optimistic rate is literally
0.000, and the labelled denominators are 1–3 pairs. **But Soto labels only recent, ≥98%-identity
segmental-duplication families**, so restricting precision to them excludes the domain-sharing bulk **by
construction**. r906 flagged this itself (*"NOT a genome-wide false-merge rate"*) and could not measure it;
measured, the unlabelled remainder is **90–97.5% of everything the mode adds**.
⭐**r906's 0.919 is not wrong, it is scoped** — and the scope was doing all the work.

## ⚠ Neither available truth can score this arm

- **protein referee** — circular by construction (§6y3): it is itself built by clustering translated CDS.
- **Soto** — SD-scoped: it calls ancient paralogues *different families*, so protein's genuine
  ancient-paralogue finds (`ABCC1~ABCC6`) are counted as false merges.

**No truth on hand separates "ancient paralogue" from "domain sharer"**, which is the distinction this arm
turns on. That is why the measurement is reported as a **bound** (optimistic/pessimistic) rather than a rate.

## Verdict and the one untested lever

⛔ **r906 stays unadopted.** Its three blockers: the prereg is now written, the third is dissolved by the
sensitive-mode framing — and **the false-merge measurement, the one that mattered, fails.**

⚠**Untested, and deliberately not tried here** (the prereg forbids re-tuning): §6ko's **0.30 coverage
floor** is exactly what lets a domain qualify. A much higher floor (0.70–0.90 of the longer protein) would
demand whole-protein homology rather than a shared domain. That is a different rule needing its own
pre-registration, and it is the only version of this idea left standing.

## Prediction scorecard

Predicted ⚠ PARTIAL with "optimistic 0.00% and the decision resting on the uncorroborated fraction; the
question is whether it is a tail or the bulk." **The mechanism was right and the magnitude was not**: it is
the bulk (90–97.5%), and the optimistic rate was not 0.00% either (0.544/0.500/1.000/0.000 — on 1–3 pair
denominators). I also committed in advance to reporting the uncorroborated fraction as *unmeasured* rather
than *false*; the symbol-root diagnostic is what let me go further and say what it actually is.
