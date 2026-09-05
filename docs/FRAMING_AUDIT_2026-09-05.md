# Framing audit — does "gene family, not duplication block" serve the objectives and answer the advisor? (2026-09-05)

Checked against `THESIS_OBJECTIVES.md` (the three objectives), `ADVISOR_QUESTIONS.md` (Q1–Q11, the Part-1
cross-examination, the Part-3 concessions), `reference_advisor_canzar` (clean combinatorial structure, provable
statements, no arbitrary thresholds, suspicious of merging without a similarity threshold, 1/k is bad), and the
chapter draft `CHAPTER_O1_O2_COMPOSITION_DRAFT.md` after §6eu.

## 1. The decision itself — holds
The chapter defines the **gene family** (loci sharing a transcribed core) and says why it did not build the
**duplication block** (SEDEF's object). This is the object O1 names ("multi-copy gene family") and the object the
advisor asked for (interest 1: "find multi-copy copies at the RNA level"). The certificate that the catalog can
tell the two apart is real and measured: under representative-only locus semantics NPIP (LCR16a, 29 loci) and the
SMG1P/LCR16u cluster (MCL7) are separate families; under attribution semantics they merge into one 48-locus block
(§6ef/§6eg). **Keep the decision. Put the 29-vs-48 certificate in §1, not only in the ledger.**

## 2. Where the framing does NOT yet meet the objectives — four gaps

**G1. Two O1 definitions coexist, and the advisor will find both.** `THESIS_OBJECTIVES.md` states O1 as
"quasi-clique in E_r; MCC = χ(H)" and `O1_O2_COMPOSITION.md` §(V, E_r) carries the γ-quasi-clique definition
with its theorem (P1) and NP-hardness note; the chapter's O1 is MCL over annotation homology + SD-core
refinement + read-supported units. Neither document says which is THE method or how they relate. For an
examiner who wants a clean combinatorial object, "Family = MCL over that graph" (chapter §2) is the weakest
sentence in the draft: MCL is a flow heuristic with an inflation parameter, so the family is defined as an
algorithm's output, not as an object with a property. **Required:** state the family as an object —
*a set of loci each of which shares a core segment (SD-linked) with at least half of the set* — and present
MCL as the pre-clustering that proposes candidates for that definition (the refinement is what certifies
membership; the drop arm's failure on patchy-SD families, §6ei, is the honest limit). Then retire or
subordinate the E_r/γ text in `THESIS_OBJECTIVES.md` and `O1_O2_COMPOSITION.md` — a thesis edit, the user's call.

**G2. "At the RNA level" is true of the node, not of the edge.** O1 says the family is defined
*topologically at the RNA level*; in the chapter, DNA (asm20 homology of annotated gene sequences, SD partner
depth) defines membership and RNA defines the unit's exon chain. That is defensible ("DNA proposes, RNA
disposes", §6dv: RNA cannot dispose membership — corroboration anti-correlates with cohesion) but the
objective's wording must change to match, or the chapter will be read as not delivering O1 as stated.
**Required:** restate O1 as "define the family's loci by shared duplicated core and its nodes by transcription",
and say in §2 that RNA was tested as a membership criterion five times and failed (register rows, §6ea).

**G3. The thresholds the advisor is suspicious of are not all justified.** Justified: inflation (2.0–4.0
invariant with a size-safe prune, §6ec), prune (size-dependence proven, §6ec), the exonic floor (1 bp is the
distribution's wall, §6dt), the core rule ("half", one constant), unit support (≥3 reads = the corroboration
floor already used). **Not justified anywhere in the MCL era:** identity ≥ 0.70, coverage-of-longer ≥ 0.30,
≥ 300 bp, min family size 3, the 50-kb giant-intron cut, and O2's PSV constants (2 reads/allele, 4 coverage,
τ 6.9, ε 0.003). `ADVISOR_QUESTIONS.md` already concedes ~25 default-reachable constants. **Required before the
chapter is shown:** a sensitivity table for 0.70/0.30/300 on the anchored families (NPIP, the tandem array,
the 3 artefacts) — the same instrument used for inflation. If the anchors are invariant over a band, say the
band; if not, the threshold is a decision and must be stated as one.

**G4. O2's certificate is conditional, and §6eu showed two ways the condition failed.** The chapter says a read
is assigned when `min_p < α/(n−1)`; the advisor will ask whether that p is an error rate. It is an error rate
*given the column set is the true PSV set*. §6eu found the column set wrong twice — the read-support filter
deleted the columns of unexpressed paralogues (11 confident wrong calls at p 6.7e-9) and a fragment unit let a
secondary record carry the decision (190 wrong calls at p ≤ 1e-133) — and both were caught only by the
placement-agreement audit, not by the certificate. **Required:** §4 must say the certificate is conditional
on the columns, name the audit that checks the condition (MAPQ-60 placement agreement, a floor and not a
score — Q8's rule stands: never claim O2 beats minimap2), and list D3 (secondary records as observations,
row 691) as open. Q5's shared-read question and Canzar's own multimapper frame land exactly on D3; it needs
its pre-registered fix before the defense, not after.

## 3. Standing questions — where the framing leaves each
| question | status under the gene-family framing |
|---|---|
| Q1 method for families | answered, once G1 states the object rather than the algorithm |
| Q2 real or artefact | answered: SEDEF corroboration 39/46, repeat library kills 3, MCL4 element-welded, CHM13 landing (§6dy–§6eh); the wrong-end concession (identity 0.83 band) stands |
| Q3 borrowing | unchanged: none; E_c splits only |
| Q4 isoforms | improved: the unit is the read-supported chain within the annotated locus (§6eu), no longer one representative intron chain; say so |
| Q5 shared reads | **open through D3** (row 691) |
| Q6 portability | unchanged: apes not reproducible; the MCL catalog has not been run on a second species |
| Q7 boundaries | answered (§6ay/§6ba); the unit-extent fix adds that boundaries now follow reads inside the annotation |
| Q8 1/k | answered; the sweep's 99.4 % is agreement with minimap2 on MAPQ-60 reads, a floor |
| Q9 NPIPA/NPIPB | **not addressed in the MCL catalog** — the chapter reports NPIP as one 29-locus family; the A/B split he asked for must be shown or its absence explained by the core rule |
| Q10/Q11 junctions, PSV | unchanged; PSV credibility now carries §6eu's two defects and the fix status |

## 4. What to do, in order
1. Restate O1 (object, not algorithm) in `THESIS_OBJECTIVES.md` and chapter §2; subordinate the E_r/γ text — **decision needed**.
2. Sensitivity table for identity/coverage/length on the anchors (G3) — one afternoon, the inflation instrument.
3. NPIPA/NPIPB under the core rule (Q9) — one run on the 29 loci with the CHM13 landing labels already in `all46_chm13_labels.json`.
4. D3 pre-registration (Q5, Canzar's frame) — score each molecule's sequence once against every copy.
5. Decide the two default flips from §6eu (`--units-follow-reads`, PSV read filter); the paired table is the evidence.
