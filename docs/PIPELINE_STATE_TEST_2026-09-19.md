# Does the pipeline, as it stands today, fulfil O1, O2 and the advisor's standing questions?

Run 2026-09-19 at `1cebbb01`. **Everything marked TESTED was executed today against the current binary**;
everything marked CARRIED is a ledger number not re-derived here, with its section named so he can chase it.
The dossier `docs/ADVISOR_QUESTIONS.md` was last revised 2026-09-16 and **predates §6o2–§6r5**, so several
of its headline sentences are stale; those are listed at the end.

## O1 — identify multi-copy gene families — ⭐ TESTED, PASSES end to end

Run from scratch: RefSeq CHM13 chr16 gene+pseudogene bodies (2,082) → `samtools faidx` → `minimap2 -x asm20
-c --eqx -P` all-vs-all (79,207 records, 69 s) → `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`
(the shipped default rule, §6kt/§6o2–§6o6).

| | |
|---|---|
| graph | 380 nodes, 627 edges (identity ≥ 0.7, cov_longer ≥ 0.3, ≥ 300 bp) |
| pairs dropped, no exonic evidence | 10,637 |
| pairs dropped, shared-exon fraction < 0.60 | 965 |
| **families (≥ 2 members)** | **90**, 271 members, largest 26 |
| size distribution | 56×2, 19×3, 6×4, 3×5, 2×6, 2×7, 1×11, 1×26 |

⭐**The advisor's own example comes back exactly: all 21 chr16 NPIP genes land in ONE cluster (MCL0),
none missing** — NPIPA1/2/5/6/7/8/9, NPIPB2–B15, NPIPB10P, NPIPB14P. MCL0's other 5 members are unnamed
`LOC` loci (LOC100505915, LOC124907807/808/834, LOC128966608), **not adjudicated here** — they may be
unnamed NPIP paralogs or genuine over-merge, and this run does not settle which.

Other clusters are biologically coherent without being tuned for: **SMG1 + SMG1P1…P7** (11), **PLA2G10EP…KP**
(7), and the **PKD1P** family (PKD1P2, PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P5-LOC105376752, PKD1P6-NPIPP1) —
note the PKD1P–NPIP readthroughs cluster **separately from NPIP**, which is the §6m1/§6q boundary story.

⚠**Two honest qualifications.**
1. This is the **guided/annotation-node mode**. The dossier's "NPIP lands in 4 clusters (25/3/2/1)" is about
   the **de novo read-node catalog** (§6hu) — both can be true, and today's run does **not** show the
   fragmentation is fixed. It shows the edge+grouping rule is sound when the nodes are right, which is
   exactly what §6kg said (annotation nodes F 0.955, de novo nodes 0.726).
2. **chr20 yields 0 families, correctly.** Its gene-body alignments have median cov_longer 0.0026 and p99
   0.286, below the 0.30 floor — chr20 simply has no gene-body-level paralog pairs at this stringency. An
   empty answer on a family-poor chromosome is the rule working.

⚠A trap worth recording: `mcl_families` keys the PAF on `chrom:start-end` (`parse_gene_key`). Feeding a PAF
whose records are named by **gene symbol** silently yields **0 nodes, 0 edges** with no error — my first run
did exactly that. Keep `samtools faidx`'s native names.

## O2 — assign copies under ambiguity — ⭐ TESTED, and the honest verdict is "abstains, correctly"

§6r4, today. O2's subject is AS-tied multimappers only. Measured on its best real substrate, the **Y
ampliconic genes** (guided catalog from the CHM13 annotation, 8 families / 30 copies, human testis A119b
chrY, 603,346 records):

| | gorilla autosomes (78 fam) | **human YAGs (8 fam)** |
|---|---|---|
| AS-tied share | 0.30% | **53.1%** |
| contested (O2's actual subject) | 21 | **3,641** |
| **assigned** | 2 | **12 (0.3%)** |
| tied (abstain) | 1 | **3,423 (94.0%)** |
| ambiguous | 18 | 206 (5.7%) |

⭐**The mechanism works and the discipline holds: 94% abstentions, never a 1/k split.** That is the direct
answer to Q8.
⚠⚠**But as a copy-resolution method on the hardest real case it resolves 0.3%.** Only DAZ (9) and RBMY (3)
yield any assignment; BPY2, CDY, HSFY, PRY, TSPY, VCY yield zero. **Do not present O2 as solving YAG copy
assignment.** The defensible claim is: it identifies the ambiguous population, refuses to guess, and says
which families are resolvable at all.
⚠**O2's decision set is bounded by catalog coverage, not aligner ambiguity** (register 876): 24,543 gorilla
molecules are AS-tied but only 91 land in ≥2 copies of a 509-copy catalog. Any O2 rate is a statement about
the catalog first.

## The standing questions — status today

| | question | status |
|---|---|---|
| Q1 | a method that identifies families? | ⭐**TESTED today**, chr16, 90 families, NPIP 21/21. ⚠dossier text is stale (describes triangle-leaders as opt-in; the shipped rule is MCL + `--min-shared-exon-frac 0.60`) |
| Q2 | real, or artifacts/overfitting? | CARRIED. Forking-paths accounting still the strongest card; ⚠the "~0.83 band" concession stands |
| Q3 | borrow information across the family? | CARRIED |
| Q4 | could two near-identical copies give the same isoforms? | ⭐**Answered by §6r4**: on YAGs 94% of contested molecules are *tied* — the copies are that similar |
| Q5 | tandem near-identical copies share a read? | CARRIED (§6aq; ⚠check junction count first — the discriminator is degenerate below ~3) |
| Q6 | ports to other families/tissues/apes? | ⭐**Strengthened**: §6q7 ran gorilla + human on the lab's own libraries; §6r4 adds chrY |
| Q7 | TSS/TES/UTR boundaries too convenient? | ⭐**Now measurable and NOT flattering**: §6r5 locus size ratio median 0.994 but **q25 0.676** — a quarter of loci are short, the §6p4 5′-truncation signature |
| Q8 | "1/k is a bad assignment" | ⭐**TESTED**: 94.0% abstention on 3,641 contested molecules, 0 × 1/k |
| Q9 | NPIPA/NPIPB distinct subfamilies | ⚠**Today's run puts A and B in one cluster (MCL0)** — that is L3; the subfamily split is L4's job (§6p5, selective at 0.995), untested in today's run |
| Q10/Q11 | non-canonical junctions; PSV credibility | CARRIED (§6m7/§6m8: the 209/249 "ceiling" IS the canonical-motif count) |

## What is stale in `docs/ADVISOR_QUESTIONS.md` (revise before the meeting)

1. **Q1's definition paragraph** — triangle-supported leaders are described as the current best and opt-in.
   The shipped DNA rule is **MCL with `--min-exonic-bp 1 --min-shared-exon-frac 0.60`** (§6o2–§6o6), and the
   RNA rule is **L3 components at `w_98 ≥ 0.985`** with **L4 at 0.995 applied selectively** (§6p1/§6p5).
2. **Q1's NPIP headline** — "4 clusters (25/3/2/1)" is the de novo catalog. State the mode; the guided mode
   returns 21/21 in one cluster.
3. **Nothing in the dossier mentions the assembler bake-off** — §6q7 (vs the lab's own StringTie/FLAIR/isoseq
   on A119b and GGO) and §6r5 (TPM ρ 0.879, locus size ratio) are new and are what he asked to run.
4. **Part 0d's IGV file** — `bench/locus_bed.py` now emits the loci BED directly (§6r5).
5. ⚠**The overall family metric to quote is F 0.628 on all 72 families, not 0.881 on 11 hand-picked ones**
   (§6o1, register in `feedback_metric_traps`).
