# Defense-Readiness Audit — will the objectives survive advisor (Canzar) scrutiny?

2026-06-27. 6 adversarial "Canzar-examiner" agents (one per objective + cross-cutting) + synthesis,
each grounded in the code/docs. Supersedes the attainment status of the 2026-06-25 loose-ends audit
(which predates the definitive O2 recompute and the O4 collapsed-copy wiring).

## Per-objective attainment

| Obj | Verdict | One-line |
|---|---|---|
| **O1** family copies | attained-with-caveats | Threshold-free conflict graph now runs GW (82 same-chrom families on disk) — but the **headline cross-chrom catalog (152–157) is a threshold pipeline** (Louvain + density-gate 0.30 + asm20 refine id≥0.80/cov≥0.50) with **no orthogonal precision number** (70.9% = annotation FLOOR; 90% DNA partly circular; real SEDEF/BISER undone). |
| **O2** assign-or-abstain | attained-with-caveats | Decision rule genuinely clean (assign/abstain/tied, no 1/k; α principled; 99.4% of "tied" carry a real min_p=1.0 K=0 certificate). BUT 75.1/24.8/0.0 is the **per-region CLI on the annotation-refined co-located subset, not the principled conflict-graph catalog** (~47/47 under its own gate). No non-circular accuracy on real reads. |
| **O3** ASJ | attained-with-caveats | Per-molecule, phasing-free linkage is real. ~~"120 genetic ASJ" conflates not-edit with not-copy; 44/120 at paralog-suspect LOC*; masquerade separator never run.~~ **RESOLVED (P4, 2026-06-28):** masquerade separator WAS run on the 18 LOC\* windows → 17/18 copy-confounded, excluded → **defensible genetic core ~77** (not 120). **CORRECTED (M1, 2026-06-28):** the PSMD2/DAXX "on the canonical dinucleotide / creates-destroys the motif" claim was genome-FALSE (anchors at donor−1 / exon boundary; 0/475 on a core dinucleotide, GT-AG intact) — retracted; mechanism = splice-REGION variant, genome-pinned in `bench/asj_motif_check.py`. Remaining caveat: allele-vs-copy at the excluded LOC\* loci still needs DNA. |
| **O4** reference-absent | **mechanism-only** | Clean architecture (same VG-path object; two-stage freeze byte-identical superset; DNA-needs first-class; fails-safe). But **ZERO real GGO copies admitted**; positive demo is synthetic tiny-intron; **gate-5 asm20 can't admit real multi-kb-intron copies**; 4 MHC = old-screen candidates in the worst het regime; raw flag 7.4% ≈ background. |
| **O5** multimappers→assembly | not-attained (by design) | Explicitly future; no claim to defend. |
| **THEORY** | attained-with-caveats | The jewel: Lemma 1 (MCC=χ(H)), Thm1 NP-hard, the K≥3 recombination obstruction, disjoint-clique self-certifier — proved AND exhaustively machine-checked (full L=3 universe), honest that Strong-Sep is sufficient-not-necessary. BUT **RECOVER/self-certifier/coloring are NOWHERE in production**; the running gate touches the theory only at the K=0 vertex. The deep theorems guarantee an algorithm you don't run. |

## Will it pass the advisor?
**Yes — with honest reframing, NOT as currently headlined.** Defensible thesis, one genuinely original
theoretical contribution, unusually disciplined self-auditing. NOT defensible if the empirical headlines are
asserted as precision/accuracy.

**Lead with (defense-grade):** (1) the identifiability THEORY (non-circular by construction); (2) the O2
decision rule as a structure (assign-or-abstain + per-read certificate + K=0→Tied certificate; frame the
24.8% as a *contribution*); (3) O3 PSMD2/DAXX exemplars (~20, needs no external catalog); (4) the honesty
discipline itself (the committee credits the self-found over-claims).

**Hardest questions (prepare, don't dodge):**
- "Which of your THREE family catalogs produces the O2 headline, and why the threshold-pipeline one not your
  principled conflict graph?" — **the killer; the build-vs-run gap moved, it did not close.**
- "What is the FDR of your cross-chrom catalog?" — no non-circular answer exists.
- "Show me one real gorilla reference-absent copy." — there are none.
- "Your deep theorems describe RECOVER; grep shows it isn't in the code — why credit Thm 2/3 as load-bearing?"

**Must be retracted/relabelled (won't survive an external check as-is):** "70.9% orthogonally confirmed"
(→ FLOOR); "89–90% DNA-confirmed" as orthogonal (→ lower bound); "silver 99.9%" as accuracy (circular —
agreement with minimap2's own placement); "0 false merges" as precision (reuses own homology rule); "75.1%
definitive O2" as genome-wide (→ co-located annotation-refined subset); "4 confirmed MHC ref-absent copies"
as O4 deliverables; "GW families WITHOUT arbitrary thresholds" (true only for the 82-family raw object, not
the headline one).

## What's missing — prioritized
**MUST DO:**
- **P0 — the one external O1 check: a SEDEF/BISER segmental-duplication map (build SEDEF from source / BISER
  off-WSL2).** Highest leverage; converts O1 from self-referential to falsifiable; the parCN standard your own
  Soto-2025 ref names as gold. **Resolves the biggest risk.**
- **P1 — reconcile the three catalogs; run O2 on the PRINCIPLED one.** Report the conflict-graph-catalog
  number as the GW headline (honestly ~47/47), or derive why the exon-sum refined catalog is the principled
  substrate. Your elegant artifact and your headline number are currently different objects — pick one story.
- **P2 — fix O4 gate-5 (`asm20`→`-x splice`) before ANY real-data O4 claim,** then attempt one real admitted
  copy (zero is itself a reportable measured result).
- **P3 — one non-circular O2 accuracy point on real reads** (pin sim5x in CI with acc|assigned≥0.99 + K0→100%
  tied; find a locus with truly external per-read labels). Reconcile the contradictory sim5x tables (20% vs
  100% resolvable for K≥2 — one is wrong).
- **P4 — run the O3 masquerade separator** (`scan_gene_copy_specific_junctions`) on the 44 LOC* calls; report
  within-gene-het vs paralog-locus separately; recompute on the corrected GGO_mm.bam.

**HONEST REFRAMING SUFFICES (relabel, no new compute):** attach FLOOR/lower-bound/co-located-subset/candidate
to every headline; delete the stale denovo_families.tsv (T_CORE=0.13 over-merge); split the min_p tail
(only ==1.0 is an impossibility *certificate*; α≤min_p<1.0 is power-limited abstention — different label);
headline O3 ~20 not 120; state O4 as "detect-and-flag, zero real copies admitted, copy-vs-allele needs DNA";
drop theory's "dichotomy fully closed" / "Thm1-3 executable in production" (certifier is Python-only).

**NICE-TO-HAVE:** a **bridge theorem** (prove min_p≥α is a sound certificate for non-Strongly-Separated /
recombination positions, so the running gate inherits a combinatorial guarantee for K≥1 not just K=0 — the
cleanest way to make the theory load-bearing without rewriting production); run RECOVER+self-certifier on
real/sim5x and report the Strong-Sep / recombination-free / rejected fractions; justify-or-remove the magic
numbers (edit_rate 0.2 wholly underived, min_clusters 3, remap 0.98, coverage 0.50, merge 0.30) with stability
sweeps.

## THE single biggest risk
**Circularity.** Every empirical accuracy headline is self-consistency, not an external check: silver =
minimap2's own placement; DNA = the span containing the building exons; 70.9% = an annotation floor rechecked
with the catalog's own homology rule; sim5x = your own generative model. The one orthogonal check your field
uses — a SEDEF/BISER segdup map — is the one not run. A committee cannot separate your method's accuracy from
minimap2's, nor falsify your catalog.

**Neutralize, in order:** (1) **build the external check (P0)** — even a partial segdup map on a subset of
chromosomes breaks the circularity charge; (2) **if P0 can't land before the defense, win by inversion** —
lead with the THEORY (non-circular by construction) and present every empirical result explicitly as
bounded/honestly-labelled/falsifiable-in-principle (floors not precisions, certificates not guesses,
detect-and-flag not fabrication). A committee forgives a *named* open external check; it does not forgive a
circular number *asserted as precision*. The deepest fix that serves both O2 and Theory: **make the principled
artifacts the load-bearing ones** (run O2 on the conflict-graph catalog; run RECOVER or prove the bridge
theorem) so your jewels and your numbers describe the same pipeline.
