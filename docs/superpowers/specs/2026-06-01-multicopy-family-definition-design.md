# Design: A formal, operational definition of a multi-copy gene family in the VG

Status: DESIGN (2026-06-01). Approved direction: **C** — state the definition rigorously
*and* make it exactly what the tool computes, validated on a labeled gallery.

Grounding: the genome-wide paralog secondary-dependence scan (GGO gorilla IsoSeq vs RefSeq GFF;
846 expressed multi-copy candidates → edit-distance taxonomy) and a validation against 14 curated known
families. Tooling: `bench/paralog_secondary_scan/`. Memory: `project_paralog_secondary_scan`, `project_vg_wiring`.

---

## 1. Problem & motivation

The thesis assembles multi-copy gene families with a variation graph (`--vg`), doing what StringTie cannot.
But "multi-copy gene family" has had no operational definition — progress was anecdotal (DAZ recovered,
GOLGA6L7 failed). The advisor's bar: a definition must be **non-trivial** (not "a single SNP separates the
copies") and **honest about the hard case** (when copies are *not* separable). We need one object that is
both a rigorous definition and a predicate the tool computes, with examples and counterexamples.

Two empirical anchors already establish the shape:
- **The scan**: of 846 expressed multi-copy candidates, an edit-distance copy-anchor finds genuine expressed
  copies recoverable only via secondary alignments — alongside a large class of *pure spillover* loci that
  naive secondary use would fabricate. Distinguishing the two is the contribution.
- **Distinguishing positions are not a given**: 14% of secondary-dependent reads carry **no** distinguishing
  position; 5 families are entirely non-identifiable (TSPY-like). The definition must make this regime
  first-class, not assume separability.

## 2. The object: the VG as a splice-graph generalization

A standard **splice graph** for one gene is a DAG `G = (V,E)` with `V` = exons, `E` = junctions; transcripts
are source→sink paths (paths = isoforms of one gene).

A **family splice graph** generalizes this to paralogous copies: `G_F = (V, E, σ, ℓ)` where
- `V` = **exon-classes**; each node `v` carries **per-copy sequences** `σ_v : C_v → Σ*` over the copies
  `C_v` that traverse it. `v` is **shared** if `|C_v| ≥ 2` (the homologous backbone), **copy-specific** if
  `|C_v| = 1`. (This is exactly `ExonClass{per_copy_sequences, copy_specific}` in `vg_hmm/family_graph.rs`.)
- `E` = junction edges (`JunctionEdge`).
- `ℓ` labels each maximal path by **(copy, isoform)** — so a path encodes *which copy* and *which isoform*;
  the copy labels partition the path set. A copy `c` is the sub-splice-graph `recover_paralog_path(c)`.

Reads multi-map across the shared nodes — that shared backbone is what makes the loci a *family* and what
makes per-copy assignment hard.

## 3. The definition

Let `H_F` be the **read-multi-mapping graph**: vertices = candidate loci, an edge `(a,b)` when ≥1 read
chain-matches both loci. (This is bottom-up from the reads, *not* from annotation names — see §6.)

> A **multi-copy gene family** (in the VG) is a connected component of `H_F` that forms a family splice graph
> satisfying:
> - **(M) multiplicity** — ≥2 copies (member loci / path-bundles);
> - **(H) homology / connectivity** — the copies share ≥1 exon-class, i.e. the component is connected by
>   multi-mapping reads. *Primary structural gate.*
> - **(X) expression** — ≥2 copies are *independently expressed*: each has ≥`k` reads **anchored** to it
>   (default `k=3`);
> - **(scope)** each copy is **spliced** (≥1 intron). Single-exon paralog loci are out of scope (§7).
>
> Each family additionally carries its **identifiability partition** (§4) and **genomic arrangement** (§5).

### Anchoring (the non-trivial core — edit distance, not a single SNP)
For a shared node `v` and copies `c,c' ∈ C_v`, a position is **(c,c')-distinguishing** if
`σ_v(c) ≠ σ_v(c')` there. For a read `r` with edit distances `NM_c(r)` to its comparable-extent placements,
`r` is **c-anchored over c'** iff `dNM = NM_{c'}(r) − NM_c(r) ≥ T` (default `T=2`). This is the
fingerprint-EM log-likelihood ratio: each distinguishing mismatch is a near-independent observation, so the
discriminant is the **raw count of distinguishing mismatches**, *not* the per-base rate (NM 14 vs 26 = 12
distinguishing positions = decisive, even though the rate gap is <1%). `|dNM| < T` ⇒ the read covers no
distinguishing position ⇒ **ambiguous** (the boundary-theorem case).

## 4. Identifiability as a graded property (the hard regime, first-class)

Distinguishing positions are **per copy-pair, per-position** — not global. Define the **identifiability
relation**: copies `c,c'` are *read-identifiable* iff some read is `c`- or `c'`-anchored over the other.
Non-identifiable pairs form **equivalence classes** (transitive closure). A family's **identifiability
partition** is these classes, summarized as:
- **full** — all copies pairwise identifiable (tie-fraction < 0.15);
- **partial** — some classes resolve, others don't (DAZ: 11% ties);
- **none** — one class; copies indistinguishable over the expressed region (TSPY: 97% ties).

For a non-identifiable class the tool **reports the combined (aggregate) abundance and flags the copies as
unresolvable, with the supporting counts**, rather than inventing a per-copy split. The honest answer there is
a *range* of per-copy splits all equally consistent with the reads — for two tied copies, any split from 0 to
the total (the non-identifiability analysis of Hjorleifsson–Pachter) — so a single number would be arbitrary. The
**boundary theorem** states when read-length × divergence yields any distinguishing position; the
**structural-linkage channel** (copy-specific junctions) is how we push the boundary — and is only available
because the copies are spliced (§7).

Identifiability is a *reported property*, not a membership gate: a fully non-identifiable but expressed,
connected cluster (TSPY) is still a **family** — labeled `FAMILY, non-identifiable`.

## 5. Genomic arrangement: segmental vs tandem

Arrangement is a property of where the copies sit, **orthogonal to membership** (the definition is
connectivity-based and location-agnostic). The scan confirms genuine families span all arrangements
(expressed families: distal 40, trans 19, tandem 19). Arrangement is **recorded per family** because it
affects two things:
- **Graph construction.** *Tandem* copies (adjacent/same locus) can fall in one assembly bundle → risk of
  merge/collapse → need intra-bundle splitting (the GOLGA6L7-type blocker). *Segmental* copies (distal/trans)
  sit in separate bundles → linked only by multi-mapping **secondary** alignments → the secondary-dependence
  the scan measures.
- **Identifiability covariate.** Recent tandem arrays skew more non-identifiable (median tie-fraction 0.10 for
  tandem vs ~0.00 for segmental) — a tendency, not a rule (DAZ is tandem+inverted yet partially identifiable;
  NBPF is segmental and identifiable).

`locus_rel ∈ {tandem, distal, trans, overlapping}` is emitted on each family.

## 6. Membership is connectivity-based, not annotation-based

Assembling families by gene **name / protein domain** is wrong: ZNF-by-name = 586 "members" (they share the
zinc-finger *domain*, not orthology); GOLGA mixes distinct golgins with the GOLGA6/8 duplicons; PRAME members
don't multi-map in the data (H≈0). The family is therefore a **connected component of the read-multi-mapping
graph** — exactly what `FamilyGraph` construction already builds (position clustering + minimizer-Jaccard
refinement in `family_graph.rs`). Annotation names/descriptions are used **only as a validation lens**, never
to define membership.

## 7. Scope & out-of-scope

**In scope:** spliced (multi-exon) multi-copy families. **Out of scope:** single-exon paralog loci.
Data justification (overlap-based exploration of 281 single-exon multi-mapping candidates,
`bench/paralog_secondary_scan/se_explore.py`): 158/281 are **pure spillover** (processed pseudogenes of
highly-expressed housekeeping genes — intronless retrocopies), 19 non-identifiable (histones, U6 snRNA), and
the 44 "expressed" are mostly retrogenes the parent dwarfs (G≈3 vs sibling 100s). Including them would import
the fabrication-risk class, not genuine families — and the structural-linkage channel needs **junctions**,
which single-exon copies lack. Genuine intronless families (e.g. olfactory receptors) are a known sequence-only
special case, noted but not handled by the core.

## 8. Examples & counterexamples gallery (validated)

| locus | M | arrangement | verdict | identifiability |
|---|---|---|---|---|
| NBPF | 23 | segmental | **FAMILY ✓** (16 expressed) | full |
| DAZ | 3 | tandem, inverted | **FAMILY ✓** | partial (11% ties) |
| RBMY | 13 | tandem (Y) | **FAMILY ✓** | partial |
| amylase, tubulin-β | — | segmental | **FAMILY ✓** | full |
| **TSPY** | 12 | tandem array | **FAMILY, non-identifiable** | **none (97% ties)** |
| **MAGEA** | 14 | segmental (Xq28) | FAMILY (2 resolvable) | none (89% ties) |
| SORD + LOC | 2 | tandem | **counterexample: spillover** (0 reads own the LOC) | — |
| ribosomal/COX processed pseudogenes | many | dispersed | counterexample: spillover (single-exon, §7) | — |
| ZNF-by-name (586) | — | dispersed | **artifact of name-assembly** → use connectivity (§6) | — |
| β-defensin, protocadherin, PRAME | — | — | **correctly excluded**: not expressed in this tissue (X fails) | — |

## 9. Operationalization (the predicate the tool computes)

`classify_family(component, reads) -> FamilyVerdict` where `component` is a connected component of `H_F`
(from `FamilyGraph`) and the result carries:
- `verdict ∈ {family, family_nonidentifiable, gene_plus_unexpressed_paralog, spillover, not_connected,
  not_expressed_here, single_exon_out_of_scope}`;
- `n_copies (M)`, `n_expressed (X)`, `connectivity (H)`;
- `identifiability_partition` (the equivalence classes) + `identifiability ∈ {full, partial, none}`;
- per non-identifiable class: `aggregate_abundance` + an explicit `unresolvable` flag with its evidence;
- `locus_rel` (arrangement).

Inputs already in hand: `FamilyGraph` (ExonClass nodes, copy paths), the per-read multi-mapping placements
with NM (`mm_placements`), the gene/exon-class spans. The formal definition (§3–§4) **is** this function's
contract; no separate notion.

## 10. Parameters

`T` (distinguishing-mismatch threshold) = 2 — robust (scan headline stable across T=2/3/5); `k` (reads to call
a copy expressed) = 3; `MIN_SHARED` for H from existing `--vg-family-min-shared`; tie-fraction cutoffs
full<0.15 / none≥0.6. All env-overridable; calibrated on the gallery.

## 11. Validation gate

1. **Known-family panel** (`validate_known_families.py`): DAZ/NBPF/RBMY/amylase/tubulin ⇒ `FAMILY ✓`;
   TSPY/MAGEA ⇒ `FAMILY, non-identifiable`; SORD+LOC ⇒ spillover; DEFB/PCDH/PRAME ⇒ not-expressed/not-connected;
   ZNF/GOLGA-by-name ⇒ flagged as name-assembly artifacts (resolved by connectivity).
2. **Scan gallery**: the edit-distance taxonomy (`expressed_real_copy` vs `pure_spillover` vs
   `ambiguous_dominated`) reproduces, with DAZ3 calibrating as weak/partial.
3. **Single-exon exclusion**: documented via `se_explore.py` (158/281 spillover).
4. **Determinism**: same inputs → same partition.

## 12. Relation to prior art

Generalizes the **splice graph** (transcript assembly) to paralogous copies; the **range of per-copy
abundances consistent with the reads** for non-identifiable classes is the non-identifiability analysis of
Hjorleifsson–Pachter (we report that range / abstain rather than pick one point); the apportionment of
ambiguous reads is the **multimapping-read assignment** problem (the advisor's domain) solved by edit-distance
likelihood + structural linkage rather than a single diagnostic SNP.

## 13. Risks / open questions

- **`H_F` construction cost** genome-wide (multi-mapping graph) — mitigated: built per candidate component, not
  globally; `FamilyGraph` already clusters.
- **Tissue dependence of X** — a real family unexpressed here is correctly "not-expressed-here"; the definition
  is about the *reconstructable* family. State this explicitly to the advisor.
- **`T=2` on indel-rich regions** — guard with comparable-extent; revisit if HiFi indels inflate dNM.
- **Within-family vs genome-wide anchor** — identifiability is computed *within* the family (vs sibling copies);
  expression (X) is vs all copies genome-wide. Keep the two anchors distinct (already the case in the tooling).
