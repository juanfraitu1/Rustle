# Multi-copy gene-family definition — canonical spec

End-to-end, reproducible. Turns the over-merged cDNA all-vs-all into a layered, trustworthy family
set + a `rustle --family-manifest`-compatible manifest. Implementation: `bench/family_def_build.py`
(`python bench/family_def_build.py [--ancient]`). Every family edge is **whole-gene / whole-protein
homology, never a shared sub-domain** (enforced + audited, see §Sub-domain safeguard).

## Pipeline

| step | criterion | effect |
|---|---|---|
| 0. cDNA homology | id≥0.90 (rep_ava.tsv) | 1216 families, max comp **238** (over-merged) |
| 1. retrocopy filter | drop one-side-intronless edges | removes 1931 processed-pseudogene edges |
| 2. strand filter | drop antisense ('-') edges | removes 482 sense/antisense + opposite-repeat artifacts |
| 3a. **recent coding families** | cDNA edge ∩ whole-protein homology (below) | 826 families |
| 3b. **immune-receptor families** | segment biotype + cDNA + same-locus (≤3 Mb) | 30 single-locus tandem arrays |
| 3c. ancient domain families (`--ancient`) | protein-graph community detection | optional super-families (ZNF/OR/GPCR) |

**Result: 856 families, max comp 238→47** (a real orthogroup), **4/4 named reals preserved.**
Manifest: `/home/juanfra/winloci_scratch/family_manifest.tsv` (6-col, rustle-compatible, `#` header)
and `family_manifest.full.tsv` (+biotype, layer).

## Whole-protein (anti-sub-domain) homology bar — step 3a

A protein edge counts only if the homology spans most of BOTH proteins, never a fractional domain:

    fident ≥ 0.30  AND  min(qcov,tcov) ≥ 0.50  AND  max(qcov,tcov) ≥ 0.70  AND  e ≤ 1e-5

- `min-cov ≥ 0.50`: the LONGER protein is ≥50% covered.
- `max-cov ≥ 0.70`: the SHORTER protein is ≥70% covered — the binding anti-sub-domain term
  (excludes two large proteins sharing one ~60% domain).

### Sub-domain safeguard — audit
mmseqs all-vs-all on 22,614 representative proteins → measured on family-bar edges:
- min-cov median **0.89**; alignment spans **97%** of the shorter protein (median); 91% span ≥80%.
- only **4% were sub-domain-suspect** (shorter protein <70% spanned) → removed by the `max-cov≥0.70`
  term (−56 edges, −9 spurious families, **0 real families lost**, max comp unchanged).

So the family edges are provably whole-protein homology. (No Pfam/InterPro/hmmer installed; the
reciprocal-coverage bars are the safeguard. Domain-architecture matching would be the gold standard
if a domain annotator becomes available.)

## Why this is principled (not threshold-fitting)
- The cDNA id≥0.90 bar keeps families **recent** → ancient domain families (ZNF 569-blob, OR, GPCR)
  are split into their recent sub-duplications, not over-merged. This is the recent-duplication /
  copy-assignment regime the thesis targets.
- Coding genes → protein-validated (the orthogonal ground truth: co-threading recovers protein
  orthogroups at F1 0.72, but protein homology IS the truth and is used directly here).
- Non-coding lncRNA "families" (the bulk of the resistant comp-238 blob) are excluded by construction
  (no protein) — they are repeat-bridges, a separate phenomenon.
- Immune segments (no protein, single V-domain, tandem) get a biotype+locus criterion instead.

## Known boundaries (explicit, not silent)
- **Partial paralogs** (GGT1~GGTLC2): protein homology sub-bar (GGTLC2 covers part of GGT1) → not in
  recent_coding. The recurring partial-duplication boundary; needs a separate truncation-aware rule.
- **Ancient families** (ZNF/OR/GPCR/HOX): only as recent sub-families unless `--ancient`.
- **lncRNA families**: out of scope for protein-coding families; on the cDNA axis if wanted.

## Files
- `family_def_build.py` — the pipeline + manifest (this spec).
- `family_def_retrocopy_filter.py`, `family_def_strand_probe.py` — steps 1–2.
- `family_def_protein_validate.py`, `family_def_immune_and_protein.py` — validation + layers B/C.
- `family_def_protein_validation.md`, `family_def_cothread_results.md`,
  `family_def_unbiased_differentiators.md` — the investigation record.

---

## OPTIONAL alternative detector — shared k-mer common-core (advisor-requested)

`bench/family_def_common_core.py` — an **opt-in** alternative to the `~B` (pairwise-alignment + clique)
detector. It finds a family by the COMMON CORE: the shared k-mer substring across the reads/copy-models
of family members, assigning members by **containment** of that core (`|cand ∩ core| / |core| ≥ cmin`).

- API: `common_core(member_seqs)` → (core, body_len, ref, body_frac); `detect_family(seed, pool)` →
  members; `detect_for_gene(gene, cmin)` for CLI.
- Run: `python bench/family_def_common_core.py GENE [cmin]` (e.g. `RFPL1 0.25` → recovers RFPL1/2/3).
- Self-validation: `python bench/family_def_common_core.py`.

**Honest limitation (reproduced, matches prior hands-on experience):** the core captures only the
CONSERVED body — measured at **7–83% of the gene length** (RFPL = 24%; LILRA = 7%) — always under the
full gene, because the boundary is set by sequence conservation, not gene structure. And membership has
**no clean threshold**: permissive (`cmin` low) OVER-merges (pulls in domain-sharers, e.g. RPL36A +2),
strict UNDER-merges (drops divergent members). This is why `~B` (whole-gene reciprocal alignment +
clique) is the DEFAULT — its boundary is the gene, not the conserved core. The common-core module is
kept available for comparison / fallback, not as the primary detector.

### Read-following variant (`--multimap`)
A second mode reproduces the *actual* prior approach: seed = one copy → take its **multimapping reads**
→ follow them to where else they map → those loci are the other copies, each giving a **separate core**
(the region the shared reads cover). `python bench/family_def_common_core.py --multimap GENE`.

This OVER-merges symmetrically to the k-mer version's under-merge: the cores span the whole SIMILAR
region the reads can't disambiguate (RABL2A → cores 116–354% of the gene; MAGEA9 → 2244%, the whole
cluster), because the boundary is the read-confusable (`~R` / collapsed) extent, not the gene. So both
variants confirm the same root cause: **the core boundary is set by sequence conservation / read-
confusability, which has no fixed relation to the gene boundary — hence over- AND under-merge.** That
is exactly why `~B` (whole-gene reciprocal alignment of copy MODELS) is the default: its boundary IS
the gene.
