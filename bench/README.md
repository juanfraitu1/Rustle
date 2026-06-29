# bench/ docs index

The analysis record, consolidated 2026-06-29 from 92 micro-docs into 21 (git history retains every original).
Themed docs (UPPER_SNAKE with `_AND_`/`_DEF`) each merge several former `*.md` verbatim under `## <name>`
sections; the canonical standalone docs below are kept whole because they are current or heavily referenced.

## Themed (consolidated) docs
- **FAMILY_DEF.md** — the family-definition body (21 docs): read-conflict/de-tie definition, formal proof,
  BAM/junction/VG signals, genome-wide PSV graph, criterion bake-offs, detection validation.
- **COPY_ASSIGNMENT_AND_GATE.md** — copy-split, the significance gate, PSV LLR-vs-votes, τ sweeps, phasing,
  identifiability bound, resolution boundary, quant, StringTie head-to-head (17 docs).
- **ASJ.md** — allele-specific junctions (3 docs).
- **DENOVO_PIPELINE.md** — de-novo family pipeline, two-pass, read-coherence, intron-chain discovery (7 docs).
- **DNA_PROTEIN_VALIDATION.md** — DNA PSV catalog, DNA/RNA overlay, Compara, transcript validation (5 docs).
- **REFERENCE_ABSENT_AND_UNMAPPED.md** — reference-absent/collapsed copies, unmapped rescue, cross-chrom (6 docs).
- **PERFORMANCE_AND_IO.md** — pipeline speedup, input formats/ties, winnowmap vs minimap2 (3 docs).
- **REVIEWS_AND_AUDITS.md** — paper-grounded review, airtight fixes, loose-ends audit, scorecards (8 docs).
- **PANELS_AND_NOTES.md** — case panels + recent lever notes (conflict-edge unification, Bayesian posterior,
  intron-retention rescue, POA-core completion, loci/copies/isoforms, SEDEF build) (10 docs).

## Canonical standalone docs
- **SIM_GROUND_TRUTH.md** — fully-simulated planted-genome benchmark (the airtight non-circular validation).
- **copy_assignment_theory.md** — the machine-checked theory (Lemmas/Theorems 1–7, MWCA capstone).
- **DEFENSE_READINESS_AUDIT.md** — defense audit (objective-by-objective status).
- **OBJECTIVES_STATUS.md** — O1–O5 partition table.
- **RELATED_WORK_METHODS.md** — methods-to-find-multi-copy-families table.
- **FLAGSHIP_CASE_STUDIES.md** — the headline case studies.
- **FAMILY_LEVELS.md** — RNA/DNA/protein three-level cross-tab.
- **COPY_ASSIGN_RECOMPUTE.md** — the definitive O2 recompute on GGO_mm.
- **IDENTIFIABILITY_LIMITS.md** — the identifiability framing.
- **F4_SCOPE.md** — facility-location capstone scope.
- **MILESTONE_reference_absent_copies.md** — reference-absent milestone.
- **UNANNOTATED_COPIES.md** — de-novo zero-annotation paralog copies.

To recover any original: `git log --follow -- bench/<old-name>.md` then `git show <rev>:bench/<old-name>.md`.
