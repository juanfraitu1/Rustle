# Family Definition (consolidated)

> Merged from 21 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [family_def_SPEC](#family-def-spec) · [family_def_readconflict](#family-def-readconflict) · [family_def_bam_signals](#family-def-bam-signals) · [family_def_cothread_results](#family-def-cothread-results) · [family_def_junction_splicing](#family-def-junction-splicing) · [family_def_newbam_validation](#family-def-newbam-validation) · [family_def_nonlinear_census_results](#family-def-nonlinear-census-results) · [family_def_protein_validation](#family-def-protein-validation) · [family_def_scaffold_wiring](#family-def-scaffold-wiring) · [family_def_unbiased_differentiators](#family-def-unbiased-differentiators) · [family_def_vg_coherence_results](#family-def-vg-coherence-results) · [family_definition](#family-definition) · [family_definition_artifact_filter](#family-definition-artifact-filter) · [family_definition_formal](#family-definition-formal) · [family_definition_note](#family-definition-note) · [family_graph_definition](#family-graph-definition) · [poa_family_definition](#poa-family-definition) · [family_to_copy_bridge](#family-to-copy-bridge) · [family_criterion_bakeoff](#family-criterion-bakeoff) · [family_detection_validation](#family-detection-validation) · [psv_graph_genomewide](#psv-graph-genomewide)


---

## family_def_SPEC

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
reciprocal-coverage bars are the safeguard. Domain-architecture matching would be the strongest available check
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


---

## family_def_readconflict

# Defining a multi-copy family at the RNA level: read-conflict vs sequence similarity

## The question
The de-novo pipeline currently defines a family by **sequence similarity** — an edge between two rep
transcripts iff `contiguous_core_coverage ≥ T_CORE = 0.13` (or, as a contrast, minimizer-Jaccard), then dense
communities of that graph. This conflates true recent paralogs with **domain-sharers** (genes sharing one
exon/domain but not paralogous), uses an arbitrary threshold, and measures a *proxy* for the thing that
actually matters downstream. Do we have a good definition, and can we do better?

## The proposed definition (operational / identifiability-grounded)
A **family is a connected component of the read-conflict graph**: two loci are linked iff reads actually
cross-map between them with tied alignment scores — a read with a primary alignment in one locus and an
**AS-tied secondary** in the other (a genuine alternative placement at a *distinct* locus, not one read
overlapping two nested gene annotations). This is Canzar's multimapping-conflict graph, and it is exactly the
unit the copy-assignment problem operates on. Properties:
- **No arbitrary similarity threshold** — the boundary is the alignment-score tie (a property of the data).
- **Provable decomposition** — by construction reads never cross-map outside their component, so the
  assignment problem separates exactly across families with no information lost.
- It distinguishes paralogs from domain-sharers by the *extent* of cross-mapping (whole-transcript vs a domain).

## The experiment (`family_def_readconflict.py`, `family_def_comparison_figure.py`)
On a Compara-labelled set — 5 confirmed recent-paralog pairs (APOBEC3D/F, RFPL1/2/3, RABL2A/B) and 7
domain-sharer pairs (ASDURF/ASNSD1, CASP8/FLACC1, …, all nested opposite-strand genes sharing one exon) — we
measure, per pair, minimizer-Jaccard (similarity) and the read-conflict count + mutual coverage from `GGO.bam`.

| pair | class | Jaccard | read-conflict | MAPQ-0 signal |
|---|---|---|---|---|
| RABL2A~RABL2B | paralog | 0.313 | **190** (mut_cov 1.00) | 170/200 secondary |
| APOBEC3D~APOBEC3F | paralog | 0.144 | **6** (mut_cov 0.92) | within-family secondaries |
| RFPL1/2/3 (×3) | paralog | 0.13–0.17 | **0** | 0 secondary, 0 MAPQ-0 |
| all 7 domain-sharers | domain-sharer | 0.10–**0.41** | **0** | 0 |

## Findings
1. **Sequence similarity does NOT separate the classes.** The domain-sharer CREB1~METTL21A (Jaccard 0.406)
   scores *higher* than every true paralog (max 0.313). The domain-sharer range (0.10–0.41) encloses the
   paralog range (0.13–0.31). No threshold cuts them — the current definition's core failure, on labelled data.
2. **Read-conflict has zero false positives on domain-sharers** (7/7 = 0). The thing the similarity definition
   gets wrong (merging domain-sharers), read-conflict gets right — because nested genes sharing one exon
   produce no *alternative placement*; a read maps to one locus and merely overlaps two annotations.
3. **Read-conflict flags exactly the families where copy-assignment is needed** (RABL2 190, APOBEC3 6) and is
   **silent on RFPL** — a true paralog whose reads place *uniquely* (0 secondary, 0 MAPQ-0). RFPL is paralogous
   evolutionarily but poses no copy-assignment problem, so excluding it from the resolution-needed set is the
   correct operational behaviour.

## Interpretation
Read-conflict is **not** a paralog classifier and should not claim to be — RNA cannot robustly recover
evolutionary paralogy (established separately; needs protein/genomic orthology). It measures
**identifiability-RELEVANCE**: the set of loci among which reads are genuinely confused. That is a strictly
better definition for this thesis than sequence similarity, because (a) it never false-alarms on
domain-sharers, (b) it has no tuned threshold, and (c) it picks out precisely the families the downstream
machinery (PSV/junction assignment, the identifiability theorem) exists to resolve.

## Secondary-cap exposure (`secondary_cap_exposure.py`)
The BAM was built with `minimap2 -ax splice:hq` (no `-N`) → default **N=5 secondaries (≤6 placements/read)**.
For a co-located array all of a read's placements fall in one region, so records-per-read = multimapping degree;
a spike at 6 = saturation (a read may have MORE unreported placements → potentially missed conflict edges).

| region | multimappers | saturated (≥6) | max degree |
|---|---|---|---|
| TSPY array (>6 copies) | 44 | **36 (82%)** | 6 |
| RBMY1 array | 127 | 9 (7%) | 6 |
| chrY amplicon | 8 | 0 | 2 |
| APOBEC3 / RABL2 / RFPL / DSFAM38 (≤7 copies) | ≤15 | **0 (0%)** | ≤2 |

Verdict: the cap bites only **>6-copy co-located arrays** (TSPY-class). Even there the family is a connected
*component* — each read still carries 6 placements, over-connecting the array — so missing edges don't fragment
it; only edge *weights* and the inclusion of a perpetually-marginal copy are at risk. For every ≤6-copy family
(all validated ones), exposure is zero. → safe to port; flag edge weights as lower bounds in mega-arrays.

## Next step (implementation)
Port the criterion into `detect_edges`: keep POA contiguous-core as a cheap homology *prefilter* (necessary to
cross-map), but **confirm and weight each edge by read cross-mapping** (the AS-tied secondary index already
exists: `tied_secondary_reads`). Families become the conflict-graph components. The portable kernel
(`read_conflict.rs`: AS-tied conflict edges + connected-component families) is TDD'd in the Rust tree; the
remaining integration is plumbing per-locus secondary placements into the detection stage.


---

## family_def_bam_signals

# More signal in the BAM: NM (and rl) beyond de — retained introns in the copy model

Two questions about whether the definition extracts the BAM fully.

## 1. Retained introns ARE summed into the copy model — by design, with one subtlety

`S(v)` = union of every read's `get_blocks()` (the aligned M/=/X segments). A normally-spliced
read contributes its **exons** (the intron is an `N` gap, excluded). An **intron-retention** read
has no `N` over the retained intron — it is aligned there — so that intron *is* a block and is
summed. This is the intended "full gene minus introns *except retained ones expressed in the RNA*."

**Subtlety (a real refinement).** The union takes *all* reads' blocks with no support threshold, so
a **single** intron-read-through read pads `S(v)` with the whole intron. Because `~B` scores
**reciprocal coverage**, padding `S(i)` with intronic sequence the properly-spliced homolog `S(j)`
lacks **lowers** `cov_i` and can push a real pair toward the `τ` edge — a plausible cause of the
validated families that sit at reciprocal coverage ~0.31. **Fix:** gate block inclusion by read
support (include a base only if ≥k reads align through it), so single-read read-through does not
enter `S(v)`; this should *raise* real families off the `τ` edge. (Not yet applied.)

## 2. de is not the only signal — NM discriminates bridges that de cannot

The new BAM (`-N 50 -p 0.1 --eqx`) carries `NM` (raw edit distance), `AS`/`ms` (DP scores),
`s1`/`cm` (chaining), and `rl` (length of query in repetitive seeds) — only `de` is used. Per
chromosome, classify `~R` cross-mappings as `~B`-copy vs `~B`-bridge and compare the reads:

**NC_073240.2 — 69 copy / 59 bridge edges, 4,560 read-placements**

| signal | real-copy (median) | bridge (median) | separation |
|---|---|---|---|
| `de` (gap-compressed divergence) | 0.0094 | 0.0083 | 0.9× — **none** |
| `NM`-rate (edit distance / aligned bp) | 0.0013 | 0.0035 | **2.6× — discriminates** |
| `rl`-fraction (repeat-seed bp / aligned bp) | 0.0041 | 0.0066 | 1.6× — weak |

**`de` does not separate bridges from copies at the read level (0.9×)** — confirming that `~B`
(copy-level homology) is doing the bridge-pruning, not the read-level tie test. But **`NM`-rate
separates them 2.6×**: a bridge read is forced onto a non-homologous locus with large indels, so
its *raw* edit distance is high while its *gap-compressed* `de` stays low. `de` compresses away
exactly the indel burden that flags a bridge; `NM` keeps it.

### Implication
- `de` is the correct **tie** criterion (it measures sequence confusability, which is genuinely
  similar for copy and bridge cross-mappings — the bake-off already preferred it over `AS`, which
  is biased on intronless retrocopies). `DE_max` as an absolute ceiling is weakly load-bearing.
- **`NM` adds an orthogonal, read-level bridge signal** that `de` misses. It could (a) tighten
  `~R` by down-weighting high-`NM` (forced-indel) placements before they form edges, or (b) serve
  as a cheap read-level pre-filter complementing `~B`'s copy-level pruning. `~B` remains the
  principled pruner (it directly tests "are these copies"); `NM` is a corroborating proxy — which
  is consistent with the conflict-criterion bake-off finding that NM *corroborates* de.
- `rl` (repeat-mediated-bridge flag) is directionally right but weak here; `AS`/`ms` are
  length/structure-biased. `NM` is the highest-value currently-unused tag.

Scope: one chromosome; `~B`-class via reciprocal coverage of exon-union models (`COV_MIN=0.30`).
Reproduce: `python bench/family_def_bam_signals.py NC_073240.2`.


---

## family_def_cothread_results

# Junction-position co-threading — genome-wide deployment (first build)

Deploys the panel-validated junction-POSITION co-threading metric
(`family_def_backbone_vgspine.py`) at genome scale to test whether it de-bridges the cDNA-homology
over-merge while keeping real (count/length-drifted) paralogs connected — the one thing intron
count/length concordance could not do.

## What was built
- `build_gene_exon_index.py` → `gene_exon_index.json` (30,453 genes, representative-transcript exon
  coordinates; needed to build spliced cDNAs).
- `family_def_cothread_genomewide.py` — aligns the two spliced cDNAs of each candidate edge
  (minimap2, batched per reference, 5-way parallel), projects junctions through the CIGAR, and emits
  per-edge position features: `frac_ref` (share of the RICHER backbone co-threaded), `frac_mem`
  (share of the smaller), `aln_cov_small`, `contig` (longest contiguous co-threaded run).
- Scope: both-spliced edges inside the large post-retrocopy components (≥20 members, where
  over-merges live) + a real-control set (panel families + Compara checkable multi-exon pairs +
  held-out GGT1~GGTLC2). 8,569 bridge-candidate + 13 real-control edges aligned.

## Result — partial win, honestly bounded

**The big components are all genuine OVER-MERGES** (heterogeneous unrelated named genes bridged by
LOC/Gnomon predictions), not real mega-families — so there is a real de-bridging target:
- comp-201: named = ACTR3B, AGAP5/6, ARPC2, BCO2, CASR, CCDC62… (unrelated); median frac_ref = **0.00**
- comp-140: named = ATG10, BCR, CPO, EPPK1, FUT3/5, IFNE… (unrelated); median frac_ref = **0.00**
- comp-238: 234/238 LOC; the 4 named (FOXD4L1, NDUFA8, RHBDF1, SLC35F5) are mutually unrelated;
  median frac_ref = **0.33** (a DENSE over-merge)

**Co-threading cleanly cuts the frac_ref≈0 over-merges** (201, 140): removing low-co-threading edges
splits them into coherent pieces (201 → 14+8+8+…; 140 → 37+19+10+…) because their genes are
homologous in SEQUENCE but their junctions do NOT align (the domain-bridge signature).

**A low global threshold de-bridges without harming real paralogs.** Filtering measured edges at
`frac_ref ≥ 0.20` (keeping all unmeasured small-component edges):

| thr | families | max comp | named reals connected |
|----:|---------:|---------:|----------------------:|
| base| 1216     | 238      | 5/5 |
| 0.15| 1251     | 221      | 5/5 |
| **0.20**| **1252** | **219** | **5/5** |
| 0.33| 1261     | 205      | 4/5 ← cuts GGT1~GGTLC2 |

At 0.20: **+36 families**, all five named real paralogs stay connected — including the held-out
partial-dup **GGT1~GGTLC2** (frac_ref 0.20: GGTLC2 shares 3 of GGT1's 15 junctions). The real
families span frac_ref 0.20–0.89 (RABL2 0.89, APOBEC3 0.67, RFPL 0.50, GGT 0.20); the pure domain
bridges sit at ~0.00. The 0.20 gap between them is what count-concordance could not find.

**The limit:** the densest over-merge (comp-238, 234 LOC genes, median 0.33) RESISTS — those genes
genuinely co-thread a shared multi-exon repeat/domain (conserved junctions despite being unrelated
families). Junction position cannot break a co-threading shared domain. Max component stays 219.
Cracking it needs community detection on the co-threading-weighted graph (weak-cut / modularity),
not a global edge threshold.

## Status
- The metric is a real, sequence-orthogonal de-bridging lever for the **junction-non-aligning**
  over-merge subtype (domain bridges proper), safe for real paralogs at frac_ref ≥ 0.20.
- It does NOT solve the **co-threading shared-domain/repeat** over-merge (the dense LOC blob).
- GGT1~GGTLC2 at exactly the 0.20 boundary is a genuine definitional question (is a partial
  duplication "the same family"?), not a metric failure.

## UPDATE — co-threading-weighted community detection (the fix for the dense blob)

`family_def_community_debridge.py`: per large component, weight edges by frac_ref and run weighted
Louvain (`networkx louvain_communities`, seed=0); small components (real families) left untouched.
This is what the global threshold could not do — it cracks the dense comp-238 repeat blob.

| config | families | max comp | real families intact | coherence |
|---|---:|---:|:--:|---:|
| base (over-merged) | 1216 | **238** | — | — |
| co-thread res=0.5 | 1266 | 76 | 5/5 | 0.67 |
| **co-thread res=1.0** | **1284** | **57** | **5/5** | 0.67 |
| co-thread res=1.5 | 1289 | 47 | 5/5 | 0.67 |
| **co-thread + strand** | 1112 | 57 | **5/5** | **0.77** |

**comp-238 (238 genes, the blob the threshold left at 219) → 17 communities (max 57).** comp-201 → 135,
comp-140 → 52. All **5/5 named real families stay together** (RABL2, APOBEC3, RFPL, and the borderline
partial-dup GGT1~GGTLC2 at frac_ref 0.20). The mutually-unrelated named anchors of each over-merge
(comp-238: NDUFA8/RHBDF1/SLC35F5) land in **distinct** communities — the bridge is cut.

**Not over-fragmentation (audited):** the components that shatter into singletons are genuine
over-merges of unrelated genes (e.g. a 25-gene component whose 10 named members are 10 different
roots — CADM2/HMGCR/GBP2/…, density 0.11, median weight 0.00). Only 2 same-root groups exist in the
big components (TMEM, ZNF — domain-name prefixes, not true families). The big components are NOT real
families, so shattering them is correct; the real families (validated by name) stay intact.

**Caveat:** the large LOC-only sub-communities of comp-238 (57/50/48) are sub-clusters of a shared
repeat/domain; whether each is a *true* family can't be settled from annotation (uncharacterized
Gnomon predictions) — needs orthogonal (protein/Compara) evidence.

## Stranded graphs (the user's idea — verified, complementary)

cDNA homology with alignment strand '-' (between two mRNA-SENSE cDNAs) = one gene's sense vs the
other's ANTISENSE = sense/antisense overlap or opposite-orientation repeat, NOT paralogy. `dna_homology`
ignored strand. `family_def_strand_probe.py`: 482/15479 (3.1%) post-retrocopy edges are antisense;
dropping them removes **435 genes / ~158 pseudo-families** that were connected ONLY by antisense edges
— spot-checked, nearly all genomically OVERLAP at id=1.00 (CELF3~SNX27, CASQ2~VANGL1, CLCC1~GPSM2 =
classic sense/antisense pairs). All 5 real families are same-strand '+' (inverted-dup paralogs are
still sense-sense in cDNA space), so it is safe. BUT strand does NOT fix the headline blob: comp-238
has only 1 antisense edge (same-strand repeat). It de-bridges comp-201 (20% antisense) and is a clean
precision axis that raises post-community coherence 0.67→0.77.

## Status / next
- DENSE over-merge SOLVED by co-threading community detection (238→57, reals intact, robust to res).
- Strand = complementary safe precision filter (removes antisense pseudo-families); combine with above.
- Remaining: validate the LOC-only sub-communities (needs orthogonal evidence); run `--all`; minor
  Louvain run-to-run jitter (seed set; ±1 family) to pin down.


---

## family_def_junction_splicing

# Can the VG / junction structure define gene families? (precision, not just resolution)

Synthesis of 4 angles (junction-architecture, genome-scale, vg-cothreading, splicing-unique),
each adversarially verified (confound/base-rate + scale/generalization lenses).
All headline numbers below were re-reproduced (this session) from the live scripts/data.

## Bottom line

**Q1 (VG for family DEFINITION): PARTIAL / mostly NO at scale.** The VG/junction-backbone is the
*correct structural object* for stating the criterion (Canzar-shaped, threshold-light), and it gives
one clean by-construction precision bit. But as a genome-wide family *caller* it does NOT separate the
real over-merge problem: no single VG/architecture operating point both cuts the 389-gene domain-bridge
mega-family AND keeps the validated real families connected.

**Q2 (a UNIQUE junction signal orthogonal to sequence id): YES, but narrow.** Exactly one clean,
threshold-free, sequence-orthogonal signal survives verification: **retrocopy intronlessness**. It flags
processed-pseudogene copies that cDNA identity provably cannot (id ~0.95-1.0). Everything broader
(continuous intron-count / intron-length / exon-count "architecture concordance") is either an arbitrary
panel-fit threshold, a gene-size proxy, or actively destroys real diverged paralogs.

cDNA id is *anti-discriminative* on this panel (AUC real-vs-homologous-counter = 0.029): the false
counters have id ~1.0, higher than the diverged real families (0.90-0.99). So an orthogonal signal is
genuinely needed — junctions provide one, but only for the retrocopy sub-problem.

## The single cleanest unique signal: retrocopy intronlessness (S1)

Definition: an edge (A,B) with cDNA id >= 0.9 where **min(n_intron(A), n_intron(B)) == 0 while the other
is spliced**. Exact from the GFF exon count; no threshold.

What it buys (that sequence id cannot):
- Flags processed pseudogenes / retrocopies. Panel: EEF1A1 (7 introns) ~ LOC retro (0 introns) id=0.98;
  CNN2 (6) ~ LOC retro (0) id=0.99. 2/2 retrocopies flagged, **0/5 false flags on real families**.
- The edges it cuts have median cDNA id ~0.947 (range 0.90-1.00) — sequence id is blind to them.
- Mechanistically clean & VG-native: a 0-intron copy **provably cannot thread a multi-exon spliced
  backbone** (no bubbles to traverse) → structural, by-construction exclusion, not a tuned filter.
- NOT a coverage proxy: retrocopies have high cDNA coverage (0.43-1.00).

Where it FAILS / its scope limits:
- It does **NOT** touch the actual over-merge. The 389-gene mega-family is built from **multi-exon ↔
  multi-exon DOMAIN BRIDGES**. The strict intronless-only filter (Cspliced) cuts only **196/1522 = 13%**
  of mega edges and shrinks the mega merely 389 → a **238 + 201 + 140-gene** triple blob. The domain
  bridge survives essentially intact. (Reproduced this session.)
- **Critical correction (verified):** "intronless ⇒ retrocopy ⇒ not a family" is biologically FALSE.
  Of the 5198 either-intronless edges only **1931 are the genuine one-side-intronless (retrocopy)
  pattern**; **3267 (63%) are BOTH-intronless** and contain real single-exon tandem families:
  **PCDHB (33 edges), IFNA (26), HNRNPCL, ELOA...** Dropping every intronless-incident edge therefore
  destroys real intronless tandem families ~10:1 over true retrocopies. **The valid signal must be
  restricted to one-side-intronless-vs-spliced-parent (1931 edges), NOT all-intronless.** The panel hid
  this because all 3 panel real families are multi-exon (zero intronless real families in the panel).

So S1 (restricted to one-side-intronless) is a real, clean retrocopy filter — orthogonal to id — but it
is a side problem, not the domain-bridge over-merge.

## The VG combinatorial criterion (Canzar-shaped) and what the empirics say

Criterion (threshold-light, provable-shaped):
> A candidate set S is a **family** iff its members are **parallel s-t paths through one variation graph
> that share its junction backbone** — i.e. there is a single bubble-chain spine whose branch nodes
> (junction boundaries) are co-threaded by all members in the sequence-aligned frame. A domain-sharer
> threads only a sub-path / lands in a separate component; a retrocopy is a single straight-through path
> with no bubbles and cannot thread the spine at all.

Operationalized (A3) as `shared_backbone_fraction` = (# reference-spine junctions co-threaded by all
members, in the cDNA-aligned frame, position-matched within TOL bp) / (# reference junctions).

What the empirics support:
- **Panel (n=9): clean.** REAL = {RABL2 0.89, APOBEC3 0.67, RFPL 0.50}; ALL domain-sharers and retros
  = 0.00. Margin REAL_min − COUNTER_max = **+0.50**, and it is **tolerance-insensitive** (margin stays
  +0.50 for TOL in {0,4,8,...,50} bp, even exact TOL=0). It beats intron-count-ratio (margin −0.167)
  and it is NOT a coverage proxy (APOBEC3-real at alnCov 0.39 beats CASP8-domain at alnCov 0.46 — the
  discriminator is junction *position*, not coverage). The mechanism is the right one: **position is the
  correct invariant; intron/exon LENGTH is the wrong invariant** (real paralog introns drift in length).

What the empirics do NOT support (verifier-downgraded from "strong" to "promising on n=9"):
- **The fractional position metric was never deployed at scale.** The genome-wide numbers substitute the
  exact intronless axis + a cDNA-cov<0.5 proxy, not the validated junction-position projection.
- **The strict genome-scale architecture form (Carch = equal intron count + ordered length concordance)
  DESTROYS all three real families** (verified this session): under Carch the mega breaks 389 → 32, but
  RABL2 keeps only RABL2A, and **APOBEC3 and RFPL have ZERO members left in the graph**. Reason: real
  paralog intron lengths/counts diverge post-duplication (RABL2A/B introns 455 vs 214, 4213 vs 2858;
  RFPL counts 3/4/2).
- **The continuous count_conc "family gate" (>=0.6) is panel-overfit and breaks a held-out positive.**
  The single Compara-validated paralog held out of the panel, **GGT1~GGTLC2 (cDNA id 0.98), has
  count_conc = 0.267** (15 introns vs 4) → its only edge is dropped → family destroyed (false negative).
  RFPL2~RFPL3 = 0.50. A real partial-duplication paralog has the *same* architecture signature
  (mismatched intron counts, asymmetric coverage) as a domain-sharer (count_conc up to 0.67), so the gate
  provably cannot separate them. The 0.6 threshold is non-eliminable and chosen only to keep RFPL
  connected — exactly the magic number Canzar distrusts.
- **exon-count concordance (S5, panel AUC 0.929) is a gene-SIZE proxy** (corr |dExon| vs max-exon-count
  = 0.71) and degrades to AUC ~0.755 genome-wide; it fails on size-matched domain-sharers (ASDURF dExon=1,
  GPR39 dExon=2 fall inside the real range).
- **Two goals are in direct tension:** the real-family-safe gate leaves a 204-238-gene blob and
  collaterally damages ~65% of comparable size-2/3 families; the blob-breaking gate annihilates the real
  families. No single operating point does both.

## Decisive vs suggestive

DECISIVE (survives both verifier lenses, reproduced):
- cDNA id is anti-discriminative for this over-merge (AUC 0.029); an orthogonal signal is needed.
- **One-side-intronless retrocopy flag** is a clean, threshold-free, sequence-orthogonal precision bit:
  1931 such edges genome-wide, median id ~0.95, 2/2 panel retros, 0 panel real. VG-native exclusion.
- The retrocopy flag does **not** solve the domain-bridge mega (only 13% of mega edges; 238+201+140 blob).
- Strict architecture concordance (Carch) destroys all 3 real families — rejected.
- The continuous count_conc gate (>=0.6) is panel-fit and produces a false negative on the held-out
  Compara paralog GGT1~GGTLC2 (cc=0.267).
- "all-intronless ⇒ retrocopy" is false; the both-intronless half (3267 edges) holds real tandem
  families (PCDHB, IFNA, ...). The signal must be the one-side-intronless restriction only.

SUGGESTIVE (clean on n=9, unproven at scale):
- The **junction-POSITION co-threading fraction** (shared_backbone_fraction, margin +0.50, tolerance-free)
  is the most promising VG-native criterion and the correct invariant (position not length), but it was
  never run at genome scale; the deployed proxy contradicts it. This is the one thing worth building.

## Honest answer

- VG helps DEFINE families only partially: it supplies the right *structure* and one clean by-construction
  bit (retrocopies can't thread a spliced spine), but it does not, at any tested operating point, cleanly
  cut the multi-exon domain-bridge over-merge while preserving real diverged paralogs.
- The unique junction signal orthogonal to sequence id is real but narrow: **retrocopy intronlessness**
  (restricted to one-side-intronless), which sequence id cannot see. The broader "architecture
  concordance" story does not survive scale.

## Recommended next step

Build and genome-wide-deploy the **junction-POSITION co-threading metric** (A3's
`shared_backbone_fraction`, position-matched in the cDNA-aligned frame, tolerance-free) as the actual
VG criterion — it is the only candidate that is both Canzar-shaped and panel-clean — and test whether a
*position*-based (not length/count) backbone-overlap fraction can cut multi-exon domain bridges
(partial-colinear / sub-path alignments) while keeping length-drifted real paralogs (RABL2/APOBEC3/RFPL,
GGT1~GGTLC2) connected. Ship the one-side-intronless retrocopy flag now as a standalone, by-construction
precision filter (1931 edges) — independent of whatever the position metric shows.

Artifacts: family_def_backbone_vgspine.py (A3, position metric, panel),
family_def_junction_genomewide.py (Cspliced/Ccount/Carch genome-wide),
family_def_intron_index.py (builds /home/juanfra/winloci_scratch/gene_intron_index.json),
family_def_splicing_signals_a4.py (S1-S5), family_def_junction_lever.py (count_conc gate),
compara_paralog_relation.json (held-out positives incl. GGT1~GGTLC2).


---

## family_def_newbam_validation

# The family definition under proper multimapper sampling (new -N 50 -p 0.1 BAM)

The definition's single biggest empirical caveat was that `~R` (read-confusability) was built on
`GGO.bam`, aligned with minimap2's **default secondary cap** — so most cross-mapping was never
emitted, and `~R` was measured on a small fraction of the real read confusion. A re-alignment
with `-N 50 -p 0.1 --secondary=yes --eqx -Y` (`GGO_mm.bam`) surfaces the suppressed secondaries
(de tags intact). This re-validates the definition on the data it should always have had.

Method (`bench/family_def_newbam_validate.py`): per chromosome, re-run the **exact** `~R` scan
(`family_def_genomewide.scan` logic, region-restricted) on OLD vs NEW, then apply `~B` (cached /
on-demand-built exon-union copy models, reciprocal coverage ≥ 0.30). `~B` copy models are
BAM-sampling-independent (built from primary alignments), so the comparison isolates the `~R`
effect.

## Result (three paralog-dense chromosomes)

| chrom | secondaries OLD→NEW | multimapper reads | `~R` candidate edges | **`~B` real-copy edges** | `~B`-pruned bridges |
|---|---|---|---|---|---|
| NC_073244.2 (RFPL) | 4,818 → 313,883 (65×) | 751 → 5,113 (6.8×) | 18 → 33 | **2 → 2** | 16 → 31 |
| NC_073248.2 (chrY-like) | 5,732 → 92,889 (16×) | 182 → 353 (1.9×) | 28 → 49 | **21 → 35** | 7 → 14 |
| NC_073240.2 | 22,309 → 258,319 (12×) | 2,487 → 10,821 (4.4×) | 68 → 128 | **58 → 69** | 10 → 59 |

## What this shows

1. **The undercount caveat was real and large.** The old BAM emitted 12–65× fewer secondaries and
   1.9–6.8× fewer multimapper reads than a proper `-N 50 -p 0.1` alignment. `~R` was genuinely
   undersampled.

2. **`~R` gains recall — it recovers hidden real families.** On the *dispersed*-paralog
   chromosomes, `~B`-validated real-copy edges **grow** (chrY 21→35, +67%; NC_073240 58→69, +19%):
   genuine paralogs whose cross-mapping the default cap had hidden, now detected. On RFPL the count
   is stable (2→2) because its copies are co-located/near-identical and already cross-mapped even
   under the old cap — exactly where the secondary cap did *not* bite.

3. **`~B` remains the precision filter, and it scales.** Every extra candidate edge the richer
   sampling surfaces beyond the real copies is a bridge, and `~B` prunes it (bridges 10→59 on
   NC_073240). Without `~B`, "candidate families" would balloon with the sampling depth (5→14 on
   RFPL); with `~B`, the validated set tracks real copies, not read counts.

**Conclusion.** The two-relation design is vindicated by the data it was missing: `~R` supplies
recall (now properly sampled, recovering hidden dispersed paralogs), `~B` supplies precision
(load-bearing, absorbing the extra bridges that richer sampling inevitably surfaces). The
definition is not merely robust to multimapper sampling — it improves with it, because `~B`
anchors family identity on copy-model homology rather than on read counts.

## Honest scope
- Three chromosomes, not genome-wide (the full `~R` scan on the 11.7 GB BAM is memory-heavy; a
  per-chromosome serial genome-wide pass is the natural next step).
- "Real copy" = `~B` pass (reciprocal coverage ≥ 0.30 of the exon-union models), the note's copy-
  level metric (validated as the correct copy-level criterion in `family_def_metric_compare.py`;
  contiguous-core is too strict on exon-union models — it fragments real copies).
- `~B` copy models reused from the old-BAM primary alignments (sampling-independent by design);
  on-demand models for new candidate loci built identically from primary reads.


---

## family_def_nonlinear_census_results

# Non-linear (tandem-repeat CNV) family census — adversarial verdicts

**Claim under test:** In 10 "K-frontier" families the copies are near-identical in sequence
(cDNA id ≥ 0.98) yet differ in **tandem-repeat unit count**, so the SNP/indel-PSV model
(a *linear* variation graph) collapses them but the tandem-repeat copy-number **structure**
(a *cyclic* vg graph) distinguishes them → structure moves the identifiability frontier.

**Method (per family, real code):** minimap2 self-alignment to count internal tandem-repeat
units independently of vg; pairwise alignment of the two most-similar members with the repeat
region masked to test whether the *unique* (non-repeat) sequence is really ≥98% identical;
CDS-vs-UTR via proteins.fa (CDS-only) length tracking; genomic arrangement from genes.bed;
repeat-array length vs a 1–10 kb FLNC read. Default to **refuted** unless numbers support it.

## Per-family table

| Family | Members / chrom | Repeat-count diff real? | Unit counts (key members) | Most-similar pair: unique-region id | Collapsed by SNP? | CDS/UTR | Arrangement | Read-spanable? | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| RCF_102 | 6 / mixed | yes | 6/4/2/2 (~2.3 kb unit) | **47.8%** (CNV-differing pair has 4.5 kb unique flank) | no | UTR | mixed (4 co-located ≤99 kb + 2 dispersed) | yes | **REFUTED** |
| RCF_432 | 12 / 5 chrom | no | only high-id pair 496/567 = equal (369 bp/13.67 u) | 99.7% but **+2 unique-region SNPs**, counts equal | no | CDS | mixed | yes | **REFUTED** |
| RCF_6 | 23 / 10 chrom | yes | 5/4/3/14 (~482 bp unit) | 99.6% core but **704 bp member-specific indel flank** | no | mixed | dispersed (near-id pair 2.98 Mb apart) | yes | **REFUTED** |
| RCF_277 | 6 / 1 chrom | yes | ZF 11/11/9/7/6/5 | 99.95% but pair has **equal ZF count (11=11)**; differs by 294 bp unique N-term | no | CDS | mixed | yes | **REFUTED** |
| RCF_53 | 13 / 1 chrom | yes | 2..7 (~800 bp unit) | 99.65% shared body but **400–620 bp member-private 5′ CDS** | no | CDS | mixed | yes | **REFUTED** |
| RCF_718 | 7 / 1 chrom | yes | 2.7–5.5 (~105–120 bp unit) | **92.1%** (CNV-differing pair); only ≥0.98 pair is byte-identical | no | CDS | co-located | yes | **REFUTED** |
| RCF_611 | 4 / 1 chrom | yes | ~9/~2 units (~30 bp unit) | whole-cDNA **≤91.8%**; unique region still has **4 callable SNPs** | no | CDS | co-located | yes | **REFUTED** |
| RCF_135 | 3 / dispersed | yes | 1/2/3 (743 bp unit) | **0%** homologous unique region (unique = fully member-specific) | no | mixed | dispersed (136.6 Mb apart) | yes | **REFUTED** |
| RCF_13 | 8 / mixed | yes | 4/3/3/2/1.. (~291 bp unit) | co-located pairs cover 31–53%; 3′UTRs **don't align**; high-id pair on diff chrom + 189 bp indel | no | UTR | mixed | yes | **REFUTED** |
| RCF_30 | 9 / 8 chrom | **no** (terminal dup, not array) | 2/2 (terminal direct dup) + 7×1 | 99.34% but **20 scattered SNPs**, repeat-CNV not real | no | mixed | dispersed (≥7.6 Mb) | yes (moot) | **REFUTED** |

## Synthesis

### Genuine K-frontier cases (genuine_kfrontier = true)
**0 of 10.** None. Every family fails the conjunction "unique region ≥98% id **AND** distinguished
only by tandem-repeat count **AND** co-located **AND** read-spanable."

### How the 10 were refuted (the recurring failure modes)
- **Repeat dominates the alignment → inflated id (the central artifact): 6 families**
  (RCF_102, RCF_6, RCF_53, RCF_277, RCF_718, RCF_13). The census ≥0.98 came from best-*local*
  chain identity or the gap-compressed `de` tag, where repeat units tile against each other.
  When the repeat is masked, every **CNV-differing** pair drops to **48–92% global id** and/or
  reveals a large **member-specific unique flank** (e.g. RCF_102's 4.5 kb 5′ flank, RCF_53's
  400–620 bp private 5′ CDS, RCF_277's 294 bp unique N-terminus, RCF_6's 704 bp indel flank).
  Those unique flanks already carry SNP/indel PSVs → **linear PSVs already resolve them**.
- **Where sequence IS collapsed, the repeat counts are EQUAL (dichotomy is fatal): 4 families**
  (RCF_432 496/567 both 13.67 u; RCF_277 both 11 ZF; RCF_718 byte-identical pair; RCF_102
  134~136 both 2 u). The one near-identical pair per family has the *same* unit count, so
  tandem-CNV structure distinguishes **nothing**; they are separated (if at all) by a handful
  of unique-region SNPs, not by structure.
- **Unique regions genuinely differ → SNPs already resolve: all 10** (the union of the two above;
  RCF_611 had 4 callable SNPs even after masking; RCF_30 had 20 scattered SNPs).
- **Dispersed → not a multimapping/copy-assignment problem: 3 families** (RCF_6 near-id pair
  2.98 Mb apart; RCF_135 136.6 Mb apart; RCF_30 9 members on 8 chromosomes, only same-chrom pair
  7.6 Mb apart). No co-location ⇒ reads don't multimap ⇒ there is no copy-assignment ambiguity to
  resolve regardless of structure.
- **Repeat-diff was not real (closest to a vg artifact): 1 family** (RCF_30 — a single head-to-tail
  terminal direct duplication in 2/9 members, periodicity ≈ 0.28 = random-DNA baseline; not a
  variable-count tandem array). RCF_432's apparent 13-vs-0 was also a self-align detection-threshold
  artifact (the high-id pair is actually equal at 13.67 u).
- **Repeat too long to span: 0 families decisively** (one array in RCF_102, ~13 kb, exceeds a 10 kb
  read, but its 4.5 kb unique flank is read-spannable anyway, so spanning is never the binding
  constraint).

### CDS vs UTR breakdown of the genuine ones
**N/A — there are no genuine cases.** For completeness, of the (refuted) families the variable
repeat was: CDS in 5 (RCF_432, RCF_277, RCF_53, RCF_718, RCF_611), UTR/VNTR in 2 (RCF_102, RCF_13),
mixed (constant copy in CDS, variable copies in UTR) in 2 (RCF_6, RCF_135), and not-a-real-array
in 1 (RCF_30). Both CDS-repeat and UTR-VNTR are read-detectable; this breakdown would matter only
had any case survived.

### Bottom line — concrete identifiability-frontier gain
**Zero families' worth.** A structure-aware ("structural PSV") VG that adds tandem-unit count as a
distinguishing variant would newly resolve **0** of these 10 currently-SNP-collapsed candidates,
because in every one either (a) the copies that differ in unit count are already ≥8–52% divergent in
their unique sequence (SNP/indel-PSVs in the unique flanks already separate them), or (b) the copies
that are truly sequence-collapsed have **identical** unit counts (structure adds no axis), or (c) they
are dispersed across chromosomes (no multimapping, so nothing to resolve). The tandem-repeat CNV is
frequently real, but it is **redundant** with — never the **sole** — distinguishing axis. The
identifiability frontier does not move.

### What a "structural PSV" is, and how copy-assignment would use it
A **structural PSV** is a copy-distinguishing variant whose state is a **tandem-unit count** (a cyclic
bubble's traversal multiplicity in the variation graph) rather than a base substitution/indel — i.e.
the allele a read carries at this site is "how many times it traversed the repeat unit." Copy-assignment
(#2) would read the count off a single FLNC molecule that spans the array and thread the read to the
copy whose unit count matches, exactly as it threads SNP/indel alleles to copies — adding one extra
column (a count, not a base) to the per-read allele profile used in the path-cover/threading.

### Recommendation
**Do NOT wire structural-PSV (tandem-unit count) into copy-assignment (#2) now.** Across all 10
hand-picked best candidates it resolves nothing the linear SNP/indel-PSV model doesn't already
resolve (the unit-count signal is redundant where present and absent where it would be needed), and
the near-identical-yet-CNV-differing cell of the dichotomy is empirically empty. Park it: revisit only
if a *co-located* family with a **read-spannable** array, **equal unique-region sequence (≥0.98 after
masking)**, **and** a genuine unit-count difference is ever found — none exists in this census.


---

## family_def_protein_validation

# Protein-level validation — the orthogonal ground truth for family definition

The cDNA-homology over-merge could never be adjudicated from the inside (0 trustworthy positives in
the big blobs). Protein homology is the orthogonal validator the project flagged as needed. Built from
`proteins.fa` (22,614 representative proteins, gene-name keyed — directly comparable to the cDNA
graph) + mmseqs all-vs-all (`prot_ava.m8`, 982,949 hits). Protein paralog edge = fident≥0.30 ∧
reciprocal-cov≥0.50 ∧ e≤1e-5 → **199,828 protein-paralog gene-pairs** = the positives we lacked.
Scripts: `family_def_protein_validate.py`, manifest at `protein_family_manifest.tsv`.

## 1. The over-merge is CONFIRMED and EXPLAINED

| cDNA component | n | coding % | non-coding | distinct protein orthogroups |
|---|--:|--:|---|--:|
| comp-238 (the resistant blob) | 238 | 38% | **62% lncRNA** | **10** (+148 lncRNA) |
| comp-201 | 201 | 71% | 29% lncRNA | 17 |
| comp-140 | 140 | 51% | 49% lncRNA | 10 |
| comp-47 | 47 | **100%** | — | **1** (a REAL family) |
| comp-41 | 41 | **100%** | — | **1** (a REAL family) |

The dense comp-238 blob that resisted every junction/topology lever is a **non-coding repeat
over-merge** — 62% lncRNA sharing a repeat element, plus 90 coding genes that form 10 *distinct*
protein families. comp-47 and comp-41 are each a single protein orthogroup = genuine large real
families. Big over-merge components are 30% lncRNA vs 16% in small families.

## 2. CORRECTION — co-threading beats cov_min (only visible with ground truth)

Recovering the protein-orthogroup partition inside the over-merged components (pairwise prec/recall):

| de-bridge weight | precision | recall | F1 (bar .30/.50) | F1 (bar .40/.60) |
|---|--:|--:|--:|--:|
| none (raw over-merge) | 0.28 | 1.00 | 0.44 | 0.42 |
| cov_min | 0.70 | 0.60 | 0.64 | 0.62 |
| **frac_ref (co-threading)** | 0.71 | 0.74 | **0.72** | **0.70** |
| product (cov_min×frac_ref) | 0.64 | 0.56 | 0.60 | 0.58 |

**co-threading recovers real families better than cov_min** — cov_min's smaller max-component (43 vs
57) was OVER-FRAGMENTATION (recall 0.60: it splits real protein families), a failure mode invisible
without ground truth. Default reverted to `frac_ref`. (This corrects the prior call to switch to
cov_min, which was based on max-comp reduction alone.)

## 3. DELIVERABLE — protein-confirmed family definition

Keep cDNA-homology edges that are ALSO protein-homologous (coding ∧ reciprocal-cov≥0.50):

- **837 families** (vs cDNA over-merged 1216), **max comp 238 → 47** (a real orthogroup), 2,944 genes.
- 58% of cDNA edges are protein-confirmed; **4/5 named real families** preserved.
- Dropped 1,775 genes with no protein-confirmed paralog: 772 lncRNA (non-coding repeat bridges,
  correctly excluded), 845 protein_coding (cDNA bridges in non-coding regions / sub-bar), and — as a
  KNOWN limitation — 103 V_segment + 13 C_region (immunoglobulin/TCR segments are short and fail the
  reciprocal-cov bar) + 33 rRNA.

This is the defensible "real protein-coding family" claim: members are coding, expressed (cDNA
homology), AND whole-protein homologous — and it breaks the over-merge to a real-family max-comp on
PRINCIPLE, not heuristic, with trustworthy positives.

## 4. Known limitations (separately handleable)
- **Immune-receptor families** (IG/TCR V/C/D/J segments) are missed — short peptides below the
  reciprocal-cov bar. Need a segment-aware bar.
- **Partial paralogs** (GGT1~GGTLC2) are missed — protein homology sub-bar (GGTLC2 covers only part
  of GGT1). The recurring partial-duplication boundary.
- **Non-coding (lncRNA) families** are a separate phenomenon (repeat-bridged) — not protein families;
  handle on the cDNA axis if of interest.

## Status / recommendation
- Protein orthogroups = the family ground truth for coding genes; use them directly.
- The cDNA de-bridge (co-threading `frac_ref`, NOT cov_min) is the proxy where protein is unavailable
  (F1 0.72 vs protein truth).
- lncRNA repeat over-merges are the bulk of the residual blob; they are not protein families.
- Next: segment-aware bar for IG/TCR; decide the partial-paralog policy; optionally a protein-level
  community detection for any protein orthogroup that itself over-merges via shared domains.

---

## ADDENDUM — segment-aware recovery (Task 1) + protein community detection (Task 2)

Script: `family_def_immune_and_protein.py`.

### Task 1 — immune-receptor families recovered (locus constraint)
The 116 IG/TCR segments (V/C) have **0 representative proteins** (not translated until VDJ
recombination), so no protein bar can reach them; and their shared V-DOMAIN over-merges segments
ACROSS loci (one naive cDNA component spanned 6 chromosomes). Fix: segment-biotype cDNA homology +
**SAME-LOCUS constraint** (same chrom, ≤3 Mb) — the inverse of the disjoint-loci rule for dispersed
paralogs. Result: **30 immune families, all 30 single-chromosome compact tandem arrays** (0.06–0.97 Mb
at the IGH/IGK/IGL/TRA/TRB loci). Recovers a whole real-family class the protein definition dropped.

### Task 2 — protein orthogroups DO domain-over-merge; the cDNA bar already handles it
The *pure* protein-orthogroup graph (fident≥0.30) over-merges ancient domain families: ZNF (569,
density 0.34), OR7D (501), GPCR, HOX, kinases — **10/78 large orthogroups** split under protein-recip-
cov-weighted Louvain. BUT the cDNA∩protein deliverable is already clean, because the strict cDNA
id≥0.90 bar splits each ancient domain blob into its RECENT sub-duplications:

| family | protein-alone | cDNA∩protein (recent paralogs) |
|---|---|---|
| ZNF | one 569-gene blob | 43 genes → families of 4,3,2,2,… |
| OR / GPCR | 501 / 231 blobs | 2 genes (the rest are cDNA-divergent ancient) |
| KRT | 78 blob | 9 genes → families of 2,2,2 |

So the two definitions capture different things and that is the right design:
- **cDNA∩protein = RECENT coding paralog families** (cDNA still ≥0.90) — clean, no domain over-merge,
  and exactly the recent-duplication / copy-assignment regime the thesis targets. 837 families, max 47.
- **protein-alone = all homology incl. ancient domain families** — over-merges via domains; needs
  community detection (the Task 2 tool) if ancient families are wanted.

## CONSOLIDATED family definition (layered)
1. **Recent coding paralog families**: cDNA-homology (id≥0.90) ∩ protein-homology (recip-cov≥0.50) →
   **837 families**, max comp 47, retrocopy- and over-merge-free.
2. **Immune-receptor families**: segment biotype + cDNA homology + same-locus → **30 families**.
3. **Ancient domain families** (ZNF/OR/GPCR/HOX): protein-level only; appear as many recent
   sub-families in (1); recoverable as super-families via protein-graph community detection if wanted.

Total **867** trustworthy families, each with an explicit, biotype-appropriate criterion — replacing
the 1,216-family / 238-blob cDNA over-merge.


---

## family_def_scaffold_wiring

# Wiring the DNA-first scaffold into the --vg pipeline — what it revealed

Goal: promote the DNA-first scaffold (family_def_dna_scaffold.py) to the production pipeline
via `--family-manifest` (template-based family assembly) and re-measure completeness end-to-end.

## Mechanically wired
- Generated the manifest from cDNA all-vs-all paralog families: `make_dna_family_manifest.py`
  -> `/home/juanfra/winloci_scratch/dna_family_manifest.tsv` (1,460 families, 5,517 loci).
- Invocation (per-chrom, OOM-safe): `rustle --vg --family-manifest <tsv> --genome-fasta <fa>
  -L <bam> -o <gtf> --vg-report <tsv>`. The flag works; bundles are matched to manifest families.

## But end-to-end recovery does NOT work cleanly — two root causes

1. **The pipeline drops manifest families by default.** `--vg` does not ingest secondary
   alignments into bundles (deliberate gate, opt-in `RUSTLE_VG_INCLUDE_SECONDARY=1`), so the
   shared multimap reads that link copies are absent -> every family is `low_shared` -> 54->0
   on the test chrom NC_073240.2. The `--vg-family-min-shared` filter is a de-novo PRECISION
   gate, redundant+harmful against a DNA-validated manifest. Even with secondaries on + filters
   relaxed, only 5/54 formed -- and with 0 shared reads (no assignment possible).

2. **The DNA scaffold is ITSELF over-merged (the decisive finding).** The cDNA all-vs-all
   paralogy (id>=0.90, recip-cov>=0.30, connected components) produces mega-"families": max
   = 389 members; 19 families >20 members hold 23% of all 5,517 member genes. The manifest
   faithfully grabs every bundle for these -> the pipeline emits an over-merged 21-copy
   "family" spanning 88 Mb. So the "DNA backwards" ground truth has the SAME bridge /
   transitive-homology problem as RNA detection -- it is not a clean external truth.

## The deep conclusion
The bridge / over-merging problem is **fundamental to homology-based family definition** --
it appears in the DNA truth (389-member components), not just RNA. "DNA backwards" does not
escape it; it relocates it. So the completeness measurements (which used these DNA families)
are inflated by the mega-families, and the scaffold cannot be cleanly wired as-is.

## Path to realize the scaffold (future work)
1. **De-bridge the DNA families** before using as a scaffold: apply the same cleaning the RNA
   side uses (contiguous-core / reciprocal-coverage / clique structure) to the cDNA all-vs-all
   components, so the manifest carries clean families, not 389-member mega-components.
2. **Add a scaffold mode to --vg**: when a manifest is provided, bypass the de-novo
   `low_shared` precision filter (the manifest is the precision gate) and use the paused
   Phase-2 secondary side-index so shared reads link copies without contaminating bundles.
Until both are done, the DNA-first scaffold is a sound prototype-level idea whose production
realization is blocked by (a) the un-cleaned DNA truth and (b) de-novo-tuned pipeline mechanics.


---

## family_def_unbiased_differentiators

# Unbiased differentiators of REAL family vs OVER-MERGE/bridge edges

**Question.** Beyond the two known levers (junction-position co-threading; same-strand filter),
are there OTHER *natural* differentiators — especially graph-topology ones (edge betweenness) —
discoverable label-free on the full 12,212-edge both-spliced cDNA-homology population? Found via
label-free intrinsic bimodality, oriented only with sparse independent labels (12 compara_pos
positives; antisense/overlap negatives), and verified against the two killer confounds
(component-size and node-degree) plus a whole-component grouped null.

All numbers below were re-measured with sklearn 1.5.2 in the sqanti3 env on
`bench/family_def_features.tsv` (scripts: `family_def_bimodality.py`, `family_def_edgebetw_probe.py`,
`family_def_confound_resolve.py`, `family_def_residual_bimodality.py`, `a2_unsupervised_structure.py`,
`a2_edgebetw_orient.py`, `a2_confound_confirm.py`, `a2_refute.py`, `family_def_refute_edgebetw*.py`,
`family_def_refute_lenmax.py`).

---

## BOTTOM LINE

**No new differentiator beats co-threading once the component-size / node-degree confound is
removed.** The two candidates that looked new both fail the orientation+confound test:

- **edge betweenness (`edge_betw`) is node DEGREE in disguise**, not an independent bridge signal.
  Its headline AUC vs compara (0.683 raw, 0.939 in an earlier unstratified run) is a near-pure
  component-size artifact: corr(edge_betw, log component-size) = **−0.71**, spearman(edge_betw,
  deg_min) = **−0.78**. Residualize on size → AUC **0.531** (chance). Inside the big over-merge blob
  the raw bridge tail *does* track artifact edges (AUC overlap **0.78**, antisense **0.664**), but
  **degree alone does it better** (deg_min AUC overlap **0.871**, antisense **0.862**) and
  degree-residualized edge_betw collapses to **0.574 / 0.282** — i.e. it adds ≈0 beyond degree.
  And it **cannot be oriented by any trustworthy positive: 0/12 compara_pos and 0/5 panel_real lie
  in any component ≥50 edges.** Verdict: confounded, strictly dominated by degree. Not worth adding.

- **`len_max` is the most intrinsically bimodal feature** (BC=0.866, silhouette=0.799,
  dBIC 1→2 = +28,529, survives size-residualization BC_resid≈0.86) and is genuinely new and
  non-circular — **but it does not orient.** It is a short-cDNA vs long-cDNA axis, not real-vs-bridge:
  AUC vs compara only **0.317** (anti-oriented vs the artifact-negative scale), and **10/12 compara
  positives sit in the LOW-len mode.** Bimodal but label-irrelevant; would be a confound if used as a
  family filter.

**The strongest CLEAN lever remains co-threading** (`co_junc`/`contig`/`frac_mem`/`frac_ref`,
AUC ≈ 0.69–0.70 vs compara), which is **size-independent** (corr(co_junc, log-size) = −0.10) and
**orientable without circularity**. Same-strand (antisense=0) stays the orthogonal second lever.

---

## Ranked verdict table

| Feature | Intrinsic bimodality | AUC vs compara (raw) | After size/degree control | Orients? | Verdict | Confound / note |
|---|---|---|---|---|---|---|
| **co_junc** (co-thread) | BC 0.59, sil 0.71 | 0.695 | 0.70 (size corr −0.10, edge_betw corr ≈0) | yes | **KNOWN clean baseline** — best clean lever | none; the bar to beat |
| **contig / frac_mem / frac_ref** (co-thread) | BC 0.53–0.88 | 0.692–0.695 | size-independent | yes | **KNOWN** baseline | circular only for community-derived labels (not used) |
| **antisense / overlap** (strand/co-loc) | n/a (binary) | independent NEG | — | NEG by def | **KNOWN** 2nd lever | orthogonal to co-threading + topology |
| **edge_betw** (topology bridge) | BC 0.94 (size-driven) | 0.683 (0.939 unstratified) | **0.531** size-resid; **0.574/0.282** degree-resid in-blob | no (0/12 pos in blobs) | **CONFOUNDED — degree in disguise** | comp-size (−0.71), node degree (spearman −0.78); deg_min alone wins (0.87/0.86) |
| **common_nbr** (topology) | BC 0.62 | 0.914 → | **0.565** size-resid | weak | **CONFOUNDED** | component size (corr +0.64) |
| **deg_max / deg_min** (topology) | BC 0.45 | 0.932 / 0.926 → | 0.728 / 0.626 size-resid | weak | **CONFOUNDED** | component size (corr +0.81 / +0.75); also anti-orients in-blob |
| **jaccard_nbr** (topology) | BC 0.65 | 0.512 (perm-p 0.87) | 0.488 | no | **NOISE** | size-driven bimodality, no label signal |
| **len_max** | **BC 0.87, sil 0.80** (most bimodal) | 0.630 / **0.317** | survives resid but still doesn't orient | **no** | **NEW bimodal but NON-discriminative** | cDNA length (spearman +0.64 len_min / −0.67 len_ratio; −0.34 seq_id) |
| **len_min / njunc_*** | BC 0.49–0.60 | 0.50–0.71 | — | mixed | bimodal, weak/none | label-irrelevant bimodality |
| **log_dist** | BC 0.28, sil 0.76 | 0.787 (same_chrom only) | — | yes | KNOWN strand-genomic | only defined on same-chrom edges |

---

## (1) NEW genuine differentiators

**None.** No feature is simultaneously (label-free bimodal) AND (orients to compara) AND
(non-circular) AND (not a pure size/degree confound). The closest structural candidate
(edge_betw within blobs) fails the degree test; the closest bimodal candidate (len_max) fails
orientation.

## (2) Confirmed-but-known (co-threading / strand)

- **Co-threading** (`co_junc`, `contig`, `frac_mem`, `frac_ref`): AUC ≈ 0.69–0.70 vs compara,
  size-independent (|corr| ≤ 0.10), essentially uncorrelated with edge_betw — the clean lever.
- **Strand / co-location** (`antisense`, `overlap`): the independent negative axis; orthogonal to
  both co-threading and topology by construction.
- `log_dist` (closer = more real) confirms but is restricted to same-chromosome edges.

## (3) Confounded / circular (downgraded by verifiers)

- **edge_betw, common_nbr, deg_min, deg_max, jaccard_nbr** — all collapse after residualizing on
  **log component-size** (the dominant confound: small components mechanically force betweenness→1.0
  and degree→low). The component-size axis is itself partly the co-threading community-detection
  output, so "small comp = real" risks the documented circularity — another reason to distrust the
  raw topology AUCs.
- **len_max** — bimodal and non-circular but confounded with **cDNA length** (and indirectly with the
  known anti-discriminative `seq_id`, whose own AUC 0.659 ≥ len_max 0.630). Does not separate
  real-vs-bridge.

---

## Edge betweenness specifically — VERDICT WITH NUMBERS

**Independent alignment-free bridge signal, or just degree in disguise? → DEGREE in disguise.
Do NOT add it as a standalone differentiator.**

| Test | Number |
|---|---|
| corr(edge_betw, log component-size) | **−0.71** |
| spearman(edge_betw, deg_min / deg_max) | **−0.78 / −0.77** |
| AUC vs compara, raw | 0.683 (0.939 in unstratified run) |
| AUC vs compara, **residualized on log size** | **0.531** (chance) |
| Whole-group permutation null (positives span ~10 groups) | **p = 0.176** (the global edge-perm p=5e-5 ignores grouping → over-optimistic) |
| compara_pos / panel_real inside components ≥50 edges | **0 / 12 and 0 / 5** (cannot be oriented by any real positive in the blob) |
| Within-blob raw AUC: overlap / antisense | 0.78 / 0.664 |
| Within-blob **degree-residualized** AUC: overlap / antisense | **0.574 / 0.282** (adds ≈0; antisense reverses) |
| Within-blob **deg_min alone** AUC: overlap / antisense | **0.871 / 0.862** (degree dominates) |
| Within-blob grouped-block permutation p: overlap / antisense | **0.001 / 0.418** (only overlap even significant) |
| Within-blob single-feature GroupKFold (raw ranking transfer): overlap / antisense | 0.86 / 0.77 (transfers as a *ranking* — but it is the degree ranking; marginal contribution beyond degree ≈ +0.00) |

Reading: the *raw* within-blob bridge tail does correlate with co-location/antisense artifact edges
and transfers across held-out components — but only because those artifact edges are the low-degree
periphery of the giant component. Once degree is controlled, edge_betw's contribution vanishes
(overlap +0.00, antisense reverses below chance). It is non-circular (independent of co-threading:
corr(edge_betw, co_junc) ≈ −0.06, corr(edge_betw, frac_ref) ≈ −0.06) but **strictly dominated by
node degree**, which is cheaper to compute. **If anything is worth a probe inside blobs, it is
`deg_min` (peripheral low-degree edges = artifact candidates), not edge betweenness** — and even that
is unoriented by any real positive and should be treated as a heuristic, not a differentiator.

---

## How a (hypothetical) new differentiator would complement the known levers

If a clean topology signal had survived, its value would be **orthogonality + cost**: topology is
computed from the homology graph alone (no per-read realignment), whereas co-threading needs the
read-to-locus threading and strand needs orientation. The measured orthogonality holds
(edge_betw ⟂ co_junc, |corr| ≈ 0.06) — the problem is not independence, it is that the only
size-independent residual signal is **degree**, which is a confound rather than a real-vs-bridge
axis, and which has **no trustworthy positive inside the blob to orient it**. So the cheap
alignment-free graph features buy nothing over co-threading here.

---

## Recommended next step

Treat topology as a confound, not a lever: keep co-threading (`co_junc`/`contig`) + same-strand as
the family-edge gate, and — if a blob-internal heuristic is wanted — try **dropping low-`deg_min`
peripheral edges inside components ≥50** as a *candidate-flagging* step only, validated against
overlap/antisense (not compara, since 0 positives live there). To find a genuinely new differentiator,
the bottleneck is **labels inside the blobs**: get a handful of trustworthy positives (e.g. Compara
paralog pairs or curated tandem arrays) that actually fall in components ≥50 edges, otherwise every
in-blob signal is orientable only by artifact-negatives and indistinguishable from degree.

---

## ADDENDUM — the de-bridge test the AUC framing missed (cov_min, jaccard_nbr)

The verification above ranked features by "AUC vs the 12 Compara positives". That test is
*uninformative inside the over-merge blobs* (0/12 Compara positives and 0/5 panel reals land in any
component ≥50 edges), so it cannot adjudicate the actual goal — de-bridging the blobs. Re-ran the
PRACTICAL test instead: plug each cheap feature in as the community-detection edge weight (same
weighted-Louvain as co-threading) and measure max-component reduction + real-family preservation.

| weight (cost) | families | max comp | reals intact | anchor sep (238/140) |
|---|---:|---:|:--:|:--:|
| base (over-merged) | 1216 | 238 | — | — |
| frac_ref — co-threading (needs minimap2) | 1286 | 57 | 5/5 | 3/3, 2/2 |
| **jaccard_nbr — graph topology (FREE)** | 1263 | 58 | 5/5 | 1/1, 1/1 |
| **cov_min — reciprocal coverage (FREE, already computed)** | 1280 | **44** | 5/5 | 3/3, 2/2 |
| common_nbr — graph topology (free) | 1252 | 68 | 5/5 | — |

**`cov_min` de-bridges BETTER than co-threading (max comp 238→44 vs 57), keeps all 5 real families,
separates the unrelated anchors (3/3, 2/2) — at ZERO extra compute** (it is just `min(cov_a,cov_b)`
from the existing all-vs-all; `dna_homology()` already loads both, we only used the `max` for the
gate). `jaccard_nbr` (pure graph topology, no alignment) matches co-threading (58 vs 57).

Mechanism (Canzar-clean): real paralogs cover each other RECIPROCALLY (both cDNAs mostly align);
domain-sharers/bridges cover ASYMMETRICALLY (the shared domain is most of the small gene but little
of the large one → low cov_min). This is the continuous form of the airtight panel's 3rd axis
(recip-cov ≥ 0.30), now used as a de-bridge weight rather than a gate.

**Honest scope:** these three agree (Spearman +0.60–0.65) — they converge on the same dense-subgraph
structure, so they are NOT an independent validation axis (the synthesizer's point stands: nothing
adjudicates whether the LOC blob sub-communities are *real* families; that needs orthogonal
protein/Compara evidence). The win is practical, not epistemic: **cov_min is a cheaper, alignment-free
de-bridge weight that is at least as good as co-threading** — and ~35–40% of its rank variance is
independent of co-threading, so multiplying the two weights is worth trying.

**Updated recommendation:** switch the community-detection weight to `cov_min` (or `cov_min × frac_ref`)
— cheaper and a tighter de-bridge — keeping same-strand as the orthogonal filter. The remaining
bottleneck is unchanged: trustworthy positives inside the big components (protein-level) to confirm
the sub-communities are real families.


---

## family_def_vg_coherence_results

# VG coherence check + non-linear structure (are we using the real VG?)

`family_def_vg_coherence.py` — per family, align every member to the longest member (minimap2 -c -N 10),
parse ALL blocks, and compute (A) member-side backbone-threading coherence and (B) the NON-LINEAR
signatures the collinear pipeline is blind to. 826 recent-coding families, 2087 member alignments.

## (A) Coherence — mostly redundant with the protein bar
Member-side forward threading <50% in 14% of families — but this is largely **UTR-driven** (cDNA
includes UTRs; the protein whole-gene bar already validated the CDS). So coherence at the cDNA level
adds little over the protein bar; it is not a clean independent bridge detector. (The ~10 true
low-density bridges from the FP audit remain the residual; protein homology already gates them.)

## (B) NON-LINEAR structure — the real finding (VERIFIED)
The whole family-definition pipeline is COLLINEAR (cDNA/protein all-vs-all, co-threading, and we even
dropped reverse-strand edges). A real variation graph models inversions (reversing traversal), tandem
repeats (CYCLES), and exon shuffling (reordered shared nodes). We use none of it. How often do these
actually occur in our families?

| non-linear case | families | % | status |
|---|--:|--:|---|
| INVERSION (≥10% reverse block within family) | 0 | 0% | none at this level (inverted gene-dups are sense-sense in cDNA) |
| TANDEM REPEAT (backbone region matched ≥2×) | 20 | 2% | **REAL — verified** |
| EXON SHUFFLE (≥20% out-of-order blocks) | 16 | 2% | real rearrangement + some repeat-driven |

**Verification (RCF_13, all-LOC tandem-repeat family):**
- Self-alignment of members shows genuine INTERNAL tandem repeats (6 and 4 off-diagonal self-repeat
  blocks) — the repeat is in the sequence, not a family-alignment artifact.
- Building the **actual `vg msga` graph** for the family yields a **CYCLIC graph** (198 nodes, 218
  edges, `cyclic`) — a cycle *is* tandem-repeat structure. The real VG represents what our collinear
  homology flattens.

## Implication
~4% of families carry non-linear (repeat / rearrangement) structure that our collinear definition
flattens to a single backbone. This matters for the thesis: a tandem-repeat COPY-NUMBER difference
between copies is a structural distinguishing feature — a PSV-like signal our SNP/indel-bubble PSV
model misses. Modeling repeats/CNV in the family VG would add copy-distinguishing power (move the
identifiability frontier), and the `vg` toolkit (installed, cyclic-graph-capable) is the right engine.

NEXT: a full non-linear census (vg graph per flagged family, classify cycle vs rearrangement vs
artifact); fold structural (repeat-CNV) variants into the PSV set for copy-assignment.


---

## family_definition

# Annotation-based RNA-level multi-copy gene family definition (completeness)

**Definition:** a multi-copy gene family = a maximal connected set of ANNOTATED genes whose representative transcripts pairwise share a POA contiguous exon core ≥ 0.13 at core identity ≥ 0.7. Roster = the 22,983 annotated gene models; grouping = the validated POA criterion; families = connected components. Annotated-only (extensible: the same criterion runs against de-novo/unannotated loci later to add unannotated copies).

- homology edges: 16,354 ; **families (components, size≥2): 1,337** ; genes in a family: 5,060
- component size distribution (top): [250, 114, 72, 47, 45, 45, 34, 32, 29, 29, 27, 26]

## Completeness — are KNOWN multi-copy families recovered?
- **curated textbook families: 8/11 recovered** (members land in one component) ; missed: CRYBG, DEFB*, SIGLEC*
- similarity-built universe families: **46/53 recovered** ; missed (sample): ASDURF, CASP8, CDPF1, CREB1, GCA, GPR39, LOC129529456

### curated families & their recovery
| family | members (annotated) | in one component? |
|---|---|---|
| APOBEC3 | 6 (APOBEC3B, APOBEC3C, APOBEC3D, APOBEC3F, APOBEC3G, APOBEC3H) | YES (5/6) |
| CRYBG | 3 (CRYBG1, CRYBG2, CRYBG3) | no |
| DAZ | 2 (DAZ1, DAZL) | YES (2/2) |
| DEFB* | 20 (DEFB1, DEFB108B, DEFB110, DEFB112, DEFB113, DEFB115…) | no |
| GGT | 6 (GGT1, GGT5, GGT6, GGT7, GGTLC1, GGTLC2) | YES (3/6) |
| MAGEA* | 5 (MAGEA1, MAGEA10, MAGEA12, MAGEA4, MAGEA9) | YES (3/5) |
| PRAMEF* | 4 (PRAMEF12, PRAMEF17, PRAMEF19, PRAMEF20) | YES (2/4) |
| RABL2 | 2 (RABL2A, RABL2B) | YES (2/2) |
| RFPL | 4 (RFPL1, RFPL2, RFPL3, RFPL4A) | YES (3/4) |
| SIGLEC* | 6 (SIGLEC1, SIGLEC10, SIGLEC12, SIGLEC15, SIGLEC5, SIGLECL1) | no |
| TAS2R* | 18 (TAS2R1, TAS2R10, TAS2R13, TAS2R14, TAS2R16, TAS2R20…) | YES (2/18) |

## Over-merge control (the precision side of 'complete')
- mega-components (size≥25, likely domain-hub chains, NOT single families): 14 (largest 250)
- these are the transitive-closure artifacts; a tighter clustering (mutual-core / community detection) would split them. Size-2/3 components are clean copy sets.

## Honest scope
- ANNOTATED genes only (this definition); unannotated/de-novo copies are added later by running the SAME POA criterion against de-novo loci (the cross-chrom discovery already does this).
- RNA-structural operational family (shared contiguous exon core), NOT the formal gene-tree/protein family — a parallel DNA tier would supply that (and the ground truth).
- universe is itself similarity-built (recovery there is partly internal consistency); the curated textbook set is the independent completeness check.


---

## family_definition_artifact_filter

# Refining the multi-copy-family definition — exclude ARTIFACTS, not non-coding (the VG does it)

**Goal correction (user):** the target is to drop **artifacts** (repeats, low-complexity, TEs,
chimeras, over-merges), NOT to require coding. A real **lncRNA multi-copy family is a real family** —
requiring an ORF would wrongly delete it. So the filter must be coding-agnostic + intrinsic.

## Two things that DON'T work (measured, `refine_family_definition.py`)
- **ORF is a bad filter.** On de-novo loci labelled by biotype: **76% of lncRNA have an ORF ≥100 aa**
  (spurious ORFs in long transcripts) and 6% of *coding* loci don't. ORF *keeps* lncRNA and *drops*
  coding — it conflates non-coding with artifact.
- **k-mer (4-mer) complexity SATURATES.** Median 0.988 across all biotypes (only 256 possible 4-mers →
  any transcript >300 bp hits the ceiling); the repeat artifacts also score 0.96. Not discriminative.

## The fix: artifacts are TOPOLOGICAL — the variation graph exposes them (no ORF, no annotation)
The same VG object that *defines* the family and *assigns* the reads also *filters* the artifacts —
the thesis's graph does triple duty, with no ad-hoc thresholds or coding assumptions.

1. **Repeat / tandem / TE → a CYCLE in the variation graph.** A real family is an *acyclic* bubble-chain
   (backbone + localised PSV bubbles); a repeat makes the graph loop back on itself. Intrinsic signal =
   **long-k-mer self-recurrence** within a transcript. Measured (`15-mer recurrence`, 2,000 de-novo tx):
   median **0.000** (real = unique k-mers, acyclic), tail to **0.57** — flags internal-repeat
   transcripts (e.g. `DN_NC_073224.2_130224868_27`, 4 kb, 57% recurrence) that 4-mer complexity (0.96)
   called "complex." **Sensitive exactly where the naive measure saturates.**
2. **Over-merge / repeat-bridge → a SPARSE, non-clique component.** A real family is a *clique* in the
   ~B backbone graph (every pair aligns reciprocally); an artifact bridge is a sparse connector through
   an articulation point (the comp-238 lncRNA blob). Intrinsic signal = graph density / articulation
   points / community detection (partly shipped: 238→57).
3. **Chimera → an INCOHERENT backbone.** A chimeric transcript's pieces map to unrelated loci → no
   single colinear spine in the VG. Intrinsic signal = backbone coherence (one dominant colinear path).

## Refined definition (coding-agnostic, intrinsic, topological)
> A multi-copy gene family = **~R ∩ ~B** AND the family's variation graph is **(a) acyclic** (no
> internal-repeat cycle — `k-mer self-recurrence < ~0.10`), **(b) a clique** (dense reciprocal ~B, not a
> sparse articulation bridge), and **(c) copy-count-bounded** (2 ≤ copies ≤ ~12; dispersed TEs make many
> copies). No ORF requirement → real lncRNA families are KEPT (they are acyclic, clean cliques);
> repeats/TEs/chimeras/over-merges are DROPPED on graph structure alone.

This keeps the definition purely topological (interest I) and unifies it with copy-assignment
(interest II): one variation graph, three jobs — define, assign, and filter — no coding model, no
arbitrary sequence thresholds.

> **Scope of the "clique" / "triple-duty" claim (L2 — now partly closed).** The clique density bar is
> enforced as a **DROP on BOTH** pipelines: the DNA cDNA-homology manifest (`make_dna_family_manifest.py`,
> `overmerge_sparse` at `n≥4 ∧ density<0.30`) AND the de-novo split path
> (`src/rustle/vg_family/family_split.rs` + `denovo_family_split.py`) — whose `web` density bar is now
> **aligned to 0.30** (was 0.15), and `Web` families are **excluded from copy-assignment**
> (`denovo_pipeline.rs`). Measured effect: the de-novo split goes 695 family / 3 web → **691 / 7**; the 4
> newly-dropped are the large multi-chromosome domain-sharing over-merges (DSFAM0 = a 164-member ZNF
> across 19 chromosomes, plus ZNF/APOL families) — the exact blobs that embarrassed "family = clique."
> `web_min_size` stays 10 on the de-novo side (vs the manifest's 4) **on purpose**: a SMALL sparse group
> can be a real divergent family (e.g. a 7-copy single-chrom MAGEB at density 0.24), so only
> *large-and-sparse* is dropped. Still partial: only the density/cliqueness leg is wired on the de-novo
> side — the k-mer-self-recurrence (cycle) leg remains cDNA-manifest-only (L5).

## Wired into the pipeline (`family_def_artifact_filter.py` → `make_dna_family_manifest.py`)

SHIPPED as a default-on filter with an opt-out (mirrors the retrocopy filter). Per family it computes
**cliqueness (density), copy-count, and k-mer self-recurrence** (genomic, capped 20 kb) and drops:
`repeat_cyclic` (recur ≥ 0.15), `overmerge_sparse` (n ≥ 4 and density < 0.30), `te_highcopy` (a high
backstop, n > 80, dense). **DENSITY — not copy-count — is the discriminator**: a real family is a dense
clique *at any size*, so a 41-copy clique (density 1.0) is KEPT and only sparse bridges are dropped.

Result on real data: **1,216 → 1,167 real families; 49 dropped** (39 repeat_cyclic + 10
overmerge_sparse — incl. the comp-238 lncRNA-repeat blob, density 0.117). Large real cliques (n =
40–47) correctly KEPT. Every family is written to the sidecar
`dna_family_artifact_class.tsv` (`family_id  class  n  density  recur  kept`).

**Opt-out / visualise:** `FAMILY_DEF_NO_ARTIFACT_FILTER=1` keeps *all* families (artifacts included) but
still tags each with its class — so you can pull a `repeat_cyclic` / `overmerge_sparse` / lncRNA family
from the sidecar and **show its graph** (a repeat family's variation graph has *cycles*; an over-merge
is a *sparse bridge*, not a clique) next to a clean `real` family. Thresholds are env-tunable
(`FAMILY_DEF_MIN_DENSITY`, `FAMILY_DEF_REPEAT_RECUR`, `FAMILY_DEF_MAX_COPIES`). Coding-agnostic
throughout — real **lncRNA** families pass (acyclic clean cliques) and are KEPT.

Artifacts: `family_def_artifact_filter.py` (filter + self-test) · `make_dna_family_manifest.py` (wired) ·
`dna_family_artifact_class.tsv` (sidecar) · `refine_family_definition.py` (the ORF/4-mer negative).


---

## family_definition_formal

> **⚠ Superseded as the O1 family *definition* (2026-06-30) — see `bench/DEFINITIONS_FORMAL.md`.** The read-conflict graph $E_c$ described below is an **ambiguity** oracle (an edge = a read that *cannot resolve* two loci), not a **homology** oracle: a multi-copy family whose copies are divergent enough to each map uniquely produces no de-tie and vanishes (measured: **~30% of multi-copy homology families / ~1/4 of copies** silently dropped). The canonical RNA family is now the **transcript-homology component $E_r$** (`bench/DEFINITIONS_FORMAL.md`), the fourth homology oracle ($E_a$—$E_b$—$E_r$—$E_p$), which *includes* the read-resolvable copies. **Everything below is correct and retained as the within-family O2 (copy-assignment) structure**: $E_c$ / `MCC=χ(H)` / the SUN 3-tier ladder are demoted to *how copies are assigned inside a fixed $E_r$ family*, with $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$. Results ledger: `bench/DEFINITIONS_FORMAL.md`.

# A formal definition of a multi-copy gene family at the RNA level — and its proof on IsoSeq reads

*Advisor interest #1: can we give a full formal definition of a multi-copy gene family and find them from
IsoSeq reads, with an honest false-positive / false-negative accounting?*

**Bottom line.** A multi-copy gene family is defined as a connected component of the **read-conflict graph**
under a divergence-tie (`de`) edge criterion. On the GGO HiFi IsoSeq substrate it reproduces a labelled
17-candidate panel exactly — **TP=7, TN=10, FP=0, FN=0 (precision = recall = 1.000)** — with every
load-bearing count independently re-derived from raw `samtools`/`pysam`, the three formal properties verified,
and the principled `de` criterion demonstrably avoiding three `AS`-tie false positives (including a 3347-read
retrocopy). Independently, **Ensembl Compara corroborates the annotated positives (RABL2, MAGEA) as recent
within-species paralogs**, and confirms the definition is a *strict subset* of recent paralogy: APOBEC3 is a Compara
recent paralog yet a correct TN here, because its copies are read-resolvable — the conflict-family is recent-paralogy
∩ read-confusability, the unit the copy-assignment problem needs. The single most important caveat: the clean ledger holds for families with **≤6 cross-mapping
copies** — `minimap2`'s default secondary cap (N=5) truncates the input and fragments a real ≥12-copy NC_086018.1 (chr23)
array; this is a false negative of the *input*, not the definition, fixed by re-aligning with uncapped
secondaries.

Artifacts: `bench/family_definition_demo.py` (+ `.tsv`/`.json`), `bench/family_definition_figure.png`,
verification workflow `wf_6aa71f9e-c3d` (4 adversarial checks + synthesis). Definition predicate is the shipped
Rust criterion (`read_conflict.rs`, de-tie, commits `2e9922e..2da06c2`).

---

## The read-conflict graph and the family object

Let the substrate be a set of long reads $R$ (GGO HiFi IsoSeq, polyA-selected) aligned to a reference with
`minimap2` (primary + secondary records retained), each alignment carrying the gap-compressed divergence tag
`de:f`.

**Vertices.** $V$ = de-novo *expressed* loci — the per-transcript intervals the upstream `detect_and_assign`
pipeline emits, not annotated gene spans. (This vertex resolution is load-bearing; see P1 and Limitations.)

**Record-attributed placement (no distance guard).** Each alignment *record* $\rho$ of a read $r$ (primary or
secondary; supplementary/chimeric records excluded as they are split-read pieces, not alternative placements) is
attributed to the single locus it best overlaps, $\mathrm{loc}(\rho)=\arg\max_i \mathrm{overlap}\big(\rho,[s_i,e_i]\big)$.
The read's *placement set* is $P(r)=\{(\mathrm{loc}(\rho),\,de(\rho)) : \rho \in \mathrm{records}(r)\}$ — at most one
entry per record, $de$ from the record's `de:f` tag. **Consequence (the principled replacement for the panel
demo's $<200$ bp guard):** two entries of $P(r)$ on *distinct* loci necessarily come from *distinct alignment
records* — a genuine multimapping — **by construction**; a single alignment spanning nested loci is attributed to
one locus and yields one entry, so it can never self-conflict. There is no distance threshold. This is the
shipped Rust formulation (`build_read_placements`), which supersedes the demo's coordinate guard.

**De-tie conflict predicate.** Fix $\Delta = 0.005$, $\mathrm{DE_{max}} = 0.05$. A read $r$ *conflicts* on the
unordered pair $(i,j)$, $i\ne j$, written $\mathrm{conf}(r,i,j)$, iff $P(r)$ contains entries $(i,de_i)$ and
$(j,de_j)$ whose divergences are **tied at the HiFi error floor**:
$$
|de_i - de_j| \le \Delta \quad\wedge\quad \max\!\big(de_i,\,de_j\big) \le \mathrm{DE_{max}}.
$$
The read fits *both* loci comparably; `minimap2` cannot decide. This uses raw divergence `de`, **not** the
aligner's composite score `AS` — the latter folds in length and is the source of the avoided false positives.

**Edges.** With $\mathrm{MIN\_READS}=3$,
$$
(i,j)\in E \iff \big|\{\,r\in R : \mathrm{conf}(r,i,j)\,\}\big| \ge 3.
$$

**Family.** A **multi-copy gene family** is a connected component $C$ of $G=(V,E)$ with $|C|\ge 2$. This is the
*identifiability-relevant* family — the maximal set of loci among which reads are genuinely confused, exactly
the unit on which the downstream copy-assignment problem (Canzar-style conflict resolution) operates. It is
explicitly **not** a claim about evolutionary paralogy: a true paralog whose reads place uniquely (RFPL1/2/3)
is correctly excluded, and a retrocopy whose reads are decidable (EEF1A1) is correctly excluded.

**Why this definition is airtight (three structural robustnesses).** The conflict-graph object is immune *by
construction* to the three artifacts that defeat a sequence-similarity family definition (POA contiguous-core,
`family_detect`):

- *Over-split fragments.* Near-identical isoform-variants the locus collapse leaves at one genomic position
  (e.g. the five PRNP 14.60–14.615 Mb fragments) cannot form a family: each record attributes to its single
  best-overlap locus, so co-positioned fragments share no conflicting read. A *similarity* definition does group
  them — measured on GGO, **~42 % of the de-tie similarity "families" (495/1,190) were one over-split locus**,
  needing a separate output-level member-merge (commit `19b348d`); the conflict definition needs none.
- *Domain-sharers.* Loci sharing only a protein domain are crossed by *uniquely-placing* reads, which produce no
  conflict edge (validated: 0 conflict on 7/7 Compara domain-sharers).
- *Retrocopies / decidable paralogs.* A read that fits one copy decisively (large $|de_i-de_j|$, or $de$ above
  $\mathrm{DE_{max}}$ on the worse copy) does not tie, so EEF1A1's 3347-read retrocopy is excluded — where the
  composite `AS` score, folding in length/gap penalties, falsely ties it (`de`-tie $\subsetneq$ `AS`-tie, the
  shipped regression invariant).

The only free parameters are the two principled de-tie constants $\Delta,\mathrm{DE_{max}}$ (disclosed below);
the former $<200$ bp coordinate guard is **eliminated** — it was an artifact of per-locus (not per-record)
placement and is structurally unnecessary.

## Verified formal properties

**(P1) Domain-sharer exclusion is by construction, not by threshold.** Because $\mathrm{place}(r,\cdot)$ is
single-valued, a read over a shared exon of two nested/overlapping genes has the *same physical alignment
record* selected as best overlap for both loci. Verified on GGO.bam: for CREB1~METTL21A all 192 shared reads,
and GCA~KCNH7 all 429 shared reads, the best-overlap record for both loci is coordinate-identical
($\Delta s = 0$; distinct placements $=0$). They contribute **zero** conflict pairs regardless of $\Delta$.
Domain/UTR sharing cannot create an edge — this is the failure mode of the sequence-similarity definition, here
excluded structurally.

**(P2) Exact decomposition.** Across the 7 recovered families, **0** reads have a de-tied placement pair
connecting loci in two *different* multi-locus components. The assignment problem separates exactly across
families with no shared-read information crossing component boundaries — the property the conflict graph requires.

**(P3) $\Delta$ and $\mathrm{DE_{max}}$ are error-model constants, not tuned similarity thresholds.**

*Derivation (Δ).* $\Delta$ is the single-read divergence-discrimination *resolution* at the HiFi error rate. A
per-read $de$ measured over aligned length $L\approx2.5\,$kb has binomial standard error
$\sqrt{\varepsilon/L}\approx0.0009$ at $\varepsilon\approx0.002$ (the within-family $de$ median); the tie statistic
$|de_i-de_j|$ then has SE $\approx\sqrt{2}\cdot0.0009\approx0.0013$, so two copies whose per-read divergence differs
by less than $\sim\!4\sigma\approx0.005$ are **statistically indistinguishable by a single read**. $\Delta=0.005$ is
therefore the read-level resolution limit set by HiFi error and IsoSeq read length — a measurement constant, not a
similarity knob. *(DE_max).* $\mathrm{DE_{max}}=0.05$ is a deliberately loose divergence ceiling separating
copy-candidate alignments — recent-paralog identity — from distinct-gene alignments; ~3% of within-family $de$
exceed it, an honest, conservative cut, not a tuned boundary.

*Empirical confirmation of the Δ valley.* Within-family per-read
divergence has median $0.0019$, p90 $0.0211$ (the HiFi floor). Within-family conflict-read $|de_i-de_j|$
medians: MAGd2 $0.0000$, MAGd3 $0.0000$, AK6 $0.0012$, MAGd1 $0.0005$, CCDC196 $0.0037$, RABL2 $0.0061$.
Resolvable decoys sit above: CNN2-retro $0.0084$, EEF1A1-retro $0.0168$, APOBEC3 $0.0241$. A $\Delta$-sweep
reproduces ground truth for all $\Delta\in[0.003,0.006]$; first failure at $\Delta=0.007$ (CNN2), then $0.009$
(EEF1A1), then $0.017$ (APOBEC3). $\Delta=0.005$ is centered in a genuine valley with two-sided margin.
*Honesty note:* RABL2's recovery is **quorum-carried, not tie-tight** — its per-read $|de\,\text{diff}|$ median
$0.0061$ exceeds $\Delta$, and only 26% (51/195) of its conflict reads are individually $\le\Delta$; RABL2 is an
edge because 51 reads $\gg 3$ clear the bar. The "cluster near 0" picture is clean for AK6/MAGd2/MAGd3; for the
recent-paralog case it is $\mathrm{MIN\_READS}=3$ that carries it.

---

## Demonstration on GGO.bam — 17-candidate panel

The ledger reproduces **exactly**:
$$\boxed{\text{TP}=7,\quad \text{TN}=10,\quad \text{FP}=0,\quad \text{FN}=0,\quad \text{precision}=\text{recall}=1.000}$$

**Seven families found** (de-tie connected components, $|C|\ge 2$): `{RABL2A, RABL2B}`, `{AK6, LOC115934278}`,
`{CCDC196, LOC129526440}`, `{MAGd1a, MAGd1b}`, `{MAGd2a, MAGd2b}`, `{MAGd3a, MAGd3b}`, and `{MAGd0a, MAGd0b}`
(displayed as `{MAGEA9, MAGd0a, MAGd0b}` — MAGEA9 is the annotation name for the same physical locus as MAGd0b,
a redundant vertex, not a third member; see note).

### Per-family evidence (raw counts re-derived from GGO.bam)

| Candidate | Regime | Truth | de-conflict reads | AS-tie reads | Verdict |
|---|---|---:|---:|---:|---|
| RABL2 (RABL2A~RABL2B) | recent paralog, separate chrom | family | **47** / 195 | 190 | TP |
| AK6 (~LOC115934278) | high-identity cross-chrom copy | family | **24** / 24 | 13 | TP |
| CCDC196 (~LOC129526440) | high-identity cross-chrom copy | family | **24** / 27 | 27 | TP |
| MAGEA_dn0 (MAGd0a~MAGd0b) | co-located array, de-novo loci | family | **103** / 104 | 104 | TP |
| MAGEA_dn1 (MAGd1a~MAGd1b) | co-located array, de-novo loci | family | **83** / 101 | 38 | TP |
| MAGEA_dn2 (MAGd2a~MAGd2b) | co-located array, de-novo loci | family | **303** / 311 | 311 | TP |
| MAGEA_dn3 (MAGd3a~MAGd3b) | co-located array, de-novo loci | family | **75** / 75 | 75 | TP |
| APOBEC3 (D~F) | diverged paralog (resolvable) | not | 0 | 6 | TN |
| RFPL (1/2/3) | paralog, reads place uniquely | not | 0 | 0 | TN |
| EEF1A1_retro (~LOC109023808) | processed-pseudogene retrocopy | not | **0** / 3410 | **3347** | TN |
| CNN2_retro (~LOC129524764) | processed-pseudogene retrocopy | not | **1** / 127 | 121 | TN |
| ASDURF~ASNSD1 | domain-sharer (nested) | not | 0 | 0 | TN |
| CASP8~FLACC1 | domain-sharer | not | 0 | 0 | TN |
| CREB1~METTL21A | domain-sharer | not | **0** / 192 | 0 | TN |
| GCA~KCNH7 | domain-sharer | not | **0** / 429 | 0 | TN |
| GPR39~LYPD1 | domain-sharer | not | 0 | 0 | TN |
| MAGEA_annot (9/4/10/1) | annotated array genes | not | 0 | 0 | TN |

### The contrast: `de` avoids three `AS`-tie false positives

The legacy `AS`-tie criterion (composite alignment-score ratio $\ge 0.9$) forms **10** components, three spurious:

- **APOBEC3** — AS links 6 reads; de = 0 (D/F genuinely resolvable, $|de\,\text{diff}|$ median $0.024$).
- **EEF1A1 retrocopy** — AS links **3347** reads; de = 0. Reads are decidable: EEF1A1 `de≈0.001` @ MAPQ 30–48
  (primary) vs LOC109023808 `de≈0.017` @ MAPQ 0 (secondary), $|de\,\text{diff}|\approx0.016 > \Delta$. `AS`
  collapses this because scores are length-driven; `de` sees the divergence gap.
- **CNN2 retrocopy** — AS links 121 reads; de = 1 (below $\mathrm{MIN\_READS}$).

This is the headline for the advisor: the criterion is principled (raw divergence at the error floor), not a
tuned score ratio — and the principled choice is exactly what kills the 3 high-count false positives.

See `bench/family_definition_figure.png`: (A) the conflict graph with families as components and the avoided
`AS` false-positive edges drawn dashed; (B) the per-candidate ledger.

### Independent confirmation (Ensembl Compara) — and the definition is STRICTER than "recent paralog"

The 17-panel truth labels are human-assigned; an orthogonal check against **Ensembl Compara** (within-species
paralogy, human-mapped, `bench/compara_cache.json`) corroborates the annotated positives and — more importantly —
shows the conflict-family is a *strict subset* of "recent paralog", separated by exactly the identifiability criterion:

- **RABL2A** — one *within-species paralog* at **Homininae** level (= RABL2B): a confirmed recent paralog pair.
- **MAGEA12** — two within-species paralogs at **Catarrhini** level: the MAGEA array is a confirmed recent tandem
  family (our de-novo MAGd0–3 loci are its expressed copies).
- **APOBEC3D/F** — two within-species paralogs at **Catarrhini** level, i.e. a *confirmed recent paralog cluster* —
  yet APOBEC3 is a **TN** in our ledger (reads decidable, 0 de-conflict). This is the crux: the definition does **not**
  recapitulate Compara recent-paralog detection — it is the **identifiability subset** of it. APOBEC3 is a recent
  paralog whose copies are RNA-resolvable, so it is correctly *not* an identifiability family. We answer a sharper
  question than Compara — not "are these paralogs?" but "are these copies confused by the reads?".
- **AK6 / CCDC196 copies (LOC115934278, LOC129526440)** — *absent from Compara* (the LOC paralogues are unannotated):
  copies the de-novo, read-level method finds that the reference annotation/Compara miss — a coverage gain, not a gap.

So the independent check (a) corroborates the annotated positives as genuine recent paralogs and (b) demonstrates the
definition's distinctive value: it is recent-paralogy **intersected with read-level confusability**, the exact unit the
downstream copy-assignment problem operates on.

---

## What the definition does NOT claim, quantified

### FN class 1 — Secondary-cap fragmentation (GENUINE false negative; the honest example)

The definition is sound but its **input is truncated** by `minimap2`'s default placement cap (N=5 → ≤6
records/read). For any array with **>6** cross-mapping copies, each read reports only its 6 nearest placements,
so edges between farther co-members are never observed and one true family **fragments**.

Concrete, reproducible instance — gorilla chr23 (NC_086018.1; mislabeled chrY in an earlier draft — real chrY is NC_073248.2) tandem array at **22.09–22.69 Mb** (~600 kb,
~35 kb period, ≥12 copy vertices, ≥8 expressed): the capped conflict graph reports this single physical array
as **1 component of 8 loci PLUS 4 singletons {V0, V3, V10, V11}**.
- 49 reads hit the 6-vertex cap exactly (1 primary `de=0.0004` + 5 secondaries `de 0.0098–0.0200`, all inside
  the de-tie window).
- Of the 28 expressed-copy pairs, **18 (64%) are never co-observed by any single read**; the two
  highest-expressed copies V0 (129 primaries) and V1 (~267), only 14 kb apart, **share zero reads**.
- Dropped expressed members: **V0, V10, V11** (129 / 26 / 32 primaries) become singletons.
- Cap-causality confirmed: removing the 49 cap-reads further shrinks the component.

The definition is correct; **mitigation = raise `minimap2 -N/-p` (uncapped secondaries) before building $G$.**
This FN is *outside* the 17-candidate panel (no panel family exceeds 6 placements), so it does not contradict
the ledger — but it bounds the demo's generality to **≤6-copy families** unless the input is re-aligned.

> **MEASURED — the re-alignment, run (2026-06-21, `bench/array_secondary_cap_fix.py`).** Earlier this caveat was
> inferred from coverage structure, not run. It has now been run: the 737 array-region reads were re-aligned
> with `minimap2 -ax splice:hq -uf -N 50 -p 0.1`, and the de-tie conflict graph was rebuilt over the same 21
> copy vertices from the capped production BAM vs. the uncapped re-map.
> - **Cap confirmed at the read level:** 119/737 reads sit at exactly 6 records (1 primary + 5 secondaries);
>   uncapping lifts them to up to 52 alignments.
> - **Fragmentation healed:** the array's cross-mapping homology core grows from **5 → 11 copies** (components
>   touching the array 17 → 11; co-observed copy-pairs 10 → 16 / 210). The cap false negative is recovered.
> - **No over-merge (FP = 0):** a far-away expressed gene (RABL2B, 48 Mb) stays an isolated singleton in both
>   graphs; and every one of the 10 remaining uncapped singletons has **0 cross-mapping reads to the core even
>   uncapped** — i.e. it maps *uniquely* and is a genuinely resolvable copy the definition correctly keeps
>   separate, not a cap FN. So after uncapping the identifiability family is exactly the 11 cross-mapping copies.
>
> Conclusion: the secondary-cap FN is an **input artifact, fully removable by uncapped re-alignment at zero FP
> cost** — not a limit of the definition. (The earlier "12 vertices / ~35 kb period" was idealized; the real
> array is a dense homology head ~4–5 kb apart plus sparser, often uniquely-mapping, distal copies.)

### FN class 2 — Unexpressed / single-transcribed-copy families (CORRECT BY DESIGN)

A family in which only one copy is natively transcribed produces no cross-mapping → no identifiability problem →
correctly not a family (the stated scope: identifiability-relevant, not evolutionary). RFPL1/2/3 are expressed
(12/4/47 primaries) but place uniquely → 0 conflict → correct TN. Sub-caveat: a copy receiving only secondary
spillover (array V3: 0 primaries / 116 secondaries) is a confusion *destination*; isolated, it is dropped — the
same cap-fragmentation mode, not an independent error.

### FN class 3 — MIN_READS=3 strictness (NO FN on this substrate)

The only near-threshold pair is the CNN2 retrocopy at $n_{de}=1$, correctly a true negative. No real co-located
family fires below quorum; best-overlap attribution drops no real secondary on the panel.

## FP-robustness (quantified, genome-wide)

An unbiased scanner (4.78 M alignments; 230,354 reads with ≥2 placements at $de\le0.05$) gave 5128 edges / 457
components; each classical confound was adjudicated at exact gene/exon resolution:

- **Heterozygosity** — excluded by construction. Within-locus de-tied *distinct* placements = **0** across ACTB
  (16,674 reads), GAPDH, SEPTIN7, CASP8, CREB1, RABL2A, EEF1A1. Two haplotypes of one gene cannot make a
  *cross-locus* edge.
- **Domain-sharers / segmental-dup spillover** — at exact intervals every cross-family bridge collapses
  (LRPAP1-like↔DOK7-like $=0$) while within-family LRPAP1-likeA↔B $=3242$. A 336-node "mega-component" was an
  artifact of the **scanner's 50 kb bins**, not the definition — it vanishes at de-novo locus resolution.
- **Apparent unrelated-gene survivors** (SEPTIN7↔OCLN, EEF1A1↔KIAA1328) — exonic test shows both-host-exon reads
  $=0$; the reads land on a retrocopy sharing the genomic window. These are **true** retrocopy identifiability
  families mislabeled by host annotation — correct positives, not FPs.
- **rRNA/tRNA/mito** — polyA selection removes rRNA/tRNA; **0** mito↔nuclear de-tied edges ≥3.

**No genuine FP mode found.** Robustness is **conditional on tight per-locus vertices** — the shipped pipeline
uses de-novo transcript-level loci, so this is a guardrail to preserve; the 336-node binning artifact proves
coarse vertices *can* break the best-overlap surrogate.

## Genome-wide reproduction — not overfit to the panel (`bench/family_def_genomewide.py`)

To rule out "it only works on the 17 hand-picked cases," the *exact* de-tie criterion was applied to **every
annotated gene** as a vertex — **34,114 independent vertices** (RefSeq/Gnomon genes, *not* our own de-novo loci, so
no circularity; and *coarser* than the production loci, so this is a conservative upper bound on FP). All 380,369
secondary placements in GGO.bam were scanned and attributed by best overlap.

- **416 families** (2,829 edges, 1,698 genes); **57 % are size-2** and **86 % (357/416) are co-located on ≤2
  chromosomes** — the structure of real recent-paralog tandem arrays / segdup pairs, not random cross-mapping
  (which would not co-locate). Co-location is the genome-wide real-signal proxy (most members are uncharacterized
  `LOC` arrays no orthology DB annotates, so genome-wide Compara enrichment is not computable; the *characterized*
  panel families RABL2/MAGEA are Compara-confirmed separately, above).
- **Δ is flat genome-wide (±7 % band), not a panel fit.** Family count is **388 / 416 / 436** at Δ = 0.003 / 0.005 /
  0.007 — a ±7 % band straddling the operating point — and only reaches 538 at the far decoy threshold Δ = 0.017.
  (The count rises monotonically with Δ, so this is an *insensitivity* band, not a literal local minimum; "valley" is
  reserved for the panel correctness plateau, where 0.005 is correct and the first error appears at 0.007.) The
  Δ=0.005 operating point sits in this flat region genome-wide, independent of the panel.
- **The only systematic FP is the documented coarse-vertex over-merge:** 14 % (59/416) are cross-chromosome bridges
  of mixed genes (e.g. a 60-gene `LOC` component spanning 4 chromosomes; a 54-gene one over 19) — *exactly* the
  best-overlap-surrogate failure the FP-robustness section names. These concentrate in a handful of mega-components.

**DEMONSTRATED — tighter (production) vertices collapse the artifact (not just asserted).** Re-running the *same*
scan over the production **de-novo transcript loci** (101,467 tight vertices via `family_def_genomewide.py denovo`):
cross-chromosome bridges fall **59 → 20** (14 % → 9 %), coherence rises **86 % → 91 %**, and the Δ-valley stays
stable (195 / 212 / 227 families at Δ = 0.003 / 0.005 / 0.007). The *worst* coarse over-merges — the 60-gene `LOC`
component and the unrelated-gene bridge (`AP3B1`/`COL24A1`/`GPHN`/`GRID2`) — **vanish entirely**. The ~20 residual
cross-chrom families at production resolution are *not* the over-merge: each is a **co-located core plus a few
cross-chromosome links** (e.g. the n=25 family's members are all at `NC_073227.2:12.08 Mb`), the signature of
genuine **dispersed paralogs / processed-pseudogene retrocopies** — the "true retrocopy identifiability families
mislabeled by host annotation" the FP-robustness section already counts as *correct positives*, not unrelated-gene
bridges.

So at genome scale the definition is **stable (Δ flat to ±7 %), structurally sensible (size-2 / co-located
dominated), and its over-merge FP is a *measured* vertex-granularity effect that the production loci cut by
two-thirds** — the criterion itself produces no genuine FP mode, and the result is the opposite of a panel-tuned one.

## Precision/recall against a DNA-sequence ground truth (`bench/family_def_dna_pr.py`)

A natural advisor question: *put a precision/recall number on the RNA family definition against the "biological"
(DNA-sequence) definition.* We can — but the honest result is that the two definitions **answer different
questions**, so the raw numbers are not error rates; they are decomposed below. Every number here was independently
re-derived from the raw alignments by an adversarial verification pass (**0 discrepancies**); the per-edge audit
table is `bench/family_def_dna_pr_edges.tsv`.

**DNA ground truth (independent of reads).** All-vs-all alignment (`minimap2 asm20`) of the longest cDNA per gene over
the *same* 34,114 gene vertices. For each unordered pair we keep the best identity and the aligned fraction in each
direction (`cov_a`, `cov_b`).
- **LOOSE** paralog edge: $\mathrm{id}\ge 0.90 \wedge \max(cov_a,cov_b)\ge 0.30$ (the pinned `config.sh`
  MIN_IDENTITY/MIN_COV_FRAC; one-directional, because real paralogs' UTRs diverge — RABL2A/B align over only 34–38 %
  of the mRNA). → **17,410 edges / 1,460 families**.
- **WHOLE-GENE** paralog edge (reciprocal): $\mathrm{id}\ge 0.90 \wedge \min(cov_a,cov_b)\ge 0.50$ — meant to exclude
  domain-sharers. → **8,698 edges / 895 families** (a strict subset of LOOSE).

**Edge-level precision/recall** of the RNA de-tie graph (2,829 edges) against this truth:

| DNA truth | TP | FP | FN | precision | recall |
|---|---|---|---|---|---|
| LOOSE, all genes | 1822 | 1007 | 15588 | 0.644 | 0.105 |
| WHOLE-GENE, all genes | 1288 | 1541 | 7410 | 0.455 | 0.148 |
| LOOSE, expressed only | 1169 | 1660 | 3088 | 0.413 | 0.275 |
| WHOLE-GENE, expressed only | 822 | 2007 | 830 | 0.291 | 0.498 |

Both axes are dominated by **definitional difference**, not error:

**Recall — 80 % of "missed" DNA paralog pairs are transcriptionally SILENT.** Decomposing the 15,588 LOOSE-DNA edges
the RNA graph does not carry: **12,500 (80.2 %)** have ≥1 *unexpressed* copy (no RNA evidence — out of scope by
definition); **2,141** have reads that place uniquely (RNA-distinguishable copies); **524** are *resolvable*
(reads cross-map but the divergence gap exceeds Δ, the APOBEC3 principle); **423** are *sub-quorum* (genuinely tied
but at only 1–2 reads, below MIN_READS=3). The honest single recall figure is therefore **recall | cross-mapping
universe = 0.658** (of expressed paralog pairs whose reads actually co-map, the fraction linked). The often-quoted
0.812 ("tied universe") is reached only by *additionally* defining the 524 resolvable pairs as out-of-scope — a
modelling choice that must be argued, not folded in silently. **Note (anti-circularity):** "genuine miss = 0" (no
expressed pair with ≥3 tied reads is unlinked) is *true but by construction* — an RNA edge **is** a pair with ≥3
codivergent reads — so it is **not** presented as an empirical validation; the operative recall loss is the
sub-quorum (423) and resolvable (524) buckets.

**Precision — the DNA bar is the threshold-dependent one, and it under-detects.** Of the 1,007 RNA-only edges,
**212 (21 %)** have real cDNA homology $\mathrm{id}\ge 0.80$ but fall *just* under the arbitrary 0.90/0.30 bar — e.g.
`TBC1D1~LOC134756953` (id 0.8975, covA 0.94, **293** tied reads) and `RABL2B~LOC134756389` (id 0.85, 37 reads): real
paralogs the bar rejects. Crediting these gives **effective precision = (1822+212)/2829 = 0.719**. Of the remaining
**757 zero-homology** RNA-only edges, only **131** are tandem (<1 Mb) — defensible local paralogs whose longest-cDNA
representative failed to align (e.g. `GSTM2~LOC115933235`, a 118 kb tandem in the real GSTM cluster); the 245
same-chrom-far and 381 cross-chrom are unvalidated. The **highest-evidence** cross-chrom edges are *not* spurious
noise: `OCLN~SEPTIN7` (**3,369** tied reads) and `BCAS4~CCDC30` (**962**) are **processed-pseudogene/retrocopy** and
**read-through/chimeric** loci.

**The definitional crux — "family" ⊋ paralogy.** The read-conflict criterion detects *read-indistinguishable locus
pairs*, a **superset** of sequence paralogy: it also fires on retrocopies, processed pseudogenes, and read-through
fusions. A cDNA-vs-cDNA DNA truth structurally cannot represent these either, so they are a **mutual blind spot**, not
an RNA error.

**Orthogonal, not nested.** On the expressed universe the RNA graph and the whole-gene DNA graph **cross**: 830
whole-gene paralog pairs the RNA graph misses (silent/resolvable) vs 1,130 RNA edges that are not whole-gene
paralogs. So the read-conflict graph is **not a proxy for a static coverage threshold** — it conditions on actual
transcribed, co-mapping evidence. And no single coverage cut separates the classes cleanly: among LOOSE-expressed
edges, RNA-confirmed pairs have median reciprocal coverage **0.75** vs **0.27** for DNA-only — shifted but
overlapping. The reciprocal-coverage bar **drops RABL2A~RABL2B** (real, recip 0.34) yet **keeps CREB1~METTL21A**
(domain-sharer, recip 0.54).

**The one clean win — domain-sharers.** All five panel domain-sharers (`CREB1~METTL21A`, `GCA~KCNH7`, `CASP8~FLACC1`,
`ASDURF~ASNSD1`, `GPR39~LYPD1`) are genomically nested/adjacent genes sharing **one** homologous domain at mapq 0;
**all five pass the LOOSE DNA bar** (and `CREB1~METTL21A` passes even the reciprocal whole-gene bar), yet the RNA
best-overlap attribution **correctly excludes all five** (`in_rna=0`). A DNA homology bar over-calls them as families;
the read-conflict criterion does not.

**Conclusion.** A precision/recall number against DNA exists and reproduces exactly, but it measures the *overlap of
two different definitions*. Read-confusability is an **orthogonal, expression-conditioned evidence axis**,
non-redundant with — and non-nested in — cDNA-homology thresholds: a clean advantage on domain-sharers, a
complementary blind spot (and a unique strength) on retrocopies/pseudogenes, and a recall ceiling set entirely by
transcriptional silence and the read-evidence quorum, not by mis-calls.

### Sharpening levers tested — only repeat-masking helps, and only modestly

Five natural ways to "sharpen" the definition were each implemented and scored against the DNA homology truth with
the same `good:bad` metric (real-paralog edges lost per junk edge lost; want < 1). **Four fail for one common reason;
only repeat-masking is net-positive** — the survey is the strongest evidence the definition is not leaving easy
precision on the table.

| lever | level | result | mechanism |
|---|---|---|---|
| better clustering object (modularity / biconnected / cut-edge / k-truss / k-core / cliques) | graph | no clean gain | over-merge is a sparse multi-edge web (vertex problem); 57 % of families are size-2 (no density) |
| high query coverage (both placements) | read | net-harmful (good:bad 1.28) | real conflicts on partially-aligned reads |
| ≥1 intron (both placements) | read | net-harmful (1.83) | real conflicts on non-spliced reads |
| multi-exon locus | de-novo vertex | **inert** (0 edges) | assembly emits ~100 % multi-exon loci |
| multi-exon gene | gene vertex | net-harmful (2.53) | many real families *are* single-exon |
| junction concordance (same exon architecture at both loci) | read structure | real signal, **steep cost** (good:bad 1.66) | separation 0.56 vs 0.41 too weak; mislabeled paralogs splice concordantly |
| **repeat-masking** (drop majority-soft-masked placements) | read/genome | **net-POSITIVE (good:bad 0.76–0.89)** | bridges are TE/satellite-driven, not genic |

**(1) A different graph object does not help** (`bench/family_def_graph_operators.py`). Replacing the connected
component with biconnected components, modularity communities, or cut-edge pruning was tested with a
*confound-free* metric — the DNA-paralog rate of CUT vs KEPT edges (granularity mechanically inflates within-family
purity, so raw purity is not trustworthy). Modularity's higher purity (0.55 vs CC 0.28) is granularity-gaming: it
cuts edges that are **53 %** real paralogs (ratio 1.27). Only cut-edge pruning cuts the *right* edges (30 % vs 65 %
paralog, ratio 2.16) but reaches just **40** of them — the over-merge mega-components are **2-edge-connected** (a read
confused among many loci via one shared element makes a *dense web* of sparse edges, not a single cut-edge), so no
clustering operator removes them.

**(2) Read-level coverage / intron filters are net-harmful** (`bench/family_def_read_filters.py`). Requiring each
conflicting read to be high-coverage and spliced at *both* placements removes genuine paralog edges faster than
artifacts (good:bad lost 1.28–1.83) for ≤ 0.024 precision: real cross-copy conflicts are often carried by reads that
do not cross a junction at both copies, while the marquee spurious bridge `OCLN~SEPTIN7` (3,369 reads) is *kept*
because its reads *are* full-length and spliced.

**(3) The vertex-level intron filter is inert or net-harmful** (`bench/family_def_denovo_intron.py`). The "retrocopy =
intronless locus" intuition is correct biologically but does not separate the classes: over the production de-novo
loci it removes **0** edges (the assembly emits 93,372/93,373 multi-exon loci; all 632 cross-chrom de-novo bridge
edges are multi-exon↔multi-exon), and over annotated-gene vertices it is net-harmful (good:bad 2.53) because **many
real paralog families are single-exon** (histones, olfactory receptors, interferons) — single-exon-ness is not
specific to the contamination.

**(4) Repeat-masking is the one net-positive lever** (`bench/family_def_repeat_mask.py`). The NCBI gorilla assembly is
already soft-masked (lowercase = TE / satellite / low-complexity; **segmental duplications are NOT masked**, so real
segdup paralogs are preserved). A conflicting-read placement is "repeat-driven" if its aligned exonic (M/=/X) blocks
are majority soft-masked; dropping those before the graph build is the **first** lever with `good:bad < 1`
(0.89 at threshold 0.5, 0.76 at 0.7), raising precision 0.644 → 0.684. The mechanism is confirmed directly: no-homology
bridge edges have mean placement repeat-fraction **0.34** vs **0.21** for real DNA paralogs — the spurious
cross-mapping is genuinely driven by shared masked repeats, not genic homology. The gain is **modest** (+0.04) because
the distributions overlap (real paralogs also carry ~21 % repeat content), so no threshold separates them cleanly; but
it is principled, mechanism-specific, and segdup-safe, and is the only structural predicate worth adopting (as a
preprocessing step, threshold ≈ 0.5).

**(5) Junction concordance — the spliced-structure lever — is a real but too-weak signal**
(`bench/family_def_junction_concordance.py`). Requiring a conflicting read to recover the *same exon architecture*
(matching exon-length signature) at both loci produces the **largest raw precision gain** (0.644 → 0.709), and the
mechanism is real: DNA paralogs splice concordantly at rate **0.56** vs **0.41** for no-homology bridges. But the
separation is too weak for a hard filter — good:bad **1.66** (~35 % recall cost) — because real paralogs read
discordantly via alternative splicing / partially-aligned diverged copies, while the bridges that are *mislabeled real
paralogs* splice perfectly concordantly and survive. So it is a legitimate soft signal, not an adoptable filter.

**Why the structural filters fail — and why repeat-masking is the exception.** This exposes the unifying reason: the
spurious contamination is **two distinct populations**. (a) *Genuine repeats / retrocopies* — caught by repeat-masking
(external genome annotation) and partly by junction-*discordance*. (b) *Mislabeled real paralogs* (coarse-vertex
attribution) — read-wise *indistinguishable from real paralogy* (full-length, spliced, multi-exon, junction-concordant,
distributed homology) because they **are** paralogy, only attributed to the wrong gene; no substring / Jaccard /
coverage / junction feature can separate them, and only **vertex resolution** can. So any predicate built from read or
locus *structure* removes real paralogs at least as fast as artifacts.
Repeat-masking is the lone exception because it adds *external* information (the genome's repeat annotation) that is
correlated with the contamination mechanism rather than with read structure — yet even it only partly separates the
classes. The dominant clean lever remains **vertex granularity** (de-novo loci, bridges 59→20), already in place; at
de-novo resolution the residual 20 cross-chrom bridges are between multi-exon loci, i.e. largely *correct* dispersed
paralogs. The precision residual is therefore architectural (vertex resolution) plus a modest repeat-mask gain, not a
missing read predicate.

### VG structural validation — the sequence-graph lever that rejects repeat-bridges (`bench/family_def_vg_validation.py`)

The read-conflict graph is a *discovery* layer (which loci are confusable); it discards the **sequence structure**. The
variation graph is the *validation* layer: a real family's copies are **long parallel paths** sharing a backbone
(PSVs = inter-copy variation), whereas a repeat-bridge shares only a short bubble. This is the one lever that reaches the
contamination the read-features could not — and it is built natively from sequence, so it sidesteps the gene-label
(vertex) problem.

**Isoform aggregation is a definitional requirement, not a detail.** A *copy* is a gene; its reads are *isoforms*. If the
graph is built from raw isoforms, alternative splicing within one copy masquerades as extra copies. So each copy must be
the **exon-union of all its isoform reads** ("the full gene minus introns"; retained introns appear where reads span
them), and the graph must separate **inter-copy variation (PSVs → copies)** from **intra-copy variation (splicing →
isoforms)**. This is also a *sensitivity* gain: the longest-single-isoform cDNA truth under-detects paralogy (the 212
sub-bar / 757 no-homology edges); the all-isoform exon-union recovers it.

Proof of concept on 8 candidate edges (all-isoform exon-union per copy → align copies → backbone):
- **True repeat-bridges are cleanly rejected** — `OCLN~SEPTIN7` (max backbone coverage **0.05**) and `BCAS4~CCDC30`
  (**0.09**): the marquee bridges that *no* read-feature could remove now show **no backbone in either direction**.
- **Real families confirmed as parallel paths** — `RABL2A~RABL2B` 0.94, `CCDC196~LOC` 0.97, `RABL2B~LOC134756389` 0.81
  (the last recovered from sub-bar by isoform aggregation).
- **Isoform aggregation works and recovers copies** — 100s–1000s of distinct isoforms collapse to **one** exon-union per
  copy (no pollution); fragmentary annotations are completed from reads (`LOC134756953` 3.5 kb → 11 kb); and
  `GSTM2~LOC115933235`, labelled *no-homology* by the longest-cDNA truth, is **recovered** (coverage 0.98) as a real
  co-located GSTM-cluster paralog the truth missed.
- **Asymmetric cases point back to vertex granularity** — `AK6`, `GSTM2`, `TBC1D1` have coarse annotation spans
  (250 kb – 1.18 Mb) that inflate the larger copy's exon-union and break reciprocal coverage; de-novo loci give clean
  exon-unions. So VG validation **composes** with vertex resolution rather than replacing it.

So the cleanest definition has three composed layers: **read-conflict graph** (discovery) over **de-novo loci** (vertex
resolution) → **VG validation** with **all-isoform exon-union** copies (reject no-backbone repeat-bridges, count copies
by PSV haplotype not isoform). VG validation is the structural lever for the repeat-bridge population that the read/locus
predicates could not separate; the `--vg` family-graph machinery already builds exactly this object.

**Genome-wide, OOM-safe, validated** (`bench/family_def_vg_reinforce.py`). The three layers were run genome-wide:
de-novo loci → read-coupled candidate families ($\sim_R$ components) → **$\sim_B$-connected core** (a member is kept iff
it shares an exon-union backbone with a family-mate, so isolated bridge members fall out and size-2 families survive iff
the pair shares a backbone; a guard rejects *only* well-modelled backbone-less loci, never on inability to validate).
Memory stayed > 15 GB free throughout (copy models built only for the ~1,300 family loci, reads capped, streamed).
Result: **212 candidate → 196 validated families**, cross-chromosome bridges (≥3-chromosome components) **20 → 12**,
**14 well-modelled backbone-less members rejected**, size-2 families **80 → 73**. The genome-wide refinement is *modest*
(de-novo candidates are already mostly coherent, so $\sim_B$ mainly clears the residual cross-chromosome bridges); the
*decisive* evidence for $\sim_B$ is the per-candidate separation (repeat/retro bridges ≤ 0.1 vs real copies ≥ 0.8). The
formal object is the $\sim_B$-connected component within a read-coupled group (components of $\sim_B \cap R^\*$), **not**
the direct edge-intersection $\sim_R\cap\sim_B$ (which over-rejects, fragmenting real arrays at any weak-backbone edge),
and **not** the gene-vertex pairwise variant (`family_def_vg_pipeline.py`, good:bad 1.33, coarse annotation inflates the
exon-union). Both ingredients are load-bearing: de-novo vertices remove annotation inflation; the within-family
backbone connectivity isolates bridges while preserving pairwise (size-2) families. The closed-form statement is in
`family_definition_note.md`.

## Residual hardcoded parameters (disclosed, not independently swept)

The two de-tie constants $\Delta,\mathrm{DE_{max}}$ are now the **only** free parameters (the former 200 bp
co-location guard is gone — record-attributed placement makes it structurally unnecessary). Both are
**error-model-derived** (P3): $\Delta$ as the single-read divergence resolution at HiFi error/length,
$\mathrm{DE_{max}}$ as a loose copy-vs-distinct-gene divergence ceiling. $\mathrm{MIN\_READS}=3$ is a
minimum-evidence quorum (a standard support floor), not a tuned constant. The honest residual: the empirical
$\Delta$-valley was confirmed against only the panel's 3 resolvable decoys (a genome-wide decoy sweep would tighten
the margin), and $\mathrm{DE_{max}}$ was set conservatively rather than swept for its own valley. Truth labels are
human-assigned priors; the *positives* (RABL2, MAGEA) are now independently corroborated as recent within-species
paralogs by Ensembl Compara (above), and the cross-mapping evidence was verified consistent with every label —
full external orthology ground truth across all 17 candidates (incl. the de-novo LOC copies) remains open.


---

## family_definition_note

# A self-contained RNA-level definition of multi-copy gene families

**Two relations on expressed loci — read-confusability and backbone-homology — and the family is where both hold.**

This note states the definition in closed form. The empirical record (panel ledger, genome-wide run, survey of
rejected alternatives) is in `DEFINITIONS_FORMAL.md`; this note is the formal object and is self-contained
given one named input (the locus set, §1).

---

## 0. The object, in one sentence

> A **multi-copy gene family** is a maximal set of expressed loci that is both **read-coupled** (its loci are mutually
> confusable to the aligner) and **backbone-connected** (its loci share a common gene backbone).

The two italicised conditions are two relations defined below; a family is a connected component of the second
relation restricted to a coupled group of the first. Neither relation alone is the right object — that is the point
(§4).

---

## 1. Primitives

Fix a set of spliced long reads aligned to a genome.

- **Vertices $V$.** The **de-novo loci** emitted by the upstream read-coherence assembly: each vertex is an expressed
  transcription unit — a maximal set of reads sharing splice structure at a genomic position. $V$ is therefore
  *read-derived and label-free* (no annotated gene is consulted), but it is an output of that named procedure, not a
  closed-form predicate in this note; the structural claims below that rely on tight per-locus vertices (P6) are
  stated as conditional on it.
- **Copy model $S(v)$.** The **all-isoform exon-union** of a locus $v$: every exonic base any read at $v$ aligns
  through (a retained intron appears where reads span it). A *copy* is a locus, and its isoforms **aggregate** into one
  $S(v)$ — so alternative splicing within a copy is one model, not several copies.
- **Divergence.** Each alignment of a read $r$ to a locus carries gap-compressed per-base divergence $d(r,v)$
  (minimap2 `de`). Each alignment **record** is attributed to its single best-overlap locus — ties broken toward the
  *smallest containing* locus, which makes the attribution single-valued and routes a record enclosed by a larger
  nesting locus to the most specific one — so a read's placement set has at most one entry per record.

---

## 2. The two relations

**Read-confusability $\sim_R$** — the aligner cannot tell the loci apart. With tie tolerance $\Delta$, ceiling
$\mathrm{DE_{max}}$, quorum $k$:
$$
i \sim_R j \iff \big|\{\, r : d(r,i),d(r,j)\le \mathrm{DE_{max}},\ |d(r,i)-d(r,j)|\le \Delta \,\}\big| \ge k .
$$
Because attribution is per-record, $r$ contributes only when it is a genuine *alternative* placement (two distinct
records), never when one locus merely nests in another. Let $R^\* $ be the **read-coupling** relation: $i\,R^\*\,j$
iff $i,j$ lie in the same connected component of $\sim_R$.

**Backbone-homology $\sim_B$** — the loci are the same gene over a real fraction of *both* copies:
$$
i \sim_B j \iff \min\!\big(\mathrm{cov}_i,\mathrm{cov}_j\big)\ge \tau ,
$$
with $\mathrm{cov}_i$ the aligned fraction of $S(i)$ (minimap2 `asm20`). Reciprocal coverage is what distinguishes a
shared backbone (both copies largely homologous) from a shared element (one copy aligns only over a short insert).
A separate identity floor is unnecessary: the aligned region is intrinsically high-identity — every pair clearing
$\tau$ has alignment identity $\ge 0.96$ (measured, $n=393$) — so identity never binds independently and is **not a
free parameter** (§6).

> **Granularity (resolves an apparent code/note discrepancy).** $\sim_B$ is a **copy-level** relation: one alignment
> of the whole exon-union models $S(i),S(j)$, scored by reciprocal coverage. The production family-graph builder
> *additionally* applies a **contiguous-core** test at the **exon level** (per shared exon node). The two are
> complementary, not competing — copy-level coverage certifies "these loci are copies"; exon-level contiguity governs
> which individual exons fuse into one shared graph node. Contiguous-core is **not** substitutable at the copy level:
> on whole exon-union models it fragments at exon boundaries and rejects real copies (measured: it drops ~20 % of clean
> two-copy families that are 99 % identical over 80 % of their length).

---

## 3. The definition

> **Definition.** A **multi-copy gene family** is a connected component, of size $\ge 2$, of the graph
> $G = \big(V,\ \sim_B \cap R^\*\big)$ — i.e. a $\sim_B$-connected component *within* a read-coupled group. Its
> **copies** are the locus models deduplicated by exon-union identity.

Equivalently and operationally: take the connected components of $\sim_R$ (read-coupled candidate families); within
each, keep the connected components of $\sim_B$ (the backbone-coherent cores). A locus that is read-coupled to a group
but shares a backbone with no member of it is in no family.

*Why these coincide (and why $R^\*$, not $\sim_R$).* $R^\*$ is an equivalence relation, so every $\sim_B\cap R^\*$ edge
lies within a single $R^\*$-class and connectivity cannot cross classes; within a class $\sim_B\cap R^\*=\sim_B$. Hence
$\mathrm{components}(V,\sim_B\cap R^\*)=\bigcup_\text{classes}\mathrm{components}(\sim_B|_\text{class})$. The closure is
load-bearing: the analogous statement with the *direct* relation $\sim_R$ is false — in a $>\!2$-copy array whose
distal copies never co-place within one read's placement set, $\sim_R$ would fragment a coupled group that $R^\*$
keeps whole so $\sim_B$ can re-knit it.

---

## 4. Why both relations — and why this is the right object

Each relation alone is a known, deficient definition; the contribution is that the family is their conjunction.

- **$\sim_R$ alone is the read-conflict graph.** Necessary — if two loci are not confusable, assigning a read between
  them is trivial and they are not one unit — but it *over-couples*: a shared transposable element or a processed
  retrocopy makes the reads of two unrelated loci cross-map, fusing them.
- **$\sim_B$ alone is sequence homology** (the "biological" definition). Necessary — copies must share gene structure
  — but it is annotation-style: it admits unexpressed copies (out of RNA scope) and says nothing about whether the
  reads are actually confusable. (It does *not*, on its own, exclude domain-sharers; those are removed by $\sim_R$, P1.)

A family is the **conjunction**: loci that are at once *coupled* ($R^\*$ — the copy-assignment problem is non-trivial)
and *the same gene* ($\sim_B$ — it is well-posed). $\sim_R$ supplies the coupling, $\sim_B$ supplies the identity; the
elegance is intrinsic — the "validation" is not a separate mechanism but **connectivity of $\sim_B$ inside a coupled
group**, so a backbone-less bridge locus is simply not in any component, and a genuine pair is a single
$\sim_B$ edge inside its coupled group (no density required).

---

## 5. Properties

**Structural — hold by construction** (P1 additionally requires tight vertices, P6).

- **P1 (domain-sharers excluded) — by construction *given* tight vertices (P6).** Single-valued per-record best-overlap
  gives nested/adjacent single-domain genes $0$ conflicting reads ⟹ not $\sim_R$ ⟹ no family. The *implication* is
  structural, but its premise — best-overlap is single-valued and faithful — holds only when vertices are tight
  transcription units; under coarse vertices best-overlap can mis-attribute and manufacture domain-sharer edges, so
  **P1 inherits P6's conditionality** (it is structural-given-P6, not unconditional). *Panel: 5/5 domain-sharers
  excluded, though they pass DNA homology.*
- **P3 (no isoform pollution).** Copies are exon-unions, so within-copy alternative splicing is an intra-copy bubble,
  never an extra copy. *100s–1000s of isoforms per locus collapse to one $S(v)$.*
- **P4 (pairs admitted).** A genuine pair is a single $\sim_B$ edge inside its coupled group — no triangle or density
  is required (which is why clique/$k$-truss objects, needing density, are wrong: 57 % of families are pairs).

**Empirical / conditional — hold on the measured substrate.**

- **P2 (repeat/retro bridges rejected).** Two loci whose reads cross-map through a shared element but which are *not
  copies of each other* fail $\sim_B$ (one side aligns only over the short shared insert) ⟹ rejected. *Decisive
  measurement: `OCLN~SEPTIN7` (3,369 confusing reads) has backbone coverage $0.05$, `BCAS4~CCDC30` $0.09$, vs real
  copies `RABL2A~RABL2B` $0.94$, `CCDC196~LOC` $0.97$.* (A retrocopy that genuinely **is** a copy of its parent shares
  a backbone and is correctly kept; $\sim_B$ separates "retro-mediated bridge between non-copies" from "retrocopy of
  parent" exactly by the backbone.)
- **P5 (resolvable paralogs excluded — by quorum).** A diverged copy is excluded iff *fewer than $k$* reads tie at the
  `de` floor; this is a quorum statement, not a per-read one — inclusion can be carried by a quorum even when the
  median per-read gap exceeds $\Delta$ (e.g. RABL2, median gap $0.006$), and `EEF1A1`'s retro-pseudogene is excluded
  because $0$ reads tie. The definition is thus *recent-paralogy ∩ read-confusability*, a strict subset of paralogy.
- **P6 (well-posed vertices).** Reciprocal coverage in $\sim_B$ is meaningful only when $S(v)$ is a tight
  transcription unit; this is conditional on the de-novo locus pipeline (§1) producing per-locus, not coarse, spans —
  the one external assumption.

---

## 6. Parameters

| symbol | role | value | basis |
|---|---|---|---|
| $\Delta$ | $\sim_R$ tie tolerance | 0.005 | single-read divergence resolution at HiFi error: per-read SE $\sqrt{\epsilon/L}\approx0.0009$, tie statistic $\sqrt{2\epsilon/L}\approx0.0013$, $\sim 4\sigma$ |
| $\mathrm{DE_{max}}$ | $\sim_R$ divergence ceiling | 0.05 | loose copy-vs-distinct-gene ceiling |
| $k$ | $\sim_R$ quorum | 3 | **quorum classifier, load-bearing** (not an inert floor): it admits true positives the per-read tie test alone misses — RABL2's median per-read $|\Delta d|=0.0061>\Delta$, so only the quorum of $k$ tied reads carries it (§5, P5) |
| $\tau$ | $\sim_B$ reciprocal coverage | 0.30 | the one data-calibrated threshold: set in the empty gap between repeat-bridges (one-sided coverage $\le 0.1$) and validated copies ($\ge 0.31$ genome-wide). $\tau=0.30$ sits at the gap's lower-copy edge — permissive (it admits partial-coverage copies); it could be centred lower in the gap at a small recall cost |
| $\mathrm{GUARD}$ | $\sim_B$ rejection guard | 20 reads | a backbone-isolated locus is rejected as a bridge only with $\ge\mathrm{GUARD}$ reads (enough to model confidently); below it the locus is held out as unmodelled, not rejected |

An identity floor $\iota$ is **omitted as measured-inert**: every pair clearing $\tau$ already has alignment identity
$\ge 0.96$ ($n=393$, min $0.963$), so identity never binds independently. $\Delta$ is a measurement constant; $k$ is a
load-bearing quorum classifier; $\tau$ is the single empirical threshold, placed in a measured gap. $\mathrm{GUARD}$
gates *rejections* only, so §3/§4's "a backbone-less locus is in no family" is exact for well-expressed loci and a
hold-out (not a rejection) for sparse ones.

---

## 7. Evidence (summary; full record and the rejected-alternatives survey in `DEFINITIONS_FORMAL.md`)

- **Panel.** 17 hand-labelled IsoSeq candidates: **TP = 7, TN = 10, FP = 0, FN = 0**.
- **$\sim_B$ is the decisive lever.** Backbone coverage separates repeat/retro bridges (one-sided coverage $\le 0.1$)
  from validated copies — the clean panel cases sit high (RABL2A~RABL2B $0.94$), and genome-wide the lowest validated
  copies reach $\ge 0.31$, so the gap is $[0.1,\,0.31]$ with $\tau=0.30$ at its upper edge. This is the population no
  read-level predicate (coverage, intron, junction concordance) could separate, because that contamination is
  read-indistinguishable from real paralogy.
- **Robust to multimapper sampling (new $-N\,50\,-p\,0.1$ BAM).** The old alignment used minimap2's default secondary
  cap, undersampling $\sim_R$. Re-aligning with $-N\,50\,-p\,0.1$ surfaces $12$–$65\times$ more cross-mapping; $\sim_R$
  then *recovers* hidden dispersed paralogs (per-chrom $\sim_B$-validated copies grow, e.g. chrY $21\to35$) while
  $\sim_B$ prunes every extra bridge (e.g. $10\to59$ on one chromosome). Family identity tracks copy-model homology,
  not read counts (`family_def_newbam_validation.md`).
- **Genome-wide, OOM-safe** (memory > 15 GB free throughout; copy models built only for the ~1,300 family loci,
  reads capped, streamed). Over de-novo loci: **212 candidate → 196 validated families**; cross-chromosome bridges
  (components spanning $\ge 3$ chromosomes) **20 → 12**; **14** well-modelled backbone-less members rejected as
  confident bridges; size-2 families 80 → 73. The genome-wide refinement is *modest* — most de-novo candidate
  families are already coherent, so $\sim_B$ mainly removes the residual cross-chromosome bridges; the *decisive*
  evidence for $\sim_B$ is the per-candidate separation above, not a large genome-wide rejection count.
- **Honest cost.** $\sim_B$ is a precision filter with a small recall cost: among the 14 rejected members the
  rejected edges are $9$ no-homology vs $12$ DNA-paralog/sub-bar, so it is conservative — calibrated to discard
  confident bridges, not to maximise rejection.

---

## 8. Relation to the copy-assignment problem

A family is, by construction, the unit on which read-to-copy assignment is at once necessary ($R^\*$ couples its reads)
and well-defined ($\sim_B$ guarantees a shared backbone, hence comparable copies). The definition hands the
copy-assignment / identifiability machinery exactly its input — mutually confusable, genuinely homologous copies — and
nothing else.

---

## 9. What it does not claim

- It defines families that are **expressed and read-confusable**; unexpressed or read-resolvable paralogs are out of
  scope **by design** (they pose no copy-assignment problem), not missed.
- "Family" here is **read-indistinguishable homologous copies** — a strict subset of sequence paralogy; genuine
  repeats and bridges between non-copies are rejected by $\sim_B$, and read-resolvable diverged paralogs by $\sim_R$.
- **Copy number is deferred.** Counting copies requires PSV-haplotype resolution across the shared backbone and is not
  computed here; the backbone supplies its substrate, and exonically identical copies collapse to one backbone (the
  identifiability floor).
- The single external dependency is the de-novo locus set (§1); the closed-form claims are conditional on it producing
  tight transcription-unit vertices.


---

## family_graph_definition

# Graph-based multi-copy family definition (per-family POA variation graph) — an option, with a formal definition

Upgrade #1 as an **option** (not a replacement for the pairwise definition): build one POA variation
graph per candidate family (all N member transcripts at once) and read the family criterion off the
graph. It yields a clean, formally-statable definition that generalizes the validated pairwise one.

## Formal definition
Members `S = {s₁,…,s_N}`. Let `G` be their partial-order (POA) alignment graph with columns
`c₁,…,c_M` in topological order. For a column `c`, let `supp(c) = #{ i : sᵢ is non-gap at c }`.
The **conserved core** is the longest contiguous run `R` of columns with
`supp(c) ≥ max(2, ⌈N/2⌉)` for every `c ∈ R` (a majority spine). Then

> `S` is a **multi-copy family at level T**  ⟺  `|R| ≥ T · median_i |sᵢ|`,
> with **family score** `= |R| / median_i |sᵢ| ∈ [0,1]`.

**Reduces to the validated pairwise core at N=2:** `⌈2/2⌉ → 2`, so a core column requires *both*
copies non-gap, and `R` is exactly the longest co-aligned block — the pairwise contiguous-core we
validated (RABL2/DAZ score = their pairwise core, 0.16).

## Why it's cleaner than pairwise + connected-components
- **One statistic from one graph** — no O(N²) pairwise comparisons, no transitive-closure family
  building. Membership is "does this copy thread the shared core?", a graph property.
- **It exposes over-merge.** Pairwise connected-components chain distinct subfamilies through domain
  hubs (mega-"families" of 250). The graph score on such a chain is LOW (no column shared by a majority
  of members) → **7 of the 14 N≥25 components score < T = NOT single families** (the other 7 are genuine
  large families with a real shared core). Pairwise CC could not make this distinction; the graph does.

## Validation (labeled set)
| class | graph score | vs T=0.045 |
|---|---|---|
| RABL2 (2), DAZ (2) | 0.162 | ✓ family |
| RFPL (4) | 0.081 | ✓ family |
| APOBEC3 (6) | 0.062 | ✓ family |
| domain-sharers (CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, …) | ≤ 0.029 | ✗ rejected |

**Clean separation: domain-sharer max 0.029 < T=0.045 < real-family min 0.062.** Perfect on the
labeled set. Genome-wide: 1,333 candidate families scored; 1,211 ≥ T; the 122 < T are domain-hub
chains / weak components the graph correctly down-weights.

## Honest caveats
- **Graded, not binary** — recent near-identical copies score high (RABL2 0.16); ancient/length-variable
  families score moderate (APOBEC3 0.062, RFPL 0.081). The score is a family-*tightness* measure, so the
  margin between diverged-real-families (≥0.062) and over-merged-chains (median 0.04) is real but thin —
  clean for recent families, tighter for ancient ones (consistent with the RNA-only tier's reach limit;
  ancient families still benefit from the protein/DNA tier).
- **T is metric-specific** (0.045 here vs 0.13 for the raw pairwise core — different normalization).
  Both are data-placed in the labeled gap.
- **Bounded inputs:** members capped at 30 (shortest) and sequences > 15 kb skipped (POA memory);
  noted, affects only the largest families.

## Verdict
The per-family POA variation graph gives a **clean, formal, N-way family definition** that reduces to
the validated pairwise core, separates real families from domain-sharers, and — uniquely — exposes the
transitive-closure over-merge. Worth keeping as the family-definition **option** for well-annotated,
multi-member families.

## Reproduce
- `MINIFORGE python bench/family_graph_define.py` (pyabpoa POA-MSA per family → core score)
- `python3 bench/family_graph_fig.py`


---

## poa_family_definition

# POA contiguous-core coverage family criterion

> *Naming note:* the separating axis is **contiguous-core coverage = largest ungapped block / shorter gene**. Because it divides by the SHORTER gene it equals the **maximum** of the two per-gene block ratios, so it is NOT a reciprocal coverage (a reciprocal = the *min* over both genes, i.e. divide by the LONGER gene). The word *reciprocal* applies only to the all-column literal metric `min(aligned/len_a, aligned/len_b)`. The two metrics use opposite denominators; only the all-column one is reciprocal.

> **POA contiguous-core coverage (largest ungapped block / shorter gene) >= 0.13 accepts 5/5 confirmed tandem dups and rejects 7/7 domain-sharers.**

## Definition under test

Two loci are in a multi-copy *transcript* family iff, aligned, the homologous region covers >= T of BOTH transcripts:

```
reciprocal_coverage(a,b) = min( aligned_len/len(a), aligned_len/len(b) )
```

i.e. the two genes must CO-ALIGN over a MAJORITY of *each* gene, not just one shared domain/exon. Both quantities (reciprocal aligned coverage, identity over aligned columns) are read off the alignment and are POA-derived. No DNA/protein domain annotation, no BLAST, no k-mer/minimizer is used as the criterion (minimizer-Jaccard is recomputed only to CONTRAST).

### POA-faithfulness note

POA generalises pairwise alignment to >2 copies by embedding each new sequence into a DAG of the previous ones; for a PAIR there is no DAG to extend, so the two-sequence POA instance IS a single pairwise global alignment. We use `Bio.Align.PairwiseAligner` in **global** mode (Needleman-Wunsch) with match=+2, mismatch=-1, gap-open=-5, gap-extend=-1 on the gene-representative sequences.

## Two POA coverage axes (the key result)

We measured coverage two ways, BOTH purely alignment-column-derived:

1. **all-column** reciprocal coverage -- the LITERAL definition above, summing every ungapped aligned column; and
2. **contiguous-core** coverage -- the largest *single* ungapped co-aligned block divided by the **shorter** gene (`biggest_block / min(len)`). This is NOT a reciprocal coverage: dividing by the shorter gene makes it the *maximum* of the two per-gene block ratios (a reciprocal would divide by the longer gene = the *min*). It is the metric that separates; the name 'reciprocal' belongs only to axis 1.

The literal all-column metric does **NOT** cleanly separate the classes: confirmed range [0.471, 0.842], domain-sharer range [0.268, 0.843] -- they OVERLAP. The reason is a global-alignment artifact: when the Needleman-Wunsch aligner is handed two NON-homologous transcripts it still pays a little gap cost to harvest scattered chance matches, fragmenting the alignment into hundreds of tiny blocks (domain-sharers here split into 89-579 blocks) whose aligned columns *sum* to a high reciprocal coverage even though the alignment is non-homologous filler (note the low all-column identity of those pairs in the per-pair table).

The **contiguous-core** metric removes that filler and is the POA-only fix. It is **SEPARABLE**: confirmed range [0.174, 0.608], domain-sharer range [0.008, 0.055], control range [0.007, 0.083]. A real copy has ONE large homologous block; a domain-sharer's largest single block covers only a few percent of the shorter gene -- indistinguishable from random cross-family controls.

## Result: the separating threshold

The principled, POA-only threshold (on the contiguous-core axis) is **T = 0.13** (midpoint between the lowest confirmed and the highest domain-sharer/control coverage -- they are strictly separated).

| category | n | accept at T (cov>=T) | reject at T (cov<T) |
|---|---|---|---|
| confirmed | 5 | **5** | 0 |
| domain-sharer | 7 | 0 | **7** |
| control (cross-family) | 60 | 0 | 60 |

Is the distribution bimodal/separable? **YES -- clean bimodal split (~3.2x vs domain-sharers, ~2.1x vs the nearest control).** (contiguous-core axis: confirmed >= 0.174 vs domain-sharer <= 0.055 and control <= 0.083; the binding barrier is the nearest control at 0.083.)

## Contrast: minimizer-Jaccard does NOT separate them as well

Recomputing the OLD criterion (canonical minimizer-Jaccard, k=15/w=10, blake2b-64 -- identical to `family_detection_validation.py` and rustle production) for the SAME labeled pairs:

| criterion | threshold | confirmed accepted | domain-sharers rejected | clean separation? |
|---|---|---|---|---|
| **POA contiguous-core coverage** | T=0.13 | 5/5 | 7/7 | **YES** |
| POA all-column coverage (literal) | T=0.91 | 0/5 | 7/7 | no (overlap) |
| minimizer-Jaccard (best-fit) | J>=0.278 | 1/5 | 6/7 | no (overlap) |
| minimizer-Jaccard (production 0.30) | J>=0.30 | 1/5 | 6/7 | - |

Even at its *best-fit* threshold the minimizer-Jaccard cannot cleanly separate the two classes -- confirmed and domain-sharer Jaccard ranges overlap, which is exactly why the similarity-only grouping conflated them. Confirmed Jaccard range [0.126, 0.322], domain-sharer Jaccard range [0.098, 0.417]. (Note RFPL1/2/3 -- genuine recent dups -- have Jaccard ~0.13-0.17, LOWER than several domain-sharers, so no Jaccard cutoff can both keep them and drop the domain-sharers.)

*Display-rounding note:* GPR39<->LYPD1's Jaccard is 0.2779, which rounds to 0.278 -- the same 3dp as the best-fit threshold J>=0.278. It sits just BELOW the threshold, so it is rejected (x < t); the apparent tie is only a rounding artifact in the table.

## Per-pair table

| pair | label | contiguous-core cov | all-column cov | aligned identity | minimizer-Jaccard | n blocks | biggest block | len A | len B | same univ. family |
|---|---|---|---|---|---|---|---|---|---|---|
| RFPL1 <-> RFPL3 | confirmed | **0.608** | 0.842 | 0.833 | 0.156 | 37 | 898 | 1476 | 1722 | yes |
| RFPL1 <-> RFPL2 | confirmed | **0.606** | 0.611 | 0.896 | 0.126 | 72 | 895 | 1476 | 2414 | yes |
| RFPL2 <-> RFPL3 | confirmed | **0.520** | 0.711 | 0.879 | 0.165 | 65 | 895 | 2414 | 1722 | yes |
| APOBEC3D <-> APOBEC3F | confirmed | **0.302** | 0.471 | 0.907 | 0.140 | 20 | 860 | 2846 | 4843 | yes |
| RABL2A <-> RABL2B | confirmed | **0.174** | 0.814 | 0.652 | 0.322 | 130 | 416 | 2715 | 2395 | yes |
| CASP8 <-> FLACC1 | domain-sharer | **0.055** | 0.613 | 0.663 | 0.200 | 259 | 163 | 2958 | 4795 | yes |
| ASDURF <-> ASNSD1 | domain-sharer | **0.031** | 0.330 | 0.765 | 0.098 | 89 | 21 | 675 | 2037 | yes |
| GPR39 <-> LYPD1 | domain-sharer | **0.020** | 0.651 | 0.646 | 0.278 | 191 | 42 | 3259 | 2142 | yes |
| CCDC188 <-> ZDHHC8 | domain-sharer | **0.014** | 0.661 | 0.655 | 0.220 | 320 | 53 | 3751 | 5638 | yes |
| CDPF1 <-> PPARA | domain-sharer | **0.012** | 0.640 | 0.632 | 0.160 | 579 | 79 | 6706 | 10397 | yes |
| CREB1 <-> METTL21A | domain-sharer | **0.010** | 0.843 | 0.577 | 0.417 | 430 | 68 | 8050 | 7022 | yes |
| GCA <-> KCNH7 | domain-sharer | **0.008** | 0.268 | 0.822 | 0.175 | 572 | 33 | 4136 | 15456 | yes |

Control pairs (random cross-family, n=60, seed=1729): contiguous-core coverage median 0.033, max 0.083; all-column coverage median 0.427; minimizer-Jaccard median 0.000.

## Interpretation

On the POA **contiguous-core** axis the distribution is **bimodal and cleanly separable** (~3.2x vs domain-sharers, ~2.1x vs the nearest control): every confirmed tandem/recent-duplicate pair has ONE large homologous co-aligned block (>= 0.174 of the shorter gene), while every domain-sharer's largest single block covers only a few percent (<= 0.055) -- indistinguishable from random cross-family controls (<= 0.083). A single principled threshold T = 0.13 accepts ALL confirmed and rejects ALL domain-sharers AND all controls.

The LITERAL all-column reciprocal coverage does NOT achieve this, because a global alignment of two non-homologous transcripts inflates coverage with gappy chance-match filler (visible as low all-column identity and hundreds of tiny blocks). Requiring the co-alignment to be a contiguous homologous CORE -- still a pure alignment statistic, no domain/BLAST/k-mer input -- is the **POA-only fix for the domain-sharing false positives.** The minimizer-Jaccard cannot do this at all: it rewards any shared subsequence, so a shared domain inflates the intersection exactly like whole-gene homology, leaving confirmed and domain-sharer Jaccard values overlapping.

## Honest caveats

- **Operational RNA-structural definition, NOT the gene family.** This defines a *POA-coalignable multi-copy transcript family* (do the mature transcripts co-align over a majority of each?). It is deliberately NOT claimed equal to the DNA/protein gene family from a gene tree -- proving that equivalence would need DNA/protein/domain evidence, which the constraint forbids. It is an RNA-level, alignment-only operational criterion.
- **Pairwise-global as the POA pairwise instance.** For two sequences the POA instance reduces exactly to a single global alignment; for >2 copies a true POA DAG would be used and could differ slightly at multi-copy column assignment. The pair tests here are the exact two-sequence case.
- **All-column vs contiguous-core coverage.** The clean separation is on the contiguous-core axis, NOT the literal all-column reciprocal coverage (which overlaps because a global aligner pads non-homologous pairs with gappy chance-match filler). Both are alignment-column statistics and thus both POA-only/within-constraint, but the operational criterion that WORKS is specifically 'largest single UNGAPPED co-aligned block >= T of the shorter gene'. NOTE this is the largest *ungapped* block, not a local-alignment span: a Smith-Waterman local alignment of these domain-sharers ALSO stitches chance matches and reports a high aligned-length (verified: domain-sharer local aligned-len/minlen reaches 0.91-0.97, overlapping the confirmed pairs) -- so it does NOT separate. It is the longest CONTIGUOUS gap-free homologous run that distinguishes a real copy (one long block) from a domain-sharer (many short blocks patched together). That largest-ungapped-block read-off is what the criterion uses.
- **Small labeled N.** Only the Compara-checkable named pairs are labeled: 5 confirmed and 7 domain-sharers. The separation is an exact count on a small sample, a directional sanity check, not a population rate.
- **Gene-representative sequences.** Coverage is computed on one representative sequence per gene (`gene_rep.fa`), not all isoforms; a different representative isoform could shift coverage at the margin.
- **RABL2A<->RABL2B label (assumption is doing real work).** Treated as confirmed-real per task spec because it is the flagship 2-copy tandem pair. Compara returned a fetch error for RABL2B, but note that RABL2A's OWN Compara query SUCCEEDED (para_status=ok, 68 paralogues) and does NOT list RABL2B's ENSG (ENSG00000079974). So the 'confirmed' label here rests on the task-spec/flagship assumption, not merely on a one-sided fetch error -- RABL2A's own successful list also omits RABL2B. Reassuringly, RABL2A<->RABL2B is the LOWEST confirmed pair (contiguous-core 0.174), yet still cleanly above T=0.13, so it does not prop up the separation.
- **Labels are external corroboration only.** Compara is used purely to *check* the POA criterion's verdicts; it is never an input to the criterion. The criterion is alignment-derived end to end.
- **Determinism.** Controls use a fixed RNG seed (1729); alignment and minimizer machinery are deterministic. Output is byte-stable.


---

## family_to_copy_bridge

# Bridge: from the family definition to copy assignment

**PSV-on-backbone *are* the columns.** One object — the family's variation graph — is both the detection output and the
resolution input; identifiability is the single thread through both.

This note connects two formal pieces:
- **Detection** — `family_definition_note.md`: a multi-copy gene family = a $\sim_B$-connected component within a
  read-coupled ($R^\*$) group of de-novo loci. Output: the family's member loci **and** their shared **backbone** (the
  $\sim_B$ alignment of the all-isoform exon-union copy models).
- **Resolution** — `copy_assignment_theory.md`: reads as partial allele-vectors over **columns** $[m]$ (PSVs); the
  **conflict graph** $H$; $\mathrm{MCC}(R) = \chi(H)$ (Lemma 1); **Strong Separation** $\Rightarrow$ unique recovery
  (Thm 2); the **K-frontier** (recombination / collapse).

Detection produces exactly what resolution consumes. The map is below.

---

## 1. The handoff

Detection hands resolution a family $C \subseteq V$ and the backbone that $\sim_B$ already computed. Resolution needs
three things from its model (§2 of the theory): the **columns** $[m]$, the **copies** as allele-vectors, and the
**reads** as partial allele-vectors. All three are read off the family's variation graph.

Build the family's graph once: align the copy models $\{S(v) : v \in C\}$ (the same exon-union sequences $\sim_B$
aligned pairwise, now taken multi-way / POA). The graph is a shared **backbone** (consensus) with **bubbles** where the
copies diverge.

---

## 2. The bridge map (three identifications)

> **Columns = bubbles.** A **PSV column** $j \in [m]$ is a backbone position at which the copy models do not all carry
> the same base — a **bubble** of the variation graph — kept only if it is **read-supported** (recurs across reads, so
> it is a paralog difference, not a one-read sequencing error; the error-vs-recurrence test). The allele alphabet
> $A_j$ is the set of bases the copies present at $j$.

> **Copies = paths.** Each copy $v \in C$ realises an allele-vector $c_v = \big((S(v) \text{ aligned})_j\big)_{j\in[m]}$
> — its **path** through the bubbles. This is the theory's "gene copy" $c_v \in \prod_j A_j$.

> **Reads = partial paths.** A read $r$ threaded through the graph observes the bubbles it spans:
> $\mathrm{obs}(r) = \{\,j : r \text{ covers column } j\,\}$, and $r(j)$ is the allele it carries there. This is the
> theory's partial allele-vector.

With these three identifications, every object of `copy_assignment_theory.md` is defined directly from the family
graph, and the theory applies **verbatim**: the conflict graph $H$ on the reads, $\mathrm{MCC} = \chi(H)$ (the copy
number), Strong Separation (Thm 2) for unique recovery, and the K-frontier for the limit.

---

## 3. The two relations, in resolution's terms

The detection relations are not discarded — they are the *preconditions* that make the columns meaningful:

- **$\sim_B$ (shared backbone) $\;\Longrightarrow\;$ the columns exist.** Two loci are $\sim_B$ iff their copy models
  align reciprocally — i.e. they share a backbone on which the bubbles $D_{ij} \subseteq [m]$ (the distinguishing
  columns of the pair, §5 of the theory) are *defined*. Without $\sim_B$ there is no common coordinate to call a
  "column," so resolution is ill-posed. $\sim_B$ is precisely the well-posedness guarantee resolution assumes.
- **$\sim_R$ (read-confusability) $\;\Longrightarrow\;$ the problem is non-trivial.** Two loci are $\sim_R$ iff reads
  cross-map at tied divergence — i.e. the conflict graph $H$ is non-empty and the reads genuinely must be assigned.
  If the loci were not $\sim_R$, the reads would already place uniquely and $H$ would be edge-free.

So **a family ($\sim_R \cap \sim_B$) is exactly the regime where copy assignment is both well-posed (columns exist) and
non-trivial (reads conflict)** — which is the entire point of detecting it first.

---

## 4. One object, two readings — and the identifiability thread

The family's variation graph is read **twice**:
- **Detection** reads its *macro* structure: do the copies form long **parallel paths** (a shared backbone)? That is
  $\sim_B$; it certifies "these loci are copies."
- **Resolution** reads its *micro* structure: do the **bubbles** distinguish the paths, and do reads **link/cover**
  them? That is Strong Separation; it certifies "these reads can be assigned to copies."

Both are governed by the same bubbles $D_{ij}$, so a single quantity — the distinguishing-column structure — controls
the whole arc. This yields the precise containment:

> **Detection $\supseteq$ resolution.** The definition can *detect* a family that resolution *cannot* resolve. The
> family needs only $D_{ij} \neq \varnothing$ for the pair to be a paralog *and* reads to cross-map; resolution needs,
> in addition, that the bubbles are **covered and linked** so no recombinant or collapsed cover competes (Strong
> Separation). The **K-frontier** names exactly the gap:
> - $D_{ij} = \varnothing$ for some true pair $\Rightarrow$ those copies are one path (**collapse**, $K=0$ between them)
>   — detectable as a family, *provably unresolvable* into individual copies;
> - $D_{ij} \neq \varnothing$ but reads do not link/cover them $\Rightarrow$ recombinant covers compete (the $K\!\ge\!3$
>   phasing failure) — resolvable only with longer reads.

So the identifiability theorem is not a separate result bolted on after detection: it is the *same* graph, read at the
bubble level. Detection establishes the backbone on which the K-frontier is even stateable.

---

## 5. Empirical instantiation (`bench/psv_graph_demo.py`)

The two regimes, on real great-ape families, are the two sides of §4:
- **RABL2 (2 copies)** — the bubbles separate the two paths and reads cover them: Strong-Separation-like. **58 reads
  thread to a single copy, 100 % agreeing** with the best-mapping copy; the rest stay on the shared backbone (cover no
  bubble) — the *coverage* face of the K-frontier.
- **RFPL4A array (5 copies)** — 18 bubbles, but copies 2–5 take the **same allele at every bubble** ($D_{ij}=\varnothing$
  within the cluster): a founder (RFPL4A) plus **4 collapsed paths**. The graph *shows* the collapse — 6 reads resolve to
  RFPL4A, 5 only to the cluster (which of the 4 is unidentifiable) — the **$K=0$ collapse** face of the K-frontier.

The variation graph makes $D_{ij}$ visible (the bubbles), so "resolvable vs not" is read directly off the picture — the
identifiability dichotomy as graph structure.

---

## 6. Exact vs. pipeline

The **column abstraction is exact**: given the family graph, the bubbles, copy-paths, and read-paths are defined with no
free parameter, and the theory's results (Lemma 1, Thm 2, the K-frontier) hold as stated. The **pipeline** supplies one
calibration — the read-support threshold that separates a paralog bubble from a one-read error (the recurrence test,
§2-bubbles) — and inherits the de-novo locus vertices from detection. Both are the same disclosed dependencies the
family definition already carries; the bridge adds none.

---

## 7. In the code — the bridge is the production pipeline

Every object above already exists in the Rust path (`src/rustle/vg_family/`), and the copy-assignment driver `run_layer2`
(`layer2.rs:537`) walks exactly this chain:

| formal object (this note) | Rust function | location |
|---|---|---|
| **$\sim_B$ family** (copies share a backbone) | `contiguous_core_coverage` $\ge$ `T_CORE` (= the family criterion) | `family_detect.rs:15,31`; `family_graph.rs:933,1046` |
| **columns = bubbles** (PSVs from the copy backbone) | `psv_columns_from_reference` (POA-aligns the copies' exon sequences, emits a column per divergent position) | `psv_linkage.rs:159` |
| (node-shared PSVs, the co-located case) | `psv_columns_for_family` | `psv_linkage.rs:53` |
| **copies = allele-vectors / paths** | `PsvColumn.per_copy` (`PsvCopyAllele`) | `psv_linkage.rs:21,31` |
| **reads = partial allele-vectors** | `genotype_family_reads` → `ReadGenotype.psv_votes` (2nd BAM pass) | `psv_linkage.rs:517,420` |
| **assignment** (thread reads → copies) | `assemble_psv_isoforms` | `psv_linkage.rs:798` |
| **K-frontier gate** (skip if $<$ `min_psv_columns`) | `family_identifiability_gated` / `psv_columns_and_identifiable` | `psv_linkage.rs:340,350,372` |
| **entry point** (end to end) | `run_layer2` | `layer2.rs:537,985` |

So the abstraction and the implementation are the *same object*: `psv_columns_from_reference` POA-aligning the copies' exon
sequences **is** the $\sim_B$ backbone yielding the bubbles; `genotype_family_reads` + `assemble_psv_isoforms` **is**
threading reads to paths; `family_identifiability` **is** the K-frontier cut. The PSV-aware variation graph is not a future
plan — it is what the `--vg` layer-2 path runs, and it is unit-tested (`family_graph.rs` tests:
`family_merge_default_merges_near_identical_paralogs`, `fragment_in_larger_exon_does_not_merge_under_jaccard_default`,
`shared_low_complexity_tract_does_not_merge_at_default`).

**The one open lever (a deliberate flip, not a gap).** The copy-assignment *family-merge* step still gates the $\sim_B$
backbone-coverage check **opt-in** — `family_min_core_coverage()` defaults to `0.0` (off), enabled by
`RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (`family_graph.rs:416,447,1066`) — for byte-identical-default discipline, even though
`family_detect.rs` already uses the same $\sim_B$ test as *the* criterion. Making it the production default is a single
threshold flip ($\tau$), but it changes output, so it is owned by a genome-wide validation pass (the byte-identical-off
discipline), not a blind change. That validation is the last step to make the formalized $\sim_B$ definition the *default*
substrate of assignment, rather than the available one.

---

## Summary

The family definition does not merely *precede* copy assignment — it **constructs its input**. The $\sim_B$ backbone is
the coordinate system; its bubbles are the columns $[m]$; the copies are paths; the reads are partial paths; and the one
identifiability structure (the distinguishing columns $D_{ij}$) governs both whether the loci are a family and whether
its reads can be assigned. Detection and resolution are two readings of one variation graph — which is why making the
graph PSV-aware is the whole pan-transcriptomic copy-resolution method in a single object.


---

## family_criterion_bakeoff

# Conflict-edge tie criterion: AS vs NM vs de (the advisor's "don't trust AS" test)

## Question
The conflict-graph family definition links two loci iff `>= MIN_READS` reads cross-map between them with a
**tied** placement. The tie was defined on minimap2's **AS** (alignment score). The advisor's objection: AS is
a length + gap-model composite — the aligner's best guess, not raw evidence — and should not be taken at face
value. Do we get a better, threshold-defensible family definition by tying on raw-er signals: **NM** (edit
distance) or **de** (minimap2's gap-compressed per-base divergence)?

## Method (`family_criterion_bakeoff.py`)
Per labelled pair, collect every read's alignments in both member spans, attribute each **record** to its
best-overlap member (the same fix that killed the nested-gene FP in the Rust pipeline), and test the read's two
placements for a tie under each criterion:
- **AS-tie**: `min(AS) >= 0.9·max(AS)` (the shipped Rust default)
- **NM-tie**: `|NM_a/alnlen_a − NM_b/alnlen_b| <= D` AND `max <= CEIL`
- **de-tie**: `|de_a − de_b| <= D` AND `max(de_a,de_b) <= CEIL`

Edge iff conflict-read count `>= MIN_READS`. Panel: 5 Compara paralog pairs + 7 domain-sharers (labelled
ground truth) + 7 cross-chromosome paralogs (5–90 % identity) + 3 co-located arrays (MAGEA/PRAMEF/XAGE).
Verified by 5 adversarial sub-investigations (workflow `wf_01631729-ab7`); script audited sound.

## Headline finding — AS has a systematic false-positive mode on retrocopies
On **diverged processed-pseudogene / retrocopy** pairs, AS gives *near-identical scores* to placements that NM
and de show are 10–20× apart in quality, because the read's AS at the spliced **parent** pays per-intron
gap-open penalties that the single-exon **retrocopy** alignment never incurs — so AS ranks the worse copy
at-or-above the true parent and fires a spurious conflict edge:

| pair (parent ~ retrocopy) | AS edges | de/NM edges | median de gap | parent favoured |
|---|---|---|---|---|
| EEF1A1 ~ LOC109023808 | **3347 (FP)** | 0 | 0.0168 (~16×) | 3347/3347 |
| CNN2 ~ LOC129524764 | **121 (FP)** | 0 | 0.0084 (~10×) | 121/121 |
| ADH5 ~ LOC101125574 | **3 (FP)** | 0 | 0.0296 | 3/3 |

AS scores the *worse* copy strictly higher in 71 % (EEF1A1) / 90 % (CNN2) of reads. These reads are cleanly
resolvable to the intron-bearing parent (parent reads carry 6–8 CIGAR `N`-gaps the copy lacks); the AS tie is
an artifact, not ambiguity. This is the advisor's thesis proven on real data — and it matters because
retrocopies/pseudogene copies are a large slice of the copy landscape (1,313 pseudogene + 665 unannotated
copies genome-wide).

## Corrected confusion matrix (operating point D=0.005, CEIL=0.05, MIN_READS=3)
Negatives (should NOT fire): 7 domain-sharers + 3 retrocopies + APOBEC3 (asymmetric, minimap2-resolved 15/15) = 11.
Positives (should fire): RABL2A~RABL2B, AK6~LOC (jac .73), CCDC196~LOC (jac .90) = 3. RFPL excluded (0 cross-mapping reads — silent by absence).

| criterion | TP | FN | TN | FP | precision | recall |
|---|---|---|---|---|---|---|
| **AS** | 3/3 | 0 | 7/11 | **4** (EEF1A1, CNN2, ADH5, APOBEC3) | **0.43** | 1.0 |
| **NM** | 3/3 | 0 | 11/11 | 0 | 1.0 | 1.0 |
| **de** | 3/3 | 0 | 11/11 | 0 | **1.0** | 1.0 |

**`de` wins** (NM corroborates, AS demoted to audit). de is clip-robust where NM is not (APOBEC3D NM-rate
inflates to 0.16–0.22 on soft-clip/poly-A tails while de stays ~0.02), and has the wider safety margin
(CCDC196 fires 24 under de vs 16 under NM). Crucially **de-tie ⊂ AS-tie** on every pair — de removes exactly
AS's 4 false positives — so `de-tie is a strict subset of AS-tie` is the portable regression check.

## Adjudications
- **APOBEC3 not-firing is correct, not a FN.** Its cross-mapping reads fit APOBEC3D 2–3× better than APOBEC3F
  (median |de gap| 0.024; 0/15 reads fit both ≤5 %; minimap2 resolves 15/15 at MAPQ60-vs-MAPQ0). Contrast the
  genuine RABL2 fire: median |de gap| 0.006, 187/195 reads fit both ≤5 %, 118/195 MAPQ-0-at-both. APOBEC3 sits
  with the resolvable cases. It only fires under de at the loosest D=0.02/CEIL≥0.07 corner — spurious.
- **Arrays:** annotated MAGEA/PRAMEF/XAGE members show 0 cross-member conflict — a **test-design artifact**, not
  a miss. The pipeline's real MAGEA conflicts are between **unannotated de-novo sub-loci** 24–207kb from any
  annotated gene, and those cross-map (103/69/311 reads). The port's conflict graph must be built on the
  pipeline's coordinate sub-loci, not annotated-gene vertices.

## Operating point & stability
`D_DE=0.005, CEIL_DE=0.05, MIN_READS=3`. The tight corner is the unique fully-clean point on the 9-config
D×CEIL grid: D is the load-bearing knob (D→0.01 re-admits CNN2; D→0.02 re-admits EEF1A1 catastrophically), and
0.005 sits below the smallest genuine-conflict de gap (RABL2/CCDC196 ~0.004–0.006) and below the smallest
retrocopy gap (CNN2 0.0084). Domain-sharers stay silent across the entire grid (that axis is easy).

## Port plan (into `read_conflict.rs` + `build_read_placements`)
1. `BamRead`/placement carries per record: **de** (f32, `de:f` tag), **nm** (u32, `NM:i`), **aln_len**
   (query-aligned length = read len − clips). Keep **AS** for audit only.
2. Tie predicate becomes **de-tie** (`|de_a−de_b| ≤ 0.005 && max ≤ 0.05`); NM-tie computed in parallel for the
   audit log; AS-tie kept only to assert `de-tie ⊆ AS-tie` on real data.
3. Edge iff de-tied count `≥ MIN_READS=3` (raise from 2 to clear the CNN2 NM noise-floor; de alone is 0-FP even at 2).
4. **Replace the 200 bp distinct-record guard** with a CIGAR/identity comparison before trusting this on
   genuinely collapsed tandems (<200 bp-separated copies) — the guard is panel-safe only because every same-chrom
   copy here is ≥7 kb apart. (Open item.)
5. Make `RUSTLE_CONFLICT_{D_DE,CEIL_DE,D_NM,CEIL_NM,MIN_READS}` env-overridable, defaulting to the operating point.

## Open questions
The 200 bp guard is unprincipled (untested on collapsed tandems); AK6/CCDC196 should-fire labels are soft
(near-error-floor de ties minimap2 itself resolves by MAPQ); de is minimap2's own tag (a re-align could shift
it); the 0-FP result is a 14-pair curated panel (genome-wide unverified); the de-tie behaviour on the
pipeline's de-novo sub-loci vertices (vs annotated-gene vertices) is not yet measured.

---

## Tag-discriminator dig (does any minimap2 tag beat or augment de?)
Full tag vocabulary in GGO.bam: NM ms AS nn ts tp cm s1 s2 de rl. Dumped the per-read tag vector for every
cross-mapping read on the panel (`tag_discriminator_dump.py` → `.tsv`; genuine=246, artifact_retro=3540,
artifact_asym=15; domain-sharers=0 rows, pre-killed by best-overlap attribution) and evaluated 12 single tags +
combinations + a confound audit (workflow `wf_09afe78c-599`, 5 analyses + synthesis).

**Verdict: `de` wins unchanged; nothing beats it.** Two failure modes, and `de` already covers the one that
survives attribution:

| discriminator | independent of de? | per-pair result | verdict |
|---|---|---|---|
| **de-tie** `\|Δde\|≤0.005 & max≤0.05` | — (the axis) | **7/7, 0 FP/0 FN** | **WINNER** |
| `\|Δnmr\|` (NM/alnlen) | no — `nm/alnlen ≈ de` | false-fires CNN2 | de-in-disguise |
| `Δms`, `Δcm`, `Δs1` | no — corr w/ de −0.63/−0.53/−0.28 AND length-confounded (corr w/ alnlen .58/.79/.92) | false-fire CNN2/APOBEC3 | circular + length artifact |
| `s2/s1` (chaining uniqueness) | yes | inverted on retros (~1.0), 98% empty on copy side | **dud** (blind to retrocopies) |
| `nintron` (CIGAR N-gaps) | yes (structural) | FN on genuine AK6 (copy is itself an intronless retrogene → same signature) | independent but **label-confounded** |
| **`mapq`** | **yes (placement confidence)** | both-mapq0 = genuine; retro reads mapq60-at-true | **independent corroborator** |
| `ts`, `rl`, `min_mapq`, `ms_raw` | — | AUC ≈ 0.5 | inert |

**Net value of the dig:** (a) no tag catches a case `de` misses — `de`'s tie clause already rejects APOBEC3
(all 15 reads `|Δde|` 0.013–0.038 > 0.005), the case we hoped `s2/s1` or `mapq` would need to catch;
(b) `mapq` is the one genuinely-independent, non-confounded axis — worth **logging** per edge (fraction of
supporting reads with both-side `mapq==0` = genuine-multimapper signature) and as an **optional default-off
AND-gate** (`mapq_a≤40 & mapq_b≤40`) that widens the safe-`D` window 0.005→0.014 (2.8×) for noisier substrates;
(c) the `max(de)≤0.05` ceiling is redundant on this panel (the tie clause carries the discrimination) but kept
as cheap insurance. `D=0.005` is a knife-edge against CNN2 (false-fires 29 reads at `D=0.008`) — do not loosen
without the mapq gate.

**Refined port delta:** BamRead carries **`de` + `aln_len` + `mapq`** per record (DROP `nm` — it is
de-in-disguise and added nothing); keep `AS` only to assert `de-tie ⊆ AS-tie`; exclude **supplementary**
(chimeric/split ≠ multimapping ambiguity; the shipped Rust path currently includes it). Named consts
`DE_TIE_DELTA=0.005, DE_MAX=0.05, MIN_CONFLICT_READS=3` (no buried thresholds). Per-edge report logs the
both-mapq0 fraction as a confidence field. Do NOT port `nintron`/`ms`/`cm`/`s1`/`s2`/`ts` into the gate.


---

## family_detection_validation

# rustle paralog-family DETECTION: criterion validation
_Deterministic. Fixed seed = **20260616**. k=15, w=10. Generated by `bench/family_detection_validation.py`._
## Headline
> At threshold T=0.30 (rustle default), the criterion recovers 24/62 truth families as PURE clusters at pairwise P=0.995 R=0.478 F1=0.646 (ARI=0.644, V=0.964) — high precision with 1 between-family false merge (a CTRL-CTRL pair LOC129531078<->LOC129531352, Jaccard 1.0000) — itself a missed paralog leaked into the control pool, NOT a genuine cross-family error. It is recall-limited: 0.30 is ABOVE the within-family Jaccard mode at the whole-gene level (median 0.217), so 38/62 families fragment. The data-optimal whole-gene threshold is 0.06 (49/62 pure, clustering F1=0.732, ARI=0.729). NOTE: the truth set is itself a sequence-similarity clustering (minimap2), so this validates one similarity criterion against another, not against curated biology.
The data-optimal threshold (max raw pair-threshold F1, before union-find) is **0.06** (F1=0.824); rustle's 0.30 vs the optimum is discussed below. (NB: the Q2 table reports POST-clustering F1, which differs from the raw pair-threshold F1 because union-find applies transitive closure.)
## What is under test
rustle (`src/rustle/vg_family/family_graph.rs`) groups loci into paralog families when their sequences have **minimizer-Jaccard >= FAMILY_MERGE_JACCARD_DEFAULT = 0.30** (minimizers k=15, w=10). This script reproduces that exact similarity gate on labeled GENE sequences and union-find clusters them, then scores the predicted clusters against the annotated paralogy truth set.
## Inputs reconciled
- Truth: `universe.tsv` = **62 families** over **195 genes** (gene_id/family_id labels).
- Sequences: `gene_rep.fa` (1 representative per gene). Header IDs are bare names; universe uses bare names, so **195/195** universe genes matched directly (no `gene-` prefix present). Unmatched: 0.
- Control set: **150** NON-UNIVERSE genes sampled (fixed seed 20260616) from the 1669 fasta genes NOT in any universe family. Each is its own singleton label so that control-control and control-paralog pairs are BETWEEN and any control that clusters is scored as an error. **These controls are NOT certified single-copy** (Finding #4): they are merely absent from `universe.tsv`, which was itself filtered to minimap2 hits at ident>=0.90 / cov>=0.50, so a paralog pair below those cutoffs falls into the control pool. Residual control-control similarity: **6 CTRL-CTRL pairs >= 0.06** including LOC129531078<->LOC129531352 at Jaccard 1.0000 (a missed paralog leaked into the controls).
- **Control-purity SCREEN (independent minimap2, ident>=0.9/cov>=0.5; controls kept, not silently dropped, so the leak stays visible).** Two presets, because the discrepancy is itself the point:
  - `-x asm20` (EXACTLY the universe-builder criterion): **0/150 impure**. This MISSES the J=1.0 leak `LOC129531078`==`LOC129531352` — those two are byte-identical 153 bp sequences, but `asm20` returns NO alignment for sequences that short (it does not seed them). The universe builder used `asm20`, which is precisely WHY this real paralog pair leaked into the controls: the truth set's own criterion is blind to short paralogs.
  - `-x sr` (sensitive short-sequence seeding): **2/150 impure** — `LOC129531078`~`LOC129531352` ident=1.000 cov=1.000; `LOC129531352`~`LOC129531078` ident=1.000 cov=1.000. The minimizer-Jaccard criterion under test (k=15/w=10) is itself MORE sensitive than the universe's `asm20` on short sequences — it scored that pair J=1.0 — so it actually CORRECTS a truth-set miss here, which is the opposite of an error.
- Working test set: **345 genes** (195 paralogs + 150 controls).
## Q1 — Threshold justification (separability)
- Pairs: within-family **406**, between **58,934**.
- **AUC (within vs between) = 0.9820** (1.0 = perfectly separable).
- within-family Jaccard percentiles [5/50/95] = [0.003 / 0.217 / 0.934].
- between Jaccard percentiles [5/50/95] = [0.000 / 0.000 / 0.000].
- BETWEEN pairs >= 0.30: **1** pair(s) (= false-merge COUNT; the fraction 0.000017 is tiny but NON-ZERO — do not read it as zero). Breakdown: 0 cross-family paralog, 0 paralog-control, 1 control-control. Fraction of WITHIN pairs < 0.30: **0.6010** (missed-merge mass).
- The between-pair(s) >= 0.30 are: `LOC129531078`<->`LOC129531352` J=1.0000 [CTRL-CTRL].
- within-family **median Jaccard = 0.217** (below 0.30).
- Local pair density in [0.25,0.35) (the 0.30 band): **0.0009** of all pairs.
**Verdict (separability):** the two classes ARE separable on the between side (AUC=0.982): between-family pairs pile up near Jaccard 0 (unrelated genes share ~no canonical minimizers), and only **1** between pair reaches 0.30 — and that one is a CONTROL-CONTROL pair (LOC129531078<->LOC129531352, J=1.0000) that is itself a missed paralog, not a genuine cross-family error. BUT the within-family distribution is NOT a tight high mode — it is broad (5th/50th/95th pctile 0.003/0.217/0.934), because whole-gene paralog pairs range from near-identical to weakly similar (divergence + differing exon content). So the data is bimodal in the SENSE that unrelated pairs are cleanly separable from related ones near 0, but it is NOT cleanly bimodal at 0.30: 0.30 cuts through the BODY of the within-family mode, not a valley.
## Q2 — Family detection at 0.30 vs data-optimal
| metric | rustle 0.30 | optimal 0.06 |
|---|---|---|
| pairwise precision | 0.995 | 0.595 |
| pairwise recall | 0.478 | 0.951 |
| pairwise F1 | 0.646 | 0.732 |
| ARI | 0.644 | 0.729 |
| V-measure | 0.964 | 0.981 |
| homogeneity | 0.999 | 0.971 |
| completeness | 0.932 | 0.991 |
| families recovered PURE | 24/62 | 49/62 |
| families split (fragmented) | 38 | 7 |
| families over-merged (impure cluster) | 0 | 6 |
| clusters fusing >1 truth PARALOG family (control-blind) | 0 | 3 |
| clusters fusing >1 truth label (CONTROL-AWARE) | 1 | 6 |
| clusters contaminated by a control gene | 1 | 3 |
| # predicted clusters (>1 gene) | 36 | 56 |

> The paralog-only over-merge count is STRUCTURALLY BLIND to any false merge that involves a control gene (Finding #3): `family_recovery()` excludes controls from the over-merge tally, so a paralog↔control or control↔control fusion reports 0 there even though pairwise precision (which counts every pair) drops below 1.0. The CONTROL-AWARE row treats each control as its own family and is the honest false-merge count.

## Where does 0.30 fall? (data-driven)
rustle's 0.30 sits ABOVE the within-family mode (within-family median Jaccard = 0.217, and 60% of within-family pairs fall BELOW 0.30 at the whole-gene level). The data-optimal threshold is 0.06 (raw pair-threshold F1=0.824 vs 0.569 at 0.30; after union-find clustering F1=0.732 at 0.06 vs 0.646 at 0.30). So for WHOLE-GENE comparison 0.30 is OFF (too high): it fragments 38/62 families and recovers only 24/62 as pure clusters, trading recall (0.48) for high precision (0.995). It is NOT a zero-false-merge regime: 1 between-family false merge(s) at raw pairs >= 0.30 (0 cross-family paralog, 0 paralog-control, 1 control-control): LOC129531078<->LOC129531352 J=1.0000 [CTRL-CTRL].

**Interpretation / reconciliation with the source.** This is a WHOLE-GENE test, whereas rustle applies 0.30 at the EXON-class level (`refine_by_minimizer_jaccard` / `merge_singletons_by_sequence`), where homologous exons are short and near-identical so their Jaccard is much higher than the whole-gene Jaccard of the two genes. The source comment explicitly tunes 0.30 to REFUSE adversarial false merges (embedded shared TE/domain J~0.22, low-complexity cores J~0.17) at the exon level, and the between-family mass here confirms 0.30 buys that precision nearly for free: only 1 between-pair >= 0.30, and that one is the CTRL-CTRL J=1.0 missed-paralog leak, not a real cross-family merge. The honest finding is therefore: for the COARSE whole-gene grouping decision tested here, 0.30 is too conservative (data-optimal 0.06); a lower bar would recover 49/62 families pure at the cost of 6 control-aware false-merge cluster(s) (3 of them paralog-paralog, 3 control-contaminated). Whether to lower it depends on the downstream cost of an over-merge vs a split — which is exactly the assembly over-merge cost the source comment calls the binding constraint, and which THIS criterion-level test does not measure.
## Honest caveats
- **TRUTH-SET CIRCULARITY (read this first).** `universe.tsv` is NOT an independent biological ground truth. It is itself built by SEQUENCE SIMILARITY: `bench/copy_recovery_eval/30_build_universe.py` runs `minimap2 -x asm20 -X --no-long-join` all-vs-all on the per-gene representative sequences, keeps pairs at **ident>=0.90 and cov>=0.50**, then takes union-find connected components (`lib_eval.build_universe`). This validation then asks whether a DIFFERENT sequence-similarity criterion (minimizer-Jaccard, k=15/w=10) recovers those same families. So the high **AUC=0.982** and the recovery numbers measure CONCORDANCE BETWEEN TWO SIMILARITY CLUSTERINGS (minimap2 alignment identity vs minimizer-set overlap), NOT agreement with curated biology. They do not validate the families against Ensembl Compara, OrthoFinder, synteny, or a curated paralog database. Any 'recovers the truth families' language must be read as 'reproduces the minimap2-defined families'; an independent cross-check against a curated paralog DB is needed before a stronger claim. The two methods share enough mechanism (both reward shared subsequence) that part of the concordance is methodological, not biological.
- **Controls are NOT certified single-copy (Finding #4).** Controls are defined ONLY as fasta genes absent from `universe.tsv`. But the universe was filtered to minimap2 hits at ident>=0.90/cov>=0.50, so a paralog pair below those cutoffs is NOT in the universe and therefore lands in the control pool. This is not hypothetical: there are **6 control-control pair(s) with Jaccard >= 0.06**, the strongest being LOC129531078<->LOC129531352 at **J=1.0000** — a genuine paralog pair the minimap2 builder missed and that leaked into the controls. It is the single between-pair >= 0.30 and the sole reason RAW-PAIR precision at 0.30 is 0.9939 (162 TP, 1 FP), not 1.0; the 0.995 in the Q2 table is the slightly different POST-clustering precision after transitive closure. To certify the controls, screen the pool with the same minimap2 step and drop any gene with any hit; until then the control set should be described as 'non-universe', not 'single-copy'.
- **Criterion-level, not the full pipeline.** This tests ONLY the minimizer-Jaccard similarity GATE on whole-gene representative sequences. The production pipeline also uses position-overlap clustering, per-exon class refinement, EM, and flow; results here bound the gate's intrinsic discriminability, not end-to-end assembly.
- **Reference sequences only.** Run on annotation-derived representative gene sequences, not on assembled reads.
- **Expressed/annotated paralogs only.** The truth set is the 62 annotated families; unannotated or unexpressed paralogs are out of scope.
- **Collapsed copies are NOT detectable this way.** Paralogs that collapse onto one locus (no distinct sequence) cannot be separated by sequence similarity; that is the separate PSV-co-segregation test.
- **Minimizer hash differs from production.** Production uses non-canonical FNV-1a; mmh3 absent, so per spec this uses a CANONICAL minimizer hashed with blake2b (64-bit). The WINDOWING is reproduced exactly. Canonicalization only adds strand robustness and does not change the bimodal-separability conclusion.
- **Determinism.** Single fixed seed = 20260616 for the control sample; no wall-clock or unseeded RNG anywhere.


---

## psv_graph_genomewide

# Genome-wide PSV resolution: where great-ape gene-family copies are read-resolvable

The PSV-aware variation graph (`psv_graph_demo.py`, the two-family demo) scaled to **every**
validated multi-copy gene family. This is the empirical, genome-wide instantiation of the
identifiability arc in `family_to_copy_bridge.md`: detection hands each family a backbone, the
bubbles on that backbone are the columns, and threading the reads through the bubbles either
assigns them to a copy or hits the **K-frontier**.

## Method

Input: the 196 backbone-reinforced families from `family_def_vg_reinforce.py` (de-novo loci,
`~R ∩ ~B`). For each family:

1. **Copy-regions.** Collapse de-novo loci that are the same genomic copy (reciprocal overlap
   > 0.5 = isoforms / read-through nests) into one region. Families whose members all collapse
   to one region are **not genomically multi-copy** and are skipped (42 of 196 — single loci
   with isoform-multiplicity, e.g. a 335 kb transcript nested over a 178 kb one).
2. **Copy model = exon-union.** Each copy's sequence is its **exon-union model** S(v) — the
   reads' aligned exon blocks unioned, introns removed (`family_def_vg_reinforce.py:build_model`,
   cached in `vg_reinforce_copies.fa`). Backbone = longest copy model; the others aligned
   `minimap2 asm20`, the reads (intron-free mature mRNA) aligned `minimap2 map-hifi`.
3. **PSV bubbles** = backbone columns where the copies carry ≥2 alleles **and** ≥3 reads support
   the column (the recurrence test separating a paralog difference from a one-read HiFi error).
4. **Copy paths** = each copy's allele string; **K** = number of distinct paths (resolvable
   haplotypes). **Thread** each read by best PSV-allele match (HiFi-tolerant).

### Why exon-union, and why only 24 families were re-run on it

A first pass aligned copies by their **genomic span** (introns included). That is cheap but
fails when copies have divergent introns — `asm20` cannot align across them, so the copies'
exonic PSVs vanish and the family is falsely called collapsed (6 `align_fail` families), and it
can also *miss* exonic PSVs that lie near intron boundaries (it over-called K=0 on, e.g.,
fam62). Aligning the **exon-union** models fixes both.

Re-running every family on exon-union is unnecessary, and the reason is a monotonicity argument:
exon-union alignment only ever surfaces **more** distinguishing columns, and adding a column
cannot lower K. So a family already resolvable under genomic-span **stays** resolvable — only
the 24 **frontier candidates** (18 genomic-K0 + 6 align_fail) can move. Control re-runs of three
resolvable families (fam18/0/4) confirmed they are unchanged (fam4: 3746/3761 reads agree).
Only the 24 candidates were re-aligned on exon-union (`psv_graph_exonunion.py`); the 121
resolvable families are carried over from the genomic-span pass (`psv_graph_combine.py`).

## Verification (`psv_graph_verify.py`, `psv_graph_combine.py`)

- **Cross-family dedup.** Nine families were the *same* genomic copies re-detected as several
  de-novo isoform-loci (e.g. 168/169/170 = one 2-copy locus counted three times). Collapsed by
  region-set match: **154 → 145 unique families**, read totals recomputed over uniques.
- **Frontier re-alignment.** Of the 24 frontier candidates, exon-union **recovered 10 to
  resolvable** (5 ex-K0: fam7/30/62/122/145; 5 ex-align_fail: fam1/14/21/75/175) and confirmed
  **14 as genuine K=0** (incl. fam34: 8 copies, 0 PSVs over 596 reads → truly identical). The
  `align_fail` bucket went to **zero** — every copy aligns once introns are removed.

## Result (corrected, 145 unique families)

| class | families | % | meaning |
|---|---|---|---|
| **fully resolvable** (K = #copies) | 123 | 84.8 | every copy read-distinguishable |
| **partial** (2 ≤ K < #copies) | 8 | 5.5 | some copies separate, some collapse |
| **genuine K = 0** | 14 | 9.7 | copies exonically identical → provably unresolvable from RNA |
| **indeterminate** | 0 | 0 | (eliminated by exon-union) |
| **read-resolvable total** | **131** | **90.3** | |

**Reads:** 64,066 threaded across the unique families — 35,329 (55 %) assigned to **one copy**,
1,344 to a collapsed group, 19,922 (31 %) cover **no PSV** (the *coverage* face of the
K-frontier — a full-length read spanning only shared exons), the rest unexplained/ambiguous.
Single-copy assignments **agree with the independent best-mapping copy 95.4 %** of the time.

## Why this matters

- **The unresolvable core is ~1 family in 10.** 9.7 % genuine K=0 from RNA read-threading is the
  same magnitude as the ~12 % K=0 from the combinatorial copy-assignment census
  (`copy_assignment_theory.md` / `resolution_improvement_bound.md`) — two unrelated methods,
  over different family sets, agree the collapsed core is order-10 %.
- **The dichotomy is structural, not anecdotal.** 85 % of families lie on the `K = #copies`
  diagonal (Fig. B): the bubbles separate every path. The K=1 row is the frontier. The demo's
  RABL2 (resolvable) and RFPL4A (collapsed) are the two ends of a genome-wide distribution, not
  cherry-picked cases.
- **Detection ⊇ resolution, measured.** All 145 families were *detected* (each is a `~R ∩ ~B`
  component); 90 % are also *resolvable*. The 9.7 % genuine-K0 gap is exactly the bridge note's
  prediction: detectable as a family, provably unresolvable into individual copies.

## Honest caveats

- **Validation is a proxy, not ground truth.** "Agreement with the best-mapping copy" treats
  minimap2's primary as truth, which is exactly what is unreliable for paralogs. 95.4 % global
  agreement is reassuring; the ~5 % disagreement concentrates in size-heterogeneous families
  (family 32: a 124 kb backbone bundled with ~30 kb paralogs) where mapping itself is arbitrary —
  there, PSV disagreement may be a *correction*, not an error.
- **Low-coverage K=0.** Of the 14 genuine-K0 families, five carry < 50 reads (fam176 has 8);
  their copies show 0 divergent columns, but at that depth a rare PSV could be unsampled. The
  well-covered K=0 calls (e.g. fam34/46 with 364–596 reads and still 0 PSVs) are firm.
- **De-novo loci are imported** from detection; the copy-region dedup (overlap > 0.5) is the one
  free parameter added here.

## Performance note

The bottleneck is the **Python pileup + per-read threading** (O(backbone × depth) per family),
not `minimap2`. High-coverage families (3,000–4,000 reads) dominate. The Rust production path
(`psv_linkage.rs::psv_columns_from_reference` + `genotype_family_reads`) does exactly this
extraction natively and would run ~10–100× faster — this Python scan is an analysis harness, not
the shipping pipeline. For the harness itself, capping reads (the classification needs only
enough depth to clear the ≥3-read PSV test) is the cheap speed-up.

Artifacts: `psv_graph_genomewide.py` · `psv_graph_verify.py` · `psv_graph_exonunion.py` ·
`psv_graph_combine.py` · `psv_graph_genomewide.png` · `psv_graph_exonunion_all.json`.

## divergence_floor — absolute reciprocal-identity floor on the E_r edge rule (SHIPPED, default-ON)

An **absolute reciprocal-identity divergence floor** on the cross-gene edge rule of the RNA-only
family definition (`bench/family_rna_refine.py`), **DEFAULT-ON with opt-out**. It scopes the
catalog to families whose members are `>= FLOOR` reciprocal whole-transcript identity, at high
precision (single-exon domain-shares and divergent paralogs excluded). Soto's SD98 move (scope to
the clean near-identical regime) with an opt-out that recovers the ambitious divergent-inclusive
catalog **byte-identically**. **FAMILY_DEF is the sole home of this floor.**

**Default FLOOR = 0.80** (`MIN_FAMILY_IDENTITY = 0.80`). Opt out with `--no-divergence-floor` /
`--min-family-identity=0` / `RUSTLE_NO_DIVERGENCE_FLOOR=1`; set a different floor with
`--min-family-identity=X` (e.g. `0.85`).

### The metric

```
recip_id_best = matches_best / max(len_a, len_b)
              = min over the two members of (aligned-AND-matching bases / member length)
```

Computed from the **SAME** spliced-exon minimap2 alignment (`-cx asm20`) that defines `aln_frac` /
`core_recip` — only the residue-match count and member lengths are added (`ri_build_recip_id.py` →
`ri_sharedlen_recip_id.tsv`, `aln_len` reproduces the universal cache byte-for-byte). **RNA-only,
NOT DNA.** It works because it is RECIPROCAL (min over both members): a domain-share where a small
gene sits inside one big exon of a large gene scores high on the small side but LOW on the large
side (`matches / len_large`), so the min cuts it. Shared-block identity or `aln_frac` alone do NOT
work — domain-shares share one exon at ~95%.

Class medians on the 1457 core+aln-passing edges (`bench/divergence_floor.tsv`): TP real cDNA
paralog n=375 median **0.615**; genuine over-merge n=619 **0.295**; truthbar divergent paralog
n=463 **0.320**. The reciprocal length penalty separates near-identical real paralogs (0.615) from
divergent/domain-share over-merges (~0.30) by design.

### P/R-vs-floor curve (diploid gold oracle; `PYTHONHASHSEED=0`, gamma 0.20, seed 0)

`R_oracle` = diploid multi-copy oracle genes recovered / 57; `P_dedup` = 1 − distinct-FP /
oracle-mapped; `E_p` = protein-family block purity; domShare = # of 4 named domain-shares excluded.

| floor | nFam | R_oracle | P_task / P_dedup | E_p | distinctFP | domShare excl | md5 |
|-------|------|----------|------------------|-----|-----------|---------------|-----|
| **off** | 573 | **0.877 (50/57)** | 0.896 / 0.917 | 0.892 | 4 | 0/4 | **dca64cbd** (opt-out) |
| 0.70 | 335 | 0.596 (34/57) | 0.879 / 0.909 | 0.961 | 3 | 3/4 | 0c433468 |
| 0.75 | 320 | 0.579 (33/57) | 0.909 / 0.939 | 0.966 | 2 | 3/4 | 4c8c88ee |
| **0.80** ⭐ | **307** | **0.561 (32/57)** | **0.903 / 0.936** | **0.967** | **2** | **4/4** | **e84dc2bc** (default) |
| 0.85 | 293 | 0.509 (29/57) | 0.893 / 0.929 | 0.976 | 2 | 4/4 | b27512ce |
| 0.90 | 284 | 0.491 (28/57) | 0.926 / 0.963 | 0.979 | 1 | 4/4 | 0b9bf57b |

`floor=off` reproduces the divergent-inclusive default catalog **dca64cbd BYTE-IDENTICAL** —
the floored build path is validated and the opt-out guaranteed exact.

### Domain-share exclusion + FP composition

The 4 named single-exon domain-shares: MOV10+RHOC 0.111, RBBP4+SYNC 0.350, DEDD+NIT1 0.533 (all
excluded at every floor), and **RHD+SDHD 0.785 = min(0.785,0.879)** — the tightest, excluded only
at **floor >= 0.80**. This is exactly why 0.80 is the minimum knee (0.75 gives 3/4). Verified at
the catalog level (each gene absent from the multi-copy catalog), not merely a direct-edge cut.

**FP residual is clean and STABLE from 0.80 on** — at 0.80 and 0.85 the distinct diploid-oracle FP
= **2**: (1) `LOC129529978 + LOC129529986` = the irreducible **MAGE-X** DNA-only array over-size
floor (no RNA metric removes it), and (2) `GSTM2 + LOC101129940` = one **GST protein-domain hub**.
E_p-impure blocks drop **79 (off) → 10 (0.80) → 7 (0.85)**. At 0.90 the GST hub is also cut
(distinctFP → 1, MAGE only) but at a cost of 4 more oracle genes.

### Chosen default = 0.80 (dominates 0.85)

0.80 is the empirically-best knee: identical distinct-FP residual (2, same two blocks), identical
4/4 domain-share exclusion, retains **3 MORE** oracle genes (32 vs 29, +10% relative recall) at
essentially equal precision (E_p 0.967 vs 0.976), and is the minimum floor excluding RHD+SDHD
(0.785). User's suggested **0.85 is validated** (4/4 out, distinctFP 2, E_p 0.976) — a defensible
marginally-conservative alternative but buys no distinct-FP reduction for a 3-gene recall cost.
Precision improves vs floor-off on the load-bearing axes: E_p 0.892 → **0.967**, P_dedup 0.917 →
**0.936**, distinctFP 4 → **2** (the task-formula P is denominator-noisy — `oracle_mapped` shrinks
48 → 31 with the catalog — so E_p / distinct-FP are the load-bearing precision signals).

**Scoped claim:** *Gene families >= 80% reciprocal whole-transcript identity: P_dedup = 0.936 (E_p
0.967; distinct diploid-oracle FP = 2 = irreducible MAGE-X DNA over-size floor + one GST domain
hub; all 4 named domain-shares excluded), R_oracle = 0.561 (32/57).*

### Opt-out (byte-identical) and honesty

`--no-divergence-floor` / `--min-family-identity=0` / `RUSTLE_NO_DIVERGENCE_FLOOR=1` all recover
the divergent-inclusive catalog **dca64cbd byte-identical** (573 families, R_oracle 0.877).
Mechanism: when floor <= 0 the per-edge `recip_ok()` returns `True` unconditionally and the
reciprocal-identity cache is never loaded — the pre-floor code path is bit-for-bit unchanged. The
floor **composes** with the other gates (`--no-repeat-gate`, `--no-split-recombinants`,
`--no-repeat-bridge-gate`, `--high-precision` gamma 0.20 → 0.40), each ablation touching only its
own gate. Honest scoping: making the floor default DROPS R_oracle 0.877 → 0.561 (18 lost oracle
genes) BY DESIGN (excludes divergent + partial/length-mismatched homologs, SD98-style); the
ambitious catalog is one flag away. RNA-only guard hard-asserts the floor feature set is exactly
`{recip_id_best}`, disjoint from any DNA/protein/genome/soft-mask column.

Files: `bench/family_rna_refine.py` (floored edge rule) · `bench/divergence_floor.py` (sweep) ·
`bench/divergence_floor.tsv` (P/R curve) · `bench/ri_build_recip_id.py` +
`bench/ri_sharedlen_recip_id.tsv` (reciprocal-identity cache) · `bench/test_family_rna_refine.py`
(6 tests: opt-out byte-identity, domain-share exclusion, scoped precision improves, env==flag,
composes, determinism).

## multi_repeat_bridge_gate — 3rd VG-native family-refinement gate (SHIPPED, default-ON)

The **multi-repeat-bridge conjunction gate**, WIRED as the third VG-native gate in
`bench/family_rna_refine.py`, **DEFAULT-ON at the conservative T=8/C=2**, opt-out
`--no-repeat-bridge-gate` / `RUSTLE_NO_REPEAT_BRIDGE_GATE=1`. It removes the DISCONNECTED
repeat-bridge over-merge class — the dominant residual (34/83 E_p-impure blocks, 41%) where the
pairwise E_r oracle merged loci sharing **no full exon** at id ≥ 0.70 (cross-component best-exon-id
median 0.586 = a sub-exonic Alu/repeat bridge, not real transcript homology).

### The gate

```
CUT family iff
  (NO full shared exon : ≥2 full-exon components AND cross-component best-exon-id < 0.70)
  AND (REPEAT-bridged  : ≥ C distinct cross-component shared VG minimizer nodes,
                         each with global multiplicity ≥ T)
  AND (same-gene guard : never separate loci of the same annotated gene)
On CUT: replace the family by its full-exon components; <2-locus components drop via multi-copy filter.
```

**Why the conjunction, not either half:** multiplicity-alone can't drop below the shipped
single-edge `min_shared_mult ≥ 20` gate without shredding GSTM2 (its internal Alu node has mult
**426**, MAGE 8–10) — but GSTM2/MAGE form a **single** full-exon component so the no-full-shared-
exon conjunct means no cross-component pair exists and the gate **structurally cannot fire**.
Exon-containment-alone was FALSIFIED (`FAMILY_DEF.md` — cut real single-gene
paralogs that fragment on UTR noise/exon-skip 1:1); the repeat conjunct + same-gene guard spare
those. Class signature (63 families): DISCONNECTED_FP `frac_no_full_shared_exon = 1.0`, sub-20
cross-component nodes; GSTM2/MAGE connected single components; real-paralog controls 82.6% share a
full exon.

### Shipped deployment (catalog-wide, T=8/C=2)

The gate is a standing predicate applied to **every** multi-copy family. At T=8/C=2 it cuts **61
DISCONNECTED families (605 → 573)**; every cut is DISCONNECTED (no full shared exon) by
construction. Named cuts: the **fam17 20-gene Alu hub** (`C9H11orf65+…+UCHL1`), PDIA3+PRKAB2, the
MPHOSPH8 and LOC134758618 collapsed-array blobs, FGD3+HIVEP1, KYAT3+RBMX, etc.

Precision (headline; improves), pre-bridge md5 **5e58378a (605)** → default **dca64cbd (573)**:

| metric | pre-bridge (605) | default gate-ON (573) |
|---|---|---|
| R_oracle (named diploid-DNA gold) | 50/57 = 0.8772 | **50/57 = 0.8772 — HELD** |
| E_p purity | 0.8694 (79 impure) | **0.8918 (62 impure) +0.0224** |
| distinct-FP blocks | 6 | **4** |
| oversize FP | 3 (MAGE-Xarray, MPHOSPH8, LOC134758618) | **1 (MAGE-Xarray only)** |

`--no-repeat-bridge-gate` recovers 5e58378a **byte-identical**. GSTM2 (fid 9/13), MAGE
(546/548/553), and the 23 real-paralog controls: **0 cut** at every (T,C) in both scopes
(structurally — single full-exon component). The only "GSTM2"-named cut is the fid-23 satellite
`GSTM2+LOC129532045+LOC134757399` (GSTM2 Alu-bridged to unrelated LOCs) — the same-gene guard
moves both GSTM2 loci together, dropping only the non-GSTM single-locus passengers.

### Sensitivity cost (disclosed) — R_oracle is blind by construction

R_oracle and E_p are BLIND to the gate's only cost (single-locus glued passengers the `<2-loci`
filter drops), because the same-gene guard structurally preserves per-gene recovery. The operative
sensitivity probe is the conservative cDNA-pair truth `kept_REAL` (`truth2_dna_loose`): **pair-
recall 0.9196 → 0.9087**, **DNA component-recovery 182/187 → 177/182** at the deployed catalog
scope (`kept_REAL −5`, NOT the roster-scope −1). Of the 5 lost REAL cDNA pairs: **4 are repeat-glue
LABEL-NOISE** (Alu-glued no-shared-exon pairs mult 231–503 the cDNA-loose oracle mislabelled REAL,
arguably correct to cut: `ASB6+NTMT1`, `CCDC152+SELENOP`, `ETFDH+PPID`, `GATC+SRSF9`), and **only 1
is near-threshold-divergent** — `LOC109024259+ZNF224` (exon-id 0.686, nicked at the strict 0.70
full-exon cutoff). At T8C2 **ZNF224 keeps all 10 KRAB-ZNF paralogs** minus only the single divergent
`LOC109024259`; **RFPL3 is UNTOUCHED** (its bridge mult 5 < T=8). The genuine paralog *cores*
survive. **Causal honesty:** the recall-safety comes from the **same-gene guard + no-full-shared-
exon conjunct, NOT the multiplicity threshold** (plain no-full-shared-exon + guard already cuts
32/35; the repeat conjunct C≥2 *narrows* to 15/35 by abstaining on the low-mult glue tail — it is a
precision knob). The name over-credits the repeat conjunct; stated plainly.

Operating points: **T=8/C=2 (shipped default)** — 15/35 roster DISCONNECTED removed, E_p +0.022,
0 genuine paralogs lost, honestly "several repeats". **T=5/C=1 (aggressive)** — 24/35 removed, E_p
+0.034, still 0 genuine paralogs lost in roster scope. Runs **after** the recombinant-split stage,
**before** allele-demote. Determinism: default md5 **dca64cbd** across runs (note md5 **e2f0a23a**
= the `--high-precision` gamma=0.40 catalog, not nondeterminism). Artifacts:
`bench/multi_repeat_bridge_gate.py` · `.tsv` (63-fam characterize, md5 5e40411a) · `.json` (gate
sweep, md5 043b22d7); tests `test_repeat_bridge_gate`, `test_repeat_bridge_r_oracle_held`,
`test_repeat_bridge_composes`.

## Catalog artifact audit (2026-07-10) — KRAB-ZNF domain bridge under --homology-primary

Ran `copy_assign` on 30 annotated multi-copy loci (2–12 paralogs, span < 900 kb from
`GGO_genomic.gff`) plus the session panel; classified every family formed (dense LOC/ZNF
mega-clusters >25 copies excluded for cost). **The default (E_c) catalog is clean of over-merges on
this panel** — 16 families formed, every one a single named gene family with span-disjoint copies
(GBP, PCDHB 5, APOBEC, APOL, TCEAL, GSTM 3, MAGEA 2, DAZ 2, TSPY 5); the 2 flagged Containment
cases (RFPL, r4) are the known coverage floor, reported-not-pruned (`DENOVO_PIPELINE.md`).

**One new artifact class, `--homology-primary`-ONLY: a KRAB-ZNF domain bridge.**
`NC_086017.1:54339728-54640140` emitted a 3-copy family of three **distinct** annotated genes —
ZNF445 / ZKSCAN7 / ZNF197. Genomic spans align at 43–68% identity over ~1% coverage (only the
shared zinc-finger exon); their **spliced exon-sum transcripts do not align at asm20 at all** (the
0.80 tier). The edge was formed by E_r's **sensitive tier (0.60 identity)**, which recovers ancient
paralogs but also bridges the zinc-finger domain shared across the whole KRAB-ZNF superfamily.
**Under the default E_c conflict oracle this cluster forms 0 families** (the three genes map
uniquely, no read is ambiguous) — so the over-merge is specific to the opt-in `--homology-primary`
path. This is the "238-blob" domain-bridge trap the genome-wide family-definition work solved with
co-threading + community de-bridging; `copy_assign`'s inline `--homology-primary` branch
(`denovo_pipeline.rs:1256–1260`) uses `homology_edges_all_reps` + a plain
`gamma_quasi_clique_partition(0.20)` — raw E_r edges with **no de-bridging**, so it recovers real
ancient paralogs *and* domain bridges indiscriminately.

**The reciprocal-coverage fix was TESTED and it FAILS — reverted (keep the lesson).** The natural
fix, reciprocal coverage `min(qcov, tcov)` instead of coverage-over-shorter, correctly shrank the
ZNF bridge (3 → 2 copies) but **destroyed four real families**: GSTM 3→0, MAGEA 2→0, PCDHB 5→0,
TCEAL 5→0, and dropped GBP 6→4. Reverted. **Reason (fundamental):** reciprocal coverage penalises
transcript-length differences, which are biologically normal (isoforms, UTR extensions, truncated
de-novo models — DAZ2 is itself a 70% model). The ZNF bridge (shorter-cov 0.55, id 0.69) sits in
the *same pairwise cell* as real sensitive-tier paralogs. **Measured lesson: no pairwise
coverage/identity threshold separates a domain bridge from a length-divergent real paralog** — the
coverage-over-shorter metric was chosen deliberately to be robust to this, and this is precisely
why domain bridges are solved by **graph structure** (co-threading + weighted-Louvain de-bridging),
not a pairwise gate. That graph fix **cannot be ported to `copy_assign`'s per-region path**: the
ZNF over-merge is a 3-node complete triangle any γ-quasi-clique accepts; de-bridging needs the
genome-wide graph where the whole KRAB-ZNF blob is visible. **Disposition:** do NOT ship a
coverage/identity threshold (breaks real length-divergent families, and is the arbitrary threshold
the advisor rejects); the KRAB-ZNF bridge is an inherent limitation of the opt-in per-region
`--homology-primary` mode (raw E_r, no genome-wide context) — for a clean real-family catalog use
the genome-wide de-bridged path (`detect_homology_catalog_genome_wide` + co-threading). Minor note:
GBP (6 annotated) fragmented as 4+2 and TCEAL (6) as 3+2 — *under*-merging (less harmful, no false
copies, but understates copy number). Related: `DENOVO_PIPELINE.md`,
`project_rna_family_homology_primary`, `project_family_def_readconflict`.

