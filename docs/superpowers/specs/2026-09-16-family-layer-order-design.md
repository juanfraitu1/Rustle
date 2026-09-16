# Family layers: an SD layer that derives Soto, and an empirical layer order — Design

Status: design approved in chat 2026-09-16 (user); this document pins the details for review before the
pre-registration is written. Read-only inventory of inputs: session scratchpad `layer_order_inputs.md` (paths are
re-verified in the implementation plan).

## Why

The user asked whether O1 has "a graph-based definition for every level that corroborates itself more the more
information we have". Two findings block a yes:
1. `docs/seeded_family_definition.md` §0★★ clause 6 ("every RNA family lies inside one DNA family") is a theorem
   for connected components only. Under clause 4 (triangle-supported leaders) it fails: a hand-checked 6-node
   counterexample, and one real violation on gorilla NPIP (0 on TBC1D3) — memory
   `project_leader_rule_breaks_nesting`.
2. The project already has more family layers than DNA vs RNA — a protein rule (§6ko) and an SD rule (shared-exon
   ≥ 0.30, §6kr–§6kt) — and they are measured NOT nested with the DNA copy graph (only 52% of DNA-family pairs are
   protein-family pairs; 25% of protein pairs are DNA pairs, §6ko).

The user asked to (a) find the adequate layer order empirically rather than assume it, and (b) have the SD layer
derive the Soto et al. 2025 (Cell) human gene families.

## ⭐ SCOPE AMENDMENT (user, 2026-09-16 13:56): NPIP + TBC1D3 ONLY, descriptive — run this first

The advisor is only interested in NPIP and TBC1D3 for now. The chromosome-wide study below (Phases A/B, held-out
chr8–11, pre-registration) is DEFERRED, kept as the follow-up for when a general order claim is needed.

Scoped study (human T2T-CHM13; descriptive; NOT pre-registered because both families were used to develop the DNA
rules, §6js — every statement is "for NPIP and TBC1D3", never a general order claim):
- **Universe U:** annotated NPIP and TBC1D3 members (RefSeq CHM13 names `NPIP*`, `TBC1D3*`, incl. pseudogenes) ∪
  every gene that ANY layer below places in the same group as a member (closure, genome-wide, cross-chromosome
  allowed). Without the closure every containment is trivially 1.
- **Layers on U:** P (§6ko tables), D (guided catalog, no SD filter), S1, S2, S3 (as defined below), EXPR (testis
  reads ≥ 3 per gene, primary `-F 2308`, overlapping exons), and **C — subfamily clades** (existing §6jp/§6jr
  results: NPIPA|NPIPB, NPIPB3-5 / B6-9 / B12-13; TBC1D3 clusters) as the finest layer.
- **Truths:** Soto family IDs containing NPIP/TBC1D3 members (for S), literature subfamilies (for C), HGNC gene
  groups (for P/D), the guided annotation truth (for D).
- **Outputs:** per-layer group membership table on U; containment matrix c(X ⊇ Y) and group-level nesting for all
  layer pairs; the resulting tournament order (reported with a cycle if one exists); per-layer truth agreement
  before/after enforcing nesting (JOIN above D, REFINE below D, EXPR per layer); a gene-level disagreement list
  (which genes P joins that D does not, which co-duplicated neighbours S splits from D, which members are
  unexpressed). SD-layer selection among S1/S2/S3 by agreement with Soto on these families (descriptive).
- **Report:** `bench/LAYER_ORDER_NPIP_TBC1D3.md` (terse tables); results under
  `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/`.

## Definitions

### Layers (each computed independently first, on the same gene universe)
- **P — protein families.** §6ko rule as shipped in `bench/protein_families.py`: longest CDS per gene, blastp
  e ≤ 1e-5, coverage ≥ 0.30 of the longer protein, MCL inflation 2.8, r2 filter (drop pseudogene biotypes and
  V(D)J segments). Defined on coding genes only.
- **D — DNA copy-graph families.** Guided mode on annotated genes, clause-2 edges (exon edge: spliced transcript
  identity ≥ 0.80 over ≥ 0.50 of the transcript; gene-body edge: identity ≥ 0.80 over ≥ 0.50 of the shorter
  body), grouped as the guided catalogs are built (`mcl_families`), with `--min-shared-exon-frac 0.0` (no SD
  filter).
- **S — SD layer candidates** (one is selected on dev, then frozen):
  - **S1** — D's graph with the shipped shared-exon rule (`mcl_families --min-shared-exon-frac 0.30`).
  - **S2** — Soto's construction without copy number: SD98 (SEDEF segmental duplications > 98% identity) region
    DNA mapped back to the genome with `minimap2 -c --end-bonus 5 --eqx -N50 -p0.5`, self-maps removed; two
    genes linked when a mapping covers ≥ 99% of an exon of each (shared exon); groups = connected components.
    Soto's famCN/parCN refinement is EXCLUDED (read-depth copy number; also circular with the truth — see
    `bench/soto/famcn_from_wssd.py` docstring).
  - **S3** — S2 computed within each D family (S2 edges only between genes in the same D family).

### Operators that make nesting hold by construction
- **JOIN upward** (for a coarser layer that disagrees both ways, e.g. protein): the coarse group = connected
  components of the graph whose vertices are the finer layer's groups and whose edges join two groups when any
  member pair is linked in the coarse layer. Every finer group lies inside one coarse group. Genes outside the
  coarse layer's universe (e.g. non-coding copies for P) ride along with their finer group.
- **REFINE downward** (for a finer layer): recompute the finer layer's grouping on the subgraph induced within
  each coarser group (edges between groups dropped). Every finer group lies inside one coarser group.
- **EXPR (RNA operator)** applied to any layer L: within each L group, the connected components of the subgraph
  induced on expressed genes (reads ≥ 3) with ≥ 2 genes.
- Properties (short proofs go in the definition doc): (T1) JOIN and REFINE yield nested partitions; (T2) if M
  refines L then EXPR(M) refines EXPR(L); (T3) for fixed L and fixed edges, adding expressed genes never moves a
  gene between EXPR(L) groups (groups only grow or merge). Not guaranteed: clause-4 leaders remain read-order
  dependent; de novo families are not nested inside guided families.

### ⚠ AUDIT AMENDMENT (2026-09-16, after the NPIP/TBC1D3 audit) — conventions later runs must follow
Found while auditing `bench/LAYER_ORDER_NPIP_TBC1D3.md` (its "Audit notes" section). Each item fixes a place where two
readings gave different numbers.
- **Layer universe.** P's universe = every gene of U with a §6ko protein, including coding genes P does not place with a
  member (they are their own group). D's universe = every gene of U that is a node of a catalog. The group tables must
  carry these rows. Outside its universe a gene is "not in layer", never a singleton.
- **Vacuous containment.** c(X ⊇ Y) over 0 pairs is not 1. A tournament comparison with 0 pairs on either side is NA, not
  "above". Report the group-level n (groups with ≥ 2 genes) next to every pair count. Call an order "transitive" only when
  all layer pairs are scored on one common gene set.
- **JOIN (amended).** Vertices are the finer layer's WHOLE groups (for a catalog: the whole catalog cluster, not the part
  inside U). U is re-closed after enforcement. Score on the coarse layer's universe and on the enforced layer's own
  universe, and say which.
- **REFINE (amended).** Where a gene of the finer layer is outside the coarser layer's universe, what REFINE does with it is
  a free choice. Report every variant: *attach* (P: to the refined part with the largest summed edge weight; C: stays
  with the largest part of its clade), *single* (its own group), and for C *together* (a clade's outside genes form one
  part). Any statement that compares REFINE with JOIN must hold under every variant, or be stated per variant.
- **EXPR.** The primary count is *any-overlap* (reads ≥ 3, primary `-F 2308`, a CIGAR block overlapping an exon), as
  defined above. "Unique" read counts are secondary, and they depend on which overlapping records exist (readthrough
  records take reads from the copies they span). Always report a threshold sweep. EXPR(L) ⊆ L holds by construction; it is
  a code check, not a result.
- **T2 needs a precondition.** For layer-specific edge sets E_L and E_M, T2 holds when, inside every L group, E_M ⊆ E_L
  (every M edge between two genes of one L group is also an L edge). Without it T2 fails. Counterexample: L = D, M = a
  REFINE of P inside D; genes a, b, c in one D group; P edge a–b; the only D path a–c–b; c unexpressed. Then EXPR(M) =
  {a, b}, but in EXPR(L) a and b are unconnected singletons and are dropped. With co-membership edges the precondition
  always holds, so T2 is guaranteed there and is not evidence.
- **C layer.** The ordered subfamily layer is the §0★★ clause-5 split system (supported SH-aLRT > 75 splits of either
  reference tree, pairwise compatible; cluster = smaller side; partitions: maximal and minimal clusters). The
  literature-anchored clades are a CIRCULAR reference column. C lies inside one family by construction (Non-goals), so
  "P/D above C" is expected. The only empirical content is whether P or D splits a clade pair.
- **D in the NPIP/TBC1D3 study** is the §6kl RefSeq E1 catalog (0.70/0.30 edges, MCL), not clause 2/4. That catalog failed
  its HGNC guard on chr16/19/20, and GENCODE E1 merges PKD1 into NPIP. Headline rows must say "E1".

## Substrates and scoring protocol

- Genome/annotation: T2T-CHM13 v2.0; RefSeq `chm13v2.0_RefSeq_full.gff.gz` (never `HSA_genomic.gff`) for D/S; the
  §6ko annotation convention for P.
- Truth for S: Soto families — `winloci_data/soto_replication/soto_gene_to_families.tsv` joined to
  `gene_v2_coords.bed` (2,332 of 2,334 genes placed; 149 genes listed in >1 family are "ambiguous").
- **Dev:** chr1, chr15, chr17 (all previously used as dev for Soto/shared-exon work).
- **Held out:** chr8, chr9, chr10, chr11 — never used as a systematic dev or held-out substrate for Soto, protein,
  guided E1 truth or shared-exon work (audit in the inventory). Unambiguous Soto families with ≥ 2 members /
  member genes: chr8 21/129, chr9 39/182, chr10 22/122, chr11 14/56. chr12 (proposed first) was swapped out: only
  4 families / 17 genes.
- **Pooled, cross-chromosome-aware scoring:** each split (dev, held-out) is scored as ONE pooled gene universe = the
  split's chromosomes. Same-family pairs count when both genes are in the split's universe, including pairs on
  two different chromosomes of the same split; pairs with a gene outside the split are ignored. A Soto family is
  restricted to its in-split members and scored only if ≥ 2 remain. Ambiguous genes are excluded from truth pairs
  and from bipartite matching (reported separately). ~20% of Soto families are interchromosomal, so this matters.
- Metrics: pairwise precision/recall; bipartite matching F (one-to-one family matching by Jaccard, as in §6ks);
  containment c(X ⊇ Y) = |pairs(Y) ∩ pairs(X)| / |pairs(Y)| on the genes present in both layers' universes (pairs
  involving P restricted to coding genes); group-level nesting = fraction of Y groups (≥ 2 genes) inside a single
  X group.

## Phases

### Phase A — development (chr1/15/17; no pre-registration needed, disclosed as development)
1. Build P, D, S1, S2, S3 on the dev universe.
2. Select the SD layer: highest bipartite F vs Soto; ties → higher pairwise precision → simpler rule
   (S1 < S2 < S3). Record all three.
3. Descriptive containment matrix for {P, D, S_sel} and the dev order (tournament: X above Y iff
   c(X ⊇ Y) > c(Y ⊇ X)).
4. Freeze S_sel, the order found on dev, and all rules; then write
   `docs/PREREG_family_layer_order_2026-09-16.md` BEFORE any held-out file is generated.

### Phase B — held out (chr8/9/10/11; pre-registered)
- **E-S (SD layer derives Soto):** S_sel's held-out bipartite F vs Soto ≥ S1's held-out bipartite F (the selected
  rule is at least as good as the shipped rule), and S_sel's held-out pairwise precision ≥ S1's. If S_sel = S1 this
  endpoint is trivially met and the claim is scored on the absolute numbers only, stated as such.
- **E-O (order replicates):** the held-out tournament on {P, D, S_sel} is transitive AND equals the dev order.
- **E-C (enforcement is cheap):** enforcing nesting in the held-out order (JOIN for layers above D that disagree
  both ways, REFINE for layers below D) lowers no layer's bipartite F vs its own truth by more than 0.02
  (P vs HGNC gene groups as in §6ko; D vs the guided annotation truth used in §6kl/§6km; S_sel vs Soto).
- **Verdict:** SUPPORTED if E-S, E-O and E-C hold; PARTIAL if E-O holds and exactly one of E-S/E-C fails; REFUTED
  otherwise.
- **Descriptive (no endpoint):** EXPR commutation per layer — agreement (pairwise F) between EXPR(L) and L rebuilt on
  expressed genes only — using per-gene read counts from `human_testis.t2t.bam` (primary, `-F 2308`, reads ≥ 3
  overlapping the gene's exons) on the held-out chromosomes. The per-gene count table does not exist yet and is
  built in Phase B.
- No rule, threshold, candidate, chromosome or scoring convention changes after any held-out number is seen.

## Outputs and code
- Python bench scripts under `bench/layer_order/` (layer builders reuse existing scripts: `bench/protein_families.py`,
  `mcl_families`, the §6ie–§6iy Soto-replication scripts minus famCN); results under
  `/mnt/linuxdisk/home/juanfraitu/layer_order/{dev,heldout}/`.
- Definition doc: `docs/seeded_family_definition.md` §0★★ gains the layer stack, the three operators, T1–T3 with
  proofs, and the corrected clause 6 (EXPR instead of the component-only nesting sentence).
- Reports: `bench/LAYER_ORDER.md` (dev + held-out, terse tables); ledger entry; negative-results register row if
  PARTIAL/REFUTED.
- No Rust changes in this project (nothing in the pipeline consumes multiple layers today).

## Non-goals
- famCN/parCN; the subfamily-clade layer (inside its family by construction; not ordered here); making clause-4
  leader order read-independent; de novo ⊆ guided nesting; gorilla (no Soto truth).

## Known limits (declared now)
- Soto families are one group's definition, built on SD98 + copy number; excluding famCN caps how exactly any
  sequence-only rule can match them (prior full replication with famCN: ARI 0.69).
- P is defined on coding genes only; containment involving P ignores non-coding copies.
- Held-out chr11 is small (14 families); results are pooled across the four chromosomes, not claimed per chromosome.
- The EXPR analysis uses one tissue (testis); expression-dependent statements are tissue-specific.
