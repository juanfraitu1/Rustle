# Family definition against Ensembl Compara's recent paralogues — design (2026-10-18)

**Status:** design approved in conversation on 2026-10-18 (sections 1-3 and the closing round). This written spec awaits the
user's review before the implementation plan is written.

## Goal

Choose the family definition, or definitions, that decide which genes form one multi-copy gene family, by their agreement
with a published, method-independent truth: **Ensembl Compara's recent within-species paralogues.** Then confirm the chosen
definition end to end with one pre-registered round of the integrated O1 -> O2 -> O3.

**The user's decisions (2026-10-18):**
- Score all candidate definitions on the same truth: the simulation's name-free rule and the thesis's guided O1.
- The truth is Ensembl Compara, "specific with definitions". Soto 2025 is not multi-copy gene families; it is reported as a
  special case only.
- Compara is used as recent-paralogue pairs, by the age of the duplication node, with the age ladder reported.
- A closing round runs the full O1 -> O2 -> O3 on families drawn by the settled definition.

## Scope

- Human only: CHM13 v2.0 with the CAT/Liftoff v2.0 annotation. No gorilla numbers, and never pooled with them.
- Protein-coding genes only, because Compara holds no pseudogenes. Non-coding members (round 19's rule) are not judged here.
- No real reads. The closing round is a simulation round, like rounds 15-19.

## Truth

- **Source:** the Ensembl REST service, release 116, `/homology/id/human/{ENSG}?type=paralogues;format=condensed`. For each gene
  it returns every paralogue with its `type` and the `taxonomy_level` of the duplication node.
  - BioMart no longer answers.
  - The FTP TSV dumps lack the duplication level.
- **A true pair:** two universe genes listed as `within_species_paralog` whose duplication level is recent.
  - **Primary:** Catarrhini or younger, the set {Homo sapiens, Homininae, Hominidae, Hominoidea, Catarrhini}.
  - **Ladder:** also scored with Simiiformes, then Haplorrhini, then Primates added, cumulatively. Reported only.
  - A pair counts if either gene lists the other. Pairs listed in one direction only are counted and reported.
- **True families:** the connected components of the true pairs. A universe gene with no true pair is a singleton.

## Universe and mapping

- **Genes:** CAT's protein-coding genes (20,048 in the full GFF). Each carries its `source_gene` (an Ensembl ID) and its
  `collapsed_gene_ids`; 727 genes absorbed two or more Ensembl genes onto one CHM13 locus.
- **In the universe:** a CAT gene with at least one Ensembl ID that Compara 116 returns. The genes left out are counted.
- **Pairs:** a Compara pair maps to a CAT pair through those IDs. A pair that lands on a single CAT gene is dropped and
  counted.
- **No span cap.** The simulation's 60 kb limit is a cost limit, not part of a definition.
  - A one-chromosome cost pilot of `mcl_families` comes first.
  - If the pilot shows the genome-wide run unaffordable, every candidate is scored on the capped universe instead. That choice
    is made and recorded before any candidate is scored.

## Candidates (fixed in advance; none tuned on Compara)

Each candidate is a partition of the universe. A gene in no family is a singleton.

| | candidate | rule |
|---|---|---|
| 1 | NF-S | name-free rule, coding genes only (rounds 17-18): canonical transcript identity >= 0.90 over >= 0.8 of it on the other locus (+-1 kb), either direction, AND one gene's exon chain carried by the other's genome |
| 2 | NF-seq | candidate 1 without the structure test |
| 3 | NF-P | the name-free rule with non-coding members (`--paralogs`, round 19), scored on its coding genes; non-coding members can bridge coding genes |
| 4 | O1-guided | the thesis's guided O1: `mcl_families --from-gtf` on a CAT coding-gene GTF with `--min-exonic-bp 1 --min-shared-exon-frac 0.60`, shipped defaults otherwise |
| 5 | Stems | the annotation-name baseline: genes grouped by name stem (the name less a final number and optional letter), groups of >= 2 |

**Secondary, never used to choose:** NF-S at identity 0.80, 0.85, 0.95 and 0.98, to read off which identity cut the external
truth supports (the §6p3 lesson: estimate the cut from an independent source).

## Scoring

The metrics mirror the native `family_score`, which is the user's standing reporting rule:
- **Pairwise:** over unordered universe gene pairs. A predicted pair is two genes in one family. Report sensitivity, precision
  and F.
- **Bipartite:** take true families and predicted families of size >= 2, with predicted families intersected with the universe.
  - Find the one-to-one matching that maximises total gene overlap (scipy `linear_sum_assignment` on the overlap matrix, the
    tie-breaking `family_score` ports).
  - **Sensitivity** = matched overlap / genes in true families.
  - **Precision** = matched overlap / genes in the matched predicted families (those with overlap > 0).
  - **F** = their harmonic mean. This family F is the primary metric.
- **Collapse:** for each predicted family, the true families (size >= 2) with a gene in it, less one; summed over predicted
  families.
- **Halves:**
  - True families are split by the SHA-1 of each family's smallest Ensembl ID: even to development, odd to held-out.
  - Singleton genes are split by their own ID's hash.
  - Each half is scored on its own genes: pairs inside the half, predicted families restricted to the half.
- **Uncertainty:** the family F's 95% CI in each half comes from a bootstrap over the overlap graph's components.
  - A component is a set of true and predicted families linked by shared genes. The matching is computed once per
    component.
  - Components are resampled with replacement, 1,000 times with seed 0, and F is recomputed from the resampled sums of
    matched overlap, true-family genes and matched predicted-family genes.
  - Resampling single true families instead would let a duplicated family compete with itself for one predicted family.
- **Soto 2025:** each candidate's concordance through `family_score --soto`, reported as a special case.

## Decision rule (pre-registered before any candidate is scored)

1. The **development winner** is the candidate with the highest family F on the development half.
2. It is **settled** if its held-out family F lies within the 95% bootstrap CI of the best held-out family F.
3. Every candidate inside that CI is reported as **tied**. That is how "one or more" definitions are settled.
4. If the development winner falls outside that CI, **no definition is settled**, and the report says so.

## Closing round (end to end)

- **The pool:** the simulation's family pool rebuilt with the settled definition. If several tie, the development winner is
  used.
  - Non-coding members join by round 19's rule, whatever the coding rule, because Compara cannot judge them and the user
    decided they are copies.
  - If the winner is O1-guided, an adapter writes its families in the simulation's format (copies with their transcripts).
- **Spent:** every earlier simulated set.
- **The draw:** the fixed rule of rounds 17-19, with seed 20.
- **The systems:** integrated v14 unchanged; the current method unchanged.
- **The worlds:** profiles v4, seeds 1001-1003, six world types under each profile.
- **The bars:** B1 and B2 in every world, with every cell required to pass, and C1 on the mainstay set, as rounds 16-19.
- **Timing:** pre-registered in its own document before any of its worlds is built.

## Implementation units

All units sit in `bench/famdef/`, on branch `famdef-compara` off `o3-fast`, in its own worktree.

- **`compara.py`:** a resumable REST fetch of every universe gene's paralogues.
  - It stays at or below 10 requests a second.
  - It retries with backoff on HTTP 429 and 5xx.
  - It writes a JSON-lines cache keyed by Ensembl ID, and the run records the release.
- **`universe.py`:** CAT coding genes with their Ensembl and collapsed IDs, from the full GFF. It counts genes missing from
  Compara.
- **`truth.py`:** true pairs at an age level, the true families (components) and the hash split into halves.
- **`candidates.py`:** builds each candidate's partition.
  - `seqfamilies.py` gains `--no-structure`, `--min-ident X` and `--max-span N` (0 = no cap). Its output is byte-identical when
    the flags are off.
  - It calls `mcl_families` and maps its clusters to CAT gene IDs.
  - It builds the stem groups.
- **`score.py`:** the pairwise, bipartite and collapse metrics, per half, with bootstrap CIs.
- **Tests, written first:**
  - the age-ladder filter;
  - ID mapping with collapsed genes;
  - the hash split;
  - the scorer on toy partitions, cross-checked once against `family_score` on a fixture;
  - the new `seqfamilies` flags, byte-identical when off.

## Order of work

1. This spec, then the implementation plan.
2. Build the units; fetch Compara; run the `mcl_families` cost pilot.
3. Pre-register (`docs/PREREG_famdef_compara_2026-10-18.md`): the candidates, the universe (capped or not), the metric and the
   decision rule, frozen before any candidate is scored.
4. Score once; report; register row; memory.
5. Pre-register, run and report the closing round.

## Risks and limits

- **Compara is protein-level and tree-derived.** An August comparison against all its paralogues bounded agreement at
  precision 0.35 (register row 316). The age restriction targets that, and the report shows precision by age level so the
  effect is visible.
- **Ensembl ID drift:** CAT's source GENCODE release against Ensembl 116. Unmatched genes are counted, not silently dropped.
- **`mcl_families` genome-wide cost is unknown;** the pilot decides the universe before any score.
- **Non-coding members are not judged** by this truth.
- **Nothing here runs on real reads.**

## Deliverables

- the pre-registration and its results;
- a report: per-candidate tables for both halves, the age ladder, the identity ladder and Soto concordance;
- register rows and a memory entry;
- the closing round's pre-registration and results.
