# Evidence-agnostic multi-copy gene families, transferable to poorly annotated relatives — design

Date: 2026-09-14 · Status: approved in conversation (sections 1–2 walked through; user: "it all looks good, please proceed")

## Goal

Define multi-copy gene families from **whatever information is available**, and use a well-annotated relative's
information for species that are poorly annotated. v1 supports two starting points:

- genome only;
- genome + partial annotation. This has two arms: a native subsample (controlled) and an annotation lifted from human
  (realistic, cross-species).

v1 is validated on T2T genomes only.

Success is measured two ways:

1. against the target species' own hidden native annotation (sensitivity, precision, bipartite matching);
2. by consistency with the human families the transferred information came from.

## Background (why this shape)

- **§6kl–§6km: gene-model families are not reproducible across annotations.** Even two expert annotations of one genome agree
  only to F ~0.8. The residual is gene-model content.
- **§6ko: protein-space families are reproducible.** Built from an independent annotation, they reproduce at F 0.977 (fresh
  hold-out chr1/2/3: 0.935 / 0.999 / 0.977).
- **§6ko: protein and DNA families are not nested** (52% / 25% pair overlap). Protein space links deep homology (OR, SLC)
  and misses protein-divergent recent duplicates (GOLGA6/8: nt 0.83, aa 0.27). DNA space sees pseudogenes and
  co-duplicated non-coding sequence.
- **§6jj–§6jk: genome-only segments carry useful signal.** Duplicated-segment atoms from SD alignments recover guided
  families genome-wide at gene-level recall ~0.57.

## 1. Architecture: one locus graph, three layers, one construction pattern

Each layer = nodes + a typed homology edge + MCL (I = 2.8, prune 1e-9). The rules never change with the input; a thinner
input only activates fewer layers.

| Layer | Nodes | Edge | Requires |
|---|---|---|---|
| **G** genomic segment | duplicated segments: SD pairs (BISER, same caller for every species) cut at all SD boundaries (`bench/dna_sd_atoms.py` atom logic) | exact-CIGAR segment alignment, identity >= 0.70, >= 300 bp, two-sided coverage >= 0.30; soft-masked repeats | genome |
| **P** protein | one protein per coding gene (longest CDS) from any model source, plus copies placed by miniprot | blastp e <= 1e-5, HSP union covers >= 0.30 of the longer protein (§6ko, r2: no pseudogene biotypes, no V(D)J segments) | any annotation with CDS |
| **E** exon / gene body | gene models: annotated, lifted, or projected onto new copies | E1: exon-to-exon, identity >= 0.70, >= 300 bp, cov_longer >= 0.30 (`mcl_families --min-exonic-bp 1`) | any gene models |

**Cross-layer map (the reconciliation).** Per locus it lists:
- P family (if coding);
- E family;
- the G segment families it overlaps;
- source (native / lifted / projected / pseudogene-attached);
- human source gene, where known.

Pseudogenes attach to P families by miniprot, as in AN-3. Layer disagreement is reported with the §6ko cause types:
deep protein-only homology, protein-divergent duplicates, and co-duplicated non-coding sequence.

**What each input activates:**

| Input | Layers |
|---|---|
| Genome only | G (no gene labels) |
| Genome + partial annotation (native subsample or human lift) | G + P + E; gene labels = annotation + discovered copies |
| Genome + full native annotation | all three; this is the truth (arm T) |

**Species:** gorilla = development; **orangutan = hold-out**; chimpanzee = report only. All are NCBI T2T v2.0 assemblies
with RefSeq annotations (GGO/PTR/PPY `_genomic.gff` in `winloci_data`) and soft-masked FASTAs (~50% lowercase). Human
(CHM13 v2.0 + RefSeq) is the source of lifted information.

## 2. Arms

**Substrate per species.** The chromosomes orthologous to human chr7, chr15, chr16 and chr17 (SD-rich; NPIP, GOLGA, TBC1D3
and KRT), about 9k genes. Orthology comes from the whole-genome Liftoff placement of human genes, not from chromosome names;
this handles the gorilla t(5;17). All homology is computed within the substrate.

**Shared preparation.**
- BISER on the substrate chromosomes (batched with `--resume`).
- Gorilla BISER vs SEDEF agreement reported.

**Arm T (truth; hidden from every other arm).** Native RefSeq through all three layers:
- E_T: E1 families;
- P_T: protein families, r2;
- G segments labelled with native genes.

**Arm 1 (genome only).** BISER pairs → segments → segment graph → MCL → segment families. No gene labels are produced. For
scoring only, each truth gene takes the segment family of the segments covering most of its exons.

**Arm 2a (native subsample, controlled).**
1. Keep 50% of native genes (seed 1).
2. P: miniprot places the kept proteins on the substrate. An alignment covering >= 0.30 of a protein and not overlapping a
   kept gene is a new coding copy; its predicted CDS is translated; blastp + MCL.
3. E: `guided_min.py` gene-body projection of the kept genes → E1 on kept + projected models.
4. G: as arm 1, with gene labels.

**Arm 2b (lifted from human, cross-species).** Liftoff places human RefSeq on the target genome, in per-human-chromosome
batches. Discovery and the three layers then run exactly as 2a, starting from the lifted genes. Lifted genes keep their human
source ID; discovered copies inherit the source ID of the protein or gene that found them.

**Outputs of every arm:** `families.<layer>.tsv` and `layers_map.tsv`. A layer without inputs is written as "not
available", never as empty families.

## 3. Evaluation

**Truth isolation.** Each arm receives an explicit input manifest: genome, BISER calls, allowed annotation file. The
orchestrator refuses to run an arm whose manifest names the native GFF, except arm T and arm 2a's subsample file (written by
arm T's prep).

**Metrics** (the §6kg/AG and §6ko conventions):
- **E layer and arm 1:** truth loci = members of E_T clusters with >= 2 members. Each is assigned the arm's family by best
  span overlap. Pairwise sensitivity/precision and bipartite R/P/F are computed among truth loci
  (`bench/rna_truth.py score`).
- **P layer:** truth genes = P_T members. Each is assigned the family of the arm's gene with the greatest CDS-base overlap
  (`bench/protein_families.py score`).
- **Cross-species consistency (arm 2b):**
  - human reference families = the same layers built on human chr7/15/16/17 from RefSeq;
  - each target family member with a human source ID is labelled with that gene's human family;
  - reported: bipartite R/P/F between target families and human families on labelled members, and the per-family
    copy-number difference (target − human).
- Species are never pooled.

**Pre-registration.** Addendum AO in `docs/PREREG_core_definition_2026-09-12.md`, written after gorilla development and
before any orangutan number exists. Proposed fixed readings (orangutan):

| id | question | reading |
|---|---|---|
| AO-1 | human-lifted protein layer (arm 2b P vs P_T) | goal met iff pair sens >= 0.90, pair prec >= 0.90, bipartite F >= 0.90 |
| AO-2 | native-subsample protein and exon layers (arm 2a P vs P_T, E vs E_T) | goal bars reported; SUPPORTED iff P meets them |
| AO-3 | genome-only segment families at gene level (arm 1 vs E_T) | reported against the bars and against §6jk's gorilla 0.567 recall |
| AO-4 | cross-species consistency (arm 2b vs human) | reported; no bar in v1 |

Chimpanzee is scored with the same readings, reported only.

## 4. Components

New (bench/, Python, same conventions as the AI–AN tools):

| file | responsibility |
|---|---|
| `ape_substrate.py` | Liftoff human → target (batched); orthologous substrate chromosomes; substrate FASTA/GFF subsets; the 50% native subsample; per-arm input manifests |
| `sd_segments.py` | BISER run (batched, `--resume`); BISER → `dna_sd_atoms.py` format adapter (validated on toy rows and on gorilla vs SEDEF); segment families (MCL) |
| `protein_projection.py` | miniprot placement of a proteome on a genome → new coding copies (node table + `.cds.tsv`) and pseudogene attachments |
| `layers_map.py` | cross-layer map per locus; cause classes of layer disagreement; cross-species consistency metrics |
| `evidence_arms.py` | orchestrator: arm → manifest → layers → outputs; truth-isolation guard |

Reused: `annotation_nodes.py`, `protein_families.py`, `node_graph_mcl.py`, `guided_min.py`, `mcl_port.py`,
`rna_truth.py`, `adjudicated_truth.translate`, `dna_sd_atoms.py`, `lit/batch_align.sh`.

**Testing.**
- Unit tests on toy inputs: BISER row parsing and CIGAR blocks, segment cutting, new-copy selection from miniprot rows,
  the truth-isolation guard, consistency metric.
- Reproduction checks:
  - `protein_families.py` on human chr1/2/3 reproduces §6ko (0.935 / 0.999 / 0.977);
  - E1 on the gorilla substrate matches the corresponding `gw_units_v3` clusters, with the difference reported;
  - `dna_sd_atoms.py` on gorilla SEDEF reproduces §6jk atoms.
- BISER vs SEDEF on gorilla is reported before BISER feeds any arm.

**Error handling.**
- Missing inputs → the layer is "not available".
- Liftoff unmapped and partially mapped genes are counted and reported.
- miniprot frameshift/stop-containing placements go to the pseudogene attachment path, not the proteome.
- BISER runs are resumable.
- Every long step is a foreground batch.

**Operational constraints (standing):**
- WSL2: foreground batches ≤ 10 min, one heavy run at a time, kill by PID only.
- Large outputs go under `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/`.
- Never pool species.
- Pre-register before looking at orangutan.
- Default flips remain the user's decision.

## Out of scope for v1

- Non-T2T / draft genomes (collapsed SDs).
- RNA reads as an input tier.
- A related-species proteome without a genome annotation (can be added as another P-layer source later).
- De novo ORF prediction in genome-only mode.
