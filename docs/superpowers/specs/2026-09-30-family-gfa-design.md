# family_gfa: a graph view of O1 families (design, 2026-09-30)

## Goal

Give every O1 family a graph representation that O1 hands over: shared exon blocks as shared nodes, copy-specific exon
sequence as private nodes, each copy a path, and fusions between families as joint graphs. It must be viewable in
Bandage and reusable for other species. It is a **re-encoding of what O1 already computed** (its own all-vs-all locus
alignments), not new inference.

**Order of use:** render the flagship families for the advisor first (NPIP, the PKD1P6-NPIPP1 fusion, one small clean
family, and TBC1D3 if a chr17 catalog is on disk), then run it on whole catalogs.

## What this is not (register constraints)

- **Not a family definer or detector.** Using a variation graph to define or detect families is circular (register
  L472, L1086). The graph is built from O1's output and changes nothing about it.
- **Wording:** "a graph view of an O1 family". Never "a family is one VG" (`docs/O1_O2_COMPOSITION.md` 25-31).
- **Not borrow strength.** Assembly-side borrowing gave 0 beneficiaries and fabricated structure (register L647-648);
  no new test is planned.
- **Not base-level.** Small differences inside a shared block become identity tags, not bubbles.
- **Not a replacement** for the existing opt-in O2 graphs (`copy_assign --phase`: `.phase.gfa`, `.exon.gfa`). Those
  merge exons only by genomic overlap (register L757), so dispersed copies never share a node; this one merges by
  homology.

## Where it lives

A new standalone binary, `src/bin/family_gfa.rs`, plus a `[[bin]]` entry in `Cargo.toml` (the repo sets
`autobins = false`). No existing source file changes, so `mcl_families` and the default pipeline stay byte-identical,
and nothing collides with work in progress on `mcl_families.rs`. Integration tests go in `tests/family_gfa.rs`. A
pipeline stage is out of scope here.

## Inputs

From one families-stage run (`mcl_families`, prefix `<out>.fam`):

| file | used for | format facts (checked on the 2026-09-25 chr16 catalog) |
|---|---|---|
| `.clusters.tsv` | family membership | `cluster_id size density frac_in corroborated chrom start end`; one row per member locus |
| `.loci.gff3` | exons per locus | `gene` span + `exon` children, 1-based inclusive |
| `.loci.fa` | genomic sequence per locus | header `chrom:start-end`, same span as the gene |
| `.loci.paf` | O1's all-vs-all alignments | minimap2 `-x asm20 -c -X -N 50 -p 0.1`; `cg:Z` CIGAR uses `M` (match or mismatch), `I`, `D` |
| `.relations.tsv` (optional) | fusions | `mcl_families --emit-relations`; `unit_ids`, `unit_loci`, `unit_families` in transcript order |

A locus is keyed by `(chrom, start, end)`; the key is identical across the first three files (353/353 clustered loci
checked). The binary stops with an error, not a guess, when a member's key is missing from any file.

## Command

```
family_gfa --prefix <out>.fam --outdir DIR [--families MCL1,MCL4] [--relations <out>.fam.relations.tsv]
```

Without `--families`, every family in `.clusters.tsv` is written.

## How one family's graph is built

1. **Vertices.** Every exonic base of every member locus (exons from `.loci.gff3`, bases from `.loci.fa`). Introns and
   flanks are not in the graph.
2. **Joins.** For every PAF record whose query and target are both members and are different loci, walk the CIGAR; each
   `M` column where both bases are exonic joins the two vertices. A `-` strand record joins with reversed orientation
   (union-find with an orientation parity bit). All records O1 used are included, as O1 includes them.
3. **Columns.** The joined sets (the transitive closure, as in seqwish) are the graph's columns. A join that would put a
   vertex with itself in opposite orientation is refused and counted (`inversion_conflicts`).
4. **Paths.** Each copy's path is its exonic bases in transcript order (5' to 3': reverse genomic order for `-` strand
   loci), so the path spells the copy's spliced transcript. Step orientations follow from the joins.
5. **Nodes.** A node is a maximal run of consecutive columns that every path through them traverses together (unitig
   compaction). Nodes are also cut wherever any traversing copy has an exon boundary, so nodes line up with exons.
6. **Node sequence and tags.** The sequence is taken from a representative copy (the member with the smallest locus key
   that passes through the node). Tags: `LN:i` length, `CN:i` number of distinct copies through it, `ID:f` mean and
   `IM:f` minimum identity of the traversals against the representative.
7. **Loops.** A column that holds two positions of the same copy (an internal repeat) produces a loop. Loops are kept and
   counted (`loops`), never hidden.
8. **Same-DNA joins.** Loci of one family can overlap in the genome (a fused locus and one of its units). A join between
   two vertices at the same chromosome position is the same DNA, not paralogy. Such joins are kept but counted separately
   (`same_dna_bp`), and nodes made only of them get the tag `SD:i:1` so figures can tell them apart.

## Outputs

Per family, in `DIR/`:

- `<family>.gfa` (GFA 1.0). The first line after `H VN:Z:1.0` is a comment: `# graph view of O1 family <id>, built by
  family_gfa from O1's own alignments; not used to define families`.
- `<family>.colours.csv` for Bandage: nodes shaded light blue (2 copies) to dark blue (all copies); private nodes (1 copy)
  orange.
- `<family>.paths.tsv`: `path_name gene_id locus strand n_exon exonic_bp`.

For the run:

- `summary.tsv`: per family `n_copies n_nodes n_shared_nodes shared_exonic_frac loops inversion_conflicts same_dna_bp`.
- With `--relations`: `fusions.tsv` (`relation_id transcript families units joint_gfa reason`) and, for every relation
  whose units fall in two or more families, `fusion_<relation_id>.gfa`:
  - the involved families' graphs, node IDs prefixed by family ID;
  - one extra path, the fused transcript: its units' paths joined in transcript order by fusion links (`L` lines tagged
    `FU:Z:<relation_id>`). A unit in no family enters as its own private exon nodes, so the fused path is complete;
  - fusion-link end nodes coloured red in its `.colours.csv`.

  Relations with all units in one family, or none, get a `reason` and no joint graph.

## Correctness checks

Test-first. Synthetic unit and integration tests:

| case | expected |
|---|---|
| two identical copies | one shared path, every node `CN:i:2`, `ID:f:1` |
| a SNP | the shared node carries `ID:f` < 1; no bubble |
| an insertion or an extra exon in one copy | a private node on that copy's path (the e3a/e3b case) |
| a reverse-strand copy | joins with reversed orientation; its path, read with step orientations, spells its transcript |
| an alignment only in introns | no join |
| an internal repeat | a loop, counted in `loops` |
| overlapping loci (same DNA) | joined, counted in `same_dna_bp`, node tagged `SD:i:1` |
| a missing locus key | an error naming the key |

Invariants asserted at run time on every family, and in the real-data test:

- each path's total length equals its copy's exonic length;
- every consecutive pair of path steps has an `L` line;
- a path, read with its step orientations, spells its copy's spliced transcript except inside shared nodes, and the
  number of differences matches the identity tags;
- `CN:i` equals the number of distinct paths through the node;
- a rerun is byte-identical.

GFA validity: every output parses with `vg view -F` (vg is installed). Bandage is not installed here; figures are
opened on the user's machine.

## Real-data smoke test

1. Rerun the chr16 families stage with current defaults plus `--emit-relations --emit-units`. The 2026-09-25 catalog
   predates the 2026-09-29 default changes, so it is used only for development.
2. Run `family_gfa` on NPIP, the PKD1P6-NPIPP1 relation, and one small clean family. Add TBC1D3 if a chr17 catalog
   exists.
3. Report per figure: copies, nodes, shared exonic fraction, loops, same-DNA bp, and the invariant results.

Machine rules: builds and runs go through `tools/rlock.sh` (a build is a heavy job), with
`CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target` and output under `/mnt/linuxdisk`.

## Out of scope

The borrow-strength test; base-level graphs; any change to family definition, O2, or pipeline defaults; a pipeline
stage and the genome-wide run (after the flagships, and after concurrent work on `mcl_families` merges); the human vs
gorilla comparison (the same exporter, run later on the gorilla catalog).

## Amendment 1 (2026-09-30): duplicon annotation

**Why.** Segmental duplications are mosaics of duplicons (ancestral duplication units). A pre-registered test on the Soto
branch (`docs/PREREG_soto_cn_vs_duplicon_2026-09-30.md`, KEY=cnduplicon) found that Soto's copy-number cut follows duplicon
boundaries (verdict FOLLOWS on the held-out and dev halves; copy-number gap vs duplicon dissimilarity, Spearman +0.58).
Fusion genes such as PKD1P6-NPIPP1 span two duplicons. A family graph should therefore show which duplicons each node sits
on, and a fusion's joint graph should show the duplicon boundary the fused transcript crosses.

**Input.** `--duplicons BED`, optional. Without it every output is byte-identical to the unamended spec. Column 4 is the
duplicon ID; overlapping records are allowed and every overlapping duplicon counts. The intended file is the DupMasker painting
of Vollger et al. 2022 (Zenodo 4726156, `data.tar` member `chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`, the source Soto
cite), fetched on 2026-09-30 to `/mnt/linuxdisk/home/juanfraitu/winloci_data/duplicons/`.

**Assembly.** The BED must be in the catalog's coordinates. The annotation is CHM13 v1.0. On 2026-09-30 gene coordinates
agreed between the CAT v4 annotation and RefSeq on CHM13 v2.0 on chr1, chr4, chr16 and chr17 (AMY1A, DUX4, NOTCH2NLA, NPIPA1,
SRGAP2C, TBC1D3, within annotation-boundary differences). The acrocentric chromosomes (13, 14, 15, 21, 22) and chrY were not
checked, and a run on them must check first. The GFA comment line records the BED's file name and sha1.

**Annotation (representation only).** Duplicons label the graph; they never change nodes, links, paths or families. Duplicon
rules were refuted as family definers (register 677, 743, 817, 997).
- Each node gets a tag `DU:Z:` listing the duplicon IDs that cover its representative copy's bases, with the base count, for
  example `DU:Z:SD9443:1179;SD9456:714`. The list is sorted by base count (descending), then by ID. A node with no duplicon
  gets no tag.
- `<family>.paths.tsv` gains a column `duplicons`: the copy's exonic duplicon composition, in the same `id:bp` form.
- `summary.tsv` gains `n_duplicons` (distinct duplicon IDs over the family's nodes) and `dominant_duplicons` (the three
  covering the most bases).
- A second colour file, `<family>.duplicons.csv`, colours each node by its dominant duplicon. The 8 most frequent duplicons in
  the family get the categorical palette in order of frequency, and the rest are grey. Bandage can then switch between the
  copy-sharing colours and the duplicon colours.
- Each fusion link in a joint graph gets `DX:Z:<left node's dominant duplicon>><right node's dominant duplicon>`.

**Tests (added to the plan as Task 7 before execution resumes).**
- A synthetic duplicon BED over the integration fixture gives the expected `DU:Z:` values, `paths.tsv` compositions and
  colour file.
- Without `--duplicons`, outputs are byte-identical to a run of the unamended binary.
- A fusion link carries the expected `DX:Z:` tag.
- A BED chromosome absent from the catalog is ignored silently, while a duplicon overlapping no node is simply not reported.
