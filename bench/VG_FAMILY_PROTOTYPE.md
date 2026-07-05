# VG-native, minimizer-free multi-copy family definition prototype

**Goal** — test whether O1 (family definition) and O2 (copy assignment) can be unified around a
variation-graph architecture instead of the shipped minimizer-based repeat catalog and k-mer
minimizer gates.  The prototype builds a pan-transcriptome exon-splice VG directly from the
BAM-derived de-novo skeletons, derives repeat hubs from VG topology, and extracts families from the
graph.

**Code**
- `bench/vg_family_prototype.py` — build the catalog.
- `bench/vg_family_prototype_eval.py` — evaluate against the same P/R truths used for the shipped
  `family_rna_refine.tsv` catalog.

Run the default O1+VG mode:

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py --threads 4
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_eval.py
```

---

## 1. What the prototype does

1. **Load de-novo skeletons** (`winloci_scratch/denovo_skeletons.tsv`) and metadata
   (`winloci_scratch/denovo_transcripts.meta.tsv`).  The skeletons are produced by
   `bench/twopass_denovo_gw_pass1.py` by grouping **primary BAM reads** (`GGO_mm.bam`) by their
   exact intron chain, so the loci are determined directly from the read alignments.
2. **Extract oriented exon paths** from `GGO.fasta` for all 110,237 loci (~1.26 M exon
   occurrences).
3. **Build a pan-exon VG.**  Each distinct canonical exon sequence is a node; consecutive exons in
   a locus are connected by a splice edge; each locus is a path through the graph.
4. **Cluster exons** (mode-dependent):
   - `exact` — no clustering; every canonical sequence is its own node (fastest, most stringent).
   - `cdhit` / `vsearch` — near-exact clustering at 95 % identity.
   - `seeded` — pure-Python k-mer seed + `edlib` verification (functional but slow on the full set).
   - `mmseqs` — fast k-mer seed + Smith-Waterman alignment via `mmseqs2 easy-cluster` (the practical
     realization of the seeded idea).
5. **Repeat-hub gate from VG topology** — any exon node used by >= `T` loci is labeled a repeat hub;
   splice edges incident to hubs are removed.  No minimizer repeat catalog is used.
6. **Family linkage** (mode-dependent):
   - `exact` / `cdhit` / `mmseqs` — graph-to-graph: two loci are linked if they share a non-repeat VG
     node or a surviving splice edge.
   - `o1vg` (default) — integrate the existing O1 transcript-homology edges
     (`bench/denovo_family_edges.tsv`, `core_recip >= 0.13`) and **gate** each edge by VG support:
     keep the O1 edge only if the two loci share a non-repeat exon node or a surviving splice edge.
     This replaces the minimizer-based repeat gate with a VG-topology support check and directly
     marries O1 homology with the exon-splice graph.
7. **Refine** raw connected components with the same γ-quasi-clique refinement
   (`genome_family_def.refine_families`, γ = 0.20) used by the shipped pipeline.
8. **Evaluate** with the shipped family-level P/R machinery (protein `E_p` purity, DNA-loose cDNA
   component recovery, diploid-DNA oracle P/R).

---

## 2. Results

| mode | multi-copy families | P(Ep) | cDNA pair recall | cDNA block-overmerge P | component recovery | P_oracle(dedup) | R_oracle |
|---|---|---|---|---|---|---|---|
| **current default** (`family_rna_refine.tsv`) | 566 | **0.8940** | **0.9043** | 0.4631 | **0.9722** (175/180) | **0.9167** | 0.8947 |
| `exact` (graph-to-graph) | 214 | 0.8879 | 0.5587 | **0.6135** | 0.9286 (104/112) | **0.9423** | 0.8596 |
| `cdhit` (95 % exon clusters) | 237 | 0.8776 | 0.6152 | 0.5902 | 0.9630 (130/135) | 0.9020 | 0.8596 |
| `mmseqs` (90 % k-mer seeded + align) | 255 | 0.8863 | 0.6630 | 0.6126 | 0.9786 (137/140) | 0.9000 | 0.8596 |
| `o1vg` + exact support | 192 | 0.8958 | 0.4935 | 0.6312 | 0.9346 (100/107) | 0.9400 | 0.8070 |
| **`o1vg` + cdhit support (default)** | 538 | 0.8866 | 0.7870 | 0.4510 | 0.9716 (171/176) | 0.9038 | **0.9123** |

Notes:
- The **default `o1vg` mode** is competitive with the shipped catalog while being fully VG-native:
  it uses the O1 homology graph for signal and the VG for repeat gating, with no minimizer
  repeat catalog.
- Pure graph-to-graph (`exact`) is extremely fast (~1–2 min) and gives the highest oracle
  precision, but it under-merges because it requires an *exact* shared exon sequence; paralogs
  separated by a few SNPs are not linked.
- Near-exact exon clustering (`cdhit`) recovers more real edges and improves recall, at a modest
  precision cost.
- A k-mer-seeded pure-Python implementation (`--mode seeded`) is also included, but it is too slow
  for the full 290 k exon set (>30 min for a fraction of the reps).  `mmseqs2` mode is the practical
  realization of the same seed-and-extend idea: it runs in ~3 min and gives better recall than
  `cdhit` while keeping oracle precision at 0.90.
- Using VG support on top of O1 edges (`o1vg`) gives the best recall (R_oracle = 0.9123) of any
  mode tested, with precision close to the shipped default.

### Residual false-positive roster (`o1vg` + cdhit default)

- **multifam** (block spans >= 2 oracle genes): 3 blocks.
  - `LOC101142904 + LOC129526550` (dl=13)
  - `FOXO1 + LOC115933254` (dl=9)
  - `GSTM2 + LOC101129940` (dl=3)
- **oversize** (RNA loci > 1.5× diploid CN): 0.
- **allele-as-copy** (RNA multi-locus, DNA CN = 1): 2.
  - `DHRSX` (dl=2, hapCN=1)
  - `LOC129530050` (dl=2, hapCN=1)
- **E_p-impure not named by DNA oracle**: 56 blocks (the usual conserved-domain residual).

---

## 3. Implementation notes / gotchas

- `cd-hit-est` clustering of ~290 k representative exons at 95 % identity takes ~10 min on 4
  threads.  A parser bug (originally parsing the cluster-line index instead of the sequence id) was
  fixed; dropped sequences are now kept as singleton clusters.
- `vsearch` was tested but `--cluster_fast` on this many short sequences was slower in practice,
  and its default `--minseqlength 32` discarded short exons.  The code now passes
  `--minseqlength 1` and guards against missing representatives.
- A pure-Python `--mode seeded` (k-mer index + `edlib` verification) is implemented as a reference,
  but it is too slow for the full 290 k exon set.  `--mode mmseqs` uses the same biological
  operation but executes it via `mmseqs2 easy-cluster` and finishes in ~3 min.
- The prototype output files are `bench/vg_family_prototype.tsv` and
  `bench/vg_family_prototype.json`.  The eval script writes
  `bench/vg_family_prototype_eval.tsv` and `bench/vg_family_prototype_eval.json`.

---

## 4. Honest take-aways

1. **VG topology can replace the minimizer repeat catalog.**  The node-multiplicity repeat-hub gate
   cuts repeat bridges using only the graph structure, and the `o1vg` mode performs on par with the
   shipped catalog.
2. **Exact graph-to-graph linkage is too strict for family definition.**  It is clean and fast, but
   whole-copy homology (the O1 signal) is still needed to link paralogs whose exons differ by SNPs.
3. **The productive integration is O1 edges + VG support.**  Keep the transcript-homology signal for
   recall, but use the exon-splice VG to decide which edges are structurally supported and which are
   repeat bridges.  This is the direction that unifies O1 and O2 around one graph object.

---

## 5. Next steps

- Tune the repeat-hub threshold `T` and the O1 edge threshold `core_recip` jointly; the current
  defaults are the first plausible values, not optimized.
- Add `aln_frac` to the O1 edge gate (currently only `core_recip` is used for the raw edge set).
- Replace the external `cd-hit-est` exon clustering with a VG-native approximate-match layer, e.g.
  bubble-aware or k-mer-based exon collapsing inside the graph, so the whole pipeline becomes a
  single VG operation.
- Feed the resulting VG directly into O2 copy assignment: the same exon-SUN/PSV graph can be
  materialized per family for read assignment.
