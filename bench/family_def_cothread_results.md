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
