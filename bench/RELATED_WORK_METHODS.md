# Methods for finding / defining multi-copy gene families (related work)

The field splits by the **level** the signal comes from. The decisive axis is the right-hand column: most
methods work for **ancient / divergent** families (divergence gives the signal) but **break on recent,
near-identical paralogs** — where the field has moved to segmental-duplication detection + read-depth copy
number + long reads. The RNA/transcriptome cell is nearly empty; that is where this thesis sits.

| Method / tool | Level | Input | Core idea | Recent / near-identical paralogs? | Ref |
|---|---|---|---|---|---|
| **All-vs-all BLAST/DIAMOND + MCL** | protein/seq homology | proteomes | pairwise similarity graph → Markov clustering | over-merges or collapses near-identical; threshold-bound | OrthoMCL; Enright 2002 |
| **OrthoFinder** | protein/seq homology | proteomes | DIAMOND all-vs-all + length/bias normalization + MCL + gene trees → orthogroups (current standard) | struggles: near-identical copies collapse to one node / over-merge | Emms & Kelly 2019 |
| **SonicParanoid / Broccoli / eggNOG** | protein/seq homology | proteomes | faster/precomputed orthogroup inference | same recent-paralog limits as OrthoFinder | various |
| **CD-HIT / MMseqs2 linclust** | seq identity | protein or NT | greedy clustering at an identity threshold | threshold-defined; arbitrary boundary on recent paralogs | Fu 2012; Steinegger 2018 |
| **Pfam / InterPro / HMMER** | domain/profile | protein | group by shared **domain architecture** (HMM profiles) | groups *domain-sharers* (FP for paralogy); domain ≠ copy | Finn 2014 |
| **Ensembl Compara (TreeBeST)** | phylogenetic | gene + species trees | gene-tree ↔ species-tree **reconciliation** → ortholog/paralog, dup/loss | hypersensitive to gene-tree error → **fails on recent near-identical** paralogs | Vilella 2009 |
| **GeneRax / Notung / ORTHOSCOPE** | phylogenetic | gene + species trees | probabilistic / parsimony reconciliation | same recent-paralog tree-error limit | Morel 2020 |
| **MCScanX / DupGen_finder / WGDdetector** | synteny | genome + annotation | **collinearity** blocks → classify WGD / tandem / proximal / transposed / dispersed duplicates | needs annotation + synteny; weak inside complex SD arrays | Wang 2012 |
| **WGAC** | DNA self-alignment | assembly | whole-genome assembly comparison → segmental duplications | classic SD detection; assembly-quality bound | Bailey 2001 |
| **⭐ SEDEF** | DNA self-alignment | assembly | **Jaccard similarity + local chaining** → segmental-duplication pairs (fast) | **strong on recent SDs** (the hard case); assembly-only, no expression | Numanagić 2018 |
| **BISER** | DNA self-alignment | assemblies | SEDEF successor: multi-genome, error/masked-aware SD detection | same SD paradigm, scaled to many genomes | Numanagić 2021 |
| **WSSD / fastCN (famCN)** | DNA copy-number | short/long reads | read-depth over multimapping windows → **family copy number** | quantifies CN; cannot phase WHICH copy | Bailey 2002; Pendleton |
| **QuicK-mer2 / fastCN (parCN)** | DNA copy-number | reads + k-mers | paralog-specific k-mers → **paralog-specific copy number** | the parCN standard (Soto 2025); needs paralog-distinguishing k-mers | Shen & Kidd 2020 |
| **CNVnator / read-depth CNV** | DNA copy-number | short reads | depth segmentation → CNV | coarse; not paralog-resolved | Abyzov 2011 |
| **Mash / minimizer-Jaccard** | alignment-free | sequences | MinHash / minimizer Jaccard distance | fast screen; **cannot separate domain-sharers** from true paralogs | Ondov 2016 |
| **⭐ SDA (Vollger)** | long-read DNA | long DNA reads | collapsed-segdup detection by depth → **PSV correlation-clustering** → per-paralog assembly | resolves collapsed segdups; **same K=0 identifiability floor** as ours | Vollger 2019 |
| **Soto 2025 (HSD pipeline)** | DNA + CN | T2T assembly + reads | SD98 self-map + shared-exon (>99% cov) + famCN MAD<1 grouping; parCN via QuicK-mer2 | the human gold-standard for recent families; **DNA, no RNA/isoform axis** | Soto 2025 (Cell) |
| **Guitart/Eichler 2024 (TBC1D3)** | long-read RNA | Iso-Seq + haplotypes | map FLNC to all haplotypes; **assign read to a paralog iff AS ≥ 10** else discard; phylo groups (graph failed) | per-copy expression assignment for ONE family; not de-novo / genome-wide | Eichler 2024 |
| **longcallR** | long-read RNA | long RNA-seq | CNN SNP caller + MEC phasing + ASE/ASJ | **uniquely-mappable only** — assigns to genes+haplotypes, never to family COPIES | Huang 2026 |
| **RSEM / Salmon / kallisto** | RNA quant (EM) | RNA-seq + transcript set | EM over multireads to apportion expression | needs a **pre-defined** transcript set; does NOT discover families | Li 2011; Bray 2016 |
| **⭐ IsoCon** | long-read RNA | targeted Iso-Seq | NN-graph cluster/correct → per-variant-position real-vs-error test → family transcripts | **closest prior art**; WITHIN a known family, targeted (RT-PCR), no de-novo / no copy-vs-allele | Sahlin 2018 |
| **⭐⭐ THIS THESIS** | long-read RNA, de-novo, genome-wide | Iso-Seq/HiFi | **read-conflict graph** (significance de-tie) defines families = identifiability components; **per-read significance gate** assigns reads to copies (no 1/k, assign-or-abstain); + allele-specific junctions | targets exactly the recent/collapsed regime (MAPQ-0); identifiability theorem (MCC=χ(H)) + per-read certificate | — |

## The dichotomy (how to read the table)

- **Ancient / divergent families** → homology clustering, domain profiles, gene-tree reconciliation. Divergence
  *is* the signal, so these work — but they **break on recent near-identical paralogs** (gene trees collapse,
  similarity grouping over-merges or can't separate domain-sharers).
- **Recent / nearly-identical paralogs (the hard case)** → the field switched to **segmental-duplication
  detection** (SEDEF/BISER) + **read-depth copy number** (famCN/parCN, QuicK-mer2) + **long reads**
  (SDA on DNA; IsoCon/Eichler on RNA). Short reads fail here — Soto 2025 measured SNV sensitivity dropping to
  **0.85% in SD98** regions, which is the long-read rationale.
- **RNA / transcriptome cell is nearly empty.** IsoCon resolves transcripts *within a known family*;
  RSEM/Salmon need a pre-defined transcript set; longcallR/Eichler stay in the uniquely-mappable or single-family
  regime. **No method does de-novo, genome-wide, RNA-level family discovery + copy-assignment** — the cell this
  thesis occupies.

## Where this thesis sits (the gap)

| | de-novo discovery | genome-wide | RNA-level | copy assignment | copy-vs-allele |
|---|:-:|:-:|:-:|:-:|:-:|
| OrthoFinder / Compara | ✓ | ✓ | ✗ (protein) | ✗ | ✗ |
| SEDEF / Soto / SDA | ✓ | ✓ | ✗ (DNA) | ~ (parCN) | ✓ (DNA) |
| IsoCon | ✗ (targeted) | ✗ | ✓ | hypothesised | ✗ |
| Eichler / longcallR | ✗ / ✓ | ✗ / ✓ | ✓ | ✓ (1 family) / ✗ | ✗ |
| **This thesis** | **✓** | **✓** | **✓** | **✓** | **✓ (ASJ; DNA for the irreducible case)** |

The three levels are kept **separate and orthogonal** (RNA = the definition; DNA = SEDEF/Soto; protein =
`protein_family_verify.py`); on real GGO, SEDEF independently confirms 73% of the RNA read-conflict families as
segmental duplications (`bench/FAMILY_LEVELS.md`).
