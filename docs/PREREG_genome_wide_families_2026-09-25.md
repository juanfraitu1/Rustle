# Pre-registration: genome-wide family definition on every sample (Figures 6 and 7)

**Written 2026-09-25 13:20, before any genome-wide family number exists.** Checked at the time of writing: no sample
has a genome-wide families run or copy catalog (`${work}/runs/*/<id>.fam.clusters.tsv` and `<id>.cat.copies.tsv`:
0 files), no genome-wide protein-homology family table exists (`truth.py protein-homology --chrom ALL`: only the
first pilot shard of each species), and no guided all-vs-all has been run on more than one chromosome. The only
genome-wide object that exists is the Ensembl Compara export itself (release 116, 3,370,724 rows, fetched
2026-09-25). User: every publication figure must be genome-wide on every sample, with apples-to-apples comparisons
and advisor-proof wording (`figures/GLOSSARY.md`).

What this document replaces:

| figure | today (development tables, kept) | genome-wide (this document) |
|---|---|---|
| Fig. 6a–c | human A119b chr16: Compara same-chromosome pairs; chr16 catalog | human A119b and human testis, whole genome: Compara pairs across chromosomes; each sample's genome-wide catalog |
| Fig. 6d | gorilla OR6737 NC_073244.2: seeding configurations vs that contig's protein-homology families | all six samples: seeding configurations vs genome-wide protein-homology families of the species |
| Fig. 7 | 5 human chromosomes, 2 gorilla contigs, per-chromosome runs of both modes | all six samples, one genome-wide run per mode; whole genome, genome minus development contigs, per-contig breakdown |

The earlier pre-registrations (`PREREG_identity_spectrum_2026-09-24.md`, `PREREG_heldout_families_2026-09-20.md`,
`PREREG_locus_read_pool_2026-09-22.md`) and their numbers stay as they are: development results for this document,
never overwritten. Every genome-wide number is a new register row.

**This figure pair decides no pipeline rule.** Nothing in `src/` or in the driver's defaults changes with any outcome
below. The outcomes decide the wording of the claims only.

## 1. Samples and exposure

| sample | species | reads | used for a family decision before? |
|---|---|---|---|
| human_A119b | human | A119b IsoSeq, CHM13 v2.0 | yes: chr16 (development), see the contig table |
| human_testis | human | public testis IsoSeq (ERR13885926), CHM13 v2.0; **mapped without `-uf`** | never |
| gorilla_OR6737 | gorilla | OR6737 testis, mGorGor1 | yes: NC_073244.2 (chr20) chose the seeding rule |
| gorilla_KB3781 | gorilla | KB3781 fibroblast cell line, mGorGor1 | never |
| chimp_PTR | chimpanzee | mPanTro3 | never |
| orangutan_PPY | orangutan | mPonPyg2 | never |

Numbers from different samples are never pooled, and numbers from different species never are either. Four of
the six samples (human testis, gorilla KB3781, chimpanzee, orangutan) have never been used for any family decision:
they are held-out samples in full. The family rules were developed on human families (NPIP, TBC1D3, AMY, the Soto
families); homologous families exist in the apes, so "never used" refers to the reads and the genome, not to the
biology.

**Contig exposure (human CHM13, both human samples).**

| class | contigs | the decision taken on it |
|---|---|---|
| development | chr16 | every early family-rule decision (NPIP), and the identity spectrum |
| threshold selection | chr5, chr7, chr21 | the 0.60 shared-exon threshold, chosen by leave-one-region-out against Soto families (register 903) |
| reused verdict set | chr2, chr8, chr10 | about 30 pre-registered guided-mode tests since 2026-09-20 (registers 903–1035) |
| scored once, no decision | chr6 | the per-chromosome Figure 7 of 2026-09-25 (no rule chosen on it) |
| never used for a family decision | every other contig | — |

**Gorilla (mGorGor1, both gorilla samples):** development NC_073244.2 (chr20: the seeding rule, registers
1059/1060/1100); scored once, no decision: NC_073234.2 (chr10, per-chromosome Figure 7). **Chimpanzee, orangutan:**
no contig was ever used.

**Substrates** (a restriction of ONE genome-wide run, never a separate per-chromosome run):

- **S0 whole genome:** every contig of the sample's annotation.
- **S1 genome minus development contigs (the headline):** human minus chr16; gorilla minus NC_073244.2; chimpanzee
  and orangutan = S0.
- **S2 genome minus every contig used for a family decision (secondary, stricter):** human minus chr16, chr5, chr7,
  chr21, chr2, chr8, chr10; gorilla minus NC_073244.2; chimpanzee and orangutan = S0.
- **S3 per-contig breakdown:** the genome-wide run restricted to one contig, for every contig; drawn as a strip with
  the classes above marked. It is not the per-chromosome run: the families were built genome-wide, so a family
  restricted to one contig can have lost members elsewhere. It shows that no contig was picked.

Restriction rule (every substrate, every truth): a reference family keeps its members on the kept contigs and is
scored when at least 2 remain; a predicted family keeps its loci on the kept contigs. Pairs (Figure 6) are kept
when both genes lie on kept contigs.

## 2. Objects

### 2.1 De novo families (Figure 7; Figure 6d, both seeding configurations)

- **Default seeding** (loci from primary alignments plus secondary alignments scoring at least 98% of the read's
  best alignment score): the run-cache stage `families` of each sample, i.e. `tools/rustle_pipeline.sh families` =
  `mcl_families --from-gtf <id>.gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60` on the genome-wide assembly.
- **Primary alignments only** (Figure 6d only): the same driver stage on `<id>.primary.gtf` (stage
  `families_primary`, requested from the run-cache owner: `families` with suffix `.primary`, needing
  `assemble_primary`).
- The all-vs-all runs through `tools/mm2_shard.sh` (`RUSTLE_MINIMAP2`), which is `cmp`-identical to one minimap2 run
  (check V1, `phase2a_shard`). Families may therefore cross chromosomes.

### 2.2 Guided families, one run per species (Figure 7)

- **Nodes:** every `gene` and `pseudogene` record of the species' RefSeq GFF on every contig
  (`chrom:start-end`, sorted unique), bodies by `samtools faidx -r`.
- **All-vs-all: the de novo flags**, `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes`, through
  `tools/mm2_shard.sh paf`. **This is a deliberate change from the recorded guided recipe** (`-x asm20 -c --eqx -P`),
  for two reasons fixed now:
  1. The two modes then differ in their node set only. The per-chromosome Figure 7 could not separate the node set
     from the flags (its caveat). `mcl_families` reads PAF columns 1–11 only (`graph_from_paf_loci`), so `--eqx`
     changes nothing it reads; `-X` = `-P -D --dual=no --no-long-join` drops self and dual (B→A) records.
  2. Scale. The recorded recipe wrote 541 MB for human chr2 alone (4,053 bodies), 6.1× the de novo PAF of the same
     chromosome (88 MB); genome-wide it projects to tens of GB, which neither the disk (178 GB free, 90% full) nor
     `mcl_families` (which reads the whole PAF into memory) can hold.
- **Bridging measurement (reported, decides nothing):** on the two development contigs (human chr16, gorilla
  NC_073244.2), per-chromosome guided families with the recorded flags vs the de novo flags, both scored against
  the contig's own protein-homology families and (human) Soto 2025. Reported as the effect of the flags.
- **Family rule:** `mcl_families --paf <paf> --gff <full GFF> --min-exonic-bp 1 --min-shared-exon-frac 0.60`,
  unchanged. Guided reads no RNA: one run per species, drawn as a reference on each sample of that species.
- `inputs` key `fig7_guided_flags=recipe` switches to the recorded flags (`tools/mm2_shard.sh guided`); the default is
  `denovo`, fixed here.

### 2.3 Genome-wide copy catalog (Figure 6a–c, human samples)

The run-cache stage `catalog` (`gw_family_catalog --piecewise`, `cmp`-identical to one run, check V2), its
`copies.tsv` and `pairs.tsv` (the shortest within-family path between two copies).

## 3. Reference families and pairs (never "truth" outside the simulation)

1. **Protein-homology families, 4 species** (`bench/truth.py protein-homology --chrom ALL`): one protein per gene
   (longest CDS), pseudogenes and immunoglobulin / T-cell-receptor segments excluded, all-vs-all BLASTP (E ≤ 1e-5),
   linked when non-overlapping HSPs cover ≥ 30% of the longer protein, MCL (I = 2.8), families of ≥ 2 genes. One
   graph over the whole genome, so families cross chromosomes; genes keyed (contig, name). Independent of both
   modes' family rule.
2. **Ensembl Compara (release 116) human paralogue pairs** (`truth.py compara`, 3,370,724 rows), cross-chromosome
   pairs kept. Figure 6a–c uses the pairs. Figure 7 uses **Compara families**, defined here:
   - nodes: genes whose symbol is a RefSeq CHM13 `gene` record with `gene_biotype=protein_coding` on the chromosome
     Compara gives for it (contig = `chr` + Compara chromosome; `MT` → `chrM`); unmatched genes are counted and
     dropped;
   - a family = a connected component of the pairs whose Compara duplication node (`subtype`) lies **at or below
     Primates** (Homo sapiens, Homininae, Hominidae, Hominoidea, Catarrhini, Simiiformes, Haplorrhini, Primates).
     Within one reconciled gene tree the relation "the two genes' duplication node is at or below level L" is
     ultrametric, hence transitive: the components are the gene-tree clades below duplication nodes of age ≤ L, not
     single-linkage chains, and no identity threshold is used. Primates is the headline because it is the youngest
     level that contains every duplication the copy-level objectives address (segmental duplications of the
     great-ape and primate lineages); it is fixed before any score exists.
   - Level sweep (table only, not the headline): cumulative cuts at Hominidae, Catarrhini, Primates, Eutheria,
     Vertebrata and all levels (Opisthokonta).
   - An unknown `subtype` value stops the build (the level order is a fixed list).
3. **Soto et al. 2025** (Table S1C, human): first family per gene, matched by exact name, expanded to every contig
   where the name has a RefSeq gene record (the `--chrom ALL` semantics of `family_score`). **Not independent:** the
   0.60 exon threshold was chosen against Soto families (register 903). Labelled so.
4. **NPIP reference set** (U2, register 990): human chr16 only, an inset, never a headline.

## 4. Figure 6 metrics and claims

### 4.1 Panels a–c: human A119b and human testis, whole genome

- **Pairs scored:** Compara pairs (both genes present) whose two genes both have a locus (representative ≥ 200 bp)
  in the sample's genome-wide assembly (`score.py spectrum --chrom ALL`), in Compara protein-identity bands
  ≥ 90, 80–90, 70–80, 60–70, 50–60, 30–50, < 30%.
- **Direct edge:** T1 (asm20, identity ≥ 0.80) or T2 (asm20 `-k11 -w5`, identity ≥ 0.60), coverage ≥ 0.50, between
  any loci of the two genes. **T3 (translated protein search) is not run genome-wide**: it is not part of Rustle's
  rule, costs hours and an estimated ≥ 13 GB of output per sample. Its chr16 values (development chromosome) are
  kept as a labelled example.
- **Same family:** both genes have copies in one family of the sample's genome-wide catalog; split by the
  catalog's own shortest within-family path (1 = direct catalog edge; ≥ 2 = only through other copies); the
  families' ceiling = both genes have a catalog copy.
- **Precision** over judgeable pairs (both genes have Compara data): direct edges by the tier's own identity band;
  catalog within-family gene pairs, all copies and multi-exon copies.
- **New strata:** same-chromosome vs cross-chromosome pairs; the largest paralogue group (connected component of the
  scored pairs) vs the rest; whole genome vs genome minus pairs with a chr16 gene.

**Claims, and what changes the wording** (per sample; both human samples must agree for a claim to be stated
without qualification):

| id | claim | stated if | otherwise |
|---|---|---|---|
| F6.1 | at ≥ 90% the families reach their ceiling | families recover ≥ 0.90 of the ≥ 90% pairs whose two genes both have catalog copies | report the fraction; "the ceiling is not reached" |
| F6.2 | at 60–90% the families and the direct edge are not separable | their Wilson 95% intervals overlap in each of the three bands | name the band and the direction |
| F6.3 | the 60–90% result is no longer one gene group | the largest paralogue group holds < 50% of the 60–90% pairs (chr16: 59 of 66) | "still dominated by <group>" |
| F6.4 | cross-chromosome pairs are recovered as well as same-chromosome pairs | at ≥ 90%, family recall on cross-chromosome pairs is within 0.10 of same-chromosome | report both |
| F6.5 | multi-exon catalog families are precise | multi-exon within-family precision ≥ 0.80, while all-copies precision is lower | report both |

### 4.2 Panel d: all six samples, both seeding configurations

- **Reference:** the species' genome-wide protein-homology families.
- **Genes scored (recall):** protein-homology genes with **≥ 2 reads whose primary alignment (not 0x100, 0x800,
  0x4) has an aligned block (M/=/X) on one of the gene's annotated exons** (`_o1.primary_counts` rule, unchanged),
  counted per contig from each sample's BAM.
- **Pair bands:** best minimap2 identity (`-x asm20 -k11 -w5 -c -X -N 100 -p 0.1 --secondary=yes`, matches ÷ block
  length) between the two genes' annotated mRNAs; the mRNA of a gene is the spliced exons of the transcript with the
  most CDS bases (the transcript whose protein entered the families), all-vs-all over every family member of the
  species, through `tools/mm2_shard.sh paf`. Bands ≥ 90, 80–90, 70–80, 60–70, < 60, none.
- **Headline:** recall of same-family pairs at **≥ 90%** mRNA identity, and pair precision over judgeable predicted
  pairs (both genes in a protein-homology family), per sample, per configuration, on S1. Also: 80–90%, the value
  without the configuration's largest predicted family (≥ 90% and precision), S0.

| id | claim | stated if | otherwise |
|---|---|---|---|
| F6.6 | seeding with secondary alignments raises ≥ 90% recall | default > primaries-only on ≥ 5 of 6 samples, none lower | list the samples where it does not |
| F6.7 | the gain is not one array | the gain survives removing each configuration's largest predicted family on ≥ 4 of 6 samples | "the gain is concentrated in single arrays on k samples" |
| F6.8 | seeding does not cost precision | default precision ≥ primaries-only − 0.05 on every sample | name the samples |

Predicted before looking: F6.6 holds (tied reads concentrate in near-identical copies); F6.7 is uncertain (on
NC_073244.2 the whole gain was one array); F6.8 holds.

## 5. Figure 7 metrics and claims

- **Scorer:** `family_score --chrom ALL --per-family --pairwise` (Rust; its per-family rows and pair counts equal
  `figures/_o1_recovery.score_arm`, checked on 16 arms, `phase2a_rust`). One-to-one bipartite sensitivity,
  precision (an upper bound: genes without a reference family are not scored) and F (pooled); pairwise
  sensitivity and precision; per reference family exact (F = 1) / partial / missed.
- **Arms:** de novo per sample (6); guided per species (4). Truths: protein-homology (4 species), Compara families
  at Primates (human), Soto 2025 (human), NPIP reference set (human chr16 inset).
- **Substrates:** S1 headline; S0, S2 in the table; S3 strip for the protein-homology families (all species) and
  the table for every truth.
- The figure is labelled **Rustle-internal: two modes of Rustle, not a comparison with other tools.**

| id | claim | stated if | otherwise |
|---|---|---|---|
| F7.1 | guided recovers the protein-homology families at least as well as de novo | guided F ≥ de novo F on S1 for every sample | name the samples where de novo is higher |
| F7.2 | the guided − de novo gap is not a development-contig artefact | the sign of the gap on S2 equals its sign on S1 for every sample | report where it flips |
| F7.3 | young families are recovered better than the whole protein-homology set | for human, F against Compara (Primates) > F against protein-homology families, in both modes, both samples | report the values |
| F7.4 | no contig is an outlier that carries the result | on the S3 strip, removing the single contig with the largest (guided − de novo) gap does not change the sign of the S1 gap, per sample | name the contig |
| F7.5 | both modes miss most protein-homology families | missed > 50% of protein-homology families in both modes, every sample (the alignability bound of Fig. 6) | report the fractions |

Predicted before looking: F7.1 holds on every sample (it held in all 12 per-chromosome comparisons); F7.3 holds;
F7.5 holds.

## 6. Checks before any figure number is quoted

- The bridging measurement (Section 2.2).
- Scoring path: `family_score` on a genes-only GFF (gene / pseudogene / ncRNA_gene records) equals the full GFF;
  `--chrom ALL` on one contig's inputs equals the per-chromosome scorer (done, 72/72, `phase2a_rust`).
- The genome-wide protein-homology families restricted to NC_073244.2 are **not** the per-contig families (one
  BLASTP database and one MCL over the genome); report both family counts, never mix them.
- `primary_counts` genome-wide on gorilla OR6737 NC_073244.2 equals the recorded per-contig counts for the genes
  present in both gene sets.
- The mRNA recipe (Section 4.2) is compared with the recorded `referee_mrna.fa` of NC_073244.2 (its builder is not
  in the repository); a difference is reported, the recipe above stands.

## 7. Scale, sharding and stop rules

Every heavy unit is resumable and runs in the foreground under `flock /mnt/linuxdisk/tmp/rustle_heavy.lock`, one at
a time, in calls of ≤ 10 minutes and ≤ 20 GB: the families and catalog stages through `make.py runs` (exit 75 = call
again); everything else through `make.py data fig6|fig7 --set figs_budget_s=540` (exit 75 = call again).

**Pilots first (the first shards of the real runs, so nothing is wasted):** one shard of each genome-wide all-vs-all
(de novo loci of human A119b and gorilla OR6737, guided bodies of human and gorilla, the spectrum's T1 and T2 on
human A119b), with wall time, peak RSS and PAF bytes read from `tools/mm2_shard.sh status` and the shard logs.

**Stop rules, fixed now:**
- A single node set projected above 24 h of all-vs-all, or above 30 GB of PAF, or a `mcl_families` /
  `score.py spectrum` projected above 20 GB of memory: stop that arm and bring it to the user. No automatic fallback
  (a per-chromosome sum is not a substitute: it cannot represent the cross-chromosome families).
- If the de novo families of both seeding configurations cannot be run on every sample, Figure 6d shows the samples
  that have both, and says which are missing; no sample is replaced by a subset of its genome.
- Disk: after a PAF is complete and consumed, the wrapper's shard directory for that key is deleted (the product and
  the Rust PAF cache remain). No genome-wide T3 output is written.

I will not change the reference definitions (including the Compara level), bands, substrates, flags or claims after
seeing any genome-wide number. Amendments go below, dated, with what had been seen when they were made.

---

## Check results (2026-09-25, before any genome-wide number; Section 6)

- `primary_counts_gw` on gorilla OR6737 NC_073244.2 equals the development counts for all 775 genes (both rules).
- The annotation cache's exon unions equal `_o1._gene_exons` for all 2,048 NC_073244.2 gene records.
- mRNA recipe vs the recorded `referee_mrna.fa` (NC_073244.2): 546 of 775 sequences identical. The recorded builder
  took each gene's longest transcript (775 of 775 match by length); the recipe of Section 4.2 (the transcript with
  the most CDS bases) stands, as fixed above. On that contig the >= 90% band holds 78 pairs under this recipe vs 83
  under the recorded one (development contig; not a genome-wide number).
- Guided all-vs-all through `tools/mm2_shard.sh paf` with the de novo flags on gorilla chrY (NC_073248.2, 1 shard):
  `cmp`-identical to a direct minimap2 run, and the `mcl_families` clusters are identical. The PAF is 2.7 MB vs 13.2 MB
  with the recorded flags (4.9x), the disk argument of Section 2.2.
- Compara families at Primates: 426 families, 1,313 genes; 15,254 Compara genes match a protein-coding CHM13 RefSeq gene
  on the chromosome Compara gives, 4,106 do not (non-coding, renamed or absent).
- The genome-scope builds of Figs 6 and 7 ran end to end on fixtures stitched from the per-chromosome development
  products (human chr16, chr2, chr6, chr8, chr10; gorilla chr20, chr10): code paths only, not genome-wide numbers,
  not recorded.
- Found while checking: `mcl_families` reads its whole `--paf` into memory (`std::fs::read_to_string`), so a
  genome-wide PAF of G GB needs more than G GB of RAM. Stop rule of Section 7 applies (20 GB); a streaming reader is
  requested (byte-identical outputs).

## Amendments

**Amendment 1 (2026-09-25 17:10, before any genome-wide family number; user decision of 16:00).**

*What had been seen when this was written.* No genome-wide family number of any kind: no sample has a families-stage
output (`${work}/runs/*/<id>.fam.*`: 0 files), no `fig6_gw_*` / `fig7_gw_*` table exists, no genome-wide
protein-homology family table exists (1 pilot shard of 21 per species), no Liftoff self-lift is merged for any
species (human 1 of 36 shards), and the two genome-wide legacy catalogs that exist (chimp_PTR, human_testis) were
never scored against any reference. Development numbers seen before this amendment (all per chromosome, none from a
genome-wide run): the Fig. 6 development tables (human A119b chr16: Compara pairs vs the chr16 legacy catalog;
gorilla OR6737 NC_073244.2: seeding configurations vs that contig's protein-homology families), the Fig. 7
development tables (human chr16, chr2, chr6, chr8, chr10; gorilla NC_073244.2, NC_073234.2; protein-homology
families, Soto 2025, NPIP set), register 1101, and the manual protein step's development tables (figS_protein,
`docs/PREREG_protein_attach_2026-09-25.md`), which include Ensembl Compara judgements of the de novo families of
human chr16 and chr6 (chr6: 14 of 14 judgeable family members correct at any duplication age).

*The decision (user, 2026-09-25 16:00).* ONE default de novo family definition, at the RNA level: reads → seeded
assembly loci → one representative per locus (its "positional exon sum": read-derived exon coordinates, genome
bases) → families = the driver's `families` stage (`mcl_families --from-gtf --min-exonic-bp 1
--min-shared-exon-frac 0.60`, MCL inflation 2.8). Copy assignment consumes the same families (the stage's copy table
`<id>.fam.copies.tsv`, `docs/PREREG_families_copy_table_2026-09-25.md`). The `gw_family_catalog` copy catalog
(`catalog` stage) becomes LEGACY: kept runnable, not the default and not a headline. Protein is not part of the
default; the manual extra-sensitive step has its own supplementary figure (figS_protein). Main figures use external
references only: Ensembl Compara (human), Soto 2025 (labelled not independent), Liftoff copies (every species).

*Changes, fixed now.*

1. **Figure 6 panels a–c (§2.3, §4.1): the default families replace the catalog.** The object scored is the
   families stage's copy table of each human sample (`samples.product(cfg, id, 'families', 'copies')`; one copy per
   member locus of every family, its representative's exons). Pairs, bands, universe (the spectrum's
   `truth_pairs.tsv`), scorer (`score.py pairs --chrom ALL --truth compara:... --universe ...`) and strata are
   unchanged. Three definitions change with the object:
   - *same family*: both genes have copies in one default family;
   - *the split of the family bar* (was: the catalog's within-family path length, from its `pairs.tsv`, which the
     families stage does not write): same family AND a direct alignment between the two genes' loci (the grey
     bar's per-pair flag: T1 ∪ T2, nucleotide) vs same family with no direct alignment;
   - *the families' ceiling* (dashed outline): both genes have a copy in some default family.
   Claims F6.1–F6.5 keep their wording and bars, now about the default families. The legacy catalog is not drawn.
2. **T3 (translated protein search) leaves Figure 6.** The chr16 `chr16_example` rows are no longer written to
   `fig6_gw_recall`; the chr16 development spectrum's tiers (T1, T1∪T2, T1∪T2∪T3, with T3's precision) move to a
   supplementary table and figure (`fig6s_protein_tiers`, `fig6s_protein`), labelled "translated protein search:
   not part of Rustle's rule". Nothing about T3 is run genome-wide (unchanged).
3. **Panel d moves to a supplementary figure (`fig6s_seeding`), rescored on external references.** Reason: it
   compares two read pools (seeding), not two family definitions; the default definition includes the seeding, and
   Fig. 3 already shows the seeding gain genome-wide against the annotation without protein; on the development
   contig the whole family-level gain was one tandem array (register 1116). Its references, in order:
   - *Compara (human samples; primary):* Compara pairs at ≥ 90% protein identity whose two genes both have ≥ 2 reads
     whose primary alignment has an aligned block on the gene's annotated exons (`_o1.primary_counts_gw`, now
     counted for every gene and pseudogene record of the annotation, the same rule); sensitivity = pairs with
     copies in one family of the configuration (`score.py pairs --universe`); precision over judgeable
     within-family pairs (both genes have Compara data; an upper bound).
   - *Liftoff copy pairs (every sample):* the Figure 8 self-lift (`docs/PREREG_liftoff_loci_2026-09-25.md`):
     pairs (a source record's annotated placement, one of its extra copies), extra copy `sequence_ID` ≥ 0.95
     (rows also at 0.98, 0.99, 1.00), both exon unions ≥ 200 bp, both loci read-supported in the sample (≥ 2 reads
     whose primary alignment, `-F 2308`, has an aligned block on the exon union: the C2 rule and table). A pair is
     recovered when a copy covering ≥ 50% of the source locus's exon union and a copy covering ≥ 50% of the extra
     copy's exon union (Liftoff's `-a`, prereg §3) share a family. Denominator: every read-supported pair (fixed,
     independent of both configurations). The conditional fraction of Figure 8's F1 (pairs whose two loci are both
     covered) is reported beside it.
   - *Protein-homology families (secondary, supplementary only):* scored only if the species' table has been built
     (it is no longer built by default; `fig6s_protein_homology=1` builds and scores it).
   Claims restated on these references: **F6.6** default sensitivity > primaries-only on Compara ≥ 90% pairs in
   both human samples, and on Liftoff pairs in ≥ 5 of 6 samples with none lower; **F6.7** the gain survives removing
   each configuration's largest predicted family (same counting); **F6.8** default Compara precision ≥
   primaries-only − 0.05 in both human samples. Predicted before looking: F6.6 holds on Liftoff pairs (the copies
   whose reads tie are near-identical by construction) and is uncertain on Compara ≥ 90% pairs; F6.7 uncertain;
   F6.8 holds. The `families_primary` stage becomes optional: without it the supplement lists the sample as not
   built, and no main figure depends on it.
4. **Figure 7 references (§3, §5).** Protein-homology families are removed from the main figure. Main references:
   Compara families at Primates (human; unchanged definition; the headline), Soto 2025 (human; not independent),
   the NPIP reference set (human chr16 inset), and **Liftoff copy pairs (every species; new)**, defined as in item 3
   (same pairs, same read-support universe per sample, same recovery rule) with the de novo families of the sample.
   The guided arm is **not scored** on Liftoff pairs: its node set is the annotation, and a Liftoff extra copy is
   unannotated by construction (it overlaps no annotated record on its strand), so the guided value is fixed by
   construction, not measured. Liftoff pairs give pairwise sensitivity only (the relation certifies pairs, it is
   not a partition: no precision, no bipartite matching). The protein-homology families move to a supplementary
   figure (`fig7s_protein_homology`, both modes, secondary), built only with `fig7_protein_homology=1`.
   - Panel d (per-contig strip): Compara families, human samples and the human guided run.
   - Panels e–g: human A119b vs Compara; human testis vs Compara; human A119b vs Soto 2025 (labelled not
     independent).
   - **F7.1** guided recovers the Compara families (Primates) at least as well as de novo: guided F ≥ de novo F on S1
     for both human samples (was: protein-homology families, every sample). **F7.2** unchanged, on Compara.
     **F7.3 replaced**: the de novo families place read-supported Liftoff copy pairs in one family — descriptive,
     per sample, no bar (the base rate is unknown; Figure 8's F1 keeps its 0.90 bar on the conditional fraction).
     **F7.4** unchanged, on the Compara strip (human). **F7.5** moves to the protein-homology supplement (secondary;
     stated only if that table is built).
   - The bridging measurement (§2.2) keeps its references (the contig's protein-homology families and Soto): it
     decides nothing and is supplementary.
5. **Checks of §6 that concern protein-homology families or the mRNA recipe** become supplementary: they are run
   only when the protein-homology tables are built.
6. **Development tables** are rebuilt under the same definitions before any genome-wide number: Figure 6 on human
   A119b chr16 (the default families of the chr16 restriction of the genome-wide assembly; the chr16 spectrum), and
   Figure 7 per chromosome, rescored against Compara families at Primates restricted to each human chromosome (a
   family keeps its members on the chromosome and is scored with ≥ 2), from the cached per-chromosome runs of
   2026-09-25 (not re-run). These are development tables. Human chr6 of A119b is also part of A119b's genome-wide
   S1, so for F7.1, F7.2 and F7.4 human testis (never used) must agree with A119b before a claim is stated without
   qualification.

Nothing in `src/` or the driver changes with any outcome (unchanged).
