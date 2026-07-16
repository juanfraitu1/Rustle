# Panels And Notes (consolidated)

> Merged from 10 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [ggo19_needy_top15.panel](#ggo19-needy-top15panel) · [ggo19_needy_top5.panel](#ggo19-needy-top5panel) · [rcf611_graph_vs_linear_results](#rcf611-graph-vs-linear-results) · [CLUSTER_FUTURE_WORK](#cluster-future-work) · [SEDEF_BUILD](#sedef-build) · [LOCI_COPIES_ISOFORMS](#loci-copies-isoforms) · [INTRON_RETENTION_RESCUE](#intron-retention-rescue) · [BAYESIAN_POSTERIOR_ZONE](#bayesian-posterior-zone) · [POA_CORE_COMPLETION](#poa-core-completion) · [CONFLICT_EDGE_UNIFICATION](#conflict-edge-unification)


---

## ggo19_needy_top15.panel

# Needy Ref Panel

- source ref buckets: `GGO_19_after_termfix_missing.ref_buckets.tsv`
- source subopt classes: `GGO_19_after_termfix_missing.subopt_by_class.tsv`
- total needy loci: `110`
- top loci written: `15`

## Global Counts

- bucket `subopt_j`: `53`
- bucket `subopt_c`: `43`
- bucket `absent_single_exon_ref`: `20`
- bucket `absent_full_junction_overlap_query`: `18`
- bucket `subopt_k`: `18`
- bucket `subopt_m`: `6`
- bucket `subopt_o`: `4`
- bucket `subopt_n`: `2`

## Dominant Classes

- class `j`: `53`
- class `c`: `43`
- class `k`: `18`
- class `m`: `6`
- class `o`: `4`
- class `n`: `2`

## Top Loci

1. `STRG.251` `NC_073243.2:36443168-36634779(+)` priority `37` needy `6` absent `0` subopt `6` classes `c:1,j:3,k:2` focus `STRG.251.1,STRG.251.2,STRG.251.4,STRG.251.5,STRG.251.8,STRG.251.9` subset `bench/ggo19_needy_top15_refs/STRG_251.gtf`
2. `STRG.151` `NC_073243.2:23033526-23064635(+)` priority `35` needy `6` absent `0` subopt `6` classes `c:1,j:4,n:1` focus `STRG.151.1,STRG.151.2,STRG.151.5,STRG.151.6,STRG.151.7,STRG.151.4` subset `bench/ggo19_needy_top15_refs/STRG_151.gtf`
3. `STRG.503` `NC_073243.2:97540381-97743862(-)` priority `29` needy `5` absent `0` subopt `5` classes `c:1,k:3,m:1` focus `STRG.503.1,STRG.503.4,STRG.503.5,STRG.503.6,STRG.503.2` subset `bench/ggo19_needy_top15_refs/STRG_503.gtf`
4. `STRG.157` `NC_073243.2:23288512-23308475(-)` priority `27` needy `5` absent `2` subopt `3` classes `j:1,k:1,m:1` focus `STRG.157.6,STRG.157.7,STRG.157.3,STRG.157.5,STRG.157.4` subset `bench/ggo19_needy_top15_refs/STRG_157.gtf`
5. `STRG.453` `NC_073243.2:82243221-82437729(+)` priority `27` needy `4` absent `0` subopt `4` classes `c:3,j:1` focus `STRG.453.2,STRG.453.3,STRG.453.4,STRG.453.5` subset `bench/ggo19_needy_top15_refs/STRG_453.gtf`
6. `STRG.442` `NC_073243.2:80297718-80330444(+)` priority `21` needy `3` absent `0` subopt `3` classes `c:3` focus `STRG.442.3,STRG.442.6,STRG.442.9` subset `bench/ggo19_needy_top15_refs/STRG_442.gtf`
7. `STRG.566` `NC_073243.2:111146297-111177575(-)` priority `20` needy `4` absent `1` subopt `3` classes `j:1,m:1,o:1` focus `STRG.566.12,STRG.566.9,STRG.566.10,STRG.566.11` subset `bench/ggo19_needy_top15_refs/STRG_566.gtf`
8. `STRG.445` `NC_073243.2:80975942-81174569(+)` priority `19` needy `3` absent `2` subopt `1` classes `c:1` focus `STRG.445.2,STRG.445.5,STRG.445.4` subset `bench/ggo19_needy_top15_refs/STRG_445.gtf`
9. `STRG.29` `NC_073243.2:17431800-17453950(+)` priority `18` needy `3` absent `0` subopt `3` classes `k:3` focus `STRG.29.2,STRG.29.3,STRG.29.4` subset `bench/ggo19_needy_top15_refs/STRG_29.gtf`
10. `STRG.440` `NC_073243.2:80273347-80288948(-)` priority `18` needy `3` absent `0` subopt `3` classes `j:3` focus `STRG.440.3,STRG.440.4,STRG.440.5` subset `bench/ggo19_needy_top15_refs/STRG_440.gtf`
11. `STRG.300` `NC_073243.2:44054004-44096254(-)` priority `17` needy `3` absent `0` subopt `3` classes `c:1,j:1,m:1` focus `STRG.300.1,STRG.300.6,STRG.300.5` subset `bench/ggo19_needy_top15_refs/STRG_300.gtf`
12. `STRG.125` `NC_073243.2:22454887-22475442(-)` priority `16` needy `3` absent `0` subopt `3` classes `j:1,k:1,m:1` focus `STRG.125.6,STRG.125.9,STRG.125.5` subset `bench/ggo19_needy_top15_refs/STRG_125.gtf`
13. `STRG.210` `NC_073243.2:30118918-30195392(+)` priority `16` needy `3` absent `2` subopt `1` classes `j:1` focus `STRG.210.4,STRG.210.5,STRG.210.1` subset `bench/ggo19_needy_top15_refs/STRG_210.gtf`
14. `STRG.52` `NC_073243.2:19033042-19077751(-)` priority `14` needy `2` absent `0` subopt `2` classes `c:2` focus `STRG.52.2,STRG.52.4` subset `bench/ggo19_needy_top15_refs/STRG_52.gtf`
15. `STRG.95` `NC_073243.2:20647910-20710598(-)` priority `14` needy `2` absent `0` subopt `2` classes `c:2` focus `STRG.95.2,STRG.95.5` subset `bench/ggo19_needy_top15_refs/STRG_95.gtf`


---

## ggo19_needy_top5.panel

# Needy Ref Panel

- source ref buckets: `GGO_19_after_termfix_missing.ref_buckets.tsv`
- source subopt classes: `GGO_19_after_termfix_missing.subopt_by_class.tsv`
- total needy loci: `110`
- top loci written: `5`

## Global Counts

- bucket `subopt_j`: `53`
- bucket `subopt_c`: `43`
- bucket `absent_single_exon_ref`: `20`
- bucket `absent_full_junction_overlap_query`: `18`
- bucket `subopt_k`: `18`
- bucket `subopt_m`: `6`
- bucket `subopt_o`: `4`
- bucket `subopt_n`: `2`

## Dominant Classes

- class `j`: `53`
- class `c`: `43`
- class `k`: `18`
- class `m`: `6`
- class `o`: `4`
- class `n`: `2`

## Top Loci

1. `STRG.251` `NC_073243.2:36443168-36634779(+)` priority `37` needy `6` absent `0` subopt `6` classes `c:1,j:3,k:2` focus `STRG.251.1,STRG.251.2,STRG.251.4,STRG.251.5,STRG.251.8,STRG.251.9` subset `bench/ggo19_needy_top5_refs/STRG_251.gtf`
2. `STRG.151` `NC_073243.2:23033526-23064635(+)` priority `35` needy `6` absent `0` subopt `6` classes `c:1,j:4,n:1` focus `STRG.151.1,STRG.151.2,STRG.151.5,STRG.151.6,STRG.151.7,STRG.151.4` subset `bench/ggo19_needy_top5_refs/STRG_151.gtf`
3. `STRG.503` `NC_073243.2:97540381-97743862(-)` priority `29` needy `5` absent `0` subopt `5` classes `c:1,k:3,m:1` focus `STRG.503.1,STRG.503.4,STRG.503.5,STRG.503.6,STRG.503.2` subset `bench/ggo19_needy_top5_refs/STRG_503.gtf`
4. `STRG.157` `NC_073243.2:23288512-23308475(-)` priority `27` needy `5` absent `2` subopt `3` classes `j:1,k:1,m:1` focus `STRG.157.6,STRG.157.7,STRG.157.3,STRG.157.5,STRG.157.4` subset `bench/ggo19_needy_top5_refs/STRG_157.gtf`
5. `STRG.453` `NC_073243.2:82243221-82437729(+)` priority `27` needy `4` absent `0` subopt `4` classes `c:3,j:1` focus `STRG.453.2,STRG.453.3,STRG.453.4,STRG.453.5` subset `bench/ggo19_needy_top5_refs/STRG_453.gtf`


---

## rcf611_graph_vs_linear_results

# Graph vs linear copy-assignment — RCF_611 simulation (decisive test)

User hypothesis: a vg graph (modeling the tandem repeat as a cycle) aligns reads to copies BETTER than
linear alignment, because it places reads correctly through the repeat instead of mis-anchoring.

Test: RCF_611 = 4 co-located paralogs (NC_073247.2, ~2 Mb span) with tandem-repeat CNV (cDNAs
599–960 bp, pairwise id 0.857–0.996). 800 KNOWN-ORIGIN reads simulated (200/copy, full-length +
fragments, ~1.2% HiFi error). Assign by LINEAR pairwise (minimap2 -N50 -p0.1, primary target) vs
GRAPH (vg giraffe to a `vg msga` graph of the 4 copies; copy = which copy's UNIQUE graph nodes the
read traverses). `bench/rcf611_graph_vs_linear.py`.

## Result — hypothesis REFUTED

| method | correct | wrong | ambiguous |
|---|--:|--:|--:|
| LINEAR-primary (pipeline) | **98.0%** | 0.9% | 1.1% |
| LINEAR-ASdecisive | 93.4% | 0.0% | 6.6% |
| GRAPH (vg giraffe) | 94.5% | 0.1% | 5.4% |

- The graph **rescues 0** reads linear got wrong; it regresses 1.
- On linear's 13 MAPQ0 reads, GRAPH gets **2/13** correct vs linear's **6/13** (graph WORSE on the hard reads).
- Linear (minimap2) is already 98% accurate.

## Why
- minimap2 already handles these tandem repeats (chaining/scoring is not fooled here) — the earlier
  "linear is mangled" stat was a SECONDARY-alignment artifact (secondaries are MAPQ0 + clipped by
  definition); primary reads are mostly MAPQ>0 and correct.
- The copy-distinguishing signal is in the UNIQUE regions, which BOTH methods align fine. The shared
  tandem repeat carries NO copy information, so collapsing it into a graph cycle adds no discriminating
  power. Reads inside shared sequence stay ambiguous either way.
- Same disjoint-levers result the project keeps finding: copy resolution lives in sequence (unique
  regions); the structural/graph representation is faithful but not a resolution lever.

## Verdict
The vg graph's non-linear (repeat) modeling is real biology but is NOT a copy-assignment improvement
here. Do not graph-align for copy resolution; linear alignment + AS-decisive copy choice is sufficient.
(One family; representative co-located tandem-CNV case with a real MAPQ0 tail. A harder regime — many
near-identical copies, very long arrays — could differ, but RCF_611 shows no benefit.)


---

## CLUSTER_FUTURE_WORK

# Deferred to the cluster (this workstation can't handle the scale)

The WSL2 workstation is memory- and CPU-bounded (≈19 GB RAM, ~5 cores) and one tool (BISER) segfaults
under WSL2. The following are validated-in-mechanism here but must be RUN AT SCALE on the HPC cluster.

## 1. Truly-orthogonal catalog precision — a real SEDEF/BISER segdup map (item #2)
The strongest external check on the family catalog (O1) is a genome-wide segmental-duplication map computed
from the assembly alone (DNA, independent of the RNA reads). Status: BISER installs (pip wheel) but its
precompiled `align` step **segfaults under WSL2** (`-11`, even with `ulimit -s unlimited`); its unrefined
putative SDs are Alu-repeat noise. The targeted genomic-DNA substitute gave 89.2% but is *partly circular*.
**On the cluster:** build SEDEF from source (or run BISER on a non-WSL2 node), produce the segdup BEDPE, then
`bench/validate_segdup.py <refined_prefix> <segdup_bedpe>` for the orthogonal precision (the script + the
within-family linkage logic are ready). Soft-masked `GGO.fasta` (≈40% lowercase) is the input.

## 2. Genome-wide copy-assignment at full scale
The definitive O2 run (74 co-located families, 43,239 reads) took ~2.5 h locally and is bounded by the
significance gate's per-read cost. A full pass (all families incl. larger arrays, the divergent-tail catalog,
abundance EM) is cluster work. Driver: `copy_assign --regions <regions> --min-copies 2` (per-region, OOM-safe);
parallelize across regions on the cluster. The local headline (75.1% assigned / 24.8% certified-tied / 99.9%
of decisive / unique-mapper agreement 99.9%) is the validated method; the cluster run scales it + the divergent-tail.

## 3. Positive demonstrations on real data (items L20 / editing)
Both the gene-conversion discriminator and the RNA-editing filter are unit-validated and conservative
(inert on the stable/diverged loci tested locally). A POSITIVE demo needs broad locus scanning:
- **Editing filter:** scan for an actively-edited locus (an A↔G PSV column with within-copy heterogeneity)
  and show the flag + the assignment change. Cluster: run `copy_assign` (default editing-on) genome-wide and
  collect families where `detect_editing_columns` fires.
- **Gene-conversion discriminator:** scan for a recently-CONVERTED (homogenized / few-PSV) array where reads
  show mid-molecule copy-switches; the `[mosaic]` log surfaces candidate-vs-confirmed per family.

## 4. Heavy installs that failed locally
BISER/SEDEF (segdup), and any tool needing a full Linux toolchain. Use the cluster's module system or build
from source there; mamba envs from this session (`biser`) are reproducible.

## Notes
- Big data lives in `/home/juanfra/winloci_scratch` (BAM `GGO_mm.bam`, the refined catalog, sim5x); mirror to
  cluster scratch. Source is committed on branch `vg/flow-capacity-apportionment`.
- Local validation that stands without the cluster: sim5x labeled-truth ladder (non-circular), the unit +
  calibration test suites (659), and the per-read certificate. The cluster items are *external* precision and
  *scale*, not new mechanism.


---

## SEDEF_BUILD

# SEDEF segdup map — builds AND runs on WSL2 (the P0 external O1 check)

2026-06-27. The "SEDEF/BISER blocked on WSL2" claim is **WRONG**: SEDEF compiles and runs end-to-end here.

## The actual diagnosis
- **Not WSL2.** WSL2 runs native Linux binaries. SEDEF (pure C++) built with the local g++ 14.3 (`-march=native`
  → SSE4.1, the CPU is an i5-8600: AVX2, **no AVX-512**) and ran the search→align→merge pipeline to a correct
  `final.bed` on a synthetic 98%-identity duplication (`fracMatch 0.985`).
- **BISER's `-11`** was a *runtime* crash inside its **Codon-compiled** `biser.exe` (`libcodonrt.so` = Codon
  runtime), NOT a load-time illegal-instruction / WSL2 wall (the binary has 0 AVX-512 insns and prints its
  usage fine). Rebuilding BISER needs the niche **Codon** toolchain. Moot — SEDEF (same group, the predecessor,
  the parCN standard in Soto-2025) works.

## Build recipe (reproducible; also portable to the cluster)
```bash
git clone --depth 1 https://github.com/vpc-ccg/sedef.git && cd sedef
mamba install -n <env> -c conda-forge boost parallel time   # Boost(ICL header-only) + wrapper deps
export CPLUS_INCLUDE_PATH=<env>/include                       # so g++ finds boost/icl
make -j4                                                      # -> ./sedef (13MB binary, SSE4.1)
```
Built binary persisted at `/home/juanfra/winloci_scratch/sedef_build/{sedef,sedef.sh}`.

## Run on the gorilla genome (the real P0 deliverable)
```bash
export PATH=/home/juanfra/miniforge3/envs/biser/bin:$PATH     # samtools, parallel, time, sedef on PATH
cd /home/juanfra/winloci_scratch/sedef_build
./sedef.sh -o GGO_sedef_out -j 2 /home/juanfra/winloci_scratch/GGO.fasta
# input MUST be soft-masked (GGO.fasta is, ~40% lowercase genome-wide) + have a .fai (it does).
# output: GGO_sedef_out/final.bed = the genome-wide segmental-duplication BEDPE.
```
Cost: heavy (3.5 GB mammalian genome; search fans over chromosome pairs, align is the long pole). Many hours on
4–5 cores; RAM-sensitive on the 19 GB box → use `-j 2` (not 4) to avoid OOM, or run on the cluster (it compiles
there identically — `module load boost` or the same mamba recipe).

## Why this matters (defense)
`final.bed` is the **non-circular external check** for O1: `bench/validate_segdup.py <refined_catalog> final.bed`
gives the cross-chrom family PRECISION against an assembly-only segdup map (independent of the RNA reads) — the
single highest-leverage item (P0) that neutralizes the circularity charge in DEFENSE_READINESS_AUDIT.md.


---

## LOCI_COPIES_ISOFORMS

# Loci, copies, and isoforms — the definitional hierarchy

A recurring source of confusion in multi-copy gene-family work is conflating three things that this framework
keeps strictly separate. This note fixes the hierarchy, the orthogonality of copies and isoforms, and the two
regimes where "copy = locus" holds or bends.

## The hierarchy

```
FAMILY   = a set of paralogous LOCI                         (read-conflict component / homology)
  └─ COPY    = one LOCUS — a genomic position (chrom:start–end)   ← the unit the DEFINITION is built on
        └─ ISOFORM  = a splice path (intron chain) of that copy    ← a WITHIN-copy axis
```

- **A family is a set of loci.** The read-conflict graph's nodes are loci (the reps produced by
  `collapse_loci`, which folds *all* the transcripts/isoforms at one genomic position into one locus). A family
  is a connected component of **loci**; a "copy" in the catalog *is* a locus.
- **`collapse_loci` collapses by genomic position, not by splice pattern.** So two loci that happen to express
  the same isoform remain two separate copies — isoforms are never what defines a copy.

## Copies and isoforms are ORTHOGONAL axes

In the variation-graph view, the family is one PSV-aware variation graph with two kinds of variation:
- **PSV bubbles** = sequence differences between copies → the **copy** axis.
- **Splice bubbles** = junction (intron-chain) choices → the **isoform** axis.

A single sequenced molecule is a **path through both**: it picks a copy (a PSV path) *and* an isoform (a splice
path). The two are independent — a copy can express many isoforms, and an isoform can be shared across copies.

| feature a read carries | what it distinguishes |
|---|---|
| a **PSV** (base differs between copies) | the **copy** (genetic) |
| a **copy-specific junction** (one copy has the exon, another lacks it) | the **copy** (genetic) — used in assignment as an extra "column" (§6 junction-as-column lift) |
| an **alternatively-spliced** junction (a copy *sometimes* uses it) | nothing about the copy — it's a *within-copy* isoform choice |
| a junction **shared** by all copies of the family | nothing — non-distinguishing |

**Rule of thumb:** a copy-specific junction distinguishes copies; a shared or alternatively-spliced isoform does
not. Isoforms are reported *per copy* (the assembled splice variants at that locus); they are output, not identity.

## "The same isoform at different loci" — three cases, all handled

This is the case that prompts the worry, and it is not a problem because the definition is locus-anchored, not
isoform-anchored:

1. **Different loci, same isoform, but the copies differ in SEQUENCE (PSVs).** The read is assigned by its PSVs;
   the shared isoform contributes nothing. Resolved — the locus is recovered from the sequence.
2. **Different loci, same isoform, AND sequence-identical over what the read observes.** This is the **K=0
   identifiability floor**: the read maps equally well to both loci (MAPQ 0), and the locus is *not recoverable
   from a single RNA read*. The method certifies it **Tied** (min_p = 1, a real impossibility certificate)
   rather than guessing. The two loci exist; one molecule simply cannot say which.
3. **Different loci, DIFFERENT isoforms (a copy-specific junction).** The copy-specific junction is a
   distinguishing feature, used alongside the PSVs to assign the read.

So "the same isoform at two loci" is, by construction, *non-confusing* — it is either irrelevant to assignment
(the PSVs do the work), or it is the honest K=0 floor (correctly abstained on), never a wrong merge or split.

## The two regimes: where "copy = locus" holds vs bends

- **Reference-present (the clean case):** each paralog occupies its own genomic position, so **copy = locus**
  exactly. The read-conflict family is a set of paralogous loci; assignment routes each read to a locus.
- **Collapsed / reference-absent:** several copies share *one* reference locus — their reads all pile there at
  MAPQ 0. Here **copy ≠ locus**: the copies are **PSV-haplotypes within a single locus**, separated by phasing,
  and the locus alone cannot distinguish them. This is the regime the whole copy-assignment apparatus exists for
  (and where reference-absent copies live).

## The one substantive boundary: copy-specific junction vs alternative isoform

There is a genuinely hard distinction the framework must respect: a junction present in some molecules and absent
in others can be **genetic** (a *copy-specific* junction — between copies) or **regulatory** (an *alternatively
spliced* isoform — within a copy). Telling them apart is the **copy-vs-allele problem**:
- The allele-specific-junctions (ASJ) machinery uses per-molecule **allele→junction linkage** (the long-read
  advantage) to separate them where a heterozygous anchor exists.
- The genuinely ambiguous case — a junction that varies *between paralog copies* that are otherwise
  indistinguishable in RNA — is the **irreducible RNA boundary** and is deferred to DNA (paralog-specific copy
  number, parCN). The framework *flags* these rather than guessing.

## Summary for the thesis

- Families and copies are defined on **loci** (genomic position + sequence), never on isoforms.
- Isoforms are an **orthogonal within-copy axis**, modelled as splice bubbles / junction-columns; they are
  reported per copy and used in assignment *only* where a junction is copy-specific.
- "Same isoform at different loci" is handled by the locus/PSV anchor; the irreducible case is the K=0 floor,
  honestly abstained on.
- The only deep question is copy-specific-junction vs alternative-isoform (copy-vs-allele), addressed by ASJ
  linkage and, for the irreducible case, DNA/parCN.


---

## INTRON_RETENTION_RESCUE

# How much could intron-retention rescue the unassignable? (measured)

The K=0 floor reads are tied because the copies are identical over the **exons** the read observes. The
intron-retention lever: a read that *retains* an intron carries extra (intronic) sequence the spliced-exon gate
never uses — potentially enough to distinguish the copies. This quantifies it on the real GGO Tied reads.

## Factor 1 — do the Tied reads even carry introns? (`bench/intron_retention_rescue.py`)

Of the **972 currently-Tied (unassignable) reads** (consensus splice sites derived from the reads themselves;
a read "retains" an intron if its reference span covers a consensus intron it did not splice):

- **215/972 = 22.1% retain ≥1 intron** (21.6% of the strict-K=0 reads with `n_decisive==0`).

So ~1 in 5 unassignable reads carries intronic sequence to exploit. (The other ~78% are cleanly spliced — for
those, intron-retention offers nothing; they remain on the proven floor.)

## Factor 2 — is the intronic sequence actually distinguishing?

All 972 Tied reads map to **one family — the GWFAM10 8-copy tandem array** (NC_073228.2 ~144.93–145.34 Mb).
Comparing its copies:

- **spliced (exon) pairwise identity: median 95.6%**
- **genomic (intron-included) pairwise identity: median 95.6%**

**The introns are NOT more divergent than the exons.** The common assumption "introns evolve faster" does **not**
hold for these recent tandem duplicates — the whole locus (exons + introns) duplicated recently and diverged
uniformly. (Note the array also contains ~40 bp-offset near-identical copy pairs — copies 0/1 and 5/6 — which are
effectively one locus and stay tied no matter what.)

## What this means — a real but limited lever

- Intron-retention does **not** give a *divergence* bonus here (introns aren't faster-evolving). Its only value
  is a **coverage** bonus: a retained intron extends a read into ~4.4%-divergent extra sequence, giving a
  read that currently spans no exonic distinguishing column more chances to span one.
- So the realistic rescue is **a fraction of the 22%** — the intron-retaining reads whose extended span happens
  to cross a position distinguishing the *specific* copies they are tied between. Reads tied between the
  near-identical (40 bp-offset) duplicates stay tied (their introns are identical too).
- The exact number requires **re-running assignment with intronic positions included** as PSV columns (the
  actual lever), not just measuring retention.

## Honest verdict on "have we exhausted the unassignable?"

- ~**78% of Tied reads are cleanly spliced** → nothing intronic to exploit → they sit on the **proven** K=0
  floor (information-theoretically unassignable from spliced RNA).
- ~**22% retain an intron** → a real but *partial* lever, worth ~a fraction of that once the intron-divergence/
  coverage factor is applied, and only realizable by adding intronic PSV columns to the gate.
- The surprise that **introns here are no more divergent than exons** means the lever is weaker than the
  textbook intuition suggests — on GGO's recent tandem duplicates the whole locus is uniformly diverged, so the
  intron is just "more of the same sequence," not a privileged faster-evolving signal.

**Bottom line:** the unassignable mass is ~78% genuinely irreducible (spliced, on the proven floor) and ~22%
*potentially* touchable by intronic sequence — but with a coverage-not-divergence mechanism, so the achievable
rescue is modest and bounded. The principled endpoint (abstain with a certificate) remains correct; intron-PSV
assignment is a bounded improvement to *try*, not a floor-breaker.

## Implemented under a flag (`RUSTLE_INTRON_PSV=1`) — and measured

To attack the unassignable with every available lever, the intron-PSV path is now implemented:
`copy_assign_pipeline::discover_intron_psvs` aligns each copy's FORWARD genomic span (exons + introns) vs copy[0]
(poasta ≤20 kb, minimap2 above), keeps the substitution columns whose position is **intronic**, and appends them
to the family profile. A read that retains an intron already carries M-aligned bases there, so the existing
per-read CIGAR sweep fills them and the significance gate uses them as extra distinguishing columns — no other
change. Threaded as `Option<&GenomeIndex>` through `assign_family[_detailed]`/`build_family_profiles`; **default
OFF = byte-identical** (env unset → `discover_intron_psvs` not called). TDD: `intron_psv_finds_a_divergent_intron_
column_when_exons_are_identical`; full lib suite 688 green.

**Measured rescue (GWFAM10 region, 8,783 reads):**
| | tied | assigned |
|---|---|---|
| OFF (exon-only) | 5267 | 3067 |
| ON (`RUSTLE_INTRON_PSV=1`) | 5261 | **3094** |

**+27 reads rescued (≈0.5% of the unassignable mass; 6 from strictly-Tied, 21 from ambiguous).** Exactly as the
introns-not-more-divergent finding predicted: the lever is real and now demonstrable, but the gain is sub-1%
because GGO's K=0 copies are uniformly diverged (the intron is "more of the same," not a faster signal). The K=0
floor is near-irreducible even with intronic sequence; abstain-with-certificate remains the right default, and
the flag exists to show the floor was attacked with the last RNA-intrinsic lever.


---

## BAYESIAN_POSTERIOR_ZONE

# Bayesian posterior + consistent zone — localizing the unassignable (2026-06-29)

The hard gate emits `assign / abstain (Tied)`. The Bayesian complement keeps the **per-copy posterior** instead
of collapsing it, which **localizes even an unassignable read to a ZONE** — the set of copies (and, for a
co-located family, the genomic sub-region) it is compatible with — with a probability distribution rather than a
bare flag. This is the soft/Bayesian reading of the same object: the frequentist gate reports the integral
support; the posterior is the fractional (LP-relaxed, Theorem 7) optimum over that support.

## What it computes

- Per read, `assign_read` now returns `posterior` = `softmax(logl)` over the candidate copies (a UNIFORM prior —
  likelihood-normalized). For an *assigned* read it is ~one-hot at `best_copy`; for a *Tied* read it spreads over
  the consistent zone.
- The binary (`copy_assign --posterior`) writes `<out>.posterior.tsv`:
  `read_name, family_id, status, n_consistent, zone_chrom, zone_start, zone_end, posterior` — where the **zone**
  is the genomic extent of the copies above a 0.01 posterior floor and `posterior` is `copy:prob,...` sorted by
  probability. An informative prior is selectable: `RUSTLE_POSTERIOR_PRIOR=abundance` weights by the EM copy
  abundance (the natural home for a DNA **parCN** prior later); default uniform.
- **Default OFF = byte-identical** (the field is populated but only emitted under `--posterior`). Full lib suite
  688 green.

## Result on GWFAM10 (8,783 reads — where the unassignable mass lives)

Zone size (number of consistent copies) for the **5,267 Tied** reads:

| consistent copies (zone) | tied reads | |
|---|---|---|
| **2** | **3,382 (64%)** | localized to a tight 2-copy zone |
| 3 | 525 | |
| 4 | 15 | |
| 5 | 186 | a ~92 kb sub-region of the array |
| 6 | 2 | |
| 7 | 28 | |
| **8** | 1,118 (21%) | fully ambiguous (whole array) |

- **~64% of "unassignable" reads are actually confined to just 2 copies** — the binary `Tied` flag was hiding
  that most of them are nearly localized; the posterior makes the residual uncertainty explicit and small.
- Example Tied row: `... tied 5 NC_073228.2 144930560 145022553 2:0.200,3:0.200,4:0.200,0:0.200,1:0.200` — the
  read is excluded from copies 5–7 (the distal array members) and spread uniformly over copies {0,1,2,3,4} in a
  92 kb zone. Assigned rows are one-hot: `... assigned 1 ... 144970572 144970572 3:1.000`.

## Why this is the right framing

- It **respects the K=0 floor**: within a truly-identical sub-zone the posterior is flat (uniform), i.e. the
  information limit shows up *as* "posterior = prior", not as a fabricated pick. It never breaks irreducibility.
- It is **not `1/k` as a hard call** (which Canzar dislikes) — it is the honest *posterior*, reported as a
  distribution over a zone, with the option of an informative (abundance / parCN) prior.
- It turns "unassignable" into "localized to zone Z with distribution π" — strictly more informative, and it is
  the principled re-entry point for the DNA copy-number prior (parCN) that would sharpen the zone.

Reproduce: `copy_assign --posterior [...]` → `<out>.posterior.tsv` (add `RUSTLE_POSTERIOR_PRIOR=abundance`).


---

## POA_CORE_COMPLETION

# POA-core completion: reaching loosely-related paralogs the read-conflict graph misses (2026-06-29)

**The two-tier family definition.** The read-conflict graph answers *where is assignment needed?* — it links only
copies a read **confuses** (empirically down to ~87% identity; below that reads resolve the copies and raise no
edge). That deliberately excludes **loosely-related paralogs** (divergent copies at other loci). The POA-core
completion is the second tier: *given that read-conflict has determined a family*, complete its roster with the
divergent homologs by sequence — so the expensive homology step runs only **when needed**, not genome-wide.

**Mechanism (`family_detect::poa_core_completion_adds`).** For each read-conflict family, attach any genome rep
whose **contiguous POA core** to a family member clears `t_core` (0.13 — a ≥13%-length contiguous identical block,
i.e. a shared conserved exon even if the flanks diverge heavily). Bounded and seeded:
- the **minimizer-LSH prefilter** (`candidate_pairs`) restricts grading to homologous candidates;
- only pairs with **exactly one endpoint in a conflict family** (the other a free rep) are graded — free-vs-free
  pairs would be NEW families, not seeded by read-conflict, so they are skipped;
- the core is confirmed with the **linear-memory longest-common-substring** metric (the documented faithful
  equivalent of the poasta "longest ungapped equal run") rather than poasta itself — O(n) per pair, so it stays
  cheap. (Using poasta here was the first cut and took ~99 min genome-wide; the LCS swap fixes it.)
- a free rep matching several families attaches to the one with the strongest core.

**Wiring.** `DenovoConfig.complete_poa_core` (default **false** → no-op, byte-identical), exposed as
`gw_family_catalog --complete-core` (cross-chrom path). Added copies are appended after the conflict core and
logged. TDD: `poa_core_completion_attaches_a_divergent_paralog_at_a_new_locus` (attaches a divergent paralog that
shares the conserved core; rejects an unrelated rep and free-vs-free pairs). Full lib suite green.

## Real-data measurement (cross-chrom, GGO_mm.bam)

- **Read-conflict cross-chrom catalog:** 265 families / ~869 copies (reach ceiling ~87% identity).
- **+ POA-core completion (`--complete-core`, t_core=0.13):** **+55 copies attached to 35 families** → 265
  families / 924 copies (baseline 869 + 55, the diff is exactly the 55).
- **Reach of the attached copies (the honest result):** measured against their family core — asm20 54/55 aligned
  median 100% min **95%**; sensitive `-k11 -w5` 55/55 aligned median 100% min **93%**. So the 55 adds are
  **NEAR-IDENTICAL (93–100%), NOT loosely-related.**

## What this means (the important finding)

The mechanism is correct (the TDD test attaches a divergent-flank paralog that shares a conserved core), but on
**real GGO the nucleotide completion does not add loosely-related paralogs** — it recovers near-identical copies
the read-conflict graph **missed for COVERAGE reasons** (too few reads to clear the `min_reads≥3` conflict-edge
floor), which is a useful completeness gain but a different thing.

**Why no loosely-related adds, and it's illuminating:** read-conflict ALREADY links divergent-flank paralogs via
their **conserved-exon reads** (a read over the shared exon ties the two loci even if the flanks diverge) — that
is exactly why read-conflict reaches ~87% *overall* identity. So the conserved-exon-sharing divergent paralogs
are *already in the families*; the only copies left for a nucleotide completion are the near-identical ones missed
for coverage. A genuinely loosely-related paralog (divergent throughout, no conserved-exon read overlap) does NOT
clear `t_core=0.13` on the exact-LCS core either (30% scattered divergence breaks any long identical run) — it is
reachable only at the **protein level** (the tier deferred in this build).

**Net:** all NUCLEOTIDE family definitions on real GGO top out at ~81–87% identity:
- read-conflict ~87% (via conserved-exon ties);
- `--refine` validates within (asm20 ~81% floor; 422 copies, 2 below 87%);
- `--complete-core` adds 55 near-identical (93–100%) missed copies — completeness, not extended reach.

**Genuinely loosely-related paralogs need the PROTEIN tier** (6-frame ORF → mmseqs, fident≥0.50; reaches ~50%
protein / ~70% genomic) — already implemented as the opt-in `--protein-tail` for `--refine`, and the natural next
lever for the completion pass if loosely-related reach is the goal. (Or DNA/SEDEF for ancient paralogs.)

## How the tiers relate (corrected)
- **read-conflict** (tier 1): identifiability — which loci need assignment; ~87% via conserved-exon ties.
- **`--complete-core`** (nucleotide, to new loci): recovers near-identical copies missed for coverage; does NOT
  extend reach on real data because the divergent-but-core-sharing copies are already in tier 1.
- **`--refine` / `--protein-tail`** (homology within / protein): the protein tier is the only one that reaches
  genuinely loosely-related (70–87%+) paralogs; nucleotide cannot.


---

## CONFLICT_EDGE_UNIFICATION

# Unifying the family-definition edge with the assignment gate (2026-06-29)

**Motivation.** The family-as-graph definition and the read→copy assignment were two principled objects that used
*different* criteria for the same underlying question ("can this read tell these two copies apart?"):

- **Family edge** (`read_conflict.rs::de_tied`): a read links two loci iff `|de_a − de_b| ≤ delta` (a fixed
  `delta = 0.005`) and `max(de_a, de_b) ≤ de_max`. The `delta` is the last hand-set constant in the definition.
- **Assignment gate** (`copy_assign.rs`): a read is *tied* between two copies iff `min_p = ε^δ ≥ α` (Theorem 4),
  with `δ` = the number of distinguishing columns the read spans — an error-model-derived, threshold-free test.

The F4 finding (raw read-conflict graphs are error-inflated, colouring ≈3× the true copy count) is the empirical
case that the conflict edge should use the *same significance test* as assignment, not a fixed tolerance.

**The unification.** A read counts as conflict evidence between two loci iff it **cannot significantly
distinguish them** under the assignment gate's own criterion. With `m_x = de_x · aln_len_x` the mismatch count to
locus `x`, the excess `δ = |m_a − m_b|` is the per-read distinguishing-column proxy, and the read ties iff

  ε^δ ≥ α            (exactly the gate's `min_p ≥ α`, Theorem 4)

keeping the `de_max` quality floor (both alignments must genuinely fit). The arbitrary `delta = 0.005` is gone;
the tie boundary is now the error model (`ε`, default `e/3 ≈ 0.001`) and the significance level (`α`, default
`1e-3`) — the *same two numbers* the assignment gate uses. One IsoCon real-vs-error criterion now governs both
the family edge and the read assignment.

**Implementation (`src/rustle/vg_family/read_conflict.rs`).**
- `Placement` gains `aln_len` (aligned-block length = #M/=/X CIGAR columns), computed at construction in
  `denovo_pipeline.rs`.
- `ConflictParams.sig: Option<(eps, alpha)>` — `None` (default) = the legacy `de_tied`, so **OFF is
  byte-identical**; `Some` = `sig_tied`. Env-gated: `RUSTLE_CONFLICT_SIG=1` (+ `RUSTLE_CONFLICT_EPS`,
  `RUSTLE_CONFLICT_ALPHA`).
- `tied()` dispatches; `conflict_edges`/`family_mapq0_support` call it. The `de-tie ⊆ AS-tie` audit invariant is
  untouched.
- Tests: `sig_criterion_ties_ambiguous_resolves_distinguishing` (the unification at unit level),
  `sig_off_default_is_byte_identical_to_de_tied` (default ships OFF). Full lib suite green (685 tests).

**Why it lands with Canzar.** It removes the last similarity-style constant from the family definition: the
boundary is now a property of the data (the per-base error rate) and a stated significance level, identical to
the assignment gate — one clean combinatorial criterion, no tuned thresholds.

## Real-data measurement (significance edge ON vs OFF; genome-wide `gw_family_catalog` on `GGO_mm.bam`)

| catalog | families | copies |
|---|---|---|
| **OFF** (legacy `de_tied`, fresh run, this binary) | **81** | **205** |
| **SIG** (`RUSTLE_CONFLICT_SIG=1`, ε=1e-3, α=1e-3) | **71** | **176** |
| existing baseline (older binary, for reference) | 82 | 207 |

- **Byte-identical OFF confirmed end-to-end:** the fresh OFF run (81/205) matches the prior baseline (82/207)
  to within ±1 family — run-to-run noise from threaded rep assembly, *not* the struct change (the new `aln_len`
  field is ignored by `de_tied`; the unit test `sig_off_default_is_byte_identical_to_de_tied` pins the criterion).
- **The significance edge is a principled REFINEMENT.** Analytically, for equal-length placements at the default
  ε/α, `sig_tied ⟹ de_tied` (the tie boundary `δ ≤ 1` is `|de_a−de_b| ≤ 1/L ≤ 0.005`), so SIG edges are a
  **subset** of de-tie edges — it can only shrink or split families, never create them.
- **Effect: a modest, defensible narrowing — 81→71 families (−12%), 205→176 copies (−14%).** The catalog drops
  the ~1/8 of families held together only by *marginally*-tied reads (within the old 0.5% tolerance but
  significantly **resolvable** under the gate). The largest families are **preserved** (9-copy GWFAM58, 7-copy
  GWFAM10 / GWFAM67), and the size distribution is essentially unchanged (mostly 2-copy paralog pairs).
- The exact per-copy coordinate overlap between the two runs (~78%) is limited by rep-assembly non-determinism
  across separate runs, not by the edge criterion — the analytic subset relationship above is the rigorous
  refinement claim; the family/copy counts are the robust empirical headline.

**Interpretation / the choice for the advisor.** OFF defines a family as loci linked by reads within ~0.5%
divergence; SIG defines it as loci a read **cannot significantly resolve** — i.e. the genuinely
collapsed/assignment-hard set, under the *same* criterion the assignment gate uses. SIG is the more principled,
threshold-free object (no `0.005`), at the cost of excluding ~12% of near-but-resolvable pairs. It ships **OFF
by default** so nothing changes silently; flip with `RUSTLE_CONFLICT_SIG=1` once the definitional choice is made.

## GSTM_REAL_COPIES — flagship worked example (real gorilla GST-Mu cluster)

A worked, corroboratable example on real gorilla Iso-Seq: the GST-Mu cluster on NC_073224.2.
Defense-slide figure: `bench/slides/gstm_real_copies.png` (regenerate with `bench/make_gstm_copies_fig.py`).

**Reproduce:**

```bash
# in /home/juanfra/winloci_scratch (GGO_mm.bam + GGO.fasta)
copy_assign --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 \
  --dump-psv --gtf --out gstm_vg
python bench/make_gstm_copies_fig.py    # reads gstm_vg.psv_{reads,copies,cols}.tsv
```

**What it corroborates (for a skeptic):**

1. **The 3 called copies ARE annotated paralogs.** copy 0/1/2 land exactly on the annotated GST-Mu genes
   **GSTM3 / GSTM5 / GSTM1** (positions from the copy tids vs `annot_intervals.tsv`) — real, textbook recent
   gene duplications anyone can check.
2. **They are genuinely similar** — grouped into ONE family only because their exon-sum identity clears the
   homology gate (≥ 80%). The clean pair GSTM5↔GSTM1 measures **83%**; GSTM3 is the divergent member (agrees
   with the others at only ~15% of PSV sites — consistent with GSTM3 being the outlier Mu paralog).
3. **The reads sort into 3 clean copy-blocks by their PSV alleles** (the SDA / Vollger read×PSV matrix, on
   real data): **411 PSV columns** distinguish the copies; the read block below each copy-consensus matches it.
4. **The assignment is verifiable.** Where minimap2 maps a read *uniquely* (MAPQ > 0), our PSV-based copy call
   agrees with it **1341/1341 = 100%** — so on the MAPQ-0 reads where the aligner gives up, the method does
   what the aligner would if it could. **2654/2673 reads assigned; 16 certified TIED** (abstain, never 1/k).
   Every copy carries hundreds of **private SUNs**, so a single read over a SUN deterministically pins its copy.

**Takeaway.** The variation graph is not a schematic: each PSV column is a bubble, each copy is a path, and the
2673 real reads visibly thread onto the 3 copy-paths — assigned when a distinguishing bubble is significant,
tied when none is.

