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
of decisive / silver 99.9%) is the validated method; the cluster run scales it + the divergent-tail.

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
