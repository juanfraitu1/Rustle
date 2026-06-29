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
