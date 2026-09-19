# The MERGED-LOCUS layer: real fusions recorded as dual membership, outside the partition

**What it is.** A node that is a genuine fused transcript belongs to TWO parent families. The lattice levels
must stay node PARTITIONS or T1 (nesting) and T1′ (laminar) fail, so dual membership is recorded as a
separate **layer** — a cover sitting on top of the partition, never a level in the nesting chain.

Artifact: `bench/merged_loci_layer.tsv`.

## Admission rule (evidence-based, not the blunt rule register 844 refused)

A node enters the layer iff all three hold:
1. **curated** as a readthrough by RefSeq (`description=` contains "readthrough"; 209 records, §6m2);
2. its **fusion junction** carries **>= 3 PRIMARY reads** (§6n8: supplementary support is 0.1% overall, so
   these are contiguous molecules, not split alignments);
3. it is **not secondary-dominated** (§6n9) — its defining junction is not a multimapping echo.

This is exactly the discrimination register 844 said was missing: it admits `PKD1P6-NPIPP1` (110 primary
reads, canonical GT-AG, §6m1) and excludes `PKD1P4-NPIPA8` (10 primary / 1,150 secondary).

## Content (NPIP/TBC1D3 region, L3 at the new 0.995 cut)

| merged locus | primary | secondary | own L3 component | families bridged | top bridges |
|---|---|---|---|---|---|
| **PDXDC2P-NPIPB14P** | 726 | 1 | **1 (singleton)** | **14** | LOC100288162(13), LOC100190986(11), CLN3(6), BOLA2(6) |
| **PKD1P6-NPIPP1** | 110 | 3 | **1 (singleton)** | **6** | LOC100288162(21), CLN3(2), MIR6511B1(1), NPIPB1P(1) |
| BOLA2-SMG1P6 | 60 | 51 | 11 | 4 | LOC100190986(4), SMG1(1), SMG1P7(1) |
| SLX1B-SULT1A4 | 120 | 4 | 8 | 1 | SULT1A1(1) |

⭐**Raising L3 to 0.995 already made the two NPIP-side fusions SINGLETON components.** The records that
held NPIP's DNA-certificate boundary for the whole of §6m1-§6m3 are now isolated by the cut itself; the
layer records what they bridge rather than letting them merge it.

## Why this is the right shape

- **The partition is untouched**, so T1, T1′ and T2 keep their proofs. §6m3 measured the alternative —
  cutting the node in two — and it was safe but harmful (fresh-arm precision 0.074 → 0.004), because both
  halves stay inside the family's component and short pieces become hubs.
- **It is small and auditable**: 4 nodes here, each with its read evidence attached.
- It answers "which families does this locus belong to" without asking the clustering to represent it.

## Caveats

- Scoped to the NPIP/TBC1D3 lattice graph; a genome-wide layer needs the same census on a full catalog.
- `BOLA2-SMG1P6` is borderline (60 primary / 51 secondary, secondary fraction 0.46) — it passes only
  because the threshold is 0.50. Worth flagging rather than trusting.
- The bridged-family counts come from L2 neighbours, so they include co-duplicated SD neighbours
  (CLN3, the LOC lncRNAs) that §6m7 showed are not true family partners. **Read the counts as "what this
  locus touches", not "what it is a member of".**
