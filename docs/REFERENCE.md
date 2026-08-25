# REFERENCE — glossary, data locations, and a worked example

Three short reference documents merged on 2026-08-25. Each was too small to justify its own file,
and all three answer *"what does this mean / where is it"* rather than recording a result.

---

## Glossary — what each term denotes

<!-- merged from REFERENCE.md on 2026-08-25 -->

The live thesis code lives in `src/rustle/vg_family/`. This note fixes the
meaning of a few overloaded words that are used in different senses across the
module and in the companion formal notes.

> **Rule of thumb.** When you see the word **"locus"** in `vg_family`, ask
> whether it means a **gene locus** (a splice-junction community of isoforms)
> or a **physical genomic span** (a `(chrom, start, end)` interval used for the
> ≥2-distinct-loci certificate). They are not the same object.

### Canonical definitions

#### Isoform

A distinct spliced transcript observed in the RNA data. Operationally it is a
specific intron chain — a ordered list of exon intervals. In code:

- `DenovoTranscript` (`family_detect.rs`) is one assembled isoform.
- `FamilyPath` (`layer2.rs`) is one recovered isoform inside a multi-copy
  family, annotated with its copy and source (`Native`, `Transferred`,
  `PsvLinked`).

An isoform is **not** a gene; one gene can produce many isoforms.

#### Gene locus

The set of isoforms that are alternative splice products of the **same gene**.
In the pipeline a gene locus is collapsed from the raw assembled isoforms by
**shared splice junctions** (`family_detect::collapse_loci`): if two isoforms
use the exact same intron `(donor, acceptor)` coordinates on the same
chromosome, they belong to the same gene locus. One representative isoform is
kept for downstream family detection.

A gene locus is therefore defined by splicing, not by span overlap. This is
intentional: dense genomes would chain many unrelated genes into one bogus
locus if we used span overlap alone.

#### Copy / member

One gene locus that belongs to a multi-copy gene family. Copies are homologous
loci. In the copy-assignment stage, reads are assigned to a specific copy.

#### Family

A set of copies (gene loci) that are homologous. In code:

- `detect_edges` (`family_detect.rs`) builds homology edges between locus
  representatives using POA contiguous-core coverage.
- `decompose_families` (`family_split.rs`) refines raw connected components
  into cohesive families (γ-quasi-clique refinement).

#### Physical span / distinct locus

The `(chrom, start, end)` interval of a gene or copy. `family_definition.rs`
uses this notion for the multi-copy certificate: a family must contain at least
**two distinct physical loci** (`distinct_loci`), where "distinct" means the
spans do not reciprocally overlap by ≥ `LOCUS_OVERLAP` (50%).

This is a genomic-coordinate check, independent of splicing. It answers the
question "do these members sit at different places in the genome?" rather than
"are these isoforms of the same gene?".

### Why two notions of locus are necessary

- **Gene-locus collapse** (shared junctions) prevents a single gene with many
  alternative splice isoforms from being mistaken for many paralog copies.
- **Physical-span distinct-loci count** prevents two nested or overlapping
  annotations of the same genomic gene from being counted as two copies of a
  family (e.g., MAGEA9 nested inside a LOC entry).

### Known limitation: isoforms with no shared junction

The shared-junction rule can be too strict: two genuine isoforms of the same
gene may have **disjoint intron sets** (for example, alternative first/last
exons that use completely different splice sites). In that case they are
currently treated as separate gene loci and may enter family detection as false
paralogs. The span-aware recovery in `family_detect::collapse_loci_span_aware`
mitigates this by additionally merging locus representatives that are
span-overlapping and either strongly contained or highly homologous.

### Mapping to the formal notes

- `DNA_FAMILY_DEFINITION_FORMAL.md` and `PROTEIN_FAMILY_DEFINITION_FORMAL.md`
  use **locus** in the physical-span sense (NCBI RefSeq gene loci).
- `family_definition_formal.md` uses **expressed locus** in the gene-locus
  sense (de-novo assembled isoforms collapsed by shared junctions).
- This note aligns the code vocabulary with the formal notes: a **gene locus**
  in code corresponds to an expressed locus in the RNA formal note, and the
  **physical span** corresponds to the DNA/protein note's locus.

---

## Where the data lives

<!-- merged from REFERENCE.md on 2026-08-25 -->

The large reference assemblies and alignment data live on a separate 1.8 TB ext4 disk. They are
not stored in the Rustle repository and must never be deleted or modified as part of a benchmark
cleanup.

### Mount or reattach after a WSL restart

The disk is Windows `PHYSICALDRIVE0`, partition 2. `wsl --mount` does not survive a WSL restart.
The commands below are copied from the guarded mount instructions in `/home/juanfra/.bashrc`.
The attach/detach syntax is also documented by
[Microsoft's WSL disk-mount guide](https://learn.microsoft.com/en-us/windows/wsl/wsl2-mount-disk).

First, in an **Administrator PowerShell** window:

```powershell
wsl --mount \\.\PHYSICALDRIVE0 --partition 2
```

Then, inside WSL:

```bash
sudo mkdir -p /mnt/linuxdisk
sudo mount --bind /mnt/wsl/PHYSICALDRIVE0p2 /mnt/linuxdisk
```

Do not hard-code `/dev/sde2` or another `/dev/sd*` name: the kernel device letter can change between
sessions. Use the stable WSL-generated `PHYSICALDRIVE0p2` path.

Verify the mount before reading data or starting a large build:

```bash
mountpoint -q /mnt/linuxdisk
findmnt -T /mnt/linuxdisk
test -r /mnt/linuxdisk/home/juanfraitu/gorilla_haps/mat.fa
```

The shell startup guard sets `CARGO_TARGET_DIR` to
`/mnt/linuxdisk/home/juanfraitu/rustle_target` only while the mount is active. If the disk is absent,
the variable is deliberately unset so a large Rust build does not silently fill the Windows C:
drive.

In a managed or sandboxed session the mount may be exposed read-only even when it is normally
writable. Check `findmnt -T /mnt/linuxdisk` and `test -w /mnt/linuxdisk` before selecting an output
directory. Never attempt to remount it from an automated task. Read reference data in place and
write small durable certificates under the Rustle `bench/` tree; use `/tmp` only for disposable
intermediates.

### Ape and human assembly inventory

The following files were verified readable on 2026-08-17:

| species/haplotype | assembly file | index/status |
|---|---|---|
| human CHM13v2.0 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` | reference used by the HSA O1 audit |
| gorilla KB3781 maternal | `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/mat.fa` | `mat.fa.fai`; 225 sequences |
| gorilla KB3781 paternal | `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/pat.fa` | `pat.fa.fai`; 24 sequences |
| chimpanzee AG18354 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna` | `.fai`; 26 sequences |
| Bornean orangutan AG05252 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` | `.fai`; 26 sequences |

The corresponding optional annotation files are `GGO_genomic.gff`, `PTR_genomic.gff`, and
`PPY_genomic.gff` under `winloci_data/`. They may be used to label results after inference, but the
O1 rooting rules must not require them.

Existing gorilla minimap2 indexes include:

```text
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.asm20.rebuild.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.asm20.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi
/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi
```

The `asm20` indexes are suitable for sensitive duplicated-block discovery. Strict `asm5` flank
mapping should use `mat.fa`/`pat.fa` directly or a separately generated `asm5` index: minimap2 index
seed parameters are fixed when the index is built.

### Detach safely

After all processes using the disk have stopped, unmount the bind mount in WSL:

```bash
sudo umount /mnt/linuxdisk
```

Then detach it in Administrator PowerShell:

```powershell
wsl --unmount \\.\PHYSICALDRIVE0
```

---

## Worked example — the DAZ locus end to end

<!-- merged from REFERENCE.md on 2026-08-25 -->

**Question this answers (for the advisor):** for a real multi-copy locus, (1) how do we use the `de` tag to
decide which copy a read belongs to, with every number shown; and (2) once a read is assigned, is it accepted
into a bundle, and is max-flow applied to assemble the copy's transcripts?

**Everything below is reproducible** from the gorilla IsoSeq BAM:
`samtools view GGO.bam NC_073248.2:42700000-43000000` (the DAZ cluster on chrY).

### The locus (real, annotated)
Two paralogous copies, ground truth known:
- **DAZ1** — NC_073248.2:42,783,133–42,859,657, **− strand** (reads align reverse).
- **DAZ3 / LOC129530216** — 42,879,918–42,945,552, **+ strand** (reads align forward).

They are near-identical inverted amplicons, so most reads **multi-map** to both copies (MAPQ 0). The job: decide,
per read, which copy it came from — then assemble each copy.

---

### Step 1 — Read-level copy assignment with `de`

minimap2 puts three relevant numbers on each alignment of a read:
- `NM:i:` — raw edit distance (substitutions **+ every indel base**).
- `de:f:` — **gap-compressed** per-base divergence (each indel counts **once**, not per base).
- `AS:i:` — alignment score.

We work in **events** = `de × aligned_length` (the gap-compressed count of distinguishing differences). The read
belongs to the copy with **fewer events**; the discriminant is `ΔEvents = events(sibling) − events(thiscopy)`.

#### Three real reads (every value straight from the BAM)

**(a) Decisive — read …46596327** (all metrics agree → DAZ1)

| | `de` | `NM` | `AS` | aligned len | **events = de·len** |
|---|---|---|---|---|---|
| DAZ1 | 0.0004 | 1 | 1686 | 2359 | **1** |
| DAZ3 | 0.0063 | 14 | 1519 | 2205 | **14** |

`ΔEvents = 14 − 1 = +13` → **DAZ1** (13 distinguishing events fewer at DAZ1).

**(b) The one that proves why `de`, not `NM` — read …53870980**

| | `de` | `NM` | `AS` | aligned len | **events = de·len** |
|---|---|---|---|---|---|
| DAZ1 | 0.0278 | **88** | 1588 | 2865 | **80** |
| DAZ3 | 0.0229 | **565** | 1381 | 2370 | **54** |

Raw `NM` says **DAZ1** (88 ≪ 565) — but DAZ3's 565 is **mostly indel slippage** (homopolymer 1-base insertions).
Gap-compress it and DAZ3 is the **better** fit: `ΔEvents = 54 − 80 = −26` → **DAZ3**. Raw edit distance would
**mis-assign this read**; `de` corrects it. *(Across this locus, `NM` is ~87% indel bases, and **8 of 177**
two-copy reads are assigned to a different copy under `de` than under `NM`.)*

**(c) Non-identifiable — read …24709554** (the honest floor)

| | `de` | `NM` | `AS` | aligned len | **events** |
|---|---|---|---|---|---|
| DAZ1 | 0.0000 | 0 | 2518 | 2518 | **0** |
| DAZ3 | 0.0000 | 0 | 2518 | 2518 | **0** |

`ΔEvents = 0`. The read covers a region **identical** between the copies — there is no evidence. We do **not**
guess; it stays split by the prior (Step 2).

---

### Step 2 — From evidence to a soft assignment (the EM)

We never hard-assign. Each read's per-copy events become a **likelihood**, and a posterior **responsibility**
`γ` (the EM E-step is `posterior = softmax(log_score + log_prior)`, `src/rustle/vg.rs:2436`). Using the event
model (each distinguishing event carries log-odds `ln((1−ε)/ε) ≈ 2.94` at ε=0.05):

| read | ΔEvents | **γ(DAZ1)** | **γ(DAZ3)** |
|---|---|---|---|
| (a) …46596327 | +13 | **1.000** | 0.000 |
| (b) …53870980 | −26 | 0.000 | **1.000** |
| (c) …24709554 | 0 | **0.500** | 0.500 |

The tie (c) splits 50/50 by the prior — apportioned, never fabricated. The EM iterates (the per-copy prior
updates from the total responsibility mass) until `γ` converges. *(On real DAZ: `[VG-FP-EM] converged in 7
iter, 179 reads adjusted`.)*

`de` is the **interpretable, by-hand** form of the same evidence the production fingerprint-EM scores from the
copies' distinguishing positions — both point the same way.

---

### Step 3 — Reads accepted into bundles, then **max-flow** assembly

This is the part after assignment.

1. **Both bundles keep the read.** A multimapper sits in **DAZ1's bundle and DAZ3's bundle** simultaneously
   (initial weight `1/NH`). Nothing is discarded.
2. **The EM rewrites the read's weight in each bundle to its responsibility γ** — in place:
   `bundles[bi].reads[ri].weight = γ` (`src/rustle/vg.rs:2480`, and `:1351` for the graph-compat EM).
   - read (a): weight **1.0** in DAZ1's bundle, **0.0** in DAZ3's.
   - read (b): weight **0.0** in DAZ1's bundle, **1.0** in DAZ3's.
   - read (c): weight **0.5** in each.
3. **Each copy's bundle becomes a splice graph**, where every node/edge **coverage = Σ of the contributing
   reads' (now-reweighted) weights**. So a copy's graph carries only the evidence the EM assigned to it; the
   spillover reads (weight ≈ 0 in the sibling) add ≈ no coverage there.
4. **Max-flow extracts the transcripts.** The assembler runs network max-flow on each copy's splice graph —
   `push_max_flow_seeded_full` / `long_max_flow_seeded_with_used_pathpat` (`src/rustle/path_extract.rs:11–13`):
   it repeatedly pushes the **heaviest source→sink path** (an isoform), subtracts that flow, and repeats, so
   each copy's isoforms are the high-flow paths through its **reweighted** graph.

**Net effect:** DAZ1's graph carries the full weight of the DAZ1-decisive reads → DAZ1's isoforms are recovered
at full abundance. DAZ3's graph carries only its genuinely-assigned reads → DAZ3 is assembled at its **honest**
(low) abundance, and the ~150 spillover reads that *align* to DAZ3 but belong to DAZ1 contribute ≈ 0 flow, so
they do **not** fabricate phantom DAZ3 isoforms.

---

### The whole chain, on read …53870980 (the hard one)

```
BAM tags         de(DAZ1)=0.0278  de(DAZ3)=0.0229   (NM 88 vs 565 — misleading: indel-inflated)
  ↓ events = de·len
events           DAZ1 = 80        DAZ3 = 54
  ↓ ΔEvents = 54 − 80 = −26  (DAZ3 fits better)
EM responsibility γ(DAZ3) = softmax(...) = 1.000
  ↓ reads[ri].weight = γ
bundle weight    DAZ3 bundle: weight 1.0   |   DAZ1 bundle: weight 0.0
  ↓ coverage = Σ weights  →  splice graph
max-flow         heaviest source→sink path through DAZ3's graph
  ↓
result           this read's evidence builds a DAZ3 isoform, not a DAZ1 one
```

Every value is checkable: the tags from `samtools view`, the EM/weight step at `vg.rs:2436/2480`, the max-flow
at `path_extract.rs:11`. The point for the advisor: **the copy decision is a concrete per-read calculation on
the `de` tag (not a black box), it is soft (ties are apportioned, not invented), and the assigned weight is
exactly what the max-flow assembler integrates into each copy's transcripts.**

