# Design: Whiteboard Explanation of Pileup + EM for Multi-Copy Gene Families

**Date:** 2026-05-28  
**Purpose:** Explain the VG mode pileup + EM method to an advisor in an upcoming meeting  
**Status:** In progress — visual prototypes built, needs final polish and implementation plan

---

## Context

The Rustle VG mode (`--vg`) assembles transcripts from multi-copy gene families by reweighting multi-mapping reads using EM. The advisor is a computer scientist who applies CS to biology. He has:
- Deep familiarity with BAM format and nucleotides
- Rejected HMM-based explanation ("how do you convert nucleotides to an HMM?")
- Challenged SNP-based explanation ("if there's a SNP there's no ambiguity — that's not the problem")
- His actual question: *"a long read maps equally to 5 locations, same alignment score, no SNP — then what do we do?"*

---

## Agreed Approach

**Format:** Whiteboard walkthrough (advisor works visually and literally)  
**Style:** Toy example + real data, every number traceable to a BAM record  
**Key principle:** No abstraction — the advisor must be able to point to any number and ask "where does this come from?"

---

## 7-Step Whiteboard Arc

### Step 1 — The problem: two nearly-identical copies
Draw DAZ1 and DAZ3 on the gorilla Y chromosome side by side.

**Real data (GGO.bam, NC_073248.2):**
- DAZ1: 42,783,133–42,859,657 (− strand), 13 exons
- DAZ3 (LOC129530216): 42,879,918–42,945,552 (+ strand), 13 exons, inverted

| Copy | MAPQ ≥ 30 (unique) | MAPQ = 0 (multi-mapper) | Without EM |
|------|-------------------|------------------------|------------|
| DAZ1 | **167 reads** | 42 reads | transcript assembled ✓ |
| DAZ3 | **0 reads** | 216 reads | **gene invisible ✗** |

DAZ3 has zero uniquely-placed reads. Standard tools discard all 216 multi-mappers. DAZ3 does not exist in the output without EM.

### Step 2 — A real ambiguous read
Read `m64076_221110_210557/50726963/ccs` (IsoSeq CCS, ~2,200 bp):
- Maps to DAZ1 at 42,766,268 with MAPQ=0 (12 junctions, forward CIGAR)
- Maps to DAZ3 at 42,903,444 with MAPQ=0 (12 junctions, mirrored CIGAR — opposite strand)
- Same alignment score at both. No SNP. The advisor's exact scenario.

### Step 3 — Signal 1: junction compatibility (always a tie in this case)
Count what fraction of the read's splice junctions appear in each copy's junction list:
- DAZ1: 12/12 matched → compat = **1.00**
- DAZ3: 12/12 matched → compat = **1.00**

Signal 1 is a tie. Move to Signal 2.

### Step 4 — Signal 2: pileup depth as prior
Before EM runs, each read contributes its fractional weight (`nreads_good += read.weight`) to the bundle's junction depth:
- Unique read (NH=1, MAPQ≥30) → weight = **1.0**
- Multi-mapper (NH=2, MAPQ=0) → weight = **0.5** per locus

```
DAZ1 depth = 167 × 1.0 + 42 × 0.5 = 188
DAZ3 depth =   0 × 1.0 + 216 × 0.5 = 108
```

DAZ1's 167 unique reads each contribute full weight = stronger prior.

### Step 5 — Score and normalize
```
score(DAZ1) = (1.00 + 0.01) × ln(188) = 1.01 × 5.24 = 5.29
score(DAZ3) = (1.00 + 0.01) × ln(108) = 1.01 × 4.68 = 4.73

weight(DAZ1) = 5.29 / (5.29 + 4.73) = 0.528  →  53%
weight(DAZ3) = 4.73 / (5.29 + 4.73) = 0.472  →  47%
```

### Step 6 — Why EM and not naive 50/50?
EM gives **53/47** in this case (nearly symmetric because DAZ3 has many multi-mappers giving it partial fractional depth).

To show when EM clearly beats naive, use hypothetical: DAZ3 has only 5 unique reads instead of 0:
```
DAZ3 depth = 5 × 1.0 + 5 × 0.5 = 7.5
score(DAZ3) = 1.01 × ln(7.5) = 2.33
weight(DAZ1) = 5.29 / (5.29 + 2.33) = 69%  ←  EM
naive:         50/50                         ←  ignores evidence
```

| | Symmetric expression (real DAZ) | Asymmetric (hypothetical) |
|---|---|---|
| Naive 50/50 | ✓ happens to be close | ✗ wrong — ignores evidence |
| Pileup EM | ✓ honest uncertainty (53/47) | ✓ correctly biased (69/31) |

### Step 7 — The three-zone picture and where IsoSeq + VG help

**Zone 1 — Diverged copies** (UTY/UTX, KDM5D/C, DDX3Y/X in GGO.bam): reads map uniquely, MAPQ=60 everywhere. No multi-mapper problem. EM not needed. (Confirmed: zero MAPQ=0 reads at these X-Y pairs.)

**Zone 2 — Similar copies, different expression** (DAZ1/DAZ3): some multi-mapping. EM + pileup prior breaks the tie using expression asymmetry. IsoSeq helps here because full-length reads reach **UTR sequences** (which diverge faster than CDS) and span more exons per read.

**Zone 3 — Truly identical copies, equal expression**: EM gives honest 50/50. Cannot do better from sequence alone. However:
- **IsoSeq + VG family graph** pushes many Zone-3-seeming cases into Zone 2 by: (a) UTR coverage revealing copy-specific SNPs, and (b) propagating expression priors from well-characterized family members to uncharacterized ones.
- Residual Zone 3 (all copies identical across entire transcript including UTRs, all equally expressed, none with unique reads) is genuinely unresolvable by sequence. The correct answer is 50/50 — honest uncertainty quantification, not a failure.

---

## Key Technical Facts (for Q&A)

- **EM implementation:** `run_pre_assembly_em_inner` in `vg.rs:2835`
- **Score formula:** `(junction_compat + 0.01) × ln(junction_depth).max(1.0)`
- **Pileup diagnostics:** `build_pileup_diagnostics` in `vg.rs:2388` — scans exonic positions only (not intronic)
- **Junction compatibility:** `junction_compatibility` in `vg.rs:2165` — fraction of read junctions matching bundle's junction list
- **Weight updates:** `nreads_good += r.weight` in `bundle.rs:2018` — multi-mappers contribute fractional weight pre-EM
- **Convergence:** scores are fixed across iterations (junction_depth is pre-computed), so heuristic EM converges in 1–2 iterations

## What Was Searched (GGO.bam)

All major X-Y homologous pairs were checked for multi-mappers:

| Pair | Y unique | X unique | MAPQ=0 at Y | MAPQ=0 at X | Multi-mapping? |
|------|----------|----------|-------------|-------------|----------------|
| DDX3Y/DDX3X | 1596 | 1358 | 0 | 0 | No |
| KDM5D/KDM5C | 442 | 221 | 0 | 0 | No |
| USP9Y/USP9X | 184 | 835 | 0 | 0 | No |
| UTY/KDM6A | 166 | 208 | 7 (→PSMA6) | 0 | Not Y↔X |
| DAZ1/DAZ3 | 167 | 0 | 42 | 216 | Yes (Y-internal) |

**Conclusion:** No X-Y pair shows inter-locus multi-mapping in this BAM. Sequences have diverged enough to map uniquely. DAZ1/DAZ3 (Y-internal tandem) is the best real example. A hypothetical with stronger expression asymmetry is needed to demonstrate EM's discriminating power over naive 50/50.

---

## Visual Prototypes Built

All HTML files are in `.superpowers/brainstorm/` on the Desktop.

1. `format-question.html` — format selection (slides/whiteboard/figure/written)
2. `approaches.html` — 3 pedagogical approaches
3. `design-arc.html` — first arc (SNP-based, superseded)
4. `design-arc-v2.html` — revised arc (junction + pileup depth, approved)
5. `toy-example.html` — toy example with made-up numbers
6. `real-daz-example.html` — real DAZ1/DAZ3 data, first version
7. `daz-em-final.html` — real data + hypothetical side by side ✓ **approved design**
8. `zone3-frameworks.html` — copy-number prior / methylation / Hi-C
9. `zone3-longread-vg.html` — IsoSeq UTR + VG family graph answer ✓ **latest**

---

## Open Questions / Next Steps

1. **Finalize the whiteboard script** — write a step-by-step spoken narrative for each whiteboard panel, anticipating the advisor's literal questions
2. **Decide on hypothetical numbers** — confirm the 167 vs 5 scenario is biologically plausible to claim
3. **Prepare Q&A answers** — "why ln()?", "why 0.01?", "what if copies have different number of junctions?", "what about Zone 3 truly unresolvable?"
4. **Optional: run Rustle in VG mode on DAZ region** to show actual EM output rather than manual calculation
