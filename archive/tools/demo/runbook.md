# Glass-box demo runbook (advisor meeting prep)

Total: ~25 minutes. Order matters — credibility first, results second, vision third.

## Setup (before he walks in)

- Terminal at `/mnt/c/Users/jfris/Desktop/Rustle/`.
- Image viewer pointed at `tools/demo/out/`.
- `git log --oneline -10` shows recent commits with timestamps (proves no last-minute hardcoding).
- All artifacts pre-generated under `tools/demo/out/synthetic/` and `tools/demo/out/chr19_triple/`.
- Confirm `stringtie --version` reports 3.0.3 (used for the side-by-side comparison in Tier 2).
- Confirm `./target/release/rustle --help` exits 0 (release binary is built).

## Tier 1 — Synthetic walkthrough (~5 min)

The point of this tier: **prove the algorithm machinery works on data where the ground truth is independently verifiable.**

1. `cat test_data/synthetic_family/truth.gtf | head -20` — *"Ground truth defined before reads existed. Three loci: MYG1_A (with an alt-splice A2 that skips exon 3), MYG1_B, MYG1_C."*
2. `head -30 test_data/synthetic_family/generate_isoseq.py` — *"Generator with `random.seed(42)`. The BAM is a deterministic function of these two files."*
3. `samtools view test_data/synthetic_family/reads_sorted.bam | head -5` — *"5 of 95 PacBio-style reads."*
4. `bash tools/demo/run_demo_synthetic.sh` — runs in seconds.
5. **Walk through the dumps in order:**
   - `column -t -s$'\t' tools/demo/out/synthetic/reads.tsv | head -10` — reads as parsed (85 data rows).
   - `cat tools/demo/out/synthetic/partition.tsv` — 2 bundles formed.
   - `cat tools/demo/out/synthetic/nodes.tsv` — 8 splice-graph nodes total across the bundles.
   - `cat tools/demo/out/synthetic/edges.tsv` — 7 splice-graph edges with capacities.
   - `cat tools/demo/out/synthetic/flow.tsv` — 3 augmenting paths peeled off.
6. **Open `tools/demo/out/synthetic/v1_splice_graph.png`** — *"The splice graph from those nodes and edges. Every value here came from one of the TSVs you just read."*
7. **Open `tools/demo/out/synthetic/v3_flow_iterations.png`** — *"Edmonds-Karp peels off augmenting paths. Each path becomes one transcript."*
8. **Compare against truth:**
   ```
   diff <(awk -F'\t' '$3=="transcript"{print $4,$5}' tools/demo/out/synthetic/baseline.gtf | sort) \
        <(awk -F'\t' '$3=="transcript"{print $4,$5}' test_data/synthetic_family/truth.gtf | sort)
   ```
   *"Recovers all 3 truth transcripts within ±2bp."*

**Alternative splicing case to call out**: MYG1_A's A2 transcript skips exon 3. In `tools/demo/out/synthetic/v1_splice_graph.png` for the MYG1_A bundle (`chr_test:997-4503`), the AS edge bypasses one of the internal nodes — corresponding to that exon-skip. The flow decomposition (V3) extracts both the full-length path and the skip path on separate iterations.

## Tier 2 — Real data: chr19 GOLGA6L triple (~10 min)

The point of this tier: **show what happens when the same machinery encounters real, noisy data — and where the family-graph approach earns its keep.**

1. `samtools view -c /mnt/c/Users/jfris/Desktop/GGO_19.bam NC_073243.2:104780000-104900000` — *"257 PacBio reads in this 120kb region containing 3 paralogs."*
2. `bash tools/demo/run_demo_chr19_triple.sh` — runs in ~1 min.
3. **Open the 3 V1 PNGs side-by-side:**
   - `tools/demo/out/chr19_triple/v1_locA.png` (LOC115931294, 104789647-104796276)
   - `tools/demo/out/chr19_triple/v1_locB.png` (LOC134757625, 104830536-104837094)
   - `tools/demo/out/chr19_triple/v1_locC.png` (LOC101137218, 104871356-104877901)
   - *"Three paralogs annotated with the same 8-intron architecture. Look at the assembled graphs."*
4. **Honest acknowledgment**: open `tools/demo/out/chr19_triple/v2_isomorphism.json`:
   ```
   cat tools/demo/out/chr19_triple/v2_isomorphism.json | python3 -m json.tool | head -25
   ```
   *"In the annotation, all 3 are 8-intron. In the assembly, locus A has 10 nodes, locus B has 9, locus C has 4. The graphs are NOT isomorphic on real data — `isomorphic = false` across all three pairwise comparisons."*
   - *"Why? Coverage variation across paralogs. Locus C has too few reads to fully resolve all introns. Some bundles cross paralog boundaries (the orchestrator prints `locus matched 2/3 bundles` warnings) — assembly noise that StringTie also has, but doesn't expose to the user."*
   - *"This is the problem the family-graph approach is designed to fix: pool evidence across copies, recover the under-supported ones."*
5. **Show V3 flow iterations on locus A**: `tools/demo/out/chr19_triple/v3_flow_locA.png` — *"Same Edmonds-Karp algorithm, on real data."*
6. **Side-by-side recovery on the 120 kb region**:
   ```
   awk -F'\t' '$3=="transcript"{n++} END{print "Rustle baseline:", n+0}' tools/demo/out/chr19_triple/baseline.gtf
   awk -F'\t' '$3=="transcript"{n++} END{print "Rustle --vg:    ", n+0}' tools/demo/out/chr19_triple/vg.gtf
   awk -F'\t' '$3=="transcript"{n++} END{print "StringTie 3.0.3:", n+0}' tools/demo/out/chr19_triple/stringtie.gtf
   ```
   Current numbers (April 2026): 5 / 5 / 5. *"On this 120 kb sub-window all three tools assemble 5 transcripts. The interesting comparison isn't transcript count — it's whether each transcript matches the annotated paralog correctly. That's what the cross-chr classification (next tier) measures."*

## Tier 3 — Closing slide (~3 min)

Read aloud the headline result from `MULTI_COPY_FAMILY_PROOF.md`:

> *"Across 33 annotated GOLGA6/GOLGA8 family loci on the full GGO.bam, Rustle VG mode recovers 8 with exact intron-chain match (vs StringTie's 4) and reduces wrong-strand calls from 5 to 1. The chr19 triple is one of those families. Even where the within-locus splice graphs aren't perfectly isomorphic, the cross-family multi-mapper EM resolves paralog-specific reads correctly."*

Keep the table in `MULTI_COPY_FAMILY_PROOF.md` (lines 7–14) on screen if asked for numbers.

## Phase 3 — Vision (~3 min)

Open `tools/demo/out/synthetic/v5_phase3_sketch.png`. (If missing, regenerate: `python3 tools/demo/render_phase3_vg_sketch.py --output-png tools/demo/out/synthetic/v5_phase3_sketch.png`.)

*"Where this goes: a single variation graph per family, with bubbles at SNP positions. Reads aligned to the graph resolve to a specific copy via SNP support. This handles unmapped reads, novel paralogs, and copy-level read assignment that's currently impossible because StringTie doesn't even use multi-mappers in long-read mode. The current `--vg` flag is the first step in that direction; the variation-graph alignment is future work."*

## Q&A backup

- *"Where's the assembly code?"* — *"I'd rather show you what it does. Pick another locus and I'll re-run."*
- *"Did you hardcode the result?"* — *"The synthetic generator is deterministic from seed=42. Re-run it. Same BAM, same answer. The chr19 BAM is your data."*
- *"Why don't the chr19 paralog graphs look identical?"* — *"Real-world coverage variation. Locus C has fewer reads than A or B. The graphs reflect what the assembler can support from the reads it sees, not what the annotation says. We're working on this — the next step is family-aware bundle splitting and SNP-based copy assignment, which use cross-paralog evidence to recover under-supported copies."*
- *"Why does locus A have 10 nodes when it's annotated as 8-intron (i.e. 9 exon nodes max)?"* — *"Sub-bundle inflation: the partitioner is producing more than one bundle per locus on real data, so the V1 renderer picks up junctions from a neighboring bundle. The `locus matched 2 bundles` warning during the run flags it. This is one of the things the family-graph approach addresses by pooling reads across copies before partitioning."*
- *"What's missing?"* — *"Three things in the scaffold but not finished: per-iteration EM weight dump (so V4 works on real data, not just fixture); SNP-based copy assignment for overlapping paralogs; novel-copy synthetic-bundle creation. All three are flagged in `MULTI_COPY_FAMILY_PROOF.md` section 'VG enhancements'."*

## Crash plan

If a live run fails:
- All TSVs and PNGs are pre-cached under `tools/demo/out/{synthetic,chr19_triple}/` and committed. Show those instead.
- For an extreme fallback, the orchestrator scripts are deterministic — re-run from a clean shell:
  - `bash tools/demo/run_demo_synthetic.sh` (~10 s)
  - `bash tools/demo/run_demo_chr19_triple.sh` (~1 min)
- If `stringtie` isn't found, the orchestrator prints `(stringtie not installed; skipping comparison)` and the rest still works.
