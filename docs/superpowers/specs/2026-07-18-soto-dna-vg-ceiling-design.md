# DNA Variation-Graph Ceiling vs Soto SDs — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-18
**Builds on:** the Soto benchmark (`bench/soto/`, `bench/SOTO_BENCHMARK.md` — RNA member sensitivity 276/362 = 76.2%), the DNA oracle prototype (`bench/soto/soto_dna_oracle_prototype.py`), and the shipped copy-graph GFA/Bandage visual language (`src/rustle/vg_family/copy_graph.rs`).
**Objective:** an advisor-facing **visual** demonstration that at the DNA level a variation graph represents **100% of the Soto segmental-duplication copies** as distinct paths (the achievable identifiability *ceiling*), while Rustle's RNA method recovers 76.2% of that ceiling — the gap being a proven K=0/expression *floor*, not a method failure.

## Goal

For a small set of flagship Soto SD families, build a base-level **variation graph** (an abpoa POA multiple-alignment → GFA; `vg msga` optional) from the family's member **genomic** sequences, emit a Bandage-renderable `<family>.gfa` + `<family>.colours.csv`, and color each copy by whether Rustle's RNA method recovered it (**green**) or not (**red** = DNA-only: K=0 exon-identity collapse / silent / coverage-limited). Every Soto member of the family appears as a path → the DNA graph is **100% complete**; the red paths are exactly the identifiability floor, made visible.

## Motivation

The objectives audit flagged O1 family-definition **circularity** (homology validating homology) as the single biggest defense risk. Soto 2025 is an *independent* external catalog — its copies/CN come from **DNA read-depth on a T2T assembly**, a different method with no shared aligner or silver standard (`SOTO_BENCHMARK.md`). So a DNA-level demonstration against Soto is *external corroboration*, and the honest "ceiling" framing turns the 24% RNA gap from an apparent weakness into a **proven identifiability floor**: at DNA level every copy is present and countable (100% achievable); RNA is bounded by what is *expressed* and by K=0 exonic identity. The variation graph is the vehicle that makes "all copies present at DNA level" concrete, visual, and tied to the thesis's VG framing (Canzar's "multimappers → assembly").

## Design

### §1 — Inputs (all on disk; per-family, no whole-genome)

- `bench/soto/80_fams.chr.bed` — Soto members, `chrN  start  end  GENE|ID_k  0  .` (362 members / 83 families).
- `soto_members.fa` = `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_members.fa` (11 MB) — member **genomic** sequences (introns included), headed `>chrN:start1-end` (1-based; BED `start` → header `start+1`).
- `bench/soto/soto_member_detection.tsv` — `family_id, gene, chrom, start, end, detected(Y/N), recovered_by` — the RNA-recovered vs not labels (the color source).
- `bench/soto/soto_floor_decomposition.tsv` — the miss-reason per undetected member (K=0 / silent / coverage), used to confirm the red flagship's copies are K=0 and to caption them.
- **`pyabpoa`** (already a dependency) — the primary graph builder (base-level POA MSA → variation-graph GFA). `vg` = `/home/juanfra/miniforge3/bin/vg` (present, v1.73) as an authentic-`vg` alternative — but its `msga` subcommand is **deprecated**. `minigraph` (present) is **rejected**: it is SV-level (default min variant ~50 bp, ignores SNPs), so it would collapse SD copies that differ only by SNPs into one path — destroying the base-level PSV/K=0 story this demonstration exists to show. Base-level resolution is required; abpoa provides it.
- Optional exon overlay: `A119b.t2t.bam` = `/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam` (indexed) — spliced Iso-Seq reads whose aligned blocks (CIGAR, `N`=intron) give exonic intervals per member locus.

### §2 — Flagship set (3 core + 1 optional)

Chosen to show the full spectrum in ≤4 pictures:
1. **SRGAP2 (`ID_462`)** — 4 members, **4/4 RNA-recovered** → all-green graph. The human-specific brain flagship; the clean *ceiling*.
2. **PMS2P (`ID_8`)** — ~9 members, all recovered → a large green fan; the *ceiling at scale*.
3. **`ID_63`** — 8 members, **4 recovered (50%)** → **green + red in one family**; the *ceiling-vs-floor picture*. (Alternates if its red members are not K=0-dominated per `soto_floor_decomposition.tsv`: `ID_14` 4/7, `ID_65` 7/9.)
4. *(optional)* **SLC9B1P1 (`ID_26`)** — RNA 0/2, DNA oracle finds 8 copies incl. 7 novel — the *DNA ceiling exceeds even Soto's annotation*. Its beyond-Soto loci come from `soto_dna_oracle_prototype.tsv` (extracted from `chm13v2.0.fa` by `samtools faidx`), a different extraction path, so it is marked optional and built only if the core three render cleanly.

### §3 — Graph construction (per family)

Per flagship family:
1. Parse `80_fams.chr.bed` for the family's members; extract each member's genomic sequence from `soto_members.fa` by its `chrom:start+1-end` header into a per-family FASTA (`<fam>.members.fa`).
2. **Build the base-level variation graph via abpoa** (primary): `pyabpoa.msa_aligner().msa(member_seqs, out_msa=True)` → gapped aligned rows (`.msa_seq`). Convert the column-MSA to a GFA:
   - Walk the alignment left→right; merge each maximal run of columns where all non-gap sequences carry the **same** base into a single node (the conserved/shared backbone segments).
   - At a column where sequences differ (SNP) or a gap run (indel) occurs, emit **one node per distinct allele** (a bubble).
   - Add `L`-lines between consecutive nodes along each member; each member's ordered node list is its `P`-line path. Deterministic; every family member is a `P`-line.
   This yields base-level SNP/indel bubbles, so K=0-similar copies still appear as distinct paths that share the backbone and split at PSVs — exactly the visual the demonstration needs.
3. **Alternative (`--builder vg`)**: `vg msga -f <fam>.members.fa` → `vg view` → GFA (a genuine `vg` graph with per-member paths, one step). Functional but `msga` is deprecated; offered for an authentic-`vg`-branded artifact if desired. The builder used is logged per family.

### §4 — Colouring (`<fam>.colours.csv` + `<fam>.legend.tsv`)

Mirror the `copy_graph.rs` convention (node-colour CSV that Bandage loads):
- Map each `P`-line (member path) to a color from `soto_member_detection.tsv`: **green `#1e8e3e`** if `detected == Y`, **red `#d93025`** if `detected == N`.
- Per node, color by which member paths traverse it: nodes traversed **only** by green paths → green; **only** by red paths → red; traversed by **both** → **grey `#9aa0a6`** (shared/conserved sequence — where copies are identical). The grey backbone with red/green arms is the visual: a red copy that shares the grey backbone with a green copy but diverges into its own red nodes is a copy DNA separates but RNA (if the shared stretch is exonic) cannot.
- `<fam>.legend.tsv` names each path: `gene, locus, detected, recovered_by, colour`.

### §5 — Optional exon overlay (K=0 made explicit)

If enabled (`--exon-overlay`), for each member locus fetch spliced A119b reads (`samtools view A119b.t2t.bam <locus>`), union their aligned exon blocks (split CIGAR on `N`) → exonic intervals; tag graph nodes overlapping an exonic interval as **exonic** (bold/darker shade) vs intronic. This makes the K=0 floor explicit: red (DNA-only) copies that share **exonic** nodes with green copies but diverge in **intronic** nodes are exactly the exon-identity collapses — DNA separates them intronically, RNA cannot. Default OFF (adds per-locus BAM fetches); the core green/red graph stands without it.

### §6 — The 100%-presence check (the ceiling claim, verified)

After building each family graph, assert that **every** Soto member of that family is present as a `P`-line path in the GFA. This is the operational "DNA = 100%" claim, checked per family (not asserted). A member that fails to appear (e.g., `vg msga` dropped a highly-divergent sequence) is logged as a graph-construction miss and excluded from the "100%" count with an honest note — never silently counted as present.

### §7 — Honesty rail (captions)

An `index.md` (or per-family caption) states, verbatim and unconditionally, under each graph:
> *DNA variation graph = the ceiling: all N Soto copies present as paths (Soto-corroborated, independent DNA-read-depth catalog). Green = RNA recovered; red = DNA-only (K=0 exon-identity / silent / coverage). The VG represents what is given; it does not "detect" families. RNA recovers 76.2% of this ceiling genome-wide; the gap is the decomposed identifiability floor, not a method failure.*

No claim that the graph *detects* anything it was not handed.

### §8 — Compute & crash-safety

- Per-family graphs only; member sequences come from the 11 MB `soto_members.fa` (a `samtools faidx`/dict lookup), never the 3 GB genome except the optional SLC9B1P1 loci. **No whole-genome VG, no `copy_assign`.**
- `vg msga` / `abpoa` run **foreground, one family at a time**, serial; outputs under `/home/juanfra/winloci_scratch/soto_vg/`. No background/nohup/pkill.
- SRGAP2 members span up to ~260 kb; if a family's graph is too large for a crisp Bandage view, log it — do not auto-truncate (honesty), but note the size.

### §9 — Testing

- **Hermetic unit test** on the colour-assignment core (`node_colour(paths_through_node, detected_map)`): green-only node → green; red-only → red; mixed → grey; matches the §4 rule. No `vg`/BAM needed.
- **Hermetic test** on the member-FASTA extraction (header `chrom:start+1-end` from a synthetic `members.fa` + a BED row) and the per-family presence check (`all members appear as P-lines`) against a tiny synthetic GFA.
- **Data-gated smoke** (Task, foreground): build ONE family (PMS2P `ID_8` — moderate member sizes ~5–15 kb) end-to-end via the abpoa primary path; assert the GFA is valid, every family member appears as a `P`-line path, and `colours.csv` is non-empty. Crash rule honored (single family, `winloci_scratch`).

## Out of scope

- The catalog-wide DNA-100%-vs-RNA-76.2% *number/table* (the user chose visual artifacts only; the number already lives in `SOTO_BENCHMARK.md`).
- Any change to the RNA detection pipeline or the Soto eval scripts.
- Rendering PNGs (we emit GFA + colours.csv; the user opens Bandage). No Bandage dependency.
- Running SEDEF/BISER genome-wide (too heavy; the member-sequence graph is the scoped stand-in).

## Global constraints

- **Visual artifacts only**: emit `<fam>.gfa` + `<fam>.colours.csv` + `<fam>.legend.tsv` (+ `index.md` captions) per flagship; no catalog-wide number is required.
- **Honest ceiling framing**: the VG *represents* all given copies (the ceiling); it does not "detect" families. Every caption carries this.
- **Presence is checked, not assumed**: the "100%" is a per-family assertion that every Soto member appears as a path.
- Per-family, foreground, serial `vg`/`abpoa`; member seqs from `soto_members.fa`; **no whole-genome, no `copy_assign`**; outputs under `winloci_scratch`.
- Colours: green `#1e8e3e` (RNA-recovered), red `#d93025` (DNA-only), grey `#9aa0a6` (shared) — consistent with `copy_graph.rs`.
- **`pyabpoa` base-level MSA→GFA is the primary builder** (deterministic, correct granularity, no deprecated tool); `vg msga` is an optional authentic-`vg` alternative (`--builder vg`, deprecated); `minigraph` is rejected (SV-level). The builder used is logged per family.
