# vg/GFA pilot on RABL2 (fam39): does the standard toolkit reproduce + extend our custom machinery?

**Question.** Adopt the real `vg`/GFA pangenome toolkit as the *standardized backbone* for our
copy / PSV / read-assignment machinery, instead of the hand-rolled Python VG in
`bench/o2_vg_visualization.py` + `bench/psv_graph_genomewide.py`. Pilot on the RABL2 flagship
family **fam39** and decide: does the standard toolkit **reproduce** our custom results
(cross-tool validation), does it **add** allele-specific expression (ASE), and is adopting it
worth the lift — or does vg **mangle the spliced-paralog structure** enough to keep the custom
machinery?

**Deliverables**
- Orchestrator: `bench/vg_gfa_pilot.py` (runs end-to-end, deterministic)
- Concordance table (machine-readable): `bench/vg_gfa_pilot.json`
- Raw vg outputs: `/home/juanfra/winloci_scratch/vg_pilot/` (`transcripts.fa`, `rabl2.vg`,
  `rabl2.gfa`, `reads.fa`, `reads.gam`, `rabl2.snarls`, `rabl2.decon.vcf`, `reads.pack.tsv`,
  `pilot_report.json`)

**Ground truth (custom o2_vg on fam39):** 5 copies (RABL2A/RABL2B + 3 LOC paralogs), K=5
fully-resolvable, all-SUN-Tier-1; **235 PSV bubbles (195 SUN)**; significance gate assigned
**191/222 reads (0 tied)**. Source: `/home/juanfra/winloci_scratch/o2vg/fam39.json`.

---

## The spliced-RNA caveat and how the transcript-level graph handled it

`vg`/`giraffe`/`GraphAligner` are **DNA-centric and do not model introns**. Building a graph on
the genomic locus (as our custom o2_vg backbone does, *with* introns) would force every read's
spliced gaps to be interpreted as large deletions. **Fix (and it works): build the graph at the
TRANSCRIPT level.** For each of the 5 copies we take the longest annotated isoform, concatenate
its exon bases (reverse-complementing minus-strand copies), and build a small pangenome of the 5
**spliced, intron-free transcripts**. Reads are IsoSeq CCS = already-spliced mRNA, so a read's
query sequence is itself exon-only — we align those directly. Every node on a copy path is exonic
by construction, so **no intron artifact is possible**. The independent hostile re-derivation
confirmed each GFA copy path reconstructs its input transcript **byte-for-byte** (copy0
2868 / copy1 2715 / copy2 3110 / copy3 3327 / copy4 2395 bp), with **no copy collapse**.

The *residual* problem is not introns but **isoform heterogeneity**: independently picking each
copy's longest isoform yields different exon/UTR content, so `vg msga` inserts large indel snarls
at UTR / exon-inclusion boundaries. Testing orthology-matched isoforms did **not** close the PSV
gap, so the gap is structural, not isoform-choice — but the hard prerequisite for production is
still a **matched-isoform / gene-level exon-union** graph (see Recommendation).

**Graph tool choice.** Used `vg msga` (progressive MSA → base-level graph), deliberately **not**
`minigraph -cxggs`: minigraph only bubbles SV-scale (≥~50 bp) events and would erase the
SNP-level PSVs entirely. `vg msga` is *deprecated* (prints a warning) but functional and
deterministic; production should migrate to `vg construct`-from-VCF or minigraph-cactus / PGGB.

---

## Concordance: standard vg toolkit vs custom o2_vg (RABL2/fam39)

| # | Metric | custom o2_vg | vg standard | Concordance |
|---|---|---|---|---|
| (a) | copies as distinct paths | 5 | 5 named paths `copy0..copy4` | **MATCH** — clean, all survive |
| (b) | PSVs: `vg snarls` | 235 bubbles | 62 snarls | vg reports snarl **SITES**, not per-column PSVs |
| (b) | PSVs: `vg deconstruct` sites | 235 | 39 (29 SNP + 10 complex) | **aggregated** into multiallelic snarls |
| (b) | PSVs: base-level SNP columns emitted | 235 | 29 | 17/39 sites keep all 5 copies; snarl-expansion yields 29 SNP columns |
| (b) | **PSVs: LITERAL per-position recall** | 235 bubbles | **22/235 (9.4%)** | vg column matched a GT bubble at the **same genomic pos (±2 bp, consistent −1 bp offset) with identical 5-copy haplotype** (minus-strand complemented) — the honest recall |
| (b) | allele-partition **precision** | 16 patterns | 9 patterns | **9/9 vg patterns ⊆ ours** — every variant vg calls is one we also call |
| (c) | read→copy assignment | 191/222 gate | 97/152 node-vote | **63.8%** (both assign) / **50.8%** (all gate-assigned reads); vg ships **no** assignment |
| (d) | ASE per-copy (`vg pack`) | none | mean ρ=0.10 / total ρ=0.50 vs gate | **NEW capability**; weak (n=5, not significant) |

### (a) Do the 5 copies survive as distinct paths? **YES.**
All five are distinct named paths. `vg stats` → **304 nodes, 366 edges, 3970 bp**. Copy-private
node counts **{copy0:6, copy1:15, copy2:12, copy3:19, copy4:20} = 72 private nodes** (private
*bases* {6,23,167,104,341}=641 — a different quantity; the two were previously conflated, now
reported separately). Strong cross-tool validation of the copy structure.

### (b) PSVs — vg captures the same variation but **aggregates** it; low per-position recall.
- **Honest headline = 22/235 (9.4%) literal per-position recall.** This replaces the earlier
  coordinate-free partition-subset proxy: each vg base-level column is mapped from copy0
  spliced-transcript coordinates back to the o2_vg genomic-backbone offset (via copy0's exon
  structure, complementing minus-strand bases) and must land on a GT bubble at the **same
  position (±2 bp, a consistent −1 bp systematic offset from the GFF-vs-backbone copy0
  definition) with a base-for-base identical 5-copy haplotype**. 22 of the 29 emitted columns
  match; the other 7 have haplotypes that exist in GT only at *other* positions.
- **Precision is high, recall is low.** All 9 distinct allele-partition patterns vg emits are
  patterns we also call (`all_subset = True`) — vg calls nothing spurious. But it recovers only
  22 of 235 per-column PSVs because **17/39 deconstruct sites drop a copy entirely** (the
  divergent/short copy4 = RABL2B, 92 SUN, shortest transcript, fails to co-traverse =
  isoform heterogeneity) and dense SNP blocks are merged into **large multiallelic snarls where
  one copy's indel masks co-located substitutions** (e.g. a 620 bp / 4-allele site, a 77 bp site
  hiding ~53 substitutions).
- **De-emphasized figure:** the "92.3% partition-pattern overlap" is **degeneracy-inflated**
  (only ≤16 allele-partitions exist for 5 copies, so ~all GT bubbles trivially share vg's 9
  patterns) and is **NOT recall** — retained in the JSON only as a clearly-labelled contrast.

### (c) Do reads assign to the same copy as our gate? **WEAKLY.**
vg ships **no** read→copy assignment. Naive private-node voting agrees with the gate on
**97/152 = 63.8%** of reads *both* tools assign, or **97/191 = 50.8%** if every gate-assigned
read is counted (vg abstentions/ties as disagreements). Systematic bias: copy1 gets 65 votes vs
the gate's 31; copy0 gets 12 vs 41 — driven by **unequal private-marker counts** (copy0 has only
6 private nodes) and the far sparser transcript-level markers. **The significance gate is not
reproduced out-of-the-box and remains necessary.** (Aside: the graph *did* place a copy on 12/14
gate-unassigned multimappers present — reference-agnostic positioning the gate abstains on,
though unverifiable.)

### (d) ASE from `vg pack` — a **NEW capability**, but weak and needs normalization.
Per-copy mean private-node coverage: copy0 15.3, copy1 25.5, copy2 28.5, copy3 19.9, copy4 14.4.
copy2 is top in both vg and the gate, but copy4 (2nd in the gate, 48 reads) is *lowest* in pack
mean coverage because its "private" mass is a large **structural isoform block, not SNPs**.
Spearman vs gate read counts: **ρ=0.10 (mean) / ρ=0.50 (total)** — a real per-copy expression
signal our machinery entirely lacks, but **ρ on n=5 copies is not statistically significant**;
production needs **SNP-only private markers + length normalization**.

---

## Determinism

Verified byte-identical across reruns: `rabl2.gfa` (md5 `cc7fdc36f0abee9746aa3f49e41a3c5f`),
sorted `reads.gam` (md5 `3e2670b9abe789f086c24be95e1f152a`, stable even at `-t 8`), and
`bench/vg_gfa_pilot.json` (dumped with `sort_keys=True` to neutralize Counter/set iteration
order). Only deprecation warnings (`vg msga`, `vg deconstruct -e`); **no vg subcommand failed**.

---

## Verdict: PARTIAL REPRODUCE + REAL EXTENSION → **ADOPT vg as substrate, KEEP custom on top**

**vg reproduces** the graph substrate (5 copies = 5 distinct paths; byte-exact path
reconstruction) and is **100%-precise** on the variation it emits (9/9 partitions ⊆ ours) —
strong cross-tool validation of the copy structure.

**vg does NOT reproduce, out-of-the-box:**
1. the 235 per-column PSV catalog — it aggregates into 39 deconstruct sites / 62 snarls, and a
   literal position+haplotype match recovers only **22/235 (9.4%)** because copies drop out of
   snarls and indels mask co-located SNPs;
2. the read→copy **significance gate** — naive node-vote agrees only 63.8% (50.8% incl. ties),
   biased; vg has no assignment primitive at all.

**vg adds** `vg pack` per-copy coverage = **allele-specific expression**, which our machinery
lacks (weak at n=5 here but genome-wide it becomes the differentiating capability).

**Recommendation: ADOPT vg/GFA as the standardized graph substrate + ASE layer, and layer the
copy-assignment theory (SUN/Strong-Sep significance gate + PSV-column machinery) on top.** vg
mangles nothing about the *copy* structure once the graph is spliced; what it does **not** do is
per-column PSV discovery or read→copy assignment — precisely our contribution. So this is a clean
division of labor, not a reason to discard the custom code.

**Hard prerequisite:** build the transcript graph from **matched orthologous isoforms or a
gene-level exon union**, not independently-chosen longest isoforms — otherwise isoform
heterogeneity drops the most divergent copy from snarls and manufactures spurious structural
bubbles (this is the single largest driver of the low PSV recall and the biased assignment/ASE).
And migrate off deprecated `vg msga` for production.

---

## Which of the 9 vg features to wire next

| # | vg feature | Pilot result | Wire next? |
|---|---|---|---|
| 1 | GFA construction (`vg construct`/`msga`) | works; msga deprecated, isoform-sensitive | **YES (P0)** — but move to `vg construct`-from-VCF or minigraph-cactus/PGGB on a **matched-isoform / exon-union** input |
| 2 | Read alignment → GAM (`GraphAligner -x vg`, or `vg giraffe`) | 205/222 reads placed, median MAPQ 60, id 0.994; splits reads into multiple GAM records (dedup per read) | **YES (P0)** — reference-agnostic per-read paths; evaluate `vg giraffe` for speed/determinism |
| 3 | `vg deconstruct` (PSV VCF) | 39 sites; 9.4% literal per-position recall (aggregates SNPs, drops copies) | **PARTIAL (P2)** — useful as a *seed* / cross-check, but needs snarl-traversal expansion; **keep the custom PSV-column caller as source of truth** |
| 4 | `vg snarls` (bubble decomposition) | 62 snarls vs our 235 bubbles | **PARTIAL (P2)** — good structural skeleton; our column machinery expands the dense snarls it merges |
| 5 | `vg stats` | clean (304 nodes/366 edges/3970 bp) | **YES (P1)** — trivial QC, adopt |
| 6 | GFA/pangenome as standardized on-disk format | valid, deterministic, interoperable | **YES (P0)** — the whole point: standardized substrate replacing hand-rolled Python VG |
| 7 | **`vg pack` (per-node/per-path coverage → ASE)** | NEW capability; per-copy expression signal (weak at n=5, needs SNP-only markers + norm) | **YES (P1)** — the genuinely additive feature; wire with SNP-only private markers + length normalization |
| 8 | `vg augment` (fold read variation into the graph; novel/reference-absent alleles) | **not run** in this pilot | **EXPLORE (P3)** — directly relevant to our O4 reference-absent-copy discovery; highest-value untested feature |
| 9 | `vg paths` (per-path support / path extraction) | used to get copy paths + owners | **YES (P1)** — needed for path-level support and per-copy attribution feeding (7) |

**Priority summary:** P0 = the substrate (GFA construct #1/#6 + GAM alignment #2). P1 = the cheap
wins to adopt now (stats #5, paths #9, and **pack #7 for ASE** — the new capability). P2 =
deconstruct #3 / snarls #4 as cross-checks/seeds but **not** as replacements for our PSV-column
caller or the significance gate. P3 = `vg augment` #8, the untested feature most worth piloting
next (reference-absent alleles / O4). The copy-assignment theory (significance gate, SUN /
Strong-Sep, per-column PSVs) stays as our layer on top of the vg substrate.
