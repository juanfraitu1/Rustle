# Genome-wide PSV resolution: where great-ape gene-family copies are read-resolvable

The PSV-aware variation graph (`psv_graph_demo.py`, the two-family demo) scaled to **every**
validated multi-copy gene family. This is the empirical, genome-wide instantiation of the
identifiability arc in `family_to_copy_bridge.md`: detection hands each family a backbone, the
bubbles on that backbone are the columns, and threading the reads through the bubbles either
assigns them to a copy or hits the **K-frontier**.

## Method

Input: the 196 backbone-reinforced families from `family_def_vg_reinforce.py` (de-novo loci,
`~R ∩ ~B`). For each family:

1. **Copy-regions.** Collapse de-novo loci that are the same genomic copy (reciprocal overlap
   > 0.5 = isoforms / read-through nests) into one region. Families whose members all collapse
   to one region are **not genomically multi-copy** and are skipped (42 of 196 — single loci
   with isoform-multiplicity, e.g. a 335 kb transcript nested over a 178 kb one).
2. **Copy model = exon-union.** Each copy's sequence is its **exon-union model** S(v) — the
   reads' aligned exon blocks unioned, introns removed (`family_def_vg_reinforce.py:build_model`,
   cached in `vg_reinforce_copies.fa`). Backbone = longest copy model; the others aligned
   `minimap2 asm20`, the reads (intron-free mature mRNA) aligned `minimap2 map-hifi`.
3. **PSV bubbles** = backbone columns where the copies carry ≥2 alleles **and** ≥3 reads support
   the column (the recurrence test separating a paralog difference from a one-read HiFi error).
4. **Copy paths** = each copy's allele string; **K** = number of distinct paths (resolvable
   haplotypes). **Thread** each read by best PSV-allele match (HiFi-tolerant).

### Why exon-union, and why only 24 families were re-run on it

A first pass aligned copies by their **genomic span** (introns included). That is cheap but
fails when copies have divergent introns — `asm20` cannot align across them, so the copies'
exonic PSVs vanish and the family is falsely called collapsed (6 `align_fail` families), and it
can also *miss* exonic PSVs that lie near intron boundaries (it over-called K=0 on, e.g.,
fam62). Aligning the **exon-union** models fixes both.

Re-running every family on exon-union is unnecessary, and the reason is a monotonicity argument:
exon-union alignment only ever surfaces **more** distinguishing columns, and adding a column
cannot lower K. So a family already resolvable under genomic-span **stays** resolvable — only
the 24 **frontier candidates** (18 genomic-K0 + 6 align_fail) can move. Control re-runs of three
resolvable families (fam18/0/4) confirmed they are unchanged (fam4: 3746/3761 reads agree).
Only the 24 candidates were re-aligned on exon-union (`psv_graph_exonunion.py`); the 121
resolvable families are carried over from the genomic-span pass (`psv_graph_combine.py`).

## Verification (`psv_graph_verify.py`, `psv_graph_combine.py`)

- **Cross-family dedup.** Nine families were the *same* genomic copies re-detected as several
  de-novo isoform-loci (e.g. 168/169/170 = one 2-copy locus counted three times). Collapsed by
  region-set match: **154 → 145 unique families**, read totals recomputed over uniques.
- **Frontier re-alignment.** Of the 24 frontier candidates, exon-union **recovered 10 to
  resolvable** (5 ex-K0: fam7/30/62/122/145; 5 ex-align_fail: fam1/14/21/75/175) and confirmed
  **14 as genuine K=0** (incl. fam34: 8 copies, 0 PSVs over 596 reads → truly identical). The
  `align_fail` bucket went to **zero** — every copy aligns once introns are removed.

## Result (corrected, 145 unique families)

| class | families | % | meaning |
|---|---|---|---|
| **fully resolvable** (K = #copies) | 123 | 84.8 | every copy read-distinguishable |
| **partial** (2 ≤ K < #copies) | 8 | 5.5 | some copies separate, some collapse |
| **genuine K = 0** | 14 | 9.7 | copies exonically identical → provably unresolvable from RNA |
| **indeterminate** | 0 | 0 | (eliminated by exon-union) |
| **read-resolvable total** | **131** | **90.3** | |

**Reads:** 64,066 threaded across the unique families — 35,329 (55 %) assigned to **one copy**,
1,344 to a collapsed group, 19,922 (31 %) cover **no PSV** (the *coverage* face of the
K-frontier — a full-length read spanning only shared exons), the rest unexplained/ambiguous.
Single-copy assignments **agree with the independent best-mapping copy 95.4 %** of the time.

## Why this matters

- **The unresolvable core is ~1 family in 10.** 9.7 % genuine K=0 from RNA read-threading is the
  same magnitude as the ~12 % K=0 from the combinatorial copy-assignment census
  (`copy_assignment_theory.md` / `resolution_improvement_bound.md`) — two unrelated methods,
  over different family sets, agree the collapsed core is order-10 %.
- **The dichotomy is structural, not anecdotal.** 85 % of families lie on the `K = #copies`
  diagonal (Fig. B): the bubbles separate every path. The K=1 row is the frontier. The demo's
  RABL2 (resolvable) and RFPL4A (collapsed) are the two ends of a genome-wide distribution, not
  cherry-picked cases.
- **Detection ⊇ resolution, measured.** All 145 families were *detected* (each is a `~R ∩ ~B`
  component); 90 % are also *resolvable*. The 9.7 % genuine-K0 gap is exactly the bridge note's
  prediction: detectable as a family, provably unresolvable into individual copies.

## Honest caveats

- **Validation is a proxy, not ground truth.** "Agreement with the best-mapping copy" treats
  minimap2's primary as truth, which is exactly what is unreliable for paralogs. 95.4 % global
  agreement is reassuring; the ~5 % disagreement concentrates in size-heterogeneous families
  (family 32: a 124 kb backbone bundled with ~30 kb paralogs) where mapping itself is arbitrary —
  there, PSV disagreement may be a *correction*, not an error.
- **Low-coverage K=0.** Of the 14 genuine-K0 families, five carry < 50 reads (fam176 has 8);
  their copies show 0 divergent columns, but at that depth a rare PSV could be unsampled. The
  well-covered K=0 calls (e.g. fam34/46 with 364–596 reads and still 0 PSVs) are firm.
- **De-novo loci are imported** from detection; the copy-region dedup (overlap > 0.5) is the one
  free parameter added here.

## Performance note

The bottleneck is the **Python pileup + per-read threading** (O(backbone × depth) per family),
not `minimap2`. High-coverage families (3,000–4,000 reads) dominate. The Rust production path
(`psv_linkage.rs::psv_columns_from_reference` + `genotype_family_reads`) does exactly this
extraction natively and would run ~10–100× faster — this Python scan is an analysis harness, not
the shipping pipeline. For the harness itself, capping reads (the classification needs only
enough depth to clear the ≥3-read PSV test) is the cheap speed-up.

Artifacts: `psv_graph_genomewide.py` · `psv_graph_verify.py` · `psv_graph_exonunion.py` ·
`psv_graph_combine.py` · `psv_graph_genomewide.png` · `psv_graph_exonunion_all.json`.
