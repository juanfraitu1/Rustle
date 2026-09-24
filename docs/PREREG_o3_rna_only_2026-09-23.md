# Pre-registration — O3 as far as RNA can carry it: flag, characterise, screen, hand to DNA

**Written 2026-09-23 (§6ze), before the tool exists and before any locus is scored.** User: *"can we polish it
so we can get as far as that before providing DNA info?"* — i.e. from one BAM, produce the divergence-pile
candidates, recognise them as the signal of a possible extra copy, take every RNA-only step that can kill
or strengthen them, and stop at the line only DNA can cross (copy number).

## What is inherited, unchanged

- **The statistic** is the held-out-validated S2 divergence mixture (08-14, `o3_excise/detector.py`,
  register rows 21/24/803): deterministic 1-D 2-means on the per-read `de` of a locus's primary reads;
  fire iff the high cluster holds a fraction `m ∈ [0.10, 0.50]` and the median gap between clusters
  `delta ≥ 0.01`, with ≥ 10 reads. TPR 0.27 overall / 0.45 when the hidden copy is ≥ 1% diverged,
  FPR 0.02 on independent negatives. Not re-fitted here.
- **The reconstruction** is §6ff's patched consensus (`bench/o3_reconstruct.py`, row 720): the reference
  patched at the sub-pile's consistent sites, aligned to the whole assembly.
- **The screens** are the ones that killed every earlier candidate: immunoglobulin / T-cell receptor
  hypermutation (`project_o3_mischain_immunoglobulin`), run exclusivity (EBV, IgG: one of four pooled runs),
  allelic magnitude (gorilla π ≈ 0.0015–0.0020; a real copy sits 30–40× above it).

## The tool — `o3_rna_flag` (Rust, `src/rustle/vg_family/o3_rna.rs` + `src/bin/o3_rna_flag.rs`)

Input: a BAM (`-F 2308` primaries, `de:f` and `--eqx` CIGARs as every project BAM has), the primary
genome FASTA, a locus set (GTF `gene_id` spans, GFF `gene` features, or BED), a minimap2 index of the
primary genome for the home search, optionally the annotation GFF for the hypermutation screen and any
number of `--confirm NAME=index.mmi` genomes (here: the two parental haplotypes).

Per locus, in this order, every rule fixed now:

1. **Pile.** Primary reads overlapping the locus, capped at 2,000 by name order (deterministic). Fewer than
   10 → `too_few`.
2. **Mixture.** 2-means on `de`; fire iff `0.10 ≤ m ≤ 0.50`, `delta ≥ 0.01`, high cluster ≥ 3 reads.
   The high cluster is the **sub-pile** (the putative hidden copy's reads), the low cluster the **host**.
3. **Consistency (new).** A hidden copy's reads share the SAME mismatches against the host (they are the
   copy's PSVs); error, editing and somatic drift do not. From the `--eqx` CIGARs: a **PSV site** is a
   reference position covered by ≥ 3 sub-pile reads with ≥ 80% of them mismatched there, and covered by
   ≥ 3 host reads with ≤ 20% mismatched. `n_psv` = number of such sites; `shared_frac` = (sub-pile
   mismatches falling on PSV sites) / (all sub-pile mismatches). **`copy_consistent` iff `n_psv ≥ 3` and
   `shared_frac ≥ 0.5`**, else `scattered`.
4. **Editing screen (new).** If `n_psv ≥ 5` and ≥ 80% of PSV substitutions are A→G or T→C
   (reference-oriented, either strand) → `rna_editing`.
5. **Patched consensus.** Template = the sub-pile read with the longest aligned reference span; its
   `=`/`X` blocks are taken from the reference and patched at PSV sites with the sub-pile majority base;
   blocks concatenated = a spliced consensus of the hidden copy over what the reads cover. Written to
   `<out>.consensus.fa`.
6. **Home search.** All consensus sequences aligned in ONE minimap2 call (`-x splice:hq -c --eqx -N 20`)
   to the primary index. Hits with query coverage ≥ 0.8 only. `host_identity` = identity of the hit
   overlapping the locus; `other_identity`/`other_locus` = the best hit elsewhere. **`unannotated_paralogue`
   iff `other_identity > host_identity`** (the consensus is closer to another reference locus than to the
   one it was flagged at — a copy the annotation or the catalog missed, not a reference-absent one).
7. **Hypermutation screen.** Locus overlaps an annotated `IG[HKL]*`/`TR[ABDG]*` gene or a gene whose
   description contains "immunoglobulin" / "T cell receptor" → `hypermutation`.
8. **Run-exclusivity screen.** Read names carry the sequencing run (prefix before the first `.` or `/`);
   only when the pile shows ≥ 2 runs. Fisher exact 2×2 (sub-pile vs host, top run vs the rest):
   **`contamination` iff p < 1e-3 and the sub-pile is ≥ 95% one run.** Testis OR6737 is a single movie
   → screen inert there, said in the output.
9. **Allelic magnitude.** `delta` in units of π (`--pi 0.002`): informational column `delta_over_pi`;
   `delta ≥ 0.01` already implies ≥ 5π.
10. **Verdict** (first that applies): `contamination` → `hypermutation` → `rna_editing` → `scattered` →
    `unannotated_paralogue` → **`reference_absent_candidate`**.
11. **What DNA would have to show** (columns, not decisions): `expected_dna_depth_ratio = 1/(1−m)` if the
    hidden copy is expressed like the host; the consensus FASTA as the probe for WGS k-mers/depth or PCR.
12. **Confirmation, when a confirm genome is given.** The same consensus batch aligned to each
    `--confirm` index; `conf_<NAME>_identity`/`_locus` = best hit at coverage ≥ 0.8.
    **CONFIRMED iff `conf_identity ≥ 0.99` and `conf_identity − host_identity ≥ delta/2`** (the consensus
    is a near-perfect match somewhere in the other genome while it is `delta`-far from the primary host).

## Validation — committed now

**Substrate with DNA truth:** the fibroblast IsoSeq of KB3781 (`fibroblasts/GCA_029281585.2_flnc_mm.bam`,
the assembly's own animal, on the GCF primary), loci = the gorilla RefSeq `gene` features
(`winloci_data/GGO_genomic.gff`), home search on `winloci_data/GGO.splice.mmi`, confirm genomes = the
parental haplotypes `mGorGor1.pat.splice.mmi` / `mGorGor1.mat.splice.mmi` (v2.0; the primary is a
16-PAT / 9-MAT mosaic, `o3_hapcnv/pri_provenance.tsv`, so copies present only on the other haplotype of a
chromosome are genuinely absent from the primary).

Reported: loci scanned / with ≥ 10 reads / fired / copy-consistent / after screens / by verdict; for the
`reference_absent_candidate` set the fraction CONFIRMED by pat or mat; the 11 expressed genes whose
probe count is higher on the non-primary haplotype (`hapcnv.tsv`, ≥ 20 fibroblast reads) and whether each
fired (information, not a recall number — those counts were never one-to-one verified).

| outcome on KB3781 fibroblast | verdict |
|---|---|
| precision of `reference_absent_candidate` against the haplotypes ≥ 0.5, on ≥ 5 candidates | ⭐ **the RNA-only flag is worth handing to DNA as is** |
| 0.2–0.5, or < 5 candidates | ⚠ **flag stands but needs a DNA step before any count is quoted** |
| < 0.2 | ⛔ **the RNA-only chain does not enrich for real copies; report the screens' kills only** |

**Second substrate, no DNA:** the OR6737 testis BAM, same loci and index, no confirm genomes. Reported:
the same funnel and the candidate rate per 1,000 expressed loci next to KB3781's. The advisor's claim
predicts a higher rate in the non-reference individual; this is reported, not judged (tissue confound).

**Predicted, before looking:** ⚠. On the matched individual the only reference-absent copies are the other
haplotype's extras (56 genes by probe count, 11 expressed at ≥ 20 reads, most with a difference of one
unit), so the confirmed set will be small (≤ 10) and the fired set will be dominated by noise the
consistency rule removes; precision 0.2–0.5. Testis rate higher than fibroblast, mostly as
`unannotated_paralogue` and `scattered`.

I will not change thresholds, the verdict order, the confirmation rule or the bar after seeing any number.

---

## Addendum 1 (2026-09-23, after the first KB3781 pass, before any foreign-genome alignment) — a cross-species screen

The first pass (446 fired → 53 `reference_absent_candidate`, 1 confirmed) showed the candidates' consensus at the
SAME identity to pat and mat as to the primary (median 0.983 / 0.984 vs 0.984), i.e. absent from the whole diploid
genome at 1.5-3% divergence, and the list is dominated by immune-cell genes (CD4, TRAF1, TCF7, SPN, THEMIS2, RASSF5,
OAS2, HERC5) in a fibroblast line. Human-gorilla transcript divergence is ~1.5-2%, so the obvious RNA-only
explanation is reads of ANOTHER SPECIES (a pooled human library), which none of the pre-registered screens tests —
the 08-13 memory dismissed human contamination only for candidates at 5.6% divergence.

**Rule, fixed before aligning anything to human:** `--foreign NAME=index` genomes (here CHM13 via
`npip_ladder/idx/target.splice.mmi`); verdict `foreign_species` iff the consensus's best foreign hit (coverage
≥ 0.8) has identity ≥ 0.995 and exceeds the primary host identity by ≥ delta/2 — the confirmation rule's shape
applied to the other species' genome. Placed after `contamination` (run-exclusive) and before `hypermutation` in
the verdict order. Chimp/orangutan are not tested (no splice index on disk); a `foreign_species` call names the
species whose genome matched, nothing more. The bar of the validation is unchanged; the candidate set it is
computed on is what remains after this screen.

---

# OUTCOME (2026-09-23)

## Positive control (simulated, chr20, 40 genes with one extra genomic copy; `bench/o3_sim_copies.py genomic`, formerly `o3/simB.py`)

| divergence | fired / 40 | `reference_absent_candidate` | confirmed against the mutated copies (`--confirm`) | median n_psv / shared_frac / host identity / confirm identity |
|---|---|---|---|---|
| 2% | **40 / 40** | 40 | **40 / 40** | 45 / 0.951 / 0.980 / 1.000 |
| 1% | 15 / 40 | 15 | **15 / 15** | 12 / 0.915 / 0.989 / 1.000 |

Every step behaves as specified on a real genomic copy: the mixture fires (1% sits on the `delta ≥ 0.01`
edge, matching the S2 detector's 0.45 stratified TPR), the PSVs are shared (0.95), the consensus finds no
other home on CHM13, and the confirmation rule recovers the copy at identity 1.000 with the `delta/2`
margin. ⚠ A first control built on `simA.py`'s reads exposed a limitation worth stating: simA mutated each
TRANSCRIPT independently, so a gene's sub-pile mixed several mutation sets and 24/40 came out `scattered`
(PSV sites need ≥ 80% of the sub-pile). The rule assumes ONE hidden copy per locus; two distinct hidden
copies at one locus will read as `scattered`. simA's names also fired the run screen (fixed: a run label
must look like an SRA/PacBio run id, `run_of`).

## KB3781 fibroblast (the assembly's own animal), 41,193 annotated genes

| step | n |
|---|---|
| loci with ≥ 10 primary reads | 18,961 |
| mixture fired | 446 (377 distinct sub-piles) |
| `scattered` (no shared sites) | 347 |
| `hypermutation` (IG/TR annotation) | 43 (18 sub-piles) |
| `rna_editing` | 2 |
| `contamination` (run-exclusive) | 1 |
| `foreign_species` (human, addendum 1) | **0** |
| **`reference_absent_candidate`** | **53** (50 sub-piles) |
| confirmed by pat or mat | **1** (KIAA0513: maternal copy at 0.9987 vs primary 0.9912, PSVs in 2 of 12 exons) |

**By the pre-registered bar: ⛔ precision 1/53 = 0.019.** Two things must be said with it. (1) The
truth the bar uses is the diploid assembly, which is complete for the OTHER HAPLOTYPE's extra copies
but by construction cannot confirm a copy the assembly lacks; the 11 expressed genes with a higher probe
count on the non-primary haplotype did not fire at all (all `no_mixture`) — those extras are below the
1% floor, as the S2 ceiling predicts. (2) The 52 unconfirmed candidates are not noise by any RNA test
available: their consensus sits at the SAME identity to pat and mat as to the primary (median 0.983 /
0.984 vs 0.984 — not in this animal's diploid genome at 1.5–3%), the PSVs spread over most exons (median
80% of blocks carry one, PSV density 1.7% of consensus length; `scattered` and `hypermutation` piles carry
PSVs in 1 of 8 blocks), and a cross-species home search gives human 0.979, chimp 0.975 (3 would pass the
foreign rule: OAS1 0.9969, LOC101145093 / LOC115935262 0.9963), orangutan 0.966 (0 pass; TGM2 0.9949
borderline). What remains is sequence closer to gorilla than to any other ape yet absent from the gorilla
assembly, expressed at 10–40% of the locus, concentrated in immune-cell genes (CD4, TRAF1, TCF7, SPN,
THEMIS2, RASSF5, OAS1/2, HERC5, SRGN, CORO1A, LAPTM5) in a fibroblast line — the same class as the three
08-13 candidates. RNA cannot take it further: **the next discriminator is DNA** — the cell line's own HiFi
DNA reads (SRA, not on disk): do the 1.7% PSV k-mers of each consensus occur in the DNA at the depth of
one copy? That is the machinery this tool hands to; its output (`fib.o3_rna.tsv`, `fib.consensus.fa` in
`/mnt/linuxdisk/tmp/gw22/o3/`) is the probe set.

## OR6737 testis (a different animal, no DNA), same loci

| step | n |
|---|---|
| loci with ≥ 10 reads | 21,454 |
| fired | 1,790 |
| `scattered` | 1,724 |
| `hypermutation` | 61 |
| `foreign_species` (human) | 0 |
| **`reference_absent_candidate`** | **5** (0.23 per 1,000 expressed loci vs 2.8 in the fibroblast line) |

The non-reference individual has FEWER candidates, not more — opposite to the direction the reference-bias
claim predicts, though tissue and library confound it (testis: single movie, run screen inert). The five:
LOC101141792, LOC134758970 (other locus at 0.963), LOC109025073, LOC101147715 (host identity 0.898, the
most divergent candidate in either library), LOC129529606.

## What is shipped

`o3_rna_flag` (`src/bin/o3_rna_flag.rs`, logic + 9 unit tests in `src/rustle/vg_family/o3_rna.rs`):
sweep-based scan (one sequential pass per contig, RecordBuf decode only for fired loci — the same 80 Mb
contig 4:20 → 0:12), `--scan-only` / `--from-scan` batching for laptops, one minimap2 call per genome,
`--confirm` and `--foreign` genomes, the pre-registered verdict, `expected_dna_depth_ratio`, the consensus
FASTA. Genome-wide gorilla: scan 3.5 min + align 4 min (three 13 GB indexes, 17 GB peak).

---

## Addendum 2 (2026-09-23, 21:14) — the exon-order rearrangement detector, ON BY DEFAULT

**What was measured first** (`bench/o3_sim_copies.py shuffled`, formerly `o3/shuf/simC.py`; 10 chr20 genes, extra copy with exons 2 and 3
swapped): the shuffled copy's reads are NOT lost — 100% map to the template as one MAPQ-60 primary with no soft
clip — but minimap2 keeps the longest colinear exon run and encodes the displaced exon as an INSERTION of exon
size (median 170 bp; 86–90% of shuffled reads carry an insertion ≥ 50 bp vs 0% of template reads), with one
intron fewer (an apparent exon skip); NM jumps 3 → 172 while gap-compressed `de` stays 0.0015 → 0.0023. So the
divergence mixture is blind to a structural copy and the assembler would call it an isoform.

**Rule, fixed before any real locus is scanned with it.** Pass 1 also records, per locus, the reads whose CIGAR
carries an insertion ≥ 50 bp. A locus with ≥ 3 such reads is decoded and each such read is tested: the inserted
sequence is searched (12-mer seeds, then direct comparison) inside every intron gap (`N`) of the SAME read; a
hit at ≥ 0.90 identity over ≥ 90% of the insertion means the read carries a reference exon it also skips —
an **exon-order rearrangement** (`rearranged`). If that matched exon is ALSO covered by the read's aligned blocks
(the exon appears twice in the read) it is a tandem/rolling-circle duplication (`duplicated_exon`: the circRNA /
back-splice class), not a rearrangement. Reads are clustered by (matched exon ± 20 bp, insertion site ± 20 bp);
**a cluster of ≥ 3 `rearranged` reads fires the locus as `structural`** (status `fired_structural`, or
`fired_both` with the mixture). The structural sub-pile is the cluster; its consensus is the template read's own
sequence (the copy's transcript in read order); the screens (run exclusivity, IG/TR, foreign genome) and the home
search apply unchanged, so an existing paralogue with that exon order elsewhere in the reference resolves to
`unannotated_paralogue`. The verdict table gains a `class` column: `divergent` / `structural` / `both`.

**Validation, committed now.** (i) `simC` at d = 0 and 2%: the shuffled genes must fire `structural` — bar
≥ 9/10 at both divergences (⭐), 6–8 (⚠), fewer (⛔) — and at d = 2% they must ALSO fire the mixture (`fired_both`).
(ii) Negative control: the simB 2% BAM (40 genes, unshuffled copies): 0 structural fires expected; any fire is
inspected. (iii) Real data: KB3781 fibroblast and OR6737 testis rescanned; the count of `structural` loci and of
`duplicated_exon` clusters is reported, with the verdict funnel, and NOT judged (no truth); the largest clusters
are listed for inspection. Divergence statistic stays `de` (the validated S2 rule is not re-fitted); the
structural detector is a second, independent trigger.

### Addendum 2 — OUTCOME (2026-09-23, 22:45)

**Simulation (`simC`, 10 shuffled genes):** ⭐ 9/10 fire `structural` at d = 0 and 9/10 `both` at d = 2%; the
miss (LINC01260) has its swapped exons under the 50 bp insertion floor (minimap2 absorbs them without an
insertion). Negative control (simB 2%, 40 unshuffled copies): 0 structural fires, 70 mixture fires as before.
On the way the matcher had to change from "the whole insertion inside one intron gap" to "a ≥ 50 bp stretch of
the insertion at ≥ 0.90 identity anywhere in the locus": minimap2 represents the displaced exon
inconsistently (partial anchor blocks, composite insertions), so only a partial, locus-wide match is robust.

**KB3781 fibroblast (real, with the haplotypes as DNA truth):** 74 loci fire structural only and 25 both
(62 distinct sub-piles over 99 annotated loci; cluster size median 7, q3 43). Verdicts of the structural class:
**64 `reference_absent_candidate`, 8 `unannotated_paralogue`, 2 hypermutation**; of the 64 candidates,
**20 are CONFIRMED by a parental haplotype** (the read-order transcript aligns at ≥ 0.99 to pat or mat while
it is 0.67–0.95 to the primary host) — against 1/53 for the divergence class. The confirmed ones are large,
haplotype-specific exon-order differences: CDK11B/SLC35E2B (cluster 448 reads, primary/pat 0.675, mat
0.991), CCN3 (90 reads, mat 0.9994), RNF168 (68, mat 0.9994), FBXL2/UBP1 (50, mat 0.998), NIPAL2 (33, mat
0.9998). These are copies whose exon order is absent from the primary and present on the other haplotype —
reference-absent at the structural level, detected from RNA and confirmed by DNA, exactly the class the
divergence mixture cannot see. The 44 unconfirmed (GPC6 77 reads, EIF2S2 52, EXT1 45, ARL15 41, CXCL13 30…)
sit at 0.64–0.94 to primary, pat and mat alike: absent from the whole diploid assembly; cDNA template
switching is the artefact class to exclude next (a read-specific junction would not cluster at one exon and
one insertion site ± 20 bp, but recurrent switching at homologous exons could), which is again a DNA question.

**OR6737 testis (different animal, no DNA):** 81 structural only + 29 both; **63 candidates (36 sub-piles)**,
2 `unannotated_paralogue`, 16 hypermutation. **The same loci recur across the two animals**: RNF168 (126
reads), SLC35E2B/CDK11B (60), ARL15 (41) — a rearrangement seen in two individuals from two tissues and, for
CDK11B/RNF168, present on KB3781's maternal haplotype: a segregating structural variant (or a primary-
assembly exon-order error), not a library artefact.

**Cost:** the structural branch decodes every locus with ≥ 3 insertion-carrying reads, so the scan is
slower on deep libraries (fibroblast batches ~8 min each, testis ~25 min); the align phase is unchanged.
The detector is on by default (no flag); `class` ∈ divergent / structural / both is a column of the table.
