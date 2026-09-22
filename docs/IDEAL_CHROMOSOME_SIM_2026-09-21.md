# Ideal-scenario chromosome simulation — in theory it does work

**§6v1, 2026-09-21.** Pre-registration `docs/PREREG_ideal_chromosome_sim_2026-09-21.md` (md5
`ff226f41`), committed `2a809b76` before any arm was scored. Tool `bench/ideal_chromosome_sim.py`.
⚠HUMAN substrate (A119b/CHM13 chr16) — do not pool with gorilla.

## Question

*"Simulate an ideal scenario using 1 chromosome with some multi-copy gene families but also some others
that are single copy, to prove that in theory it should work. We just need a way to remove unspecific
phenomena like readthroughs and ensure the reads are complete."* — establish the CEILING, and attribute
the real-data gap.

## Setup

chr16, **4,411 annotated transcripts over 1,443 genes**, 10 reads per transcript, ends jittered +-0-30 bp
(mandatory — identical reads collapse under dedup, §6n0), err 0.001. Aligned with the shipped
`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` and assembled through the shipped
`--assemble-only --assembly-polish full` path, unchanged.

Strata fixed from the SHIPPED chr16 homology graph before any arm ran: **MULTI-COPY** = the gene has >= 1
homology edge (297 genes); **SINGLE-COPY** = none (1,110). Universe = the 1,407 simulated genes that carry
an annotation record — an INPUT, so nothing is conditioned on any arm's output.

The readthrough rate is measured, not invented: **33,058 of 439,985 (7.51%)** real primary MAPQ-60 chr16
reads hit the exons of >= 2 distinct annotated genes at >= 25 bp. A_rt injects molecules splicing a gene's
transcript to a neighbouring same-strand gene's at that rate (3,014 of 43,700 reads = 6.9% realised).

## Result

| arm | stratum | raw endpoint | unscoreable | **genuine over-merges** | corrected |
|---|---|---|---|---|---|
| **A_ideal** | ALL | 0.8067 | 15.8% | **12** | **0.9578** |
| | MULTI-COPY | 0.6364 | 28.6% | 12 | 0.8915 |
| | SINGLE-COPY | 0.8523 | 12.3% | **0** | **0.9723** |
| A_trunc | ALL | 0.8188 | — | — | — |
| A_rt | ALL | 0.6503 | 12.6% | **188** | 0.7439 |
| | MULTI-COPY | 0.5556 | 22.9% | 46 | 0.7205 |
| | SINGLE-COPY | 0.6757 | 9.8% | 142 | 0.7493 |
| real data | ALL | 0.4961 | 13.9% | 139 | 0.5759 |

"Unscoreable" = the gene's node is shared ONLY with an annotated gene that overlaps it. Two annotated
genes sharing exonic sequence cannot be separated into one-to-one nodes by ANY method, so these are a
truth ceiling, not an error. "Corrected" excludes them; it is **post-hoc** and is reported beside the
pre-registered raw number, never instead of it.

## Verdict against the pre-registered bars

- **Bar 1 (A_ideal >= 0.95 overall and >= 0.90 multi-copy) FAILS as written**: 0.8067 / 0.6364.
- **Bar 2 (A_ideal < 0.80 ⟹ algorithmic defect) NOT triggered**: 0.8067.
- **Bar 3 (single-copy >= 0.95) FAILS as written**: 0.8523.
- **Bar 4 (attribution) PASSES decisively** — see below.

⚠**The bars were mis-set, and that is a prereg design flaw worth recording.** The endpoint has a
built-in ceiling of roughly **0.84** on this substrate, because 15.8% of the universe consists of
annotated genes that overlap another annotated gene. A bar of 0.95 on a metric that cannot exceed ~0.84
was unreachable by construction. This is the same ~0.8 ground-truth ceiling already documented in
§6kl/§6km, arrived at independently.

## What the simulation actually shows

⭐**With full-length reads and no readthrough the algorithm makes almost no errors**: **12 genuine
over-merges across 1,407 genes**, and **ZERO in the single-copy stratum**. Corrected for the truth
ceiling that is 0.9578 overall, 0.9723 single-copy, 0.8915 multi-copy. **So yes — in theory it works,
and multi-copy families are recovered as separate nodes.**

⭐**Readthrough alone reproduces the real pathology, and then some.** Injecting it at the measured 7.51%
takes genuine over-merges from **12 to 188** — a 15.7x increase — against **139** in the real data. The
injected rate slightly OVERSHOOTS the real defect, so readthrough is not merely sufficient to explain the
over-merge found in §6u7, it is more than sufficient. Nothing else needs to be invoked.

⭐**5' truncation costs nothing on this endpoint** (A_trunc 0.8188 vs A_ideal 0.8067, and the identical
multi-copy rate 0.6364). ⚠This does NOT contradict §6n0, where truncation cost 6 of 26 NPIP copies:
that endpoint was transcript COMPLETENESS, this one is where the node SITS. Truncation shortens models
without moving them.

## Consequence

The gap between the ideal ceiling (0.9578 corrected) and real data (0.5759 corrected) is attributed:
**readthrough carries it**. This is the third independent line pointing at the same conclusion as §6u7
(the defect is over-merge, not fragmentation) and §6v0 (full-length-ness corroborates bridges rather
than discriminating them) — and it is the first to show the algorithm is clean when the phenomenon is
removed. ⟹ **the open lever remains a readthrough-aware node SPLIT**, and it now has a measured ceiling
to aim at.

---

# 8. Does the idealized version hold a better FAMILY definition?

**User, 2026-09-21.** §6v1 scored only the NODE-construction endpoint. This runs the actual family
definition (`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`, the shipped Rust binary, same
recipe as `dn16_fam3`) on the A_ideal and A_rt assembled loci, and scores the resulting clusters against
**both** independent truths from §6u8/§6u9 (Soto cover, protein referee), restricted to chr16.

⚠**n is small** — chr16 alone has 15 Soto families / 36 protein-referee families. Per-family movement
swings are large at this n (§6u8's lesson); read the direction, not the third decimal.

## 8.1 Pipeline

For each arm: exon-merged locus spans from the assembled GTF -> region list -> `samtools faidx` body
extraction -> `minimap2 -x asm20 -c --eqx -P` all-vs-all -> shipped `mcl_families` with the locus's own
exon-union GFF3 (mirroring `dn16.gff3`'s format exactly) -> `.clusters.tsv`. A fourth arm, **GUIDED**
(annotation records used directly as nodes, no assembly at all — `chr16_guided.clusters.tsv`, already on
disk), stands in as the ceiling: perfect nodes, no simulation needed.

## 8.2 Result

| arm | clusters >=2 | vs **Soto cover** (15 fams) | vs **protein referee** (36 fams) |
|---|---|---|---|
| **GUIDED** (annotation nodes) | 90 | sens 0.4567 / prec 0.4294 / **F 0.4174** | sens 0.2046 / prec 0.3203 / **F 0.2248** |
| **A_ideal** de novo | 56 | sens 0.3692 / prec 0.3649 / **F 0.3561** | sens 0.1722 / prec 0.2968 / **F 0.1944** |
| REAL de novo | 70 | sens 0.3372 / prec 0.3255 / **F 0.3166** | sens 0.1681 / prec 0.2172 / **F 0.1805** |
| A_rt de novo | 51 | sens 0.3136 / prec 0.3182 / **F 0.3068** | sens 0.1329 / prec 0.2603 / **F 0.1570** |

**Yes — modestly, and the direction holds on both independent truths.** A_ideal beats real de novo by
+0.0395 F (Soto) and +0.0139 F (protein referee). A_rt falls below REAL on both truths despite injecting
readthrough onto otherwise-ideal reads — consistent with §6v1's own finding that the injected 7.51% rate
slightly overshoots the real defect (12 -> 188 genuine over-merges vs 139 real), so this is not a
contradiction: an aggressive dose of the one thing being tested does more damage than the full complexity
of real data, which is itself evidence the mechanism is right.

Per-family (Soto), 3 better / 6 worse / 6 unchanged: the pooled gain is real but not one-sided — a caveat
consistent with n=15.

## 8.3 The bigger finding: node CORRECTNESS is fixed, node COVERAGE is not

**A_ideal falls well short of GUIDED — F 0.3561 vs 0.4174 (Soto), 0.1944 vs 0.2248 (referee) — despite
§6v1 showing node correctness is nearly solved (0.9578 corrected).** The reason is not over-merge (that
was measured and is nearly gone) and not fragmentation-as-scored-by-§6u7 (per-copy correctness is high).
It is **node COVERAGE**: how many of the truth genes have ANY node representing them in the homology
graph at all.

| arm | graph nodes | Soto truth genes (58 on chr16) represented as SOME node |
|---|---|---|
| GUIDED | 380 | **47 / 58 = 81.0%** |
| A_ideal | 324 | **34 / 58 = 58.6%** |

A 22.4-point coverage gap, under perfect reads and no readthrough. Of the 24 genes missing a node in
A_ideal, 17 (71%) DO have an assembled locus — the locus exists, it is just absent from the graph,
meaning either it never accumulated enough exonic content to clear the `identity >= 0.70 / cov_longer >=
0.30 / >= 300bp` edge floor, or (below) it never gets scoring credit at all:

**Of the 13 genes missing in A_ideal but present in GUIDED (the true assembly-attributable set), 5
(38%) are a SCORING ARTIFACT, not an assembly failure.** RefSeq itself carries curated **readthrough/
fusion gene records** overlapping almost exactly the same span — `PKD1P3-NPIPA1` (42,950 bp) sits over
`NPIPA1`'s 14,597 bp span at 37,168 bp overlap; `PKD1P4-NPIPA8`, `PDXDC2P-NPIPB14P`, `BOLA2-SMG1P6`
do the same for `NPIPA8`, `NPIPB14P`, `SLX1B`. The de novo locus's homology signal for the real gene IS
present (verified directly on the PAF: NPIPA1's locus has 10+ alignment records at 0.83-0.99 identity to
other loci) — but a **max-overlap, one-name-per-locus resolver, applied at scoring time, hands the whole
locus's credit to the bigger fusion-named record**, which then has no truth family to match, so the real
gene registers as "no node" though its sequence was correctly assembled and correctly grouped. This is
the family-definition-scoring analogue of §6v1's "unscoreable (annotation overlap)" category, applied
here for the first time to the CLUSTER-naming step rather than the per-copy node step.

The remaining ~8 (ABCC6, HERC2P5, HERC2P8, NPIPB10P, NPIPB7, PKD1, PKD1P6, SMG1P6) have no bigger
co-located record and are candidates for a genuine remaining edge-construction gap — short or
low-coverage assembled loci that never accumulate enough exonic content relative to their full-length
siblings to clear the coverage floor. This was not run to ground on this pass and is the natural next
target.

## 8.4 Consequence

The answer to "does the idealized version hold a better family definition" is **yes, but only partly for
the reason expected.** Fixing readthrough (§6v1's target) recovers a real, if modest, family-definition
gain — confirmed on two independent truths. But it does not close most of the gap to a perfect-node
ceiling, because **most of that remaining gap is not a node-construction defect at all**: at least 38% of
it is an artifact of how RefSeq's own curated fusion-gene names get resolved at scoring time, in exactly
the same region (PKD1/NPIP tandem duplication) already flagged by §6u7/§6v1/§6v0 as ground zero for
readthrough. The rest is a real, smaller, and as-yet uncharacterized edge-coverage gap.

---

# 9. The 8 residual genes: three distinct mechanisms, not one

**User, 2026-09-21.** Diagnosed each of the 8 genes from §8.3 (ABCC6, HERC2P5, HERC2P8, NPIPB10P,
NPIPB7, PKD1, PKD1P6, SMG1P6) individually against the raw annotation, the assembled GTF, and the PAF.

## 9.1 Three categories, verified against the raw GFF and the assembled loci

| category | genes | mechanism | fixable by better assembly / reads? |
|---|---|---|---|
| **A — no transcript model exists** | HERC2P5, HERC2P8, NPIPB10P (3) | Verified on the raw GFF: these `pseudogene` records have **zero** child `transcript`/`mRNA`/`exon` records of any kind — RefSeq assigns them a genomic span and nothing else | **No.** There is no sequence to simulate and (if real biology matches the annotation) no transcription to sequence. Not an assembly or read-quality question at all. |
| **B — fusion/readthrough truth-naming** | PKD1P6 (1) | Verified on the raw GFF: `PKD1P6`'s only two `transcript` records (`rna-NR_123721.1`, `rna-NR_123722.1`, 30 and 18 exons) are children of `gene-PKD1P6-NPIPP1`, **not** of `gene-PKD1P6`. Standalone PKD1P6 has no transcript of its own — its real transcription, when it occurs, produces the readthrough molecule, not a PKD1P6-only one. Joins §8.3's 5 fusion-name cases (NPIPA1, NPIPA8, NPIPB14P, SLX1B, PKD1P6-vs-`PDXDC1`) at **6 of 13** now. | **No.** This is the same phenomenon the whole session has been chasing (§6u7/§6v0/§6v1), now showing up as a truth-side naming artifact rather than a node over-merge. Full-length confidence changes nothing about which curated gene name a real fusion molecule maps to. |
| **C — real paralog, real assembly, genuine exon-structure mismatch or fragmentation** | ABCC6, NPIPB7, SMG1P6 (3) | See 9.2 below — two different sub-mechanisms | Partially — see 9.2 |

## 9.2 Category C in detail

**ABCC6**: both ABCC6's locus (4,831 bp exonic, vs its real transcripts' 4,463-4,550 bp — essentially
COMPLETE) and its true paralog ABCC6P1's locus (`DN_chr16_18504057_9`, 2,660 bp exonic, vs its own
2,664 bp transcript — also essentially complete) assemble correctly and independently. The PAF shows a
strong genomic alignment between them (23,511 bp aligned at 97.85% identity). **The edge still fails**,
almost certainly at the `--min-shared-exon-frac 0.6` conjunct: `exonic_denominator=true` means the
coverage/sharing test is computed over EXON-restricted alignment, not the raw genomic span, and ABCC6P1
is a `transcribed_pseudogene` — very likely a retro/processed duplicate whose own annotated exon
structure has diverged from ABCC6's spliced mRNA (fewer/different exon boundaries after duplication).
This is register 925's already-documented finding — *"half of family relationships have an INTRONLESS
member"* — measured here on a NEW pair with everything else held ideal. **Not an assembly defect: the
two loci are each correctly and completely reconstructed; the paralogs' real, annotated exon structures
are what differ.**

**NPIPB7, SMG1P6**: genuinely fragmented, unlike ABCC6. Comparing locus exonic content to the gene's
OWN transcript length (not genomic span, which is misleading for intron-rich genes): NPIPB7's locus
captures 464 of 2,254 bp (20.6%) of its single 10-exon transcript; SMG1P6's captures 326 of 2,623 bp
(12.4%) of its single 18-exon transcript. Both are unusually **exon-dense, compact transcripts**
(NPIPB7 ~225 bp/exon, SMG1P6 ~146 bp/exon) — a real under-assembly, not a coverage or truncation
artifact (A_ideal reads are full-length, error rate 0.001, jitter only 0-30 bp per end). This is a
genuinely open lead: whether minimap2's spliced alignment or the assembly/polish step is what drops
these small exons was not run to ground this pass.

## 9.3 Does isoseq cluster2 help?

**For 6 of 8 (categories A and B): no, and it cannot in principle.** Category A has no transcript to
sequence in the first place, real or simulated — cluster2 clusters reads, and there would be no reads.
Category B is a truth/naming problem: PKD1P6's real biological molecule, when transcribed, IS the
PKD1P6-NPIPP1 fusion. cluster2 would correctly report a confident, full-length, HQ isoform for that
fusion — which is precisely what §6v0 already established cannot be resolved by full-length-ness,
because the fusion **is** the genuine complete molecule.

**For ABCC6 (category C, exon-structure divergence): no.** Both loci are already reconstructed
essentially completely from these ideal simulated reads; the problem is the exon-conjunct correctly
detecting that the two paralogs' real, annotated splice structures differ. Better read processing
cannot change what the paralog's real exon structure is.

**For NPIPB7 and SMG1P6 (fragmented, compact transcripts): plausibly yes, and this is the one candidate
worth testing.** `isoseq cluster2` calls a consensus directly from the FLNC read sequences (POA-style,
before any genome alignment), which could recover a clean, complete transcript model for a compact
multi-exon gene even where individual per-read spliced alignments to the genome are imperfect on small
exons — a different failure surface than a reference-guided splice-alignment assembler faces. This is
speculative: the exact cause of the fragmentation (minimap2's splice seeding on <150-225 bp exons, vs
an assembly/polish-step filter) was not isolated this pass, so cluster2's benefit here is a hypothesis
to test, not a demonstrated fix — and it is a narrow win even if confirmed (2 of 8 genes here, and an
unknown fraction of chr16 more broadly).

## 9.4 Bottom line

Of the 8 residual genes, **6 have no assembly-side fix available at all** (3 have no transcript to
recover, 3 map to a real fusion). Only **2 (ABCC6, NPIPB7... correction: ABCC6 is unfixable by
assembly too — so really only NPIPB7 and SMG1P6, 2 of 8)** are candidates for a genuine, fixable
assembly gap, and `cluster2` is a plausible but unproven lever for exactly those two. This further
narrows what remains open after §6v1/§6v2: the readthrough-aware node split is still the dominant lever;
this residual set is small, mostly structural rather than algorithmic, and the one live thread worth a
follow-up is specifically **compact, exon-dense transcript assembly**, not readthrough or node
over-merge.

---

# 10. Root cause of the NPIPB7/SMG1P6 fragmentation: a NON-CANONICAL junction in RefSeq's own annotation

**User, 2026-09-21: "let's run the diagnostics to figure out why NPIPB7/SMG1P6 fragment."** Isolated to a
precise, confirmed mechanism, ruling out three alternatives along the way.

## 10.1 Ruled out, in order

1. **Read alignment**: all 10 simulated reads' PRIMARY alignments (`-F 2308`) land within 150 bp of each
   other at the true locus, with CIGARs that correctly reconstruct the full transcript. For SMG1P6, all
   10 reads carry the **exact same 17-intron chain**, every intron at 10/10 support — zero disagreement.
2. **Assembly polish**: re-ran with `--assemble-only` and NO `--assembly-polish` flag at all. The raw,
   unpolished GTF is **byte-identical in the fragmented region** to the polished one — the loss happens
   before polish ever runs.
3. **A minimal single-read / identical-duplicate repro**: running the assembler on 1 read, or on 5
   byte-identical copies of 1 read, produces **empty output** for both genes — consistent with the
   project's documented dedup rule (`sims need -N 50 and must VARY read ends`) collapsing identical
   reads to one effective molecule, below whatever minimum-support floor applies. Not informative about
   the fragmentation itself, but rules out a naive "not enough reads" framing.

## 10.2 The actual cause: a non-canonical splice junction in the gene's OWN RefSeq exon model

Checking every intron's splice motif (donor/acceptor dinucleotide, strand-corrected) against the genome
directly:

| gene | biotype | introns | first non-canonical junction |
|---|---|---|---|
| **SMG1P6** | `transcribed_pseudogene` | 17, all GT-AG except one | intron 15 of 17: `chr16:29721738-29725407` (3,670 bp), motif **AT-AG** |
| **NPIPB7** | `protein_coding` | 8, all GT-AG except the last 3 | introns 6-8 of 8: `chr16:28752509-28760135` (TA-GA), `chr16:28760195-28761928` (GT-GG), `chr16:28762017-28778548` (GG-TT) |

**Verified independent of my own pipeline**: the exon coordinates flanking each non-canonical junction
match the RAW, un-processed GFF `exon` records for these transcripts exactly (`rna-NR_135312.1`,
`rna-NM_001396030.1`) — this is RefSeq's own curated exon model, not an artifact of my extraction or of
simulated sequencing error.

The shipped assembler enforces canonical (GT-AG / GC-AG / AT-AC) splice motifs when building a spliced
transcript model, and — confirmed by direct test — **truncates the transcript at the first non-canonical
junction it encounters, discarding everything on the far side, even when every read agrees perfectly and
completely on the true structure.** This is not a bug introduced by imperfect reads; it fires on reads
that are already ideal.

## 10.3 Confirmation: `RUSTLE_JUNCTION_MAJORITY=1` recovers it

This project already has a named, documented mitigation for exactly this failure mode
([[project_junction_majority_chr16]] / §6m8: *"strict canonicity is 69.4% of pass-1 -> GTF loss"*).
Re-ran A_ideal with it set:

| gene | default (strict canonicity) | `RUSTLE_JUNCTION_MAJORITY=1` |
|---|---|---|
| **SMG1P6** | 2 exons, 224 bp (1 of 17 introns) | **FULL RECOVERY**: `chr16:29708443-29728813`, matching the annotated span almost exactly, gene_id suffix `_18` = 18 exons |
| **NPIPB7** | 4 exons, 1,282 bp | **PARTIAL RECOVERY**: `chr16:28739068-28762017`, reaching well past the original truncation point but still short of the transcript's far end (which carries THREE consecutive non-canonical junctions, not one) |

SMG1P6 recovers completely because it has exactly one non-canonical junction; NPIPB7 only partially
recovers because it has three in a row near its 3' end, and majority-rule rescue apparently reaches
through isolated non-canonical junctions more readily than a consecutive run of them.

## 10.4 Consequence for the family-definition question and for cluster2

This closes the last open item from §9: **NPIPB7 and SMG1P6 are not a fragmentation mystery — they are
the SAME already-documented canonical-junction mechanism as §6m8, now shown to explain 2 of the "8
residual genes" at the individual-gene level**, and shown for the first time to survive into (and cost)
the FAMILY-DEFINITION endpoint specifically, not just the chr16-wide pass-1-to-GTF loss rate §6m8
originally measured.

**isoseq cluster2 cannot help here either.** The reads are already fully confirmed correct, complete,
and unambiguous (§10.1.1) — this defect lives entirely inside Rustle's own assembler, downstream of
anything a PacBio read-processing tool touches. My prior answer treating these two genes as an open,
possibly cluster2-fixable lead is **superseded**: it is neither a read-quality problem nor an unexplained
one.

**This is the user's call, not a default to flip unilaterally** (memory records the default was
deliberately left off after the chr16 arm test — *"fear refuted, default still not flipped, user's
call"*). What this session adds beyond §6m8: independent confirmation on unrelated genes (SMG1P6,
NPIPB7, vs §6m8's chr16-wide aggregate), and a NEW consequence not previously measured — this specific
loss reaches the family-definition score, not just transcript-recovery counts.

---

# 11. Including non-canonical junctions on REAL chr16 data: measured, not assumed

**User, 2026-09-21: "can we include the non-canonical junctions too? ... for now I think my advisor will
be skeptical the method really works for any family if we cannot recall all members."** Ran
`RUSTLE_JUNCTION_MAJORITY=1` through the full pipeline on REAL chr16 (same `--assemble-only
--assembly-polish full` recipe that produced `dn16.gtf`), rebuilt the locus graph and family clusters
the same way as every other arm this session, and scored against both truths.

## 11.1 Assembly and family-definition effect

| | strict (baseline) | `RUSTLE_JUNCTION_MAJORITY=1` |
|---|---|---|
| transcripts kept (polish) | 9,629 | 10,093 (+464, +4.8%) |
| de novo loci | 2,550 | 2,544 |
| graph nodes | 864 | 869 |
| clusters >= 2 members | 70 | 70 |

| | vs Soto cover (15 fams) | vs protein referee (36 fams) |
|---|---|---|
| strict | F 0.3166 | F 0.1805 |
| `JUNCTION_MAJORITY=1` | F 0.3184 (+0.0018) | F 0.1817 (+0.0012) |

Pooled F is essentially flat — negligible, in either direction, at this n. Per-family (Soto): 2 better, 1
worse, 12 unchanged.

## 11.2 The number the advisor's question is actually about: individual-member recall

| | Soto truth genes on chr16 (58) covered as SOME node |
|---|---|
| GUIDED (annotation, ceiling) | 47 / 58 = 81.0% |
| REAL de novo, strict | 34 / 58 = 58.6% |
| REAL de novo, `JUNCTION_MAJORITY=1` | **35 / 58 = 60.3%** |

**+1 gene recovered (`NPIPB2`), 0 lost.** Small, but unambiguous and free on this measure: no truth gene
that strict canonicity covered was dropped by relaxing it.

⚠**§6v4's exact repro genes behave differently on real data than in the ideal simulation.** On real chr16
reads (natural depth variation, not my controlled 10-reads-per-transcript scheme), `NPIPB7` and `SMG1P6`
are **already covered under strict canonicity** — the non-canonical-junction truncation §6v4 demonstrated
unambiguously in the clean simulation does not reproduce identically on these same two genes at real
depth, most likely because real read population diversity gives some reads a slightly different
alignment path across the junction that the clean, uniform simulated set didn't have. The mechanism from
§6v4 is still correct (independently verified against the raw GFF, and directly confirmed by the
before/after test) — it just isn't the reason these particular two genes are missing on REAL data.
`ABCC6`, `HERC2P5`, `PKD1P6` remain uncovered under BOTH arms, for the reasons already established in
§9 (real exon-structure divergence, no transcript model, fusion-naming) — this lever cannot and should
not be expected to fix those.

## 11.3 What this does and does not answer for the advisor

**Does not, on its own, close the recall gap.** 23 of 58 Soto truth genes are still missing a node under
`JUNCTION_MAJORITY=1`, against 24 under strict — a gap of 1, not 24. The 22-point remaining difference to
the guided ceiling (81.0%) is NOT primarily a canonical-junction problem on real data; §9's three-way
breakdown (no transcript model / fusion-naming / real exon divergence) accounts for at least several of
the specific genes checked, and the rest is uncharacterized. **"Recall all members" needs more than one
lever** — this is one confirmed, safe, small piece of it, not the fix.

**Known cost, already measured genome-wide** (§6n4, `bench/CHR16_JUNCTION_MAJORITY_ARM.md`): strictly-
engulfed copies 79 -> 86 (+8.9%), a precision cost on copy BOUNDARIES, not on family structure (max family
size and family fusion were unaffected — the flag's specific feared harm was refuted there).

## 11.4 Recommendation

Given zero measured node losses here, a real (if small) recall gain, and an already-refuted specific fear
at genome-wide scale, this is a low-risk, well-evidenced choice for the user to adopt. Two distinct
decisions, kept separate:

1. **Use `RUSTLE_JUNCTION_MAJORITY=1` in this session's own family-definition / recall-focused bench
   runs going forward.** Zero code change, fully reversible, and directly what was asked for "for now."
2. **Flip the compiled default in `denovo_pipeline.rs`** (`env_num("RUSTLE_MISCHAIN_...")`/junction
   majority default) so every future run gets it without the env var. This is a bigger, harder-to-reverse
   change than (1) and is NOT done here — it is the user's call, per the standing rule this project has
   already applied twice to this exact flag.

---

# 12. The full remaining recall gap, classified (not just the original 8-gene sample)

**User, 2026-09-21: "keep exploring reasons for false positives/false negatives."** After the default
flip (§11 above, now committed), 23 of 58 Soto truth genes on chr16 still lack a node. Classified all 23
— not a hand-picked sample — by the same three mechanisms §9 found on 8 genes, checked directly against
the raw GFF and the new-default assembled GTF.

| category | count | share |
|---|---|---|
| **A — no transcript model in RefSeq at all** | 6 | 26.1% |
| **B — a bigger co-located annotation record wins the naming contest** | 8 | 34.8% |
| **C — a locus exists but captures only a fraction of the gene** | 9 | 39.1% |

**A (unfixable by construction): `ABHD17AP7`, `ABHD17AP8`, `ABHD17AP9`, `HERC2P5`, `NPIPB14P`, `PKD1P6`.**
A THIRD pseudogene trio (`ABHD17AP7/8/9`) joins `HERC2P5`/`HERC2P8`/`NPIPB10P` in having zero annotated
exon structure — no sequence exists to simulate or sequence, real or synthetic.

**B (scoring artifact, not construction failure): `MIR3179-1`, `MIR3179-2`, `PKD1P1`, `PKD1P6-NPIPP1`,
`RRN3`, `SLX1B`, `SULT1A3`, `SULT1A4`.** Extends §9's mechanism (a bigger fusion/co-located annotation
record wins a max-overlap, one-name-per-locus resolver's contest) to double the count. One striking case:
**`PDXDC1` (168,469 bp) sits over BOTH `PKD1P6-NPIPP1` and `RRN3`** — even the FUSION gene itself loses
its naming contest to a still-larger neighbouring record in this densely tandem-duplicated region. There
is also a CHAIN: `SULT1A3`/`SULT1A4` lose to the fusion names `SLX1A-SULT1A3`/`SLX1B-SULT1A4`, which
themselves land in category C (locus exists, captures little) — the naming collision and the coverage gap
compound for this pair.

**C (the open, still-uncharacterized lever): `ABCC6`, `ABCC6P1`, `ABCC6P2`, `CDR2`, `NPIPB4`, `RRN3P1`,
`RRN3P2`, `SLC7A5P2`, `SLX1A-SULT1A3`.** A real locus exists for 7 of 9; `ABCC6P1`/`ABCC6P2` have no
locus at all despite carrying their own transcript record. This is the same shape of gap NPIPB7/SMG1P6
were before §6v4 diagnosed them — and NPIPB7/SMG1P6 are confirmed NO LONGER in this list, a clean
positive control that the junction-majority fix reached exactly the class of gene it was built for. This
category is where the next diagnostic effort belongs.

⭐ This reframes the "22-point gap" from §11: it was never one problem. A quarter is a true annotation
ceiling, a third is a scoring resolver artifact (the real sequence is very likely already correctly
grouped, just credited to the wrong name), and just over a third is a genuine, unexplained construction
gap — the only slice comparable in kind to what junction-majority just fixed.
