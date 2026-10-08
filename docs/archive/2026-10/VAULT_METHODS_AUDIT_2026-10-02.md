# ThesisVault methods audit, addendum — 2026-10-02

**Question:** the four papers added to `ThesisVault/` on 2026-10-01 (Shaw 2025 devider, Hosseini 2025 pHapCompass,
Chaisson 2017 polyploid phasing, Bosch 2007 FAM90A), plus the O3 chain built since the 09-28 audit
(`docs/archive/2026-09/VAULT_METHODS_AUDIT_2026-09-28.md`): is there a method we have not applied that would move O1/O2/O3?

**Answer:** one. devider's positional de Bruijn graph (PDBG) over (PSV site, allele) letters is a published,
threshold-light replacement for the read-clustering + consensus + merge core of the `o3_candidates` stage
(`docs/superpowers/specs/2026-10-02-o3-candidates-design.md` §5.3–5.5), and it has not been tried
(no hit for `de bruijn|pdbg|devider` in `src/`, `bench/`, `docs/`, the register, or memory). Everything else in the
four papers is already in the code, already refuted, or a framing point. The 09-28 verdict stands for the other 20 papers.
Method: read of the four summary notes and the devider Methods section; grep of `src/ bench/ docs/`, the register and memory
for each method; no experiment run.

## 1. devider (Shaw, Boucher, Yu, Noyes, Li; Genome Research 2025) — NOT APPLIED, worth a pre-registered head-to-head

What it does: reads aligned to ONE reference sequence + a VCF of sites become strings of (site index, allele); a PDBG over those
strings (k-mers collapse only when sites and alleles match, so the graph is a DAG); unitigs are aligned back to the graph by exact
DP and dropped by a one-sided binomial test against their best alternative path (error classes del / ref→alt / alt→ref at fixed
rates .35/.15/.10); reads are aligned to the filtered graph, a read with two equally good paths is assigned to NONE; haplotypes =
well-supported paths, each with abundance, assigned reads and a base-level consensus; haplotype count is NOT preset. Rust,
`hi-fi` preset (M = 100, α = 500, ρ = .001). GitHub bluenote-1577/devider.

Why it maps onto the O3 chain exactly:

| `o3_candidates` step (spec §5) | devider equivalent | difference that matters |
|---|---|---|
| §5.3 greedy first-fit template clustering of WHOLE reads at `de` ≤ δ, length-sorted | PDBG over PSV sites; a read contributes only at the sites it covers | 5′-truncated / sub-chain reads (the real-NPIP IsoCon failure: "reads too truncated/sparse", Amendment 11 cause) no longer need to clear a 50 %-coverage pairwise test; order-independent |
| §5.4 majority vote on the template | base-level majority consensus of the reads assigned to a path | same |
| §5.5 pairwise binomial merge of clusters (eps^k) | binomial test of a unitig against its BEST ALTERNATIVE PATH found by DP, then ρ-dedup of haplotypes on unambiguous sites | the "alternative" is the optimum over the whole graph, not one pairwise partner |
| abstention | read with ≥ 2 equally good paths assigned to none; score < 0 unassigned | the assign-or-abstain shape O2 uses, inside the candidate stage |
| isoform structure | an exon absent from a read is a run of missing sites = `del`, UNPENALISED in path finding and ignored by ρ-dedup | copies split on sequence (ref→alt / alt→ref), not on structure: the Pal 2026 structure/sequence decoupling for free (untested on RNA; must be measured) |

What it needs from us: per family, (a) one reference sequence to align the net to (the family's exon-union of reference copies, or
the longest copy as Pal 2026 did), (b) a VCF of sites: LoFreq via its wrapper, or our PSV catalog written as VCF, (c) the family read
net of spec §5.1 (the thing Amendment 11 showed is indispensable). The flag / link (δ) / merge / ≥ 2-cluster floor / exon-union
representative of §5.6–5.7 stay unchanged downstream.

Known limits (from the paper, not from a run): DNA / cDNA amplicons and metagenomes only, no RNA, no splicing; cannot "recover new
sequences de novo" (a copy that differs from the rep by an insertion is representable only at shared sites, same as our link rule);
fixed error constants and a 2/−3/−5/−1 score are heuristics (say so to Canzar: borrow the DAG + exact DP + abstention, not the
constants); copies without private sites merge, which is the abstention the chain already accepts; built for 20 haplotypes at
19,500×, per-copy RNA depth is lower (its floor is 5 reads, ours 3).

Proposed test (write `docs/PREREG_devider_candidates_*.md` first): same 53-family Aug-excision panel, no-deletion control, the real
reference-absent set (GWFAM175_B0) and both YAG panels, same metrics (D-right %, false moves, C1 false-flag family rate, M2 one
candidate per missing copy), arms = IsoCon (done) / in-house §5.3–5.5 (being built) / devider. BAMs and truth exist under
`/mnt/linuxdisk/tmp/rna_allele/`; `cargo` is installed. If devider matches IsoCon's 74 % D-right with ≤ IsoCon's false moves, the
in-house clustering can be re-specified as a PDBG; if not, the row goes to the register and the spec stands.

## 2. The other three papers — nothing new to apply

| Method | Source | Status |
|---|---|---|
| Discrete matrix completion (A·B factorisation, MEC, EM refinement) | Chaisson 2017 | = `copy_assign --em` (Vollger soft relaxation); changed 0/3,081 decisions (register 579, 617). Ploidy must be given, so no O3 use. |
| Correlation clustering of a PSV graph, cluster number as output | Chaisson 2017, SDA | Applied (`vg_family`, `psv_graph`). Register 740 rejected its clusters as "within-locus PSV haplotypes, not a copy count" — under the 10-01 framing (count haplotypes, then link at δ to separate alleles from copies) that object is the intended intermediate, so `psv_graph` over the family read net is a second, zero-cost candidate generator for the §1 head-to-head. Not a new method. |
| "Unmatched PSV cluster = O3 flag" | Chaisson 2017 (vault note's inference) | = the chain's flag + link step (Amendments 7–8). |
| Genotype-count constraint pruning of phasings | pHapCompass | Needs K and genotypes; copy number is the unknown. Not transferable. |
| Distribution over phasings / uncertainty quantification | pHapCompass | Probabilistic weighting within the candidate set is what O2 withholds by design; Canzar prefers the combinatorial rule. Not adopted. |
| Partial-phase metrics: charge unphased blocks at chance + penalty per break | pHapCompass | Already our trap T1 ("never condition the denominator on the prediction", `feedback_metric_traps`); `copy_assignment_definition.md` §9 names the same defect. A citation for the scorer, not a method. |
| Subfamilies sharing CDS but differing in first exon; CNV; gene-conversion mosaics | Bosch 2007 | Biology, no method. L4 (0.995) is the nested level; gene-conversion mosaics are the K ≥ 3 recombination obstruction already machine-checked. |

## 3. Items re-read in the light of the O3 chain (not new papers)

- **Cross-pipeline agreement as the candidate filter** (LRGASP: calls made by > 50 % of pipelines validated at 100 %). The chain's
  no-deletion control FAILED C1 (30.2 % vs 28 %, LR 2.75) and was rescued post hoc by a ≥ 2-cluster floor (LR 8.2). Requiring a
  candidate from ≥ 2 of {IsoCon, in-house, devider} is the LRGASP form of that floor and costs nothing once §1 runs. Cheap, untested.
- **Allelic divergence as the independent source for the O1 cut** (Guitart's 1.5 × intra-allelic rule). §6p3 refuted every
  self-tuned L3 cut and asked for an independent source; the chain now measures exactly that quantity (δ = p99 allelic divergence of
  single-copy genes on KB3781's own haplotypes = .00958, identity .990, between L3 .985 and L4 .995). §6kq (09-14, throwaway) found
  the allelic-scaled grouping recovers TBC1D3-CDKL. Not a metric lever (grouping is saturated, §6o8) but it answers the advisor's
  "arbitrary threshold" objection with a measured constant. Would need a prereg; low priority.
- **isONclust** (Tomaszkiewicz 2023's family-level read clustering) is minimizer-greedy — the advisor dislikes minimizers, and the
  net already attributes reads by k-mer sharing (§5.2). Not adopted.

## Caveats

- Absence of a grep hit is "not applied" only after checking the register; a method under another name could have been missed.
- Nothing here was run. The devider claims about spliced reads (`del` tolerance, ρ-dedup across isoforms) are read off the Methods
  and are the first thing the pre-registration must test.

## RESULT (2026-10-03, branch `devider-arm`): the §1 head-to-head was run and devider is REFUTED on all six pre-registered clauses

Prereg `docs/PREREG_devider_arm_2026-10-02.md`; results `docs/DEVIDER_ARM_deletion_2026-10-02.md`, `_control_2026-10-03.md`,
`_refabsent_2026-10-03.md`; register drafts `docs/REGISTER_DRAFTS_devider.md` (all on the branch). Deletion held-out: D-right 23.5 %
vs IsoCon 74.0 % (bar 59.2 %), false moves 0.72 % vs 0.06 %, M2 5/9. Control: C1 45.3 % vs 30.2 %, C2 0.49 % vs 0.05 %. Real
reference-absent: D1 fails (the B0 copy is recovered as one contig at identity .9939, under the .999 class rule and the >= 2-transcript
floor). Cause: the method is reference-bound; in 30/53 deletion families 0 % of the missing copy's reads align to any surviving copy,
exactly the families Amendment 11 found. Where the premise holds it separates perfectly (GWFAM175 499/500). The `o3_candidates`
spec stands; the Python port (phase 2) was not started. Still untested: a PDBG as the within-cluster phasing step of a reference-free net.
