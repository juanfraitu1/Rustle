# Collapse gate: admitting a single-rep locus as a multi-copy family — Design

**Date:** 2026-07-09. **Substrate:** gorilla (GGO) HiFi Iso-Seq, `denovo_pipeline::colocated_families` (O2).

## The defect

`colocated_families` requires `>= min_copies` **assembled reps**. A collapsed locus has one rep, so it is
dropped — even when its reads plainly carry several copies. On chrY this is not hypothetical:

- **DAZ**: 2 annotated copies. After the readthrough filter the window leaves **one** rep (DAZ1). DAZ2's reads
  exist — 20 primaries (19 at MAPQ 0) and 1119 secondary placements over its span — piled onto DAZ1 because the
  copies are near-identical. Result: **0 families, 0 assignments.**
- **BPY2**: same picture. Both copies produce nothing.

`recover_collapsed_copies` already reports **3 PSV-distinguishable collapsed copies** at the DAZ rep. Its return
value is `eprintln!`'d and discarded.

## Why the obvious fix fails

Counting PSV-distinguishable haplotypes into the family gate — "copies, not reps" — was measured and refuted
(`bench/COLLAPSED_COPY_GATE.md`). Gorilla is diploid, and the haplotype detector cannot tell two het alleles from
two copies:

| locus | collapsed candidates | truth |
|---|---|---|
| DAZ | 3 | 2-copy, collapsed |
| **TSPYL1** | **12** | **single-copy** |
| DERPC | 5 | single-copy |

A single-copy gene scores four times higher than the true collapsed pair. The haplotype count is not noisy here;
it is anti-correlated with truth.

## The papers: collapse first, variants second

**SDA (Vollger et al., Nat Methods 2019)** — the closest prior art, and the paper the advisor sent. It solves
exactly this, in two stages, *both by depth*. Verbatim, Methods:

> "We defined collapsed regions as those with a mean sequence coverage > 3 s.d. beyond the mean sequence coverage
> of the de novo assembly and that were at least 9,000 bp in length."

> "Three thresholds were applied to determine whether an SNV was also a PSV. First, the total depth at the given
> position had to be at least the mean coverage plus 3 s.d. Second, the frequency of the second-most frequent base
> had to be less than the mean coverage. Finally, the frequency of the second-most frequent base had to be greater
> than the mean coverage minus 3 s.d. or half the mean coverage, whichever was greater. **This process favors the
> selection of PSVs over allelic variants.**"

And the main text states the ordering explicitly:

> "We begin by identifying all collapsed duplications within each assembly based on an excess of sequencing read
> depth. … Next, we define PSVs corresponding to each collapsed segment. … requiring sequence coverages consistent
> with a single-copy locus in order to distinguish PSVs from allelic variants."

**We ran stage 2 without stage 1.** That is the entire bug.

**Bailey et al. (Science 2002)** — depth as the copy-number instrument, the WSSD ancestor of our `depth_cn`:

> "…a simple statistical test to determine whether a given stretch of sequence is duplicated based on its
> overrepresentation and average sequence identity…"

**Clair3-RNA (Zheng et al., Nat Commun 2025)** — why SDA's instrument cannot be transplanted verbatim:

> "the coverage is uneven across genomic regions in RNA-seq, and the variable coverage poses challenges for
> accurate variant calling, particularly in regions characterized by inadequate read support or excessive read
> coverage…"

and, on the allelic-balance signal that would otherwise separate a het allele from a PSV:

> "…zygosity flipping can happen."

In DNA, depth **is** copy number. In RNA, depth is copy number × expression, and allele-specific expression
destroys allelic balance. So we keep SDA's *structure* (collapse first, variants second) and replace its
*instrument* with one matched to the modality.

**IsoCon (Sahlin et al., Nat Commun 2018)** contributes the statistical style already adopted elsewhere in this
codebase: a per-variant-position real-vs-error hypothesis test rather than a tuned threshold. The collapse test
below is written in that shape.

## The instrument: ambiguity, not depth

In DNA a collapse shows up as **excess depth**. In RNA it shows up as **reads the aligner cannot place uniquely**.
Same phenomenon, different observable. Measured (primary records only, `-F 2308`):

| locus | reads | MAPQ-0 | fraction | truth |
|---|---|---|---|---|
| DAZ2 | 20 | 19 | **0.950** | 2-copy, collapsed |
| TSPY | 34 | 30 | **0.882** | 6-copy |
| RBMY (distal) | 87 | 12 | 0.138 | 6-copy |
| DAZ1 | 200 | 22 | 0.110 | 2-copy, collapsed |
| GSTM1 | 1074 | 0 | 0.000 | 3-copy family, **copies resolvable** |
| TSPYL1 | 2151 | 0 | 0.000 | single-copy |
| DERPC | 1512 | 0 | 0.000 | single-copy |
| ATXN7L3B | 1585 | 0 | 0.000 | single-copy |
| GSPT2 | 685 | 0 | 0.000 | single-copy (retrogene) |
| EEF1A1 | 3516 | 0 | 0.000 | single-copy (has a retrocopy) |

**Zero MAPQ-0 primaries across 9,449 control reads.** The statistic is expression-invariant: TSPYL1 has 2151 reads
and no ambiguity; DAZ2 has 20 reads and 95%. Depth would have ranked TSPYL1 far above DAZ2 — the exact confound
Clair3-RNA names, and the reason the haplotype-count fix failed.

GSTM1 reads 0.000 despite being a genuine 3-copy family, because its copies are divergent enough to map uniquely.
That is correct behaviour, not a miss: the signal detects **collapse** (unresolvable placement), not
multi-copy-ness, and GSTM has three assembled reps so the gate is never consulted.

## Design

One new decision point; no new pipeline. The gate runs **only on loci that fail the existing rep test**, so the
multi-rep path (GSTM, MAGEA, RBMY) is untouched and byte-identical.

```
colocated_families:  >= min_copies REPS ?
      yes ─────────────────────────────────────► family (unchanged)
      no  ──► collapse gate
                 leg 1: is this locus collapsed?      (ambiguity test)
                 leg 2: how many haplotypes?          (split_locus_copies, chi_h)
                 └─► family with n_copies = chi(H), reads certified TIED
                     (--admit-collapsed-copies: materialise copies + assign)
```

### Leg 1 — the collapse test (runs first)

A unique locus places its reads unambiguously. The statistic is `k` = the count of ambiguously-placed **primary**
reads at the locus (MAPQ 0, or carrying an AS-tied secondary placement) out of `n` primary reads. Under the null
"this locus is unique", `k ~ Binomial(n, eps_amb)` and the gate fires when the upper-tail p-value clears the
family-wide Bonferroni level already used by the assignment gate. A **p-value, not a cutoff** — the shape the
shipped significance gate uses, and the shape IsoCon uses.

**`eps_amb` must not be estimated as zero.** The controls give 0 ambiguous reads in 9,449, whose MLE is 0, under
which a *single* stray MAPQ-0 read would be infinitely significant. Estimate `eps_amb` per region from the reads
of that region's **uniquely-mappable reps** (reps that are not gate candidates), with a Jeffreys prior —
`eps_amb = (k_bg + 1/2) / (n_bg + 1)`. This keeps the rate strictly positive, adapts to sample-specific
mismapping, and reduces to a background of roughly `1/(2 n_bg)` when a region is clean. If a region has **no**
uniquely-mappable rep to estimate from, the gate must **abstain** rather than fire on an unbounded statistic.

Primary records only (`-F 2308`): secondary alignments are ambiguous by construction and would make every locus
look collapsed. Supplementary records are excluded for the same reason they are excluded from conflict edges — a
chimeric segment is adjacency, not ambiguity.

### Leg 2 — the haplotype count (conditional on leg 1)

Reuse `copy_split::split_locus_copies` and count `identifiable` haplotypes; obtain χ(H) from the existing
`readonly_copy_number::chi_h`. Nothing new is built.

### Output

**Default:** the locus becomes a family with `n_copies = χ(H)`, its reads **certified tied**, and
`depth_cn` / `famcn_readonly` populated. No per-copy consensus is materialised, so **no assignment can be wrong**.
This is the recorded division of labour: `count(χ(H))` is O1, `assign(Strong-Sep)` is O2.

`min_copies` applies to the **copy count, not the rep count**: a gated locus forms a family iff
`χ(H) >= min_copies`. A collapsed locus resolving to a single haplotype is not a family, and is dropped exactly as
today. The emitted `n_copies` is χ(H); `rescued_copies` and `collapsed_copies` record how it was obtained, so a
gated family is never mistaken for an assembled one in the output.

**`--admit-collapsed-copies` (opt-in):** additionally build per-copy consensus from the identifiable haplotypes and
assign reads under the usual significance gate. Off until the simulation below says it can be trusted.

## Validation, in order

1. **Mechanism (planted sim, non-circular).** `bench/sim_genome.py::plant_collapsed` already plants this regime.
   Score: does leg 1 fire on planted collapses and stay silent on planted unique loci (false-positive rate); does
   χ(H) equal the planted copy number.
2. **The discriminating experiment.** Extend the sim to plant **2 copies × 2 het alleles** and **4 copies** with
   matched read counts and matched PSV columns. Ask whether the gate separates them. **If it cannot, that is the
   finding**, and `--admit-collapsed-copies` stays off permanently. This is the copy-vs-allele wall, planted.
3. **Real data.** `asm_hapCN` on the 70 oracle-matched families (the existing copy-number oracle). Sanity only:
   DAZ, BPY2, RBMY, TSPY must produce families; TSPYL1, DERPC, GSPT2, ATXN7L3B, EEF1A1 must not.

## Non-goals

- **Reference-absent copies.** A copy missing from the reference produces no second placement, hence no ambiguity
  signal. Detecting it requires the depth-excess channel, which needs an RNA expression model we do not have, and
  which memory already records as DNA-bound. Out of scope; the gate must not claim to see it.
- No change to the multi-rep path, to PSV discovery, to assignment, or to the EM.
- No attempt to *estimate* copy number. χ(H) is a **lower bound**, and the theory note already says so:
  `max(n_loci, χ(H)) <= true`.

## Known failure mode, stated up front

Two copies × two alleles yields four haplotypes — which is what DAZ shows (3 emitted candidates = 4 identifiable).
Leg 1 establishes that a collapse is present; it cannot say **how many** copies collapsed. Therefore χ(H) is
reported as a lower bound, reads are tied by default, and copy materialisation is gated behind a flag whose licence
to exist is experiment 2 above.

## Files

- **Modify** `src/rustle/vg_family/denovo_pipeline.rs`: the collapse gate at the `colocated_families` rep test;
  stop discarding `recover_collapsed_copies`.
- **Modify** `src/bin/copy_assign.rs`: `--admit-collapsed-copies`.
- **Modify** `bench/sim_genome.py`: plant het alleles alongside collapsed copies (experiment 2).
- **Create** `bench/COLLAPSE_GATE_VALIDATION.md`.

## Reproduce the motivating numbers

```
samtools view -c -F 2308 GGO_mm.bam NC_073248.2:42879918-42945552   # DAZ2: 20 primaries, 19 at MAPQ 0
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --region NC_073248.2:42778133-42950552 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --recover-copies --out DAZ
# -> "collapsed copies recovered past family gate (PSV-distinguishable): 3", then "0 families"
```

Related: `bench/COLLAPSED_COPY_GATE.md` (the refutation), `bench/YAG_CHECK.md`, `bench/PRIOR_ART_LEVERS.md`,
`reference_sda_vollger`, `reference_bailey_2002_segdups`, `reference_clair3_rna`, `reference_isocon_sahlin`.
