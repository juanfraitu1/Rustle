# Pre-registration — do the gorilla individuals' own NPIP copy sequences and/or 5′-truncated reads reproduce the real NPIP loss?

**Written 2026-09-29 00:55 PDT, before any variant is called and before any read is generated.** GORILLA ONLY
(OR6737 testis `GGO_mm.bam`, KB3781 fibroblast `GCA_029281585.2_flnc_mm.bam`; the two samples are never pooled).
Dev diagnostic: nothing in `src/` is edited, no default changes, nothing is committed. Scratch
`/mnt/linuxdisk/tmp/rustle_figures_dev/ggo_npip_varsim/`.

## What is already known (not predictions; `ggo_npip_loss` 2026-09-29, `ggo_npip_sim` 2026-09-28)

- **Real signature, per sample (OR / KB).** Truth = the 25-copy T_member (`ggo_npip_sim/ann/truth.json`).
  - Presence (a same-strand shipped transcript with exonic overlap on the copy territory): 11/25 / 8/25, so
    14 / 17 absences.
  - First losing step of the absences: pass-1 floor 3 / 5; gate 8 / 8 (strict junctions + '+' single-exon
    placeholder 4 / 5, placeholder only 3 / 3, strict junctions only 1 / 0); polish mono floor 3 / 4.
  - Complete copies (one transcript carrying the whole primary chain): 0 / 0.
  - Own spliced primaries at the copies carrying ≥ 1 non-canonical junction: **0.698 / 0.725** (control 0.029 /
    0.034). Median `de`: 0.0098 / 0.0071 (reads with a non-canonical junction), 0.0059 / 0.0032 (canonical-only).
  - **Recurrent non-canonical set R** (`recur.py`, a non-canonical junction carried by ≥ 1 read of a ≥ 2-read
    skeleton on a copy territory in BOTH libraries): 12 junctions on NC_073242.2 — 1 × `gt..gc` 168 bp (NPIPB8,
    inside an annotated exon), 4 × `TT..AC` 2,958-2,965 bp, 7 × `tt..ga` 7,384-8,434 bp (+-strand motifs; the
    copies are '-' strand, so on the transcript strand these read GT..AA and TC..AA).
- **Annotation check done while designing (not an outcome):** no junction of R is an annotated intron of the raw
  gorilla GFF. Each `TT..AC` junction shares its right end (the '-'-strand donor) exactly with an annotated intron
  whose left end lies 1,119-1,137 bp further left; each `tt..ga` gap lies inside an annotated exon (the long
  3′-terminal exon); the 168-bp gap lies inside an annotated exon.
- **Prior sim.** Error-free full-length reads (FL, 30 per transcript) give 25/25 complete; 5′-truncated reads (TR,
  3′-anchored, 5′ start uniform) give 12/25 complete. Neither was at real depth, and neither carried the
  individuals' sequence.

## Question

Is the real loss (presence, first losing step, completeness) and its junction signature (non-canonical fraction,
the recurrent set R) reproduced by (a) the individuals' own copy sequences (variants vs the reference copy), (b)
the real 5′-truncated read-length distribution, (c) only their combination — or by neither (a structural cause)?

## Frozen instruments

- Binaries: `/mnt/linuxdisk/tmp/rustle_figures/cc_bin_frozen/` (= HEAD `main@a9797aee` + source.diff, the binary the
  diagnostic used and validated byte-for-byte): `copy_assign` 57168f3f0ffb4870e2d28b30be8e5752762c39d5, `as_table`
  f2abaede055052ec35ace9f98bcc7b5a5082b1ba, `mcl_families` 748ea0bdfcb989c473b694c6cd69a1cf6e9ddee3.
- Mapping: `minimap2 2.30-r1287 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4` on
  `winloci_data/GGO.splice.mmi` (whole genome, 13,649,916,704 bytes, 2026-07-16). ⚠ The real BAMs were made by
  minimap2 2.31-r1302 with the same options on the same genome; the version difference is a stated limit.
- Genome `GGO.fasta` (registry `gorilla_fasta`); annotation = the prior sim's `ann/tx.tsv` FL rows (131 transcripts,
  25 copies; each copy's native gorilla record clipped to the territory; non-canonical annotated introns merged, as
  in the prior sim).
- Scoring code: copies of the diagnostic's `trace.py`, `polish_port.py`, `arms_score.py`, `attribute.py`,
  `noncanon.py`, `dediv.py`, `recur.py`, with **path-only edits** (arm directory, BAM, as_table from environment
  variables). The diff against the originals is recorded in the Outcome. **Instrument check (must pass before any
  sim arm is scored):** the copies, pointed at the real BAMs and the diagnostic's own GTFs, reproduce the
  diagnostic's 31 absences with the same first-losing-step label per copy.

## Step 1 — the individuals' copy sequences (variant calling rule, fixed now)

- **Reads used, per sample and copy:** the copy's *own primaries* exactly as the diagnostic defines them — primary,
  not supplementary, mapped, orientation (FLAG 0x10 under `-uf`) = copy strand, ≥ 1 aligned base (`cigar_exons`
  blocks) on the copy territory. Multi-copy handling: **only reads whose primary is at that copy are used**; a
  read's secondaries never contribute; a read counts for one copy only (territories do not overlap). No MAPQ
  filter (MAPQ 60 is only 40-65% at the copies; a MAPQ-60 filter would remove the very reads under test). The
  count of these reads is the copy's **real primary depth**.
- **Pileup:** over every reference position those reads align (exons, and introns/flanks wherever a read is
  unspliced there); `=`/`X`/`M` give a base, `D` gives a deletion; `N` (a splice gap) and soft clips give no
  coverage. Coverage at a position = reads with a base or a deletion there.
- **Call rule (majority allele at a pre-stated depth):** at coverage ≥ 4, a non-reference allele is called when it
  is carried by ≥ 3 reads AND by > 50% of covering reads. Alleles: a substituted base (SNV); a deletion with an
  identical start and length (`D` op, counted at its first base); an insertion with an identical inserted sequence
  after the same reference base (counted over reads aligned on both flanking bases). A deletion masks SNV calls
  inside its span. Heterozygous sites at ≤ 50% are therefore not called (the V copy is the majority haplotype).
- **Not called (by design):** `N`-encoded gaps. minimap2's splice mode writes any long deletion as `N`, so a genomic
  deletion ≥ ~30 bp in the individual looks exactly like a non-canonical intron. Calling such gaps as deletions and
  putting them into V would make V reproduce them by construction (circular). They are the phenomenon under test;
  if V cannot produce them, the verdict is "structural" and the reads at those gaps are characterized (below).
- **Reported per copy and sample:** number of SNVs / insertions / deletions; how many are in simulated exons,
  in introns inside the copy span, or in flanks; how many sit near a splice site (in an annotated splice
  dinucleotide, in the splice region = 3 exonic / 8 intronic bp, or ≤ 20 bp from an exon boundary); the fraction
  of the copy's primary-transcript exonic bases that are callable (coverage ≥ 4).
- **R explainability:** for every junction of R, the called alleles and coverage at its four motif bases and within
  ± 10 bp of each end. "Explained by a called variant" = substituting the called alleles makes the motif canonical
  (GT-AG / GC-AG / AT-AC) on the copy strand; "covered, not explained" = coverage ≥ 4 at the motif bases and no such
  call; "unobservable" = coverage < 4 there (the motif bases are intronic and only unspliced reads cover them).

## Step 2 — the four arms (2 × 2), per sample

For each copy, **n = its real primary depth** in that sample (copies with 0 primaries get 0 reads). Read *i* of copy
*c* is paired across the four arms: the same transcript (drawn uniformly at random, seeded by sample/copy/i, from the
copy's FL transcripts), the same end jitter and the same length.

| arm | copy sequence | read length |
|---|---|---|
| **C** | reference | full length, both ends jittered 0-30 bp (the prior FL rule: only when the body stays > 100 nt) |
| **T** | reference | the real length distribution: read *i* takes the aligned query length (`query_alignment_length`, soft clips excluded) of the copy's *i*-th real own primary (a seeded permutation), 3′-anchored at the transcript end minus a 0-30 bp jitter; if that length ≥ the transcript, full length |
| **V** | variant copy (reference + the sample's called variants at that copy, applied to the exons through the indel offsets) | full length, as C |
| **VT** | variant copy | real length distribution, as T |

- Error model: HiFi-like, `bench/sim.py simulate_reads` with sub 0.001 / indel 0.0003 (the house HiFi rate of the
  tandem and copies simulators), applied after cutting. ⚠ The prior `ggo_npip_sim` FL/TR reads were error-free, so
  the prior FL arm is NOT reused (its generator and depth also differ: 30 reads per transcript).
- Alignment: the mapping above, in batches under `tools/rlock.sh heavy`, each batch sized to finish inside the
  600 s lock timeout.

## Step 3 — assembly in the real context (spike-in)

The losing steps depend on context the copy's own reads do not carry: pass-1 skeletons pool every kept record in
the window, and the polish mono floor is a quantile over the whole contig (OR 12/13/15, KB 16/16/17). A BAM of
simulated reads alone would change the mono floor by construction. So each arm is assembled on a **spike-in BAM**:

- background = the sample's real records on NC_073241.2 / NC_073242.2 / NC_073244.2 minus every record of every
  read that is an own primary at a copy (exactly the reads whose count set the depth);
- plus every record of the arm's simulated reads on those three contigs;
- as_table = the real genome-wide table (the diagnostic's `runs/{S}.molecules.tsv`) followed by the as_table of the
  arm's simulated BAM (whole genome, `cc_bin_frozen/as_table`);
- **arm E (control)** = the background alone (no reads at the copies) — what the context alone gives.

Assembly = the diagnostic's `arm.sh` commands with the BAM and as_table swapped: `ship` (the shipped command:
seeding 0.98, strict junctions, full polish, identical to the driver's `assemble` on these contigs), `nopol` and
`maj` (the two extra arms `attribute.py` needs).

## Step 4 — scoring (same code as the diagnostic)

Per arm and sample: presence (`arms_score.py`); the first losing step per absent copy (`trace.py` →
`polish_port.py` → `attribute.py`); completeness; the non-canonical fraction of own spliced primaries at the copies
and the control (`noncanon.py`); per-read `de` (`dediv.py`); reads carrying the junctions of R (exact coordinates)
and the two classes `TT..AC` 2,900-3,050 bp and `TT..GA` 7,300-8,500 bp (case-insensitive, + strand motif) at any
copy (new `recurrent.py`); source placement (a simulated read's primary on its source copy's territory). Family
membership (driver `families` = `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units` on the
ship GTF, then `npf_audit` placement) is SECONDARY and run only if one run on the real 3-contig ship GTF finishes
inside one 600 s heavy-lock slot; the real 3-contig GTF is then its baseline.

## Pre-stated decision rules

Components, each required in BOTH samples:

- **J — junction signature.** J1: the fraction of own spliced primaries with ≥ 1 non-canonical junction is within
  real ± 0.15 (OR 0.548-0.848, KB 0.575-0.875). J2: ≥ 2 reads carry a `TT..AC` 2,900-3,050 bp junction AND ≥ 2
  reads carry a `TT..GA` 7,300-8,500 bp junction at the copies. J = J1 ∧ J2.
- **L — loss signature.** L1: |absent − real absent| ≤ 3 (real 14 / 17). L2: with the step categories {no reads /
  seeding / dedupe, pass-1 floor, gate (any label with junctions or placeholder), polish, other}, half the L1
  distance between the arm's and the real category counts ≤ 3. L = L1 ∧ L2.

Verdict mapping:

| verdict | condition |
|---|---|
| **sequence sufficient** | V or VT reproduces J and L |
| **truncation sufficient** | T reproduces L and V does not (T carries reference sequence and annotated canonical models, so it cannot produce J; J is then judged separately below) |
| **interaction** | only VT reproduces J and L |
| **structural / not explained** | neither V nor VT satisfies J2; the reads at the R junctions are then characterized: block sizes on each side, whether the gap removes annotated exon sequence (missing segment) or adds unannotated sequence (extra exon / intron retention), and where the read ends lie |

The L verdict and the J verdict are reported separately when they differ (e.g. "truncation sufficient for the loss;
junction signature structural"). Descriptive, not decisive: per-copy present/absent concordance with real, complete
copies, `de`, placement, family membership, arm E.

## Predictions (committed now)

- **C:** non-canonical fraction ≤ 0.05; J2 absent; 3-7 absences per sample (copies whose few reads spread over
  several isoforms fall to the pass-1 floor), so L1 fails; complete ≥ 10/25.
- **T:** non-canonical fraction ≤ 0.05; J2 absent; 9-16 absences, via pass-1 floor, '+' placeholder and mono floor;
  L1 passes in at least one sample; complete ≤ 5/25.
- **V:** non-canonical fraction < 0.15 (J1 fails); J2 absent; loss close to C.
- **VT:** close to T; J1 and J2 fail.
- **Predicted verdict:** loss "truncation sufficient" (possibly only partly: T's gate losses come from the
  placeholder, not from strict junctions); junction signature "structural / not explained".

## Falsifiers

- V's non-canonical fraction ≥ 0.55 in both samples → the junctions are an aligner response to point/small-indel
  divergence, and the "structural" prediction is wrong.
- V or VT produces J2 → the recurrent pair is sequence-driven.
- C already reproduces L → the loss is depth and isoform spread, neither truncation nor sequence.
- Arm E gives presence at ≥ 5 copies → the background alone decides presence at those copies, and the per-copy
  comparison must be read net of E.

## Hostile self-review

1. **Variants come from the reads under test.** Reads mis-assigned from another copy (or from an unreferenced copy)
   make the "individual's copy" a majority blend; hets are dropped. This biases V toward the dominant source, which
   is what the aligner sees anyway.
2. **Coverage is 3′-heavy.** Variants can be called only where reads align, so V ≈ C over the 5′ parts. The callable
   fraction is reported per copy, and a V failure is only informative where the copy is callable.
3. **`N`-encoded deletions are excluded by design.** V can then produce R only through the aligner's response to
   point/small-indel sequence, not by copying the gaps; this is the point of the test, but it means "structural" is
   also the verdict for a simple ≥ 30 bp deletion in the individual.
4. **Annotated isoforms only.** Unannotated exons, intron retention and alternative 3′ ends in the individuals
   cannot appear in any arm; the characterization of R is where they are looked for.
5. **T is 3′-anchored at the annotated 3′ end.** Real reads can end upstream (OR NPIPB8: all 66 primaries in the first
   5 kb of the 8.75 kb terminal exon), so T reads may be more often unspliced than real. The real 3′-end offsets are
   reported next to the T result.
6. **The spike-in keeps real non-own reads** (antisense, secondaries of reads placed elsewhere); arm E measures what
   they alone produce.
7. **Thresholds are mine** (± 0.15, ≤ 3 copies, ≥ 2 reads); n = 2 samples; post-hoc dev contigs.
8. **Sufficiency, not necessity.** An arm that reproduces the signature shows that cause is sufficient in this
   model; it does not show the real data arose that way.
9. **minimap2 2.30 vs 2.31** for the real BAMs.

## Frozen code and instrument check (recorded 00:58, before any variant is called or any read generated)

Code in `ggo_npip_varsim/code/` (sha1): `callvar.py` 592846d1, `build_reads.py` d000f3a8, `common.py` 52a2b6fa,
`recurrent.py` 382864df, `map_batch.sh` be09bd0f, `spike.sh` 8d459c39, `arm.sh` c6ed4cbb; the diagnostic's scoring
copies (path-only diff: `D` and the BAM paths read from `VS_D` / `VS_BAM_OR` / `VS_BAM_KB`): `trace.py` a556b01f,
`polish_port.py` 854c4d78, `arms_score.py` 970356be, `attribute.py` d0b6b5f7, `noncanon.py` 4b7d3eff, `dediv.py`
48dae80d.

**Instrument check PASSED:** the copies, run on the real BAMs with the diagnostic's GTFs (`arms/R`), give
`trace.{OR,KB}.json` identical to the diagnostic's, polish port "ALL IDENTICAL" on both samples, presence 11 / 8,
and the same first-losing-step and completeness label for 50/50 copy-samples. Real J2 baseline (`recurrent.py R`,
own spliced primaries): `TT..AC` 2.9-3.05 kb 7 (OR) / 38 (KB) reads, `TT..GA` 7.3-8.5 kb 5 / 23 reads; R at exact
coordinates 29 / 51 reads.

---

# OUTCOME (2026-09-29, run exactly as registered; scratch `ggo_npip_varsim/`)

## Answer first

- **Pre-registered verdict: STRUCTURAL / NOT EXPLAINED.** Neither V nor VT produces the recurrent pair (J2: 0 reads
  of either class in every arm and sample, against 7 / 5 (OR) and 38 / 23 (KB) real reads), and no arm puts the
  non-canonical fraction in the real band (J1).
- **The loss is not reproduced by any arm** (L fails in at least one sample for every arm). Truncation accounts for
  most of the excess absences but through the wrong steps; the individuals' point/small-indel sequence adds 0-2.
- **No falsifier fired.** V's non-canonical fraction is 0.29 / 0.30 (bar 0.55); V/VT produce no J2 read; C does not
  reproduce L; arm E (background only) gives 0/25 present copies in both samples.

| arm | sample | absent (real) | steps early/pass-1/gate/polish | L1 | L2 (½L1) | P/A concordance | complete | spliced own primaries | non-canonical fraction (J1) | TT..AC / TT..GA / R reads (J2) |
|---|---|---|---|---|---|---|---|---|---|---|
| real | OR | 14 | 0/3/8/3 | — | — | — | 0 | 364 | 0.698 | 7 / 5 / 29 |
| real | KB | 17 | 0/5/8/4 | — | — | — | 0 | 313 | 0.725 | 38 / 23 / 51 |
| C | OR | 4 | 1/2/1/0 | ✗ | ✗ (6.0) | 15/25 | 17 | 584 | 0.267 ✗ | 0 / 0 / 0 ✗ |
| C | KB | 3 | 1/2/0/0 | ✗ | ✗ (8.0) | 11/25 | 17 | 757 | 0.262 ✗ | 0 / 0 / 0 ✗ |
| T | OR | 11 | 1/1/9/0 | ✓ | ✗ (3.5) | 14/25 | 2 | 40 | 0.100 ✗ | 0 / 0 / 0 ✗ |
| T | KB | 12 | 1/1/7/3 | ✗ | ✗ (3.5) | 20/25 | 2 | 81 | 0.358 ✗ | 0 / 0 / 0 ✗ |
| V | OR | 4 | 1/2/1/0 | ✗ | ✗ (6.0) | 15/25 | 17 | 584 | 0.288 ✗ | 0 / 0 / 0 ✗ |
| V | KB | 5 | 2/2/1/0 | ✗ | ✗ (8.0) | 13/25 | 16 | 758 | 0.299 ✗ | 0 / 0 / 0 ✗ |
| VT | OR | 11 | 1/1/9/0 | ✓ | ✗ (3.5) | 14/25 | 2 | 84 | 0.417 ✗ | 0 / 0 / 0 ✗ |
| VT | KB | 13 | 1/1/8/3 | ✗ | ✓ (3.0) | 21/25 | 2 | 82 | 0.402 ✗ | 0 / 0 / 0 ✗ |
| E | OR / KB | 25 / 25 | 22/3/0/0 | ✗ | ✗ | 14 / 17 | 0 | — | — | — |

"early" = no reads / seeding / dedupe; the one early copy in C-VT is the 0-primary copy (OR NPIPA7, KB NPIPA9), absent
by construction. Median `de` of own spliced primaries (non-canonical / canonical-only): real 0.0098 / 0.0059 (OR),
0.0071 / 0.0032 (KB); C 0.0034 / 0.0016, 0.0035 / 0.0016; V 0.0062 / 0.0024, 0.0055 / 0.0017; VT 0.0046 / 0.0024,
0.0036 / 0.0024. Simulated primaries on their source copy: 94-98% in every arm (MAPQ 60: 69-90%).

## Step 1 — the individuals' copies

- Own primaries 589 (OR) / 759 (KB), identical to the diagnostic's depth; every read counted for one copy.
- Calls: OR 268 (256 SNV, 9 ins, 3 del; 237 exonic, 22 intronic, 9 flank), KB 173 (167 / 3 / 3; 141 exonic, 22
  intronic, 10 flank). Copies with ≥ 1 call: 20 / 11. Median callable fraction of the primary transcript
  (coverage ≥ 4): 0.28 / 0.38 — coverage is 3′-heavy, as feared (hostile review 2).
- Calls concentrate at a few copies: NPIPB4 81 / 78 (76 identical in both animals), NPIPB2 38 (OR only),
  NPIPB11 34 (OR), LOC128966608 31 (OR), NPIPB7 9 / 40, NPIPB15 18 / 25. 95 calls are shared between the animals,
  76 of them at NPIPB4 — the reads placed on NPIPB4 carry the same foreign haplotype in both gorillas.
- Near a splice site (dinucleotide, splice region or ≤ 20 bp from an annotated exon boundary): 17 (OR) / 18 (KB),
  17 of them the same NPIPB4 cluster in both animals; 1 per sample in an annotated splice dinucleotide (NPIPB4
  21087852, an alternative-isoform boundary).
- **R explainability: 0 / 12 junctions explained by a called variant in either sample.** Covered, not explained:
  OR 1 (the NPIPB8 168-bp gap), KB 3 (that gap; NPIPB8 `tt..ga` 8,434; NPIPB2 `tt..ga` 8,413). All four `TT..AC`
  junctions and the other `tt..ga` junctions are unobservable (motif coverage < 4: their motif bases are intronic).
- V transcripts carry 1,496 / 61 / 20 (OR) and 1,192 / 25 / 16 (KB) applied SNV / ins / del over the simulated
  isoforms; V reads are still 2-3× less divergent than the real reads (canonical-only `de` 0.0024 / 0.0017 vs
  0.0059 / 0.0032).

## What the real reads show at R (the registered characterization; `charR.py`, `ncjunc.py`, `termgap.py`)

- **`tt..ga` 7,384-8,434 bp = a missing segment of the 3′-terminal exon.** At NPIPB4 (26 reads: OR 5, KB 21;
  LOC128966608 and NPIPB3 1 each) every read keeps exactly 54 bp at the annotated 3′ end (ending 17 bp past it) and
  resumes at 21,081,752, the terminal exon's 5′ part: a 7.5-kb segment, 94-96% repeat-masked, is absent from these
  molecules. No alignment-equivalent shift is canonical. Other reads end at the breakpoint itself (21,081,752-6)
  without the 54-bp tail. MAPQ 1-24.
- **`TT..AC` 2,958-2,965 bp = an extra exon.** The reads (NPIPB4 38: OR 6, KB 32; LOC128966608 2, NPIPB3 2)
  carry an unannotated 108-bp exon [21,084,640, 21,084,748) inside the annotated 4.1-kb intron; its donor-side
  junction reuses the annotated acceptor, and the `TT..AC` junction reuses the annotated donor. Only its acceptor is
  non-canonical: AA in the reference, one base from AG. That base is intronic and cannot be observed from RNA.
- **`gt..gc` 168 bp at NPIPB8 = a deletion inside a repeat.** 34 reads (OR 19, KB 15), all MAPQ 60; the gap is 100%
  repeat-masked, inside the 3′ UTR, and its motif bases are covered (42-71 reads) with no variant.
- **The real non-canonical signal is mostly gaps inside annotated exons:** 329/355 (OR) and 277/373 (KB)
  non-canonical instances. 288/364 (OR) and 167/313 (KB) spliced own primaries are "spliced" ONLY by gaps inside
  annotated exons (308 / 235 carry a gap inside the 3′-terminal exon; ≥ 1 kb: 39/589 and 107/759 reads, 13 / 20
  copies).
- **Real reads are 3′-anchored** (86% / 84% end within 100 bp of an annotated 3′ end), so T's anchoring was right.
  T's reads stay unspliced 3′ fragments (7% / 11% spliced) because the reference exon has no gaps; the real reads'
  "splicing" is these structural gaps.

## Predictions vs outcome

| prediction | outcome |
|---|---|
| C: NC ≤ 0.05; J2 absent; 3-7 absences; complete ≥ 10 | ✗ NC 0.267 / 0.262 (see below); ✓ J2 0; ✓ 4 / 3; ✓ 17 / 17 |
| T: NC ≤ 0.05; J2 absent; 9-16 absences; L1 in ≥ 1 sample; complete ≤ 5 | ✗ 0.100 / 0.358 (of 40 / 81 spliced reads); ✓; ✓ 11 / 12; ✓ (OR); ✓ 2 / 2 |
| V: NC < 0.15, J1 and J2 fail, loss ≈ C | ✗ 0.288 / 0.299, but ✓ J1 fails; ✓ J2 0; ✓ 4 / 5 |
| VT ≈ T, J fails | ✓ 11 / 13; ✓ |
| verdict: loss "truncation sufficient", junctions "structural" | loss: **not met** (T fails L in KB, 12 vs 17 absences; L2 3.5 in both); junctions: ✓ structural |

- **Why C is non-canonical at 26%.** Reference-identical full-length NPIP reads are 10-24 kb long.
  - 217/218 (OR) and 296/296 (KB) of C's non-canonical instances are gaps inside annotated exons, i.e.
    chaining artifacts of very long reads over repeat-rich terminal exons. Example: a 22.7-kb NPIPB7 read carries
    3,429 inserted bp balanced by in-exon `N` gaps.
  - Only 1/162 (OR) and 14/173 (KB) of C's distinct junctions occur in the real data.
  - So the long-read regime has an aligner-on-reference floor, but it is not the real signal.
- **Truncation moves absences 4 → 11 (OR) and 3 → 12 (KB); sequence moves them 0 to +2** (C → V, T → VT).
  T's losses are the wrong kind:
  - 9 / 7 are gate placeholder losses (unspliced single-exon 3′ fragments on '-' copies), with 0 junction-clause
    losses against real 5 / 5.
  - T keeps copies that real loses: the pass-1-floor copies NPIPB10P and NPIPA5, and OR's mono-floor copies
    NPIPB14P, NPIPB12 and NPIPB5.
  - T loses copies that real keeps: OR NPIPB4, NPIPA2, NPIPB3 and NPIPB7, which real keeps because its reads there
    are "spliced" by the structural gaps.

## Families (secondary, as registered)

The real OR 3-contig `families` run finished in 504 s, inside one 600 s slot, so every arm was run (KB runs 67-98 s,
OR 481-526 s). Placement of the 25 copies (`famscore.py` = `npf_audit/audit.py run()` unchanged, config override):

| arm | OR: present / in NPIP family / other family / singleton | KB: present / NPIP / other / singleton |
|---|---|---|
| real | 11 / 10 / 1 / 0 | 8 / 8 / 0 / 0 |
| C | 21 / 13 / 7 / 1 | 22 / 20 / 2 / 0 |
| T | 14 / 7 / 7 / 0 | 13 / 7 / 5 / 1 |
| V | 21 / 13 / 7 / 1 | 20 / 19 / 1 / 0 |
| VT | 14 / 7 / 7 / 0 | 12 / 8 / 4 / 0 |

Arm E was not run: it has 0 present copies. Every simulated arm places more copies OUTSIDE the NPIP family than real
does (OR 7 vs 1; KB 1-5 vs 0), because the extra copies it keeps join other clusters. Sequence (C → V, T → VT)
changes placement by at most 1 copy.

## Deviations and notes

- `charR.py`, `ncjunc.py`, `termgap.py`, `summarize.py`, `fam.sh` and `famscore.py` were written after the runs. They implement the registered
  characterization and decision table; they are not new decision rules.
- The real J2 baseline counts own spliced primaries (`recurrent.py R`); R itself was defined from ≥ 2-read skeletons
  that include seeded secondaries, which is why R has 12 junctions but some copies show 0 primary reads carrying theirs.
- Sufficiency, not necessity; n = 2 samples; minimap2 2.30 vs the real BAMs' 2.31.
