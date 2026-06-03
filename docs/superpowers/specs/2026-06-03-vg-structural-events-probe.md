# Can the VG detect inversions / segmental duplications / complex copy events? (2026-06-03)

Grounded probe + controlled-synthetic proof. Verdict: the VG has genuine DETECTION/
CLASSIFICATION primitives for every structural-event type, at varying maturity. Detection
is real (and VG-unique); RESOLUTION (assembling/flagging the rearranged copy) is mostly unbuilt.

## Inversions — DETECTION PROVEN, resolution unbuilt
Controlled synthetic (gen_synthetic.py --invert-copies): copies 0,1 forward + copy2 placed as
a true inverted duplication (revcomp genome; its 30 reads align to the **- strand** at 401kb).
- **PROVEN:** VG discovers all 3 copies as ONE family ("3 linked bundles, 12 multi-map reads").
  The gene's mRNA matches copy0/1 on + and copy2 on - , so multimap reads BRIDGE the inversion —
  a single-strand assembler would treat the inverted copy as unrelated. Mechanism: canonical
  k-mer hashing (min of fwd/revcomp, vg.rs:450) + the multimap mRNA. **This is genuinely VG-unique.**
- **NOT resolved:** the inverted copy assembles to 0 transcripts. The fingerprint-EM is
  `decisive=1 uncertain=11` (reads ambiguous across the inversion), and partition_family_by_strand
  splits the - strand copy off; it does not reach output. The inversion relationship is also not
  surfaced (no inversion tag). RC-aware exon UNIFICATION across the inversion is the unbuilt piece
  (audit theme C).

## Segmental duplications / tandem arrays / translocations — CLASSIFIED (real-data confounded)
`locus_relationship` (vg.rs:2020) already classifies every family's inter-copy spatial relation,
emitted as `family_locus_rel`: **Trans** (different chromosomes → translocation-like), **Distal**
(gap ≥ 1 Mb → dispersed segmental duplication), **Tandem** (gap < 1 Mb → tandem array/expansion),
**Overlapping** (same coords). On the clean synthetic this is correct (200 kb copies → Tandem).
- **Caveat (real data):** the emitted class is dominated by **Overlapping** (chrY 55/64, GOLGA8
  262/272) — an ARTIFACT of per-strand bundling: each locus has a + and a - bundle at the same
  coords (strand mirrors) that "overlap." So the real-data spatial class is confounded and would
  need de-mirroring before it's trustworthy. A segmental duplication's EXTENT (gene + flanking)
  vs a bare paralog is also not distinguished — only the gene's exons are modeled.
- Copy NUMBER (duplication multiplicity) = family copy count; depth-based estimate is audit #3
  (expression-confounded, diagnostic only).

## Gene conversion / mosaic recombination — BUILT
The mosaic-read detector (commit fc6f087) detects recombinant reads whose per-site copy pattern
switches — gene conversion observed directly in one read. (Coverage limited to family-linked
reads pending the assembly-stage hook.)

## Honest bottom line
YES — the VG detects these structural events, and I can prove inversion detection on controlled
ground truth + the spatial classification exists. The boundary is DETECTION vs RESOLUTION:
- Detected today: inverted-paralog membership, spatial class (tandem/distal/trans), gene-conversion
  reads, copy count.
- NOT yet built: assembling/flagging the inverted copy (RC-aware unification), de-mirroring the
  real-data spatial class, breakpoint calling, segdup-flanking extent. These are the realizable
  next builds if structural resolution (not just detection) is the goal.

## UPDATE: inversion RESOLUTION fixed (2026-06-03, commit da1f24d)
The inverted copy now ASSEMBLES (was: detected-only). Grounding corrected the framing — the
blocker was NOT missing exon unification (the inverted copy assembles cleanly de-novo). It was
the **joint-strand EM** (the DAZ1−/DAZ3+ same-locus handler) being wrongly engaged on a
DISPERSED inverted duplication: mixed-strand families have no joint graph, so the joint EM
apportions the inverted copy's reads onto the forward copies with a neutral fingerprint and
starves it. Fix: use the joint family only for single-strand + OVERLAPPING inverted pairs;
use the un-jointed (per-strand) partitions for dispersed inversions (new helper
`mixed_strand_copies_overlap`; opt-out RUSTLE_VG_INVERSION_JOINT=1). Validated: inversion
synthetic 100/100 (inverted − strand copy assembles by default); chrY unchanged (DAZ
overlapping → joint kept); GOLGA8 +3 (real dispersed inverted copies recovered, incl. − strand
@45.6M); suite 204/0. So of the structural-event boundary, INVERSION moves from "detected" to
"detected + resolved." Remaining: de-mirror the real-data spatial class; breakpoint/segdup extent.

## UPDATE 2: spatial classification DE-MIRRORED (2026-06-03, commit d0182ee)
The "Overlapping"-confound is fixed. `family_locus_rel` was dominated by Overlapping (GOLGA8
262/272, chrY 55/64) because per-strand bundling puts a +/- mirror pair at each locus and the
old check returned Overlapping if ANY two bundles overlapped — one mirror pair flipped a whole
13 Mb dispersed segdup to "Overlapping". Fix: `classify_physical_loci` merges overlapping bundle
intervals into DISTINCT PHYSICAL loci before classifying (genuine same-locus inverted pairs still
read Overlapping). Result (annotation only): GOLGA8 262 overlapping -> 259 distal (+4 tandem);
chrY 55 overlapping -> 62 distal. The classification now correctly reads GOLGA8/chrY paralog
families as DISPERSED segmental duplications. 4 unit tests; suite 208/0.

Structural-event map now: inversion = detected+resolved; spatial class (segdup/tandem/trans) =
classified + de-confounded; gene conversion = detector built. Remaining: breakpoint calling and
segdup-flanking extent (distinguishing a segdup's duplicated flanks from a bare paralog).

## UPDATE 3: segdup extent / breakpoint calling BUILT (2026-06-03, commit a8113c7)
Distinguishes true segmental duplication (gene+flanks) from bare paralog and calls the
duplication breakpoints (where cross-copy GENOME homology ends; mRNA doesn't cover flanks).
Validated: planted flank 600→600bp / 400→399-401bp exact; bare→~0/not-segdup; real GOLGA8
5 segdups (multi-kb flanks) vs 4 bare paralogs. Opt-in RUSTLE_VG_SEGDUP_EXTENT, default-off,
suite 213/0. Doc: 2026-06-03-vg-segdup-extent-breakpoint.md.

STRUCTURAL-EVENT MAP — FINAL: inversion detected+resolved; spatial class classified+de-mirrored;
gene-conversion mosaic detector built; segdup extent/breakpoints built. The probe's full
checklist (inversions, segdups, complex copy events) is now addressed end-to-end.
