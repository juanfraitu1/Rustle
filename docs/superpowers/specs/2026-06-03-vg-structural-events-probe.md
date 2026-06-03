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
