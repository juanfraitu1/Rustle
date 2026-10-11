# Direct-transfrag assembly (no network flow), 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** User (2026-10-11), after `--compat-support` failed
(register 1264): build, from scratch, an assembly that uses transfrags directly — no network flow — the way an earlier (now deleted)
StringTie-era mode reproduced StringTie except for flow-only transcripts.

**Why it can differ from register 1264.** Compatible support admitted a chain only when EVERY junction had >= 2 compatible reads; the
missed junctions sit in chains whose terminal junction one read carries (dev: 864 of 927). StringTie decides on a path's ABUNDANCE
(its reads plus the reads of the partial transfrags it explains), not a per-junction minimum. This rule does the same without flow.

## The rule (`copy_assign --direct-transfrags`, default off = byte-identical; exclusive with `--compat-support`)

Pass-1 groups (one per exact `(chrom, intron chain)`, with exact reads n, outermost start and end) are the transfrags. Transfrag T is
CONTAINED in C (C != T) when T's chain is a contiguous sub-chain `C[i .. i+m)` and T's reads stay inside C's exons around it (smallest
start >= acceptor of C's intron `i-1` when i > 0; largest end <= donor of C's intron `i+m` when i+m < len C) — the containment of
`--compat-support`. A transfrag contained in no other is MAXIMAL. Every non-maximal T transfers its n reads DIRECTLY to the maximal
transfrags that contain it, split in proportion to their own exact reads. Each maximal transfrag's abundance = its reads + what it
received; it is admitted when the abundance reaches the existing read floor (2) and its skeleton carries the rounded abundance.
Non-maximal transfrags are not emitted (their reads are in their containers). Unspliced reads, the gate, polish and every other step
are unchanged. No new parameter.

## Test (as `docs/PREREG_assembler_compatible_support_2026-10-11.md`)

Arms OFF / ON in the `primary` and `good` pools; gffcompare against the figure 1 references. Dev: gorilla OR6737 NC_073244.2, no
tuning, at most ONE amendment before any held-out number. Held-out (decides): gorilla OR6737 minus NC_073244.2; human A119b genome-wide.
**B1** intron sensitivity +1.0 point or more; **B2** intron-chain precision down by <= 0.5 point; **B3** intron-chain sensitivity not lower.
Reported: transcripts, transcript and locus levels, the share of the StringTie intron gap closed. **Verdict:** B1-B3 in both samples and
both pools => ADOPT as the default; primary pool only => default for the benchmark arm; else NOT ADOPTED, opt-in with the numbers.
