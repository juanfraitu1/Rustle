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

## Result (2026-10-11): NOT ADOPTED as the default (B2 fails); kept opt-in

**Dev (NC_073244.2, primary):** default byte-identical; ON 3,862 -> 3,999 transcripts; intron 54.8 / 83.3 -> 58.5 / 81.2 (StringTie 58.9 /
83.9); intron chain 28.3 / 40.8 -> 27.5 / 38.2 (StringTie 24.7 / 37.0). No amendment.
**Held-out** (gffcompare sensitivity / precision; gorilla without NC_073244.2, StringTie and FLAIR scored on the same contigs; human on
`eval_annotated/ref.gtf`):

| arm | transcripts | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla primary OFF | 68,828 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 |
| gorilla primary ON | 72,607 | 63.2 / 80.0 | 27.0 / 33.6 | 24.2 / 33.5 | 30.2 / 70.1 |
| gorilla good OFF | 70,134 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 |
| gorilla good ON | 73,996 | 63.6 / 79.7 | 27.1 / 33.3 | 24.4 / 33.1 | 30.5 / 68.4 |
| gorilla StringTie | — | 64.2 / 83.1 | 24.4 / 34.4 | 21.9 / 34.1 | 30.7 / 64.1 |
| gorilla FLAIR | — | 62.5 / 61.3 | 25.5 / 17.3 | 23.0 / 15.7 | 28.4 / 48.4 |
| human primary OFF | 238,406 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human primary ON | 251,710 | 64.4 / 59.0 | 21.9 / 17.1 | 20.4 / 15.3 | 34.0 / 24.1 |
| human good OFF | 252,217 | 61.2 / 62.7 | 22.0 / 17.7 | 20.5 / 15.4 | 34.8 / 22.0 |
| human good ON | 263,515 | 64.7 / 57.7 | 22.0 / 16.5 | 20.5 / 14.7 | 34.2 / 23.4 |
| human StringTie | — | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.8 | 33.8 / 20.3 |
| human FLAIR | — | 67.7 / 30.2 | 22.4 / 6.3 | 21.1 / 3.6 | 31.9 / 6.3 |

**B1 holds** (intron sensitivity +3.6 / +3.6 / +3.6 / +3.5); **B2 fails** (intron-chain precision -1.7 / -1.7 / -1.5 / -1.2, bar <= 0.5);
B3 holds (+0.2 / +0.1 / 0.0 / 0.0). Not the default. Against the tools on the same contigs (primary pool): human — direct transfrags is
above StringTie on every intron and chain number (intron 64.4 / 59.0 vs 63.5 / 57.8; chain 21.9 / 17.1 vs 18.1 / 15.0); gorilla —
intron sensitivity 1.0 below StringTie and precision 3.1 below, chain sensitivity 2.6 above and precision 0.8 below. FLAIR has higher
intron sensitivity in human (67.7) at half the precision. The flag stays opt-in (`--direct-transfrags`, `RUSTLE_DIRECT_TRANSFRAGS=1`);
whether the benchmark arm uses it is a figure decision, not this test's.
