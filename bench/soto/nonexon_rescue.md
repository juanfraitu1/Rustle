# Non-exon-signal rescue POC — 24 K=0-bearing Soto families

## Headline

**0 of 34 K=0 members across 20 families are genuinely rescuable via non-exon signal**, after adversarial verification excluded all **26 artifact "rescues"** claimed by the first-pass analysis.

- Families examined: **24**; analyzed with a multi-copy reference: **20** (the other 4 are single-copy families — no sibling to align against, so rescue does not apply).
- K=0 members in scope: **34**.
- Rescues **claimed** by the distinguishability analysis: **26**.
- Rescues that **survived** the three-lens adversarial check: **0**.

The information-theoretic finding still stands: on the 20 multi-copy families, full-genomic (UTR + intron + flanking) sequence makes the large majority of reads distinguishable (per-family `pct_distinguishable` mostly 80–100%). What the verification kills is the *causal attribution* — none of the 26 K=0 members that the analysis flagged as "rescued by non-exon signal" is actually a distinct, previously-undetectable copy lifted above the floor by that signal.

## Mechanism under test: soft-clip flank

The POC re-aligned family reads against **full-genomic copy sequences** (not exon-only references) and counted a read as "distinguishable" when it resolved uniquely to one copy. The proposed driver was **non-exon signal**: UTRs, introns, and — most concretely — **>=20 bp soft-clip flank**, i.e. read bases that hang off the exon body into copy-specific flanking sequence. In most "rescued" families 90%+ of the distinguishable reads carried such a flank (e.g. ID_332 94%, ID_431 353/354), which is what made the mechanism look real on paper.

Adversarial verification applied three lenses to every claimed rescue:

1. **Expressed-copy** — does the "K=0" copy already carry primary reads above the seed floor (3) in the standard genome-wide alignment? If so it was never dark; there is nothing to rescue.
2. **Flank / coordinate-window** — are the counted reads actually *aligning into the copy body*, or merely spanning the copy's coordinate window via long introns / large soft-clips while anchored elsewhere?
3. **Leakage** — is any read's genome-wide **best (primary) hit** actually this copy, or are the counts fed by MAPQ-0 secondary multimappers whose true home is an expressed sibling, force-binned onto a reduced N-copy reference?

Every one of the 26 claimed rescues failed at least one lens. The two recurring failure modes:

- **Expressed-copy relabeling** (e.g. ID_146, ID_240, ID_411, ID_78): the "rescued" member already had tens of unique MAPQ-60 primaries — it was the dominant expressed copy, fully detectable *without* any non-exon signal. Counting its own reads is not a rescue.
- **Multimapper leakage into a coordinate window** (e.g. ID_215, ID_391, ID_65, ID_431): the copy has ~0 native primaries, and the claimed reads are expressed-sibling secondaries (or long spliced pass-through reads) that the reduced reference forces onto this locus. The soft-clip "flank" is then a symptom of a poor secondary alignment overhanging the extracted window, not copy-specific divergence.

Several families combined both (all three lenses fired): ID_14, ID_391, ID_407, ID_65.

## Per-family table

| family | genes (claimed) | copies | reads | distinguishable | %dist | K=0 | claimed | verified | dominant artifact cause |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ID_131 | AC105272.1 | 6 | 11381 | 2369 | 65.7 | 1 | 1 | 0 | leakage + coordinate-window (spliced pass-through) |
| ID_14 | AC005562.2 | 17 | 41178 | 30380 | 96.1 | 1 | 1 | 0 | all 3: expressed upstream gene + coordinate + sibling leakage |
| ID_146 | AC233280.19 | 2 | 6712 | 6068 | 100 | 1 | 1 | 0 | expressed-copy, already detectable (not truly K=0) |
| ID_167 | AC133548.1; CHEK2P2 | 3 | 604 | 268 | 100 | 2 | 2 | 0 | leakage/flank (AC133548.1) + expressed-copy (CHEK2P2) |
| ID_213 | - | 4 | 1056 | 69 | 80.2 | 1 | 0 | 0 | not rescued (below threshold) — no claim to verify |
| ID_215 | AC244669.2 | 7 | 10768 | 3646 | 82.8 | 1 | 1 | 0 | genome-identical to sibling + oversized window + leakage |
| ID_222 | AC243829.6 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | n/a single-copy family (no multi-copy reference) |
| ID_226 | H3P4; H3C13 | 5 | 1046 | 394 | 50.6 | 2 | 2 | 0 | flank/coordinate + upstream readthrough (intronless histones) |
| ID_240 | SYT15; AL157893.1 | 3 | 218 | 123 | 98.4 | 2 | 2 | 0 | expressed-copy + leakage (both members) |
| ID_26 | AC006453.4 | 11 | 4486 | 1954 | 98.7 | 1 | 1 | 0 | expressed-copy + sibling-secondary leakage |
| ID_261 | CNTNAP3P5 | 5 | 2154 | 1176 | 98.4 | 1 | 1 | 0 | expressed-copy + leakage (missing true source locus) |
| ID_302 | BOLA2B | 1 | 0 | 0 | 0 | 1 | 0 | 0 | n/a single-copy family (no multi-copy reference) |
| ID_313 | CDH12P3; CDH12P1 | 2 | 1301 | 358 | 83.8 | 2 | 2 | 0 | expressed-copy + coordinate pass-through + leakage |
| ID_332 | DDX11L16 | 6 | 2043 | 1932 | 97.2 | 1 | 1 | 0 | leakage + flank/coordinate (edge fragments) |
| ID_334 | DEFB104B | 1 | 0 | 0 | 0 | 1 | 0 | 0 | n/a single-copy family (no multi-copy reference) |
| ID_348 | DUX4L50 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | n/a single-copy family (no multi-copy reference) |
| ID_391 | WASH6P | 5 | 2483 | 2358 | 98.3 | 1 | 1 | 0 | all 3: expressed siblings + subtelomeric flank + leakage |
| ID_402 | NCF1C | 3 | 568 | 20 | 83.3 | 2 | 1 | 0 | exonic-not-flank + near-tie leakage from expressed NCF1 |
| ID_407 | TRIM74 | 4 | 3912 | 1974 | 83.4 | 1 | 1 | 0 | all 3: co-located expressed NSUN5P2 read-through, double-counted |
| ID_411 | OR11H13P | 3 | 430 | 117 | 92.1 | 1 | 1 | 0 | expressed-copy (the wrong, dominant member was "rescued") |
| ID_431 | PPIAL4A; PPIAL4C; PPIAL4F | 5 | 942 | 354 | 57.5 | 5 | 3 | 0 | leakage (A, C) + below-threshold/identical-copy (F) |
| ID_65 | CR383656.5 | 9 | 14680 | 12892 | 96.9 | 1 | 1 | 0 | all 3: expressed FAR2P1/chr2 cluster casts secondaries onto chr14 |
| ID_78 | AC243562.1 | 5 | 17651 | 14769 | 99.8 | 1 | 1 | 0 | expressed-copy, 66 MAPQ-60 uniques (already detectable) |
| ID_92 | AC243829.7; AC053481.4 | 2 | 741 | 111 | 100 | 2 | 2 | 0 | leakage/coordinate (AC243829.7) + expressed-copy (AC053481.4) |

Totals: 20 families analyzed, 34 K=0 members in scope, **26 rescues claimed, 0 verified**.

## HONESTY RAILS

**(a) This proves the INFORMATION exists, not that the pipeline uses it.** The distinguishability numbers show that non-exon sequence (UTR/intron/flank) *carries* copy-separating signal in these families. It does **not** show that Rustle's `copy_assign` currently exploits it — wiring non-exon columns into `copy_assign`'s PSV search is a separate production follow-up, and the 0-verified result here says nothing about that follow-up's ceiling. This POC is an information audit, not a method evaluation.

**(b) All-copies test.** "Distinguishable" and "rescued" are defined against the requirement that a read resolve uniquely among **ALL** siblings of the family, not merely against one competitor. A read that separates copy A from copy B but stays tied with copy C is *not* counted as distinguishing A.

**(c) Verified rescues exclude expressed-copy / flank / leakage artifacts.** A claimed rescue survives only if the K=0 copy is genuinely below the seed floor in the genome-wide alignment (not the already-expressed copy), the reads align into the copy body (not merely span its coordinate window via long introns or large soft-clips), and at least the required number of reads have their **best full-genomic (primary) hit** at this copy (not sibling multimappers leaking onto a reduced reference). Under that bar, 0 of 26 survived.

**Note — some "rescues" were the expressed copy and were correctly excluded.** In several families (ID_146, ID_240, ID_411, ID_78, and the CHEK2P2 / AC053481.4 members) the "rescued" member already carried tens of unique MAPQ-60 primary reads in the standard alignment. It was never a dark K=0 copy; the analysis was re-detecting an already-detectable copy and relabeling its own reads as a non-exon rescue. Excluding these is not the verification being conservative — it is the verification refusing to double-count copies the pipeline already resolves. (In ID_411 the *truly* under-detected members, OR11H12 and OR11H1, were not the ones the analysis rescued — the K=0 label there is a family-collapse artifact.)
