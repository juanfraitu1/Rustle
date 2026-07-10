# YAG copy recovery vs IsoCon, and a full-panel regression after the DAZ2/harmony work

**Date:** 2026-07-10. Current binary (commit `7ffec03`), `GGO_mm.bam`, `--min-copies 2 --skip-poa-diagnostic
--homology-primary`, each region run alone, foreground.

## The comparison is not apples-to-apples — IsoCon used targeted amplicons

IsoCon (Sahlin 2018) analyzed the nine human YAG families from **targeted RT-PCR amplicons** — one primer pair
per family, PacBio Iso-Seq — in **two testis samples**, and reports transcript "groups" (its copy-number proxy)
that a transcript may belong to more than one of, so a group count is an *upper bound* that folds in isoforms.

Its per-family read depth is 10–30× ours because the amplicon targets one family:

| family | IsoCon reads (sample 1) | IsoCon "groups" (copy estimate) | our primary reads in window | our copies |
|---|---|---|---|---|
| DAZ | 495 | **3** | 220 | **2** |
| RBMY | 6615 | **18** (14 coding) | 87 distal + 77 proximal | **6 distal + 3 proximal = 9** |
| TSPY | 2121 | **20** (14 coding) | 34 | **5** |
| HSFY | 933 | **0 shared between samples** | 10 | **0** |
| BPY | 36 | **1** | see note | — |
| CDY / PRY / VCY / XKRY | 36–1766 | 1–2 each | ~0 expressed here | — |

So **no, we do not report as many "copies" as IsoCon** — but the numbers are not measuring the same thing on the
same input:

- **Depth.** IsoCon had 495 DAZ reads from an amplicon; we had 220 from the whole transcriptome. TSPY: 2121 vs
  34. We are coverage-limited exactly where IsoCon is coverage-rich.
- **Groups ≠ copies.** IsoCon's RBMY "18 groups / 61 transcripts" counts isoforms; the annotation has 13 RBMY
  copies in this gorilla assembly, and we recover **9 of them** across the array (6 in the distal cluster, exact;
  3 of 6 in the proximal), from 1/40th the reads.
- **Concordance where it matters.** HSFY: IsoCon found **zero transcripts shared between its two samples**; we
  find **zero families**. Both methods agree there is nothing callable there. That is agreement, not a miss.
- **DAZ is the honest headline.** IsoCon 3 groups from 495 targeted reads; we recover **2 copies from 220
  whole-transcriptome reads**, where before this session we recovered **0**. The gap to IsoCon's 3rd group is
  depth, not method — our DAZ2 model is even 5′-truncated (70% of span) for lack of reads at the 5′ end.

**One concrete gain from the pooling fix:** RBMY proximal went from **2 copies to 3** — junction-incidence
pooling assembled a copy whose reads had fragmented below the per-isoform gate.

The like-for-like statement: to match IsoCon we would need IsoCon's input — targeted amplification of one family
at 10–30× depth. On shared whole-transcriptome Iso-Seq the method recovers the copies the depth supports, and the
DAZ2 fix is what turned DAZ from 0 to 2.

## Does the full panel still hold? Mostly yes; two flagged, one still off-target

| region | annotated | copies | result | status |
|---|---|---|---|---|
| **GSTM** | 3 | 3 | 1 fam, 2673 asg | ✅ clean (GSTM3+GSTM5+GSTM1) |
| **MAGEA** | 2 | 2 | 1 fam, 931 | ✅ clean (MAGEA4+MAGEA10) |
| **DAZ** | 2 | 2 | 1 fam, 2213 asg | ✅ recovered this session |
| **RBMY distal** | 6 | 6 | 1 fam, 735 asg | ✅ exact |
| **r1** | — | 2 | 1 fam, 665 | ✅ unchanged |
| **TSPYL1 / EEF1A1 / DERPC** | 1 | 0 | 0 fam | ✅ controls hold |
| **CEACAM / PRAMEF / HSFY** | 6/3/5 | 0 | 0 fam | ✅ correct silence — not expressed here (16/3/10 reads) |
| **r4** | — | 6 | 2 fam, 818 | ⚠ one `Containment` pair flagged |
| **RFPL** | 2 | 4 | 2 fam, 814 | ⚠ off-target + 2 `Containment` pairs flagged |

The controls hold and every clean family is intact. Two caveats, both real:

**RFPL regressed in cleanliness.** It is a 104-read region. This session's primary-only junction support (the
correct fix — the old count included secondaries) makes the readthrough filter *less* aggressive, so a borderline
transcript that used to be dropped now survives. The result is two families that are **partly off-target**: RFPL2
(`30208643`) forms no family, and the emitted copies sit in the SLC5A4→RFPL3 interval with two `Containment`
overlaps. The `catalog_overlaps` warning fires on both pairs — the safety net works — but the families themselves
are not clean. This is the known, unfixed **Containment** defect: a spliced long transcript enclosing a shorter
one is admitted as two copies. It is the exact analogue of the chimeric-bridge case, and wants the same
treatment (drop the contained fragment), which was out of scope this session.

**BPY2 cannot be windowed.** Its two copies sit at `42645394` and `43077068` with the DAZ locus *between* them, so
any window over both BPY2 copies contains DAZ — the "BPY2" family that forms is DAZ. BPY2 itself has ~0 primary
reads (all-secondary). Unrecoverable from this RNA, as recorded.

## Verdict

- **vs IsoCon:** not the same count, and it cannot be — IsoCon targeted each family at 10–30× our depth and
  counts isoform-inflated groups. On shared whole-transcriptome input the method recovers what the depth allows,
  and this session took DAZ from 0 to 2 and RBMY-proximal from 2 to 3.
- **Panel:** clean families and controls all hold; not-expressed families correctly stay silent; RFPL and r4
  produce `Containment`-flagged families — the warning catches them, but RFPL is genuinely off-target and the
  Containment defect should be fixed next.

Related: `bench/YAG_CHECK.md`, `bench/DAZ2_RECOVERED.md`, `reference_isocon_sahlin`.
