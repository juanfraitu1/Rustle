# Per-copy recovery of NPIP and TBC1D3: our assembler and the StringTie, FLAIR and isoseq baselines (2026-09-29)

This document consolidates `docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md` for the user and the advisor.
- **Pre-registration.** Frozen on 2026-09-30 (sha1 40b581b8, recorded in `FROZEN.sha1`), before any tool's model or any
  read at a family copy was opened. The Outcome was appended later that night; the text above it is byte-identical to
  the frozen version.
- **Status.** One run, scored by the session that ran it. **Not yet independently recomputed.** The one recheck done
  since is in §4 (the no-polish arm), with its own code.
- **Species** are never pooled. StringTie, FLAIR and isoseq collapse are trusted baselines; the question is where each
  method loses a copy, not who wins.
- **"Ours"** is the shipped default at commit 3007c3d4 (strict junctions, pass-1 floor 2, good-secondary seeding), so
  it predates the 2026-09-29 default flips. `--bridge-regroup f1v2` leaves intron chains unchanged
  (`docs/archive/2026-09/PREREG_f1v2_readshare_2026-09-29.md`), so the chain counts below are not expected to move; they were not rerun.

## 1. The question

For every annotated copy of NPIP and TBC1D3, on the same alignments, which method emits:

| measure | definition |
|---|---|
| LOCUS | at least one model overlapping the copy's exons on the copy strand |
| COMPLETE | a model whose intron chain equals an annotated chain of that copy (gffcompare `=`, cross-checked by an exact-chain test: 0 disagreements) |
| UNIQUE / FUSED | UNIQUE = LOCUS and not FUSED; FUSED = the model also overlaps another copy or an annotated partner gene |

- **Data.** Human A119b (NPIP 26 copies, TBC1D3 16 records) and gorilla OR6737 (NPIP 25, TBC1D3 14). The baselines were
  run on the same two BAMs (StringTie 3.0.1 `-L`, FLAIR 3.0.1, isoseq collapse); their recipes were confirmed one by one.
- **Truth** comes from annotation only: CHM13 RefSeq (Dishuck-checked NPIP copies) and the gorilla RefSeq records.
- **Read classes** (primary reads, `-F 2308`): E2 = copies with at least 2 reads carrying an annotated chain exactly;
  E1 = exactly 1; R2 = at least 2 primary reads at the copy; R0 = none.

## 2. Answer first

1. **Where the reads carry the whole chain (E2), ours is COMPLETE at least as often as every baseline, in all three
   testable cells** (StringTie / FLAIR / isoseq in brackets): human NPIP 11 of 13 (10 / 6 / 9); human TBC1D3 9 of 9
   (9 / 9 / 8); gorilla TBC1D3 3 of 3 (3 / 3 / 2). r1065's genome-wide result holds per copy.
2. **Over all copies StringTie is one ahead at human NPIP (12 of 26, ours 11, isoseq 9, FLAIR 6).** Its one-read floor
   recovers the two copies whose only annotated-chain read is single (NPIPB12, NPIPB13), and the gorilla copy
   LOC129533792. Ours is 0 there by construction (pass-1 floor 2), isoseq 0 of 3. FLAIR reaches the gorilla one by
   folding two 5′-truncated sub-chain reads into one full read, the fold-in r1071 refuted for us. No admission rule for
   these is learnable (r1067-r1073).
3. **Gorilla NPIP is out of reach for every method.** No copy has a read carrying its full chain, and every arm is 0 of
   25 COMPLETE. **Pre-registered falsifier F4 fired here:** LOCUS on the 23 copies with at least 2 primary reads is ours
   11, StringTie 16, FLAIR 18, isoseq 23. After seeing this, a qualifier was added (post hoc): a model with at least one
   annotated junction exists at only 3 of StringTie's 16, 3 of FLAIR's 19 and 10 of isoseq's 24 present copies (ours
   0 of 11), and isoseq's best model is single-exon at 19 of 24. The extra presence is fragments; no method
   reconstructs a gorilla NPIP isoform.
4. **Seeding separates us from the baselines at low-read copies.** Six of our COMPLETE copies have fewer than 2
   exact-chain primary reads and rest on good secondary alignments: gorilla LOC129533806 and LOC129533813 (0 primary
   reads), LOC129533808 (1), LOC129533797, and human TBC1D3B and TBC1D3I. No baseline has a locus at any copy with 0
   primary reads, and our primaries-only arm loses exactly these (COMPLETE 7 → 3 gorilla TBC1D3, 11 → 9 human TBC1D3).
   **Caveat, pre-registered:** these rest on reads whose best placement is another copy. They show the good-secondary
   rule places reads there, not that the copy is expressed; deciding that is the copy-assignment question (O2).
5. **Our two human NPIP E2 misses are polish steps, not assembly.** NPIPB2: the isoform-fraction rule (0.02) inside a
   locus fused with GSPT1 (948 reads, the copy's own 23). NPIPB6: the 2-read chain is an ISM sub-chain of a 67-read
   model. Closed rows r1063 and r1064; not reopened.

## 3. Per-cell counts

Human A119b, **NPIP** (26 copies; 13 with at least 2 exact-chain reads, 2 with exactly 1):

| arm | LOCUS | COMPLETE /26 | on E2 /13 | on E1 /2 | copies with a FUSED model | models per copy, median |
|---|---|---|---|---|---|---|
| ours (default) | 26 | **11** | **11** | 0 | 13 | 16 |
| ours, primaries only | 26 | 10 | 10 | 0 | 12 | 7 |
| StringTie | 26 | 12 | 10 | **2** | 14 | 13 |
| FLAIR | 26 | 6 | 6 | 0 | 17 | 20 |
| isoseq collapse | 26 | 9 | 9 | 0 | 20 | 74 |

Human A119b, **TBC1D3** (16 records; 11 scorable for COMPLETE, 9 with at least 2 exact-chain reads):

| arm | LOCUS /16 | COMPLETE /11 | on E2 /9 | copies with a FUSED model | models per copy, median |
|---|---|---|---|---|---|
| ours (default) | 13 | **11** | 9 | 5 | 12 |
| ours, primaries only | 13 | 9 | 9 | 6 | 6 |
| StringTie | 14 | 9 | 9 | 5 | 12 |
| FLAIR | 13 | 9 | 9 | 9 | 14 |
| isoseq collapse | 16 | 8 | 8 | 8 | 44 |

Gorilla OR6737, **NPIP** (25 copies; 23 with at least 2 primary reads; **no copy has an exact-chain read**):

| arm | LOCUS on R2 /23 | COMPLETE | present copies with a model carrying an annotated junction | best model multi-exon | copies with a FUSED model |
|---|---|---|---|---|---|
| ours (default) | 11 | 0 | 0 / 11 | 10 / 11 | 1 |
| ours, primaries only | 9 | 0 | 0 / 9 | 8 / 9 | 1 |
| StringTie | 16 | 0 | 3 / 16 | 10 / 16 | 4 |
| FLAIR | 18 | 0 | 3 / 19 | 11 / 19 | 3 |
| isoseq collapse | 23 | 0 | 10 / 24 | 5 / 24 | 11 |

(The junction and multi-exon columns are computed over each arm's present copies, so FLAIR and isoseq show all 19 and 24.)

Gorilla OR6737, **TBC1D3** (14 records; 7 with at least 2 primary reads, 5 with none; 3 with at least 2 exact-chain reads):

| arm | LOCUS /14 | at the 5 zero-read copies | COMPLETE /14 | on E2 /3 | on E1 /1 | copies with a FUSED model |
|---|---|---|---|---|---|---|
| ours (default) | 7 | **2** | **7** | 3 | 0 | 1 |
| ours, primaries only | 4 | 0 | 3 | 3 | 0 | 1 |
| StringTie | 7 | 0 | 4 | 3 | 1 | 1 |
| FLAIR | 5 | 0 | 4 | 3 | 1 | 1 |
| isoseq collapse | 9 | 0 | 2 | 2 | 0 | 2 |

## 4. Where each method loses a copy

- **Ours, gorilla NPIP (14 absent copies):** pass-1 floor 3, gate 8 (strict junctions and the single-exon `+`
  placeholder), polish mono floor 3. StringTie, FLAIR and isoseq have a locus at 6, 10 and 13 of these, but with an
  annotated junction only at NPIPB14P (all three), NPIPA5 (FLAIR, isoseq), NPIPB5 and NPIPB15 (isoseq).
- **Ours, gorilla TBC1D3 (7 absent):** pass-1 floor 4 (1-6 primary reads, each a distinct chain; the baselines are
  present there with one-read or folded models) and seeding 3 (0 primary reads; every baseline is absent too).
- **Ours, human (3 absent, all TBC1D3 pseudogene spans):** a single-exon cluster dropped as a mono shadow (TBC1D3P4),
  the gate (strict junctions plus the single-exon `+` placeholder on a `-` span; TBC1D3P3), the mono floor (TBC1D3P7).
- **Recheck of the two polish misses (own code, `runs/hsa.nopol.gtf` and `hsa.ours.gtf` against the truth GTF):** the
  no-polish arm holds exact-chain models at NPIPB2 (2) and NPIPB6 (1); the shipped default holds 0 at both.
- **Baselines.** FLAIR misses 12 copies that have reads (8 have no chain shared by 3 reads, FLAIR's support floor);
  StringTie misses 10 (9 of them at 9 primary reads or fewer); isoseq misses no copy that has a primary read, so its
  loss is completeness, not presence.

## 5. Fusion and over-splitting

- **All five arms are FUSED at the four RefSeq readthrough copies** (NPIPA1, NPIPA6, NPIPA9, NPIPB14P): fusion is in
  the reads. Human NPIP copies with a FUSED model: ours 13 (11 without the two readthrough-defined copies), StringTie
  14, FLAIR 17, isoseq 20.
- **Gorilla NPIP:** isoseq is FUSED at 11 copies (one-read fragments joined to neighbouring loci), StringTie 4, FLAIR
  3, ours 1.
- **Over-splitting** (models per present human NPIP copy, median): isoseq 74 (max 321), FLAIR 20, ours 16, StringTie
  13, ours-primaries 7. Human TBC1D3: 44 / 14 / 12 / 12 / 6 in the same order.

## 6. What it means, and what is open

- **For the assembler.** On these two families it is as complete as the baselines wherever the reads allow it, has
  models at copies the primary-only baselines cannot see, and its remaining gap is the read-completeness wall that
  stops every method (gorilla NPIP; `docs/archive/2026-09/PREREG_ggo_npip_variant_sim_2026-09-29.md`).
- **Prediction record.** Of 9 predictions, 6 held in full (P1, P2, P3, P5, P7, P8), 2 partly (P4: StringTie, not
  isoseq, recovers single-read chains; P6: our fused count 13 against a bar of 8-12, and gorilla isoseq 11), and 1
  failed or is vacuous (P9). One falsifier fired (F4, §2.3).
- **Open, untested.**
  - The NPIPB2 miss is the isoform-fraction denominator being set by the dominant gene of a fused locus. Computing the
    fraction inside the bridge-free group would need its own pre-registration and a check against r1063.
  - Rerun on the new defaults, and an independent recompute of the counts in §3.
  - Hostile readings kept from the pre-registration: LOCUS at 1 bp is lenient (§2.3 is why); n is at most 26 per
    cell, so every claim is a count; the human default GTF was checked against the frozen binary regionally on
    chr16 and chr17 (byte-identical), not genome-wide.

## 7. Statement for the advisor

> On the same reads and alignments, for each annotated copy of NPIP and TBC1D3 with at least two reads carrying its full
> intron chain, our assembler reconstructs that chain at least as often as StringTie, FLAIR and isoseq collapse (human
> NPIP 11 of 13 copies against 10, 6 and 9; human TBC1D3 9 of 9 against 9, 9 and 8; gorilla TBC1D3 3 of 3 against 3, 3
> and 2). Where the only supporting read is a single read, StringTie's one-read floor recovers copies that ours does
> not. Where no read carries the full chain (gorilla NPIP, 0 of 25 copies), no method reconstructs it, and the extra loci
> the other methods report there are mostly single-exon or unannotated fragments. Our assembler also has models at
> copies whose reads all place better on another copy, which primary-only methods cannot.

## 8. Sources

- **Pre-registration and Outcome (O1-O8):** `docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md`. Register rows 1189-1193.
- **Instruments and products (scratch, no backup):** `/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools/`.
  `code/` holds `build_truth.py`, `reads.py`, `models.py`, `gc.sh` (gffcompare), `report.py`, `trace_port.py`,
  `arm_hsa.sh`, `polish_port.py` and `posthoc.py`; the frozen sha1s are in prereg §8 and `FROZEN.sha1`. `score/report.md`
  has the per-copy tables and `score/posthoc.json` the post-hoc qualifiers. The instruments are not archived in the repo.
- **Prior per-copy work reused:** `docs/archive/2026-09/PREREG_ggo_npip_variant_sim_2026-09-29.md` and
  `docs/archive/2026-09/PREREG_ggo_tbc1d3_holdout_2026-09-29.md` (first losing step for ours at gorilla copies).
