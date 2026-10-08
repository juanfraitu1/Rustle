# Pre-registration — `copy_assign` collapses distinct molecules with identical coordinates; fix it and price it

**Written 2026-09-22 (§6z6), before the fixed binary is run on any scored substrate.** Follows
`docs/PREREG_external_tool_bakeoff_2026-09-22.md` §"The second question".

## The defect, established

`src/bin/copy_assign.rs:2962` de-duplicates the region's PRIMARY reads by
`(chrom, ref_start, ref_end, intron_chain)` — no read name — before pass-1. The comment says it exists
because "PrimaryRead has no name; its (chrom, span, intron chain) identifies the placement", i.e. it was
written to drop a record fetched twice from two OVERLAPPING copy windows. But two different molecules with
the same alignment coordinates collide on that key too, and in a deep FLNC library that is common:

| substrate | primary records | dropped by the key | spliced chains with ≥2 reads | pushed below the pass-1 floor of 2 |
|---|---|---|---|---|
| human A119b chr21 | 301,903 | **24.6%** | 13,201 | **3,107 (23.5%)** |
| human A119b chr20 | 414,711 | **31.8%** | 16,166 | **3,732 (23.1%)** |
| gorilla NC_073244.2 | 158,328 | **38.7%** | 7,096 | 155 (2.2%) |

Traced case (ZBTB21, `XM_054324526.1`): two MAPQ-60 all-GT-AG reads carry the exact reference chain, both
at 40380999–40400569; the key keeps one, pass-1 sees n = 1, the chain never exists. Of the **75** human
chr20-22 reference transcripts that FLAIR/isoseq match, our raw output lacks, and ≥ 2 exact-chain reads
support, **75 (100%) have their exact reads collapse to < 2 distinct spans under this key.**

The other callers of `reads_in_region` (`gw_family_catalog`, `mcl_families`, `collapse_enumerate`) do NOT
apply this key, so the family-catalog node path is unaffected; the defect is confined to `copy_assign`
(`--assemble-only` and O2). StringTie, FLAIR and isoseq collapse all count such reads.

## The fix

A record is a double-fetch iff it overlaps a window fetched EARLIER for the same region (both readers return
every record overlapping `[lo, hi)`: the scan path tests `start >= hi || end <= lo`, the indexed path is
noodles' overlap query). So: keep a primary read from window *i* unless it overlaps a window *j < i* on the
same chromosome. With one window — the entire `--assemble-only` path — nothing is dropped. The
`bam_reads` side already de-duplicates by `(name, chrom, ref_start)` and is untouched.

⚠ NOT byte-identical for O2 either: multi-window regions now keep identical-coordinate molecules. That is
the correct behaviour, and its O2 effect is reported, not hidden.

## Substrates and arms

- **Development: human A119b chr20/21/22**, RefSeq CHM13 (10,851 reference transcripts on the three
  contigs). Arms: raw (`--assembly-polish none`) and shipped polish, old binary vs fixed binary. The lab's
  StringTie 3.0.1 / FLAIR 3.0.1 / isoseq collapse on the same contigs are the comparators.
- **Held out: gorilla `GGO_mm.bam`, all 26 contigs**, RefSeq GCF_029281585.2 (106,508 transcripts), the
  shipped polish only, fixed binary vs the old binary's run already on disk, scored against the same three
  lab tools with the 15-cell rule of the parent prereg. No tuning on gorilla.

## Metrics and bar — committed now

gffcompare v0.12.10, query = ours, reference = annotation: matching intron chains, intron-chain SN/PR,
transcript SN/PR; plus the 15 tool cells on gorilla.

| held-out gorilla outcome (fixed vs old, shipped polish) | verdict |
|---|---|
| matching chains **up**, intron-chain precision down by **< 1.0 pt**, 15-cell count **not lower** | ⭐ **ADOPT** |
| chains up but precision down ≥ 1.0 pt, or the cell count falls | ⚠ **TRADE** — report both, do not flip the default without the user |
| chains not up | ⛔ **NO** — the dedupe was not what lost them |

**Predicted, before running:** ⭐ on human (the 75 traced cases become reachable, and 23% of 2-read chains
return), and ⭐ or ⚠ on gorilla, where only 2.2% of 2-read chains were below the floor, so the gain will be
small — a few hundred chains at most — and precision should barely move because the polish is self-tuned
on the run's own support distribution. If gorilla shows NO chain gain, the human effect is library-specific
(A119b's duplicate rate), and the fix is still correct but not a bakeoff lever.

I will not change the arms, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO by the bar: the raw assembly gains chains, the shipped polish then loses more.**

Fixed binary = the window rule (`--keep-coordinate-duplicates`); old = the coordinate key. gffcompare
v0.12.10, query = ours, reference = RefSeq. The ZBTB21 chain is emitted by the fixed binary (31 vs 27
transcripts at the locus), so the mechanism is confirmed.

## Development — human A119b chr20/21/22 (10,851 reference transcripts)

| arm | transcripts | intron SN/PR | **intron chain SN/PR** | transcript SN/PR | **matching chains** | novel loci |
|---|---|---|---|---|---|---|
| raw, old | 42,121 | 62.5 / 51.0 | **26.4 / 8.7** | 24.7 / 6.4 | **2,597** | 81.9% |
| raw, fixed | 50,238 | 66.0 / 40.7 | **27.2 / 7.1** | 25.5 / 5.5 | **2,678 (+81)** | 82.7% |
| polished, old | 15,876 | 61.0 / 56.3 | **23.3 / 16.0** | 21.8 / 14.9 | **2,295** | 55.2% |
| polished, fixed | 17,296 | 62.0 / 49.1 | **21.8 / 13.7** | 20.5 / 12.8 | **2,148 (−147)** | 60.8% |
| StringTie 3.0.1 (lab) | 14,169 | 62.7 / 55.9 | 18.9 / 14.8 | 17.4 / 13.3 | 1,863 | 54.5% |
| FLAIR 3.0.1 (lab) | 54,765 | 66.5 / 28.2 | 23.1 / 6.3 | 21.3 / 4.2 | 2,268 | 79.7% |
| isoseq collapse (lab) | 169,422 | 76.9 / 15.0 | 27.7 / 2.4 | 26.3 / 1.7 | 2,723 | 88.1% |

15-cell score on these contigs: polished old **12/15**, polished fixed **7/15**.

## Held out — gorilla `GGO_mm.bam`, 26 contigs, shipped polish (106,508 reference transcripts)

| arm | transcripts | intron SN/PR | **intron chain SN/PR** | transcript SN/PR | **matching chains** | novel loci |
|---|---|---|---|---|---|---|
| polished, old | 78,602 | 59.6 / 81.7 | **27.0 / 33.1** | 24.3 / 32.9 | **25,829** | 9.4% |
| polished, fixed | 68,555 | 59.6 / **84.3** | 25.9 / **36.3** | 23.3 / **36.1** | **24,733 (−1,096)** | 9.9% |
| raw, old | 192,388 | 60.0 / 76.0 | 28.6 / 17.2 | 25.8 / 14.3 | 27,357 | 57.1% |

15-cell score: old **10/15**, fixed **12/15** (fixed now beats StringTie 5/5: 36.3 vs 34.6 chain precision).

## Verdict and why

⛔ **Matching chains are NOT up on either substrate's polished arm** (the only arm the bar names), so the
fix is not adopted. The prediction was wrong about the polish, not about the assembler: the raw arm does
recover chains (+81 human, i.e. the 75 traced cases), but every polish dial — the isoform fraction, the
mono quantile, the ISM ratio — was fitted (§6p8-§6q4) on read counts that had ALREADY been de-duplicated by
this key, and fed the true counts they cut deeper: containers gain more duplicate reads than their
fragments, so the ISM ratio drops more real short isoforms, and the locus maximum rises, so the 2% fraction
drops more minor isoforms. **The coordinate key is a de facto duplicate filter the shipped polish depends
on.** On gorilla the trade is −4.2% chains for +3.2 pts precision and two more cells; that is a different
operating point, not a better assembler, and per the bar it stays the user's call.

**Shipped:** `copy_assign --keep-coordinate-duplicates` (default off, byte-identical unset — verified on
the ZBTB21 region against the old binary and on full chr21 against this morning's GTF). Re-fitting the
polish on true counts is a separate, pre-registrable experiment; not done here.
