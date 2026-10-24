# Pre-registration: the Pareto set against a null arm, on the near-tie band alone, fresh one-genome simulations (2026-10-20)

**User (2026-10-20):** after register 1348's audit (a PSV-free null passes the same bars): "ok, lets try both, in order" —
first this test, then the FCGR2A lift check.

## Why

Register 1348 showed that the bars of `PREREG_o2_window_set` certify the 30-point window, not the set: a rule using no read
bases (minimap2's primary when no rival is within 10 points, every copy otherwise) passes them too. The set can show value
only where that null gives nothing: the **near-tie band**, reads whose closest rival is within 10 intron-free points. There
the null keeps every copy, so the question is whether the PSV columns narrow the pool, and at what cost in dropped truths.

## The rule (unchanged from `PREREG_o2_window_set`, code b4594ea2; `bench/o3_graph/o2window.py`)

- Window: alignments within 30 points of the best `ms` on two or more copies of one family, none off every copy, the
  primary among them; stitched alignments excluded.
- **Band under test:** reads first tied within 10 points (`window_copies` at W = 10, i.e. `ms > best − 10`).
- The set: `pareto_set` over the substitution columns of the common window.
- Single pick (fewest edits): descriptive.

## Arms, all scored on the band's reads

1. **The set.**
2. **Null, keep-all:** every copy in the 30-point window. It holds the truth whenever the truth is in the window and
   never narrows.
3. **Null, primary:** minimap2's primary alone (a one-copy set on every read).

## Substrate: round 26, fresh one-genome simulations

- `families26.json` (set `r26`, sha1 1826986addb5, 90 copies): 30 families drawn with `random.seed(2026)` from the 130 qualifying non-PAR families of the
  pool (`r19seq/families19cal_full2_all.json`) that no earlier round (families2-25) used. 2-9 copies each.
- Worlds INTACT, CROSS and SILENT under `profiles/{A119b,public}_testis_v4_ref.json` (no donor alleles), seeds 1601-1603:
  6 cells, 18 worlds, built by `fastset.py run` with round 25's `$V17` (`fast/heldout25/run_sim.sh`).
- Truth: each read's simulated copy. Every copy is in the reference.
- Scored read: determined, truth among the window's copies, and in the band. Pooled over the 18 worlds, reads keyed by
  world.

## Bars (the band only)

- **P:** at least 300 scored reads from at least 11 families.
- **H:** the set holds the truth on at least 99% of scored reads.
- **N (beats keep-all):** the set is one copy for at least 25% of scored reads. (Gorilla development: 42%; S1's band,
  seen after the fact: 38.6%. The keep-all null scores 0.)
- **Q (beats the primary):** of the set's singletons, at least 98% are right, and more are right than minimap2's primary
  on those same reads.

**Verdict:** PASS = P, H, N and Q hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails.

**A PASS would claim:** in one genome, among reads whose rival copy is within 10 intron-free points, the PSV columns
narrow at least a quarter of them to a single copy that is right at least 98% of the time and more often than minimap2's
coin toss, while keeping the true copy on at least 99% of all such reads. It would not claim anything about donors whose
copies differ from the reference (register 1348's S2), copies absent from the reference, or exact ties.

## Reported alongside, no bars

- the three arms on the whole 30-point window and on the 10-20 and 20-30 bands;
- the single pick's fixes and breaks;
- mean pool and set sizes;
- per family: scored, held, singletons, singletons right.

## Rules

- Round 26 is built and scored once with the code at this file's commit; bugs go into an amendment with the rerun.
- Heavy jobs foreground, one at a time, under 10 minutes each (one `fastset.py run` per command).
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Human numbers only; nothing pooled with gorilla.

## Amendments

### Amendment 1 (2026-10-20): substrate G, gorilla KB3781's held-out families — written before any gorilla band was scored

User: "I thought that for this particular experiment we would use gorilla and not use the paternal and maternal
haplotypes till O3, can we check if it holds there too?"
- **Reads and reference:** his fibroblast library on GCF_029281585.2 (`personal/sub.bam`), the 143 held-out families of
  `PREREG_o2_personal_copies` (`qual.json`, half `heldout`). No arm uses his haplotypes.
- **Answer key only:** the intron-free key (`o2window.py truth-ms`: the copy with the best `ms` in his two assemblies, by
  at least 5 points, through `orthologs.tsv`), from alignments of the band's reads to `mGorGor1.{mat,pat}.splice.mmi`.
- **Same band, arms and bars** (P, H, N, Q). Verdict reported separately from round 26; gorilla and human never pooled.
- **Caveats fixed in advance:**
  - These families are not wholly fresh: their 10-AS core (896 reads) was scored once for the single pick in register
    1347, and that run's SM392 chimeras motivated the stitched rule. No set or window result was ever seen on them.
  - A read whose copy equals the reference cannot be decided by the key when the reference itself has the near-tie
    (both copies score alike on his identical haplotype); such reads will be undetermined and unscored. The scored band
    will therefore be weighted to reads carrying the other parent's copy, as register 1348's S2 was. Reported by origin.

## Results (code 7a7be665; round 26 built and scored once; `personal/window/r26/band.json`)

**The near-tie band (first tied within 10 points): PASS** — 8,955 scored reads from 22 families.

| arm | holds the truth | one copy | singletons right | minimap2's primary right on those reads |
|---|---|---|---|---|
| **the set** | 8,953 (99.98%) | 5,045 (56.3%) | 5,043 (99.96%) | 5,041 |
| null, keep-all | 8,955 (100%) | 0 | – | – |
| null, primary | 6,972 (77.9%) | 8,955 | 6,972 (77.9%) | – |

P, H, N and Q hold. **Q holds by two reads** (5,043 against 5,041): the set almost never corrects minimap2.

**What the set does, read by read:**
- On the 5,045 reads it narrows to one copy, minimap2's primary was already right on 5,041 (99.9%).
- On the 3,910 reads it keeps two or more copies, minimap2's primary is right on 1,931 (49.4%): a coin toss. Those are
  reads identical to their copies over the aligned bases (pool 2 → set 2: 1,341 reads; 3 → 3: 401; 4 → 4: 455).
- So in one genome the set **separates minimap2's safe picks from its coin tosses** rather than fixing the coin tosses:
  it says which tied reads are decidable from their bases (and decides them, 99.96% right) and which are not.
- Single pick: carries no bar; dropped truths 2 (both S19PRAMEF7).

**Other bands (no bars):** 10-20: set holds 6,098 / 6,098, one copy on 6,077, all right; minimap2's primary right on all
6,098. 20-30: 5,024 / 5,024, singletons 5,014 all right; primary right on all 5,024. Whole window: 20,077 scored, held
20,075, singletons 16,136 (right 16,134), primary right 18,094.

## Results, substrate G (gorilla KB3781, 143 held-out families; code 92f54122; `personal/window/g/band.json`)

Answer key: the intron-free key from his two assemblies (statuses over the 2,755 window-30 reads: refcons 1,514, altonly
1,212, undetermined 26, off-copy 3). O2 never saw the assemblies.

**The near-tie band (first tied within 10 points): UNDERPOWERED** — 175 scored reads from 12 families (P needs 300 / 11).
Descriptive: set holds 173, one copy for 84 (right 82), minimap2's primary right on those 77 (5 fixes, all SM392, 0 breaks);
on the 91 reads kept ambiguous minimap2 is right on 80.

**Other bands (no bars; the whole window passes P, H, N, Q descriptively):**

| band | scored | set holds | one copy (right) | minimap2 right on those | set fixes / breaks | kept > 1: minimap2 right |
|---|---|---|---|---|---|---|
| 10-20 | 1,336 | 1,335 | 802 (801) | 598 | 203 / 0 | 527 of 534 |
| 20-30 | 1,200 | 1,199 | 787 (786) | 770 | 16 / 0 | 412 of 413 |
| window 30, all | 2,711 | 2,707 | 1,673 (1,669) | 1,445 | 224 / 0 | 1,019 of 1,038 |

- **Every one of the 224 fixes is SM392**, the family whose minimap2 primaries are two-copy chimeras (register 1347's
  audit): the stitched primary is credited to SM392.0; O2 takes no stitched alignment as a placement and the PSV columns of
  the plain alignments decide SM392.1, right 224 / 224 by his paternal assembly.
- Outside SM392 the set decides only reads minimap2 already had right (1,445 of 1,449), and keeps ambiguous reads on which
  minimap2 is right 98% (unlike the one-genome simulations, where those were a coin toss): on gorilla the key cannot break
  a near-tie the reference itself has, so the band's refcons reads are mostly undetermined and unscored.
- Human and gorilla numbers are not pooled. Audit pending.

### Audit of substrate G (metric-auditor, 2026-10-20): holds with caveats — the mechanism stated above is wrong

- Every number reproduces. UNDERPOWERED is the right verdict; the band also fails H (98.86%) and Q (97.6%), so with power it
  would have been FAIL.
- **No PSV column decided any of the 224 SM392 reads.** Each has exactly one non-stitched alignment in the window
  (`candidates` = SM392.1 alone), so the set is a singleton by elimination of the chimeric primary; `pareto_set` never ran.
  SM392.0 entered the window only through the chimera that was then excluded. This is an alignment-defect repair (the plain
  SM392.1 secondary already has the higher AS on 217 / 225), not a copy decision, and the stitched rule was designed on
  this family in this half: a named diagnostic, not a held-out finding.
- **Without SM392: fixes 0 / breaks 0** (1,445 = 1,445), the same picture as register 1348's null.
- **The narrowing is self-scored.** Singletons are 1,442 / 1,503 refcons reads (keyed by the reference's own alignment)
  and 7 / 983 altonly reads outside SM392; on those altonly reads minimap2 is right 98.0%.
- The SM392 truth itself is sound (paternal SM392.1 exon 3, margin >= 19; both copies present on the paternal haplotype in
  order).
- Wording not to use: "the PSV columns decide SM392.1"; "O2 fixes 224 reads"; "held-out replication" for SM392; "the set
  separates safe picks from coin tosses" on gorilla.
