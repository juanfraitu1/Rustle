# Win-band paralog candidates — copies StringTie misses, rustle-VG should recover

> ## ⚠ VALIDATION OUTCOME (read first) — the `dec_secondary` signal is largely artifactual
> Running the worked example end-to-end (`LOC115934290`, sister DOK7) showed the scan's
> ranking signal does **not** identify genuine rustle-VG wins:
> - StringTie assembles nothing at `LOC115934290` ✓ … but **rustle-VG also assembles
>   nothing there**.
> - Its 28 "decisive secondary" reads have their genome-wide **primary** at a *third*
>   locus (`NC_073236.2:24,822,889`, AS 2434 ≫ 2028 at LOC115934290), which StringTie
>   **does** assemble. So the reads genuinely belong elsewhere; `LOC115934290` is
>   **spillover**, and both tools are correct to skip it.
>
> **Why the signal is flawed:** a read flagged *secondary* at copy C has, by definition,
> a better placement elsewhere (its primary). If it truly belonged to C, minimap2 would
> flag it *primary* at C and StringTie would see it. So "decisive secondary" is
> self-contradictory — genuine secondary-only support exists **only for TIED reads**
> (equal AS, primary chosen by coin-flip), i.e. the high-identity DAZ/RABL2 regime that
> is *hard to prove*, not a cleaner ≤90–95 % band. The scan also used **windowed,
> annotation-defined families** (often incomplete), so "decisive vs sisters" missed the
> reads' real owner at an un-listed paralog — a third false-positive mode on top of the
> two below.
>
> **Corrected conclusion:** rustle-VG's secondary-alignment advantage is confined to
> **tied/near-tied, co-located** copies (RABL2/DAZ tandem families), where tied reads
> supply coverage *and* a minority of decisive PSV reads anchor each copy. The
> cross-chromosome, lower-identity "win band" hypothesised below does **not** yield
> secondary-rescue wins. A correct scan must use a **genome-wide primary check** (not a
> windowed family) and look for **co-located families where ST emits fewer transcripts
> than copies**. The decisive-evidence + three-bucket guard (now in
> `50_authenticity_guard.py`) correctly classifies `LOC115934290` as phantom.
>
> The candidate table below is retained as a record of the (superseded) approach.

Systematic genome-wide scan (2026-06-08) for the **win band**: paralog copies that
are (a) similar enough that StringTie collapses/drops them, yet (b) divergent enough
that decisive read evidence proves they are independently expressed. Motivated by the
DAZ (99.97 %, too hard) vs amylase (98 %, too easy — ST already gets it) contrast.

## Method (all reproducible from `/tmp` scripts; data in `~/winloci_scratch`)
1. **Reps** — one longest mRNA per protein-coding gene from `GGO_genomic.gff`
   → 22,622 reps (`/tmp/pick_reps.py`).
2. **Paralog families** — `minimap2 -cx asm10 -p0 -N50` all-vs-all of rep transcripts;
   keep cross-gene pairs with **identity 0.90–0.99** and aligned-cov ≥ 0.5; union-find
   → **816 families** (`/tmp/scan_band.py`).
3. **StringTie miss** — a copy is "missed" if no `genome_st.gtf` transcript overlaps it
   → **451 families** have ≥1 covered + ≥1 missed copy.
4. **Decisive support** — for each missed copy, count **spliced** reads whose best
   alignment is strictly better on that copy than on any sister (or unique to it),
   split by whether that alignment is **primary** (StringTie can also see it) or
   **secondary** (StringTie discards it) (`/tmp/final_scan.py`). Families > 25 copies
   (large messy LOC clusters, the precision-crater regime) are excluded from the BAM step.

**The win signal = `dec_secondary`** (decisive *secondary* reads). High `dec_secondary`
with `ST_tx = 0` = a copy StringTie cannot assemble but rustle-VG can, with proof.
High `dec_primary` instead means StringTie can see the reads too → not a VG win
(its apparent "miss" is usually a partial/threshold artifact).

## Result
- 908 missed copies scored → 46 raw clean wins (`ST_tx = 0`, `dec_secondary ≥ 5`).
- **⚠ Two false-positive modes were then found and filtered (see below): 19 of the 46
  are intronless retrocopies. The trustworthy set is the 27 GENUINE multi-intron
  wins.** Full table (with `real_introns` column): `winband_candidates.tsv`.

### Two false-positive modes (validated, now filtered)
1. **Single-exon / intronless retrocopies.** A processed retrocopy (e.g. `LOC101134109`
   of UAP1, `LOC101143340` of POLRMT) is annotated with one or several *abutting* exon
   segments (0 intron between them). Reads from a neighbouring multi-exon host gene
   span *through* it with huge introns and get counted as "spliced decisive". Filter:
   require the missed copy to have ≥1 **real intron** (gap > 50 bp between exons) —
   exon *count* is not enough (retrocopies show 2–3 abutting exons).
2. **Span-through reads.** Even for multi-exon copies, reads whose alignment extends far
   beyond the gene (giant introns) must be excluded. Filter: decisive reads must be
   **contained** within the gene span (± a few kb), not spanning beyond.

Both my original headline picks (`LOC101134109`, `LOC101143340`) were mode-1 retrocopies
— caught by direct CIGAR inspection before any claim was made.

### Headline candidates — GENUINE multi-intron clean 2–4-copy wins
| missed copy | chrom:pos | sister(s) | identity | real introns | dec_secondary | ST |
|---|---|---|---|---|---|---|
| **LOC115934290** | NC_073246.2:11,244,143 | DOK7 | 0.968 | **6** | 25 | 0 |
| LOC115933253 | NC_073248.2:7,003,704 | (n4) | 0.965 | 12 | 20 | 0 |
| LOC129523555 | NC_073242.2:29,864,793 | LOC115934567 | 0.900 | 6 | 19 | 0 |
| LOC129527454 | NC_073242.2 | SLC9B1 (n3) | 0.965 | 8 | 13 | 0 |
| LOC101131829 | NC_073229.2:104,200,054 | ARHGAP21 | 0.980 | 2 | 12 | 0 |
| **WASHC2A** | NC_073232.2 | WASHC2C (n3) | 0.985 | **22** | 11 | 0 |
| LOC109027880 | NC_073231.2 | SBDS | 0.969 | 5 | 9 | 0 |

`LOC115934290` (sister **DOK7**, 0.968 identity, 6 real introns, 8 exons, 25 decisive
secondaries all contained, ST assembles nothing) is the worked example — genuine
multi-intron structure, clean 2-copy family, dead-center win band. `WASHC2A` is a
*named* gene (22 introns) StringTie misses entirely — a striking non-LOC example.

### Larger families (more copies, messier — higher recall, watch precision)
`LOC115931247` (FTH1 family, n=6, 1118 dec_secondary) tops raw support but is a
ferritin retrocopy cluster (mode-1); `GGTLC2`/GGT clusters (n=15–16) recur. These match
the known precision-crater regime — recover with care.

## Caveats / next step
- This is a **screen**: `dec_secondary` predicts recoverability; the proof is to run
  `rustle --vg` vs StringTie on the locus and confirm rustle assembles the copy and
  SQANTI3 calls it FSM-to-self with primary/decisive authenticity. Recommend running
  the `copy_recovery_eval` protocol on the top ~10 (or the chromosomes containing them).
- Requiring **spliced** decisive reads filters out processed (intronless) pseudogenes:
  a spliced read can only align decisively to a copy that has the matching intron
  structure.
- Identity here is mRNA-level (most-conserved); per-read PSV count is what actually
  drives decisiveness — which is why 0.98 copies can still be clean wins.

---

## CORRECTED SCAN (co-located tandem + locality filter) — final result, 2026-06-08

Rebuilt the scan around the *genuine* win pattern (RABL2): **co-located** paralog
families (same chrom, < 10 Mb) where StringTie emits **fewer transcripts than copies**,
with a **primary-locality filter** (reads genuinely primary in the locus, not spillover).

- 400 tandem families → 76 with a provably-distinct multi-intron copy ST collapses.
- **Locality filter is decisive:** ~all 76 are only **3–7 % primary-in-locus** = still
  spillover (the LOC115934290 trap, now exposed systematically). **Exactly one** family
  had genuine local primary support (34 %): the **ZDHHC11/ZDHHC11B palmitoyltransferase
  tandem array** (NC_073243.2, 3 copies ~0.97 id, 244 kb), where ST emits 1 partial
  transcript for 3 copies.
- **Ran rustle-VG vs StringTie on ZDHHC11:** rustle also emits **just 1 transcript at
  the middle copy** (same as ST) — it does **not** split the array. Per-copy support
  explains why: middle copy has 10 spliced primary reads (both tools get it); the two
  flanking copies have only **2 spliced primary reads each** (lowly expressed). Their
  secondary reads (9, 21) are **not** converted by rustle's `--vg` flow.

### Bottom line
The corrected, locality-filtered, retrocopy-free scan finds **no new genome-wide
rustle-VG wins** beyond the known RABL2 portfolio. This independently reproduces the
earlier conclusions: genuine wins are a tiny handful of high-identity **co-located
tandem arrays** where tied reads (coverage) + a minority of decisive PSV reads (anchor)
+ co-bundling coincide; and the multimapper "treasure" is otherwise bounded by an
**annotation/coverage floor** (lowly-expressed flanking copies that neither tool
assembles and rustle's flow does not rescue from secondaries).

### Concrete improvement opportunity (scoped, small, precision-risky)
ZDHHC11 flanking copies show the one realizable lever: copies with a *few* decisive
primary anchors (≥2 spliced) **plus** a pile of locus-local secondary reads (9–21) that
rustle currently drops on the floor. Seeding bundle/flow from secondary piles anchored
by ≥k decisive primaries could recover them — but the evidence is thin and at 0.97
identity the secondaries are near-tied, so expect a small Sn gain at real precision risk.
