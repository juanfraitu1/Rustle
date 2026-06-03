# Are the tied DAZ reads breakable? — a tag-by-tag analysis

**Question (advisor, 2026-06-02).** For the DAZ multimapping reads that look "tied"
between DAZ1 and DAZ3, are they tied *literally everywhere*, or can we build a
discriminator from the other BAM tags to break the ties (or at least shrink the
ambiguous set)?

**Short answer.** The fully-tied reads are byte-identical on **every** alignment
tag — there is no per-read signal to exploit, and that is *provable* rather than a
failure to look hard enough. But "tied vs resolved" is not binary: most reads near
the diagonal are only *weakly* ambiguous and can be sharpened by using the tags
jointly. Only a hard floor is truly unbreakable, and the anchored prior contains it.

---

## Method

- Input: `/tmp/daz_aln.sam` = `samtools view GGO.bam NC_073248.2:42,783,133-42,945,552`
  (all reads over DAZ1 + DAZ3, primary and secondary alignments).
- For each read, keep its best alignment to **each** copy (DAZ1 and DAZ3). 215 reads
  place on both copies.
- "Fully tied" = identical `AS`, `NM`, and `de` between the two placements (39 reads).
- For those 39, compare **all** per-alignment tags between DAZ1 and DAZ3:
  `AS, ms, NM, de, s1, cm, nn, rl`, and the splice-junction count from the CIGAR.
- Script: `bench/paralog_secondary_scan/` → data `/tmp/daz_tied_data.py`
  (writes `/tmp/daz_tied.json`), figure `daz_tied_breakability.py`.

## Result — the 39 fully-tied reads are identical on all 9 tags

| tag | meaning | # of 39 where DAZ1 ≠ DAZ3 |
|-----|---------|---------------------------|
| `AS` | alignment score | **0** |
| `ms` | max-scoring-segment score | **0** |
| `NM` | edit distance | **0** |
| `de` | gap-compressed divergence | **0** |
| `s1` | chaining score | **0** |
| `cm` | chained minimizers | **0** |
| `nn` | ambiguous bases | **0** |
| `rl` | repeat-seed length | **0** |
| junctions | # splice junctions (CIGAR `N`) | **0** |

Worked example — read `m64076_221110_210557/50726963/ccs`:

```
tag   DAZ1    DAZ3
AS    1882    1882
ms    2202    2202
NM    0       0
de    0.0000  0.0000
s1    1995    1995
cm    652     652
nn    0       0
rl    231     231
junc  13      13
```

MAPQ is 0 on both placements — minimap2 itself flags them ambiguous.

## Why no tag can ever break them

These reads fall entirely on the **shared backbone** — sequence that is identical
between the two 99.97%-identical, inverted copies. When the read's bases are the
same in both copies, minimap2 computes the *same alignment* against each, so every
derived tag is identical *by construction*. This is the mathematical signature of
true non-identifiability — a property of the molecule, not a missing feature.

Phasing/linkage cannot rescue them either: they cover **zero** copy-distinguishing
positions, so there is nothing to link on.

## What *can* reduce the ambiguity (it is a spectrum, not yes/no)

By the AS-ratio of the two placements (worse ÷ better; always present on every
record; matches minimap2 `s2/s1` to ~0.004):

- **39 fully tied** (AS-ratio ≈ 1) — **apportion.** Identical on every tag; the EM
  does not guess. The **anchored prior** sends their weight to the copy that other,
  resolvable reads anchor (DAZ1), so they add ~nothing to DAZ3 (no fabrication).
- **84 low-confidence** (AS-ratio 0.95–0.999) — **sharpen.** NOT truly tied: **all
  84/84 differ on ≥1 tag** (small `AS`/`NM`/`de` gaps). No single tag decides them,
  but combining the tags into one likelihood (`AS` + `de` + `ms` + junction
  compatibility) firms them up — exactly what the EM E-step does.
- **92 confident** (AS-ratio < 0.95) — assigned directly from the tags.

**For other paralog families** (not DAZ's 39): a read covering even one PSV, plus
cross-read linkage, breaks ties. DAZ is the worst case — inverted, identical
backbone, all markers intronic — so its 39 are the floor that *proves the
apportionment is necessary*, not a place a new tag would help.

## Bottom line

The 39 are provably non-identifiable; everything near them is sharpenable by using
the tags jointly; and the anchored prior keeps the unbreakable core from doing harm
(the DAZ3 cov 163 → ~4 correction).

Figure: `daz_tied_breakability.png`.
