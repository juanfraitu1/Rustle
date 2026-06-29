# DAZ — answers to advisor's questions (meeting brief, 2026-06-08)

Region: **NC_073248.2** (chrY). **DAZ1** 42,783,133–42,859,657 (− strand) and
**LOC129530216 = "DAZ protein 3" (DAZ3)** 42,879,918–42,945,552 (+ strand).
Inverted, ~99.97 %-identical copies. All numbers below are freshly re-derived from
`GGO.bam` with correct CIGAR parsing (see caveat at bottom).

## TL;DR for the meeting
1. **The unresolvable reads DO span splice junctions** (median 5, up to 15) — the
   advisor's premise is not borne out. They are unresolvable because their junctions
   are *identical on both copies*, not because they lack junctions. The circled
   unspliced read is the 5 % minority (12/220 reads).
2. **DAZ1 is independently provable; DAZ3 is not.** 3 spliced reads map **uniquely**
   to DAZ1 (MAPQ 60), two of them **full-length 31-junction** transcripts. **Zero**
   spliced reads map uniquely to DAZ3. The advisor is right that DAZ3 has no genuine
   spliced read (its 13 "primary" reads are coin-flip multimappers, MAPQ 0).
3. **Phantom definition is the real issue — but counting is the wrong axis.** Switch
   the guard from "primary-read *count*" to "**decisive** own-copy evidence", and
   report a third bucket (`unresolvable`) so genuinely-unidentifiable copies aren't
   silently dumped as phantoms. This is what addresses "are we losing true matches?".

---

## Q1 — "the 39 unresolvable do not span any splice junction"
**Not supported by the data.** Of 220 reads in the DAZ region, only **12 (5 %) are
unspliced**; 95 % span ≥1 junction.

The "unresolvable" set = **41 fully-tied reads** (identical `AS`, `NM`, `de` on both
copies — minimap2 itself sets MAPQ 0):

| | count |
|---|---|
| fully tied (unresolvable) | **41** |
| …spliced (≥1 junction) | **29** (median 5, max 15) |
| …unspliced | 12 |

So most unresolvable reads *are* spliced. They are uninformative **for copy
assignment** because the junctions they show are shared by both copies — the
copy-distinguishing positions (PSVs) sit in introns/UTRs these reads don't
discriminate. They remain informative **for isoform structure** (they prove a
multi-exon DAZ transcript exists). The read the advisor circled is one of the 12
unspliced ones — real, but the exception, not the rule.

> Honest caveat to raise: the "39" we quoted earlier came from
> `daz_tied_breakability.md`, which (correctly) reported these reads as heavily
> spliced (13-junction worked example). The "39" is a count of tied reads, **not** a
> count of unspliced reads — those are different things and we should not let the
> advisor conflate them.

## Q2 — "do reads over DAZ1-specific exons map uniquely or to both copies?"
**Almost everything multimaps; the tiny unique signal is all on DAZ1.**

| spliced reads | count |
|---|---|
| map to **both** copies | 205 |
| **unique to DAZ1** | **3** (MAPQ 60) |
| unique to DAZ3 | **0** |

The 3 DAZ1-unique reads:
- `…/37881717/ccs` — **31 junctions**, AS 3140, MAPQ 60, spans 42,784,177–42,859,171 (≈full locus)
- `…/146867202/ccs` — **31 junctions**, AS 2108, MAPQ 60, 42,784,917–42,858,753 (≈full locus)
- `…/46531156/ccs` — 3 junctions, AS 3455, MAPQ 60, 42,837,308–42,841,787

Decisiveness of the *multimapping* spliced reads (best alignment to each copy):

| | by AS | by strict NM |
|---|---|---|
| decisive for **DAZ1** (better on DAZ1) | 174 | 173 |
| decisive for **DAZ3** | 2 | 3 |
| tied | 29 | 29 |
| unique to DAZ1 / DAZ3 | 3 / 0 | — |

**Why the advisor "cannot see spliced reads in DAZ3":** DAZ3 carries 13
primary-flagged spliced alignments, but **12 of 13 are exact ties** (AS equal to the
DAZ1 placement) — minimap2's primary flag among MAPQ-0 multimappers is arbitrary.
Only 1–2 reads are genuinely decisive for DAZ3. Filter IGV to MAPQ > 0 and DAZ3 goes
empty, while the 3 MAPQ-60 reads light up DAZ1. DAZ3's apparent coverage is DAZ1
spillover — consistent with our earlier "DAZ3 cov 113 → ~4" correction and the
decision to retire DAZ3 as a showcase.

---

## The phantom worry ("too broad → losing true matches")
DAZ shows the guard problem cuts **both ways**, and the current
`50_authenticity_guard.py` (primary-read **count** ≥ k) is wrong on both:

- **Too lenient:** DAZ3 has 13 primary spliced alignments → count-guard calls it
  **authentic**, though it is a phantom (12/13 are coin-flips, ≤2 decisive).
- **Too strict on the wrong axis (the real "losing true matches" risk):** a genuinely
  expressed copy whose reads happen to *all* be tied (no PSV coverage) can win **zero**
  primaries if the coin-flips land on its sister → it would be flagged phantom even
  though it is real.

**Recommended fix:** make authenticity = **decisive own-copy evidence**, not a primary
count. A read counts for copy C iff it aligns *strictly better* to C (lower NM /
higher AS) **or** maps **uniquely** to C (MAPQ > 0, no sister placement). Then report
**three** buckets, never two:

| bucket | criterion | DAZ1 | DAZ3 |
|---|---|---|---|
| **authentic** (recovery win) | ≥ k decisive reads | ✅ ~174 | ❌ |
| **phantom** | coverage only from sister spillover, ≤ small decisive | — | ✅ (≤2) |
| **unresolvable / unidentifiable** | supported only by tied reads | (its 29 tied) | (its 29 tied) |

The `unresolvable` bucket is the key: a copy supported *only* by tied reads is not a
proven win **and** not provably fake — reporting it separately is what stops us from
"losing true matches" by lumping them into phantoms. It is also exactly honest about
what the data can and cannot show.

## Advisor's "start from copies ≤ 90 % identity"
Agree, and we already have the knob: `MIN_IDENTITY` in `config.sh`. At lower identity
there are more PSVs → more reads become **decisive** → the `unresolvable` bucket
shrinks and recovery claims get clean. Proposal: present the headline recovery wins on
**≤ 90–95 %** families (e.g. RABL2, which is the clean win: rustle-VG 9 vs ST 4), and
keep DAZ (99.97 %) + RBMY (>99 %) as the **identifiability-floor illustration** —
showing *why* decisive evidence, not coverage, is the right currency.

---

### Caveat on provenance of these numbers
Two intermediate scratch passes today miscounted because pysam `cigartuples` are
`(op, length)` and were unpacked as `(length, op)`; that swap made spliced reads look
unspliced. The numbers in **this** document use the corrected parser
(`/tmp/daz_fixed.py`, `/tmp/daz_decisive.py`, `/tmp/daz_unique.py`) and agree with the
original `daz_tied_breakability.md` (heavily-spliced tied reads). If anyone re-runs,
use the corrected scripts.
