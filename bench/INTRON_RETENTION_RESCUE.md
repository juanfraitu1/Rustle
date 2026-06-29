# How much could intron-retention rescue the unassignable? (measured)

The K=0 floor reads are tied because the copies are identical over the **exons** the read observes. The
intron-retention lever: a read that *retains* an intron carries extra (intronic) sequence the spliced-exon gate
never uses — potentially enough to distinguish the copies. This quantifies it on the real GGO Tied reads.

## Factor 1 — do the Tied reads even carry introns? (`bench/intron_retention_rescue.py`)

Of the **972 currently-Tied (unassignable) reads** (consensus splice sites derived from the reads themselves;
a read "retains" an intron if its reference span covers a consensus intron it did not splice):

- **215/972 = 22.1% retain ≥1 intron** (21.6% of the strict-K=0 reads with `n_decisive==0`).

So ~1 in 5 unassignable reads carries intronic sequence to exploit. (The other ~78% are cleanly spliced — for
those, intron-retention offers nothing; they remain on the proven floor.)

## Factor 2 — is the intronic sequence actually distinguishing?

All 972 Tied reads map to **one family — the GWFAM10 8-copy tandem array** (NC_073228.2 ~144.93–145.34 Mb).
Comparing its copies:

- **spliced (exon) pairwise identity: median 95.6%**
- **genomic (intron-included) pairwise identity: median 95.6%**

**The introns are NOT more divergent than the exons.** The common assumption "introns evolve faster" does **not**
hold for these recent tandem duplicates — the whole locus (exons + introns) duplicated recently and diverged
uniformly. (Note the array also contains ~40 bp-offset near-identical copy pairs — copies 0/1 and 5/6 — which are
effectively one locus and stay tied no matter what.)

## What this means — a real but limited lever

- Intron-retention does **not** give a *divergence* bonus here (introns aren't faster-evolving). Its only value
  is a **coverage** bonus: a retained intron extends a read into ~4.4%-divergent extra sequence, giving a
  read that currently spans no exonic distinguishing column more chances to span one.
- So the realistic rescue is **a fraction of the 22%** — the intron-retaining reads whose extended span happens
  to cross a position distinguishing the *specific* copies they are tied between. Reads tied between the
  near-identical (40 bp-offset) duplicates stay tied (their introns are identical too).
- The exact number requires **re-running assignment with intronic positions included** as PSV columns (the
  actual lever), not just measuring retention.

## Honest verdict on "have we exhausted the unassignable?"

- ~**78% of Tied reads are cleanly spliced** → nothing intronic to exploit → they sit on the **proven** K=0
  floor (information-theoretically unassignable from spliced RNA).
- ~**22% retain an intron** → a real but *partial* lever, worth ~a fraction of that once the intron-divergence/
  coverage factor is applied, and only realizable by adding intronic PSV columns to the gate.
- The surprise that **introns here are no more divergent than exons** means the lever is weaker than the
  textbook intuition suggests — on GGO's recent tandem duplicates the whole locus is uniformly diverged, so the
  intron is just "more of the same sequence," not a privileged faster-evolving signal.

**Bottom line:** the unassignable mass is ~78% genuinely irreducible (spliced, on the proven floor) and ~22%
*potentially* touchable by intronic sequence — but with a coverage-not-divergence mechanism, so the achievable
rescue is modest and bounded. The principled endpoint (abstain with a certificate) remains correct; intron-PSV
assignment is a bounded improvement to *try*, not a floor-breaker.

## Implemented under a flag (`RUSTLE_INTRON_PSV=1`) — and measured

To attack the unassignable with every available lever, the intron-PSV path is now implemented:
`copy_assign_pipeline::discover_intron_psvs` aligns each copy's FORWARD genomic span (exons + introns) vs copy[0]
(poasta ≤20 kb, minimap2 above), keeps the substitution columns whose position is **intronic**, and appends them
to the family profile. A read that retains an intron already carries M-aligned bases there, so the existing
per-read CIGAR sweep fills them and the significance gate uses them as extra distinguishing columns — no other
change. Threaded as `Option<&GenomeIndex>` through `assign_family[_detailed]`/`build_family_profiles`; **default
OFF = byte-identical** (env unset → `discover_intron_psvs` not called). TDD: `intron_psv_finds_a_divergent_intron_
column_when_exons_are_identical`; full lib suite 688 green.

**Measured rescue (GWFAM10 region, 8,783 reads):**
| | tied | assigned |
|---|---|---|
| OFF (exon-only) | 5267 | 3067 |
| ON (`RUSTLE_INTRON_PSV=1`) | 5261 | **3094** |

**+27 reads rescued (≈0.5% of the unassignable mass; 6 from strictly-Tied, 21 from ambiguous).** Exactly as the
introns-not-more-divergent finding predicted: the lever is real and now demonstrable, but the gain is sub-1%
because GGO's K=0 copies are uniformly diverged (the intron is "more of the same," not a faster signal). The K=0
floor is near-irreducible even with intronic sequence; abstain-with-certificate remains the right default, and
the flag exists to show the floor was attacked with the last RNA-intrinsic lever.
