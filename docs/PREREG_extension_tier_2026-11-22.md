# Row 1892: extension as a tier of the ranked list only, never in `final.gtf` (pre-registered 2026-11-22)

User on the lost StringTie comparison (Sumatran orangutan testis): "Keep looking". Row 1883: transfrag extension, switched on,
REPLACES chains, and on complete truth it replaces true ones (mouse simulated: 461 true lost, 10 gained), so it stays off.
Row 1882: on Sumatran testis it turns the comparison at StringTie's size. This row asks whether the chains extension would add
can sit in the ranked list as a tier of their own, with `final.gtf` and its tier 1 untouched: then nothing is lost anywhere,
and the question is only whether an added tier of mostly unproven chains keeps the list ahead of every rival on every dataset.

**Design, fixed now.** Tier E = the multi-exon chains of the extension run's `final.gtf` (`RUSTLE_EXTEND_TRANSFRAGS=1`, row
1883's runs) that the default `final.gtf` lacks, by `reads`, most first. List = tier 1 (default `final.gtf` by reads), tier E,
then the one-read tier by `oneread_score` (chains already in tier E left out). One alternative position is not tried: after the
one-read tier it could not change any comparison this row is about.

## Measure (existing runs of rows 1880 and 1883, no new run; `r1892_extension_tier/etier.py`)
Per dataset and rival: correct chains of the list at the rival's size, without and with tier E; tier E's size and precision.

## Reading, fixed now
- **ADOPTABLE** = with tier E the list is AHEAD of StringTie on Sumatran testis, and no comparison that row 1880 found AHEAD
  turns to BEHIND on any of the twenty datasets (the three simulated sets included, where rivals are smaller than tier 1 and
  nothing can move).
- **HARMLESS, NOT ENOUGH** = no AHEAD turns, Sumatran testis stays BEHIND.
- **HARMFUL** = an AHEAD turns. Then the idea is closed.
- In every case this is the last design tried on this sample: tier position, tier content and ranking were each registered
  once (rows 1890, 1892); a further variant would be a search on one library.
- ADOPTABLE would make tier E a candidate for an opt-in switch like the third tier of row 1891 (it costs one more assembly);
  the user decides.

## Limits
- Tier E on complete truth is expected to be almost entirely wrong (row 1883); it cannot hurt there only because no rival is
  larger than tier 1 on the simulated sets. That is a property of these three datasets, not a guarantee.
- The truth that rewards tier E on the non-human apes is the kind row 1883 found suspect.

## Result (2026-11-22): HARMLESS, NOT ENOUGH; Sumatran testis -686 -> -60, still behind

`r1892_extension_tier/e.NAME.txt` (`etier.py`), existing runs of rows 1880 and 1883.
- **No AHEAD turns:** all 65 comparisons that row 1880 found AHEAD stay AHEAD with tier E.
- **Sumatran testis against StringTie:** 16,295 (-686) -> 16,921 (**-60**) correct chains at StringTie's size. Still BEHIND.
- FLAIR on HG002 (two tiers): -6,171 -> -5,384, still BEHIND (row 1891's opt-in third tier is what closes that one).
- Tier E by dataset, chains / annotated share: simulated 31-642 chains at 0.02-0.03; mouse libraries 309-767 at 0.04-0.20; human
  630-15,795 at 0.10-0.24; non-human apes 1,879-9,255 at 0.13-0.28 (Sumatran testis 5,254 at 0.164).
- **Verdict as registered: HARMLESS, NOT ENOUGH.** Tier E is not built. As fixed above, this is the last design tried on this
  sample: extension switched on (rows 1882 / 1883), the seeding pool of `all` (1889), the order of the one-read tier (1890) and
  extension as a tier of its own (this row) were each registered once; none puts the list ahead of StringTie on Sumatran
  orangutan testis against RefSeq. The remaining gap is 60 chains of 16,981 with tier E and 582-686 without it.
