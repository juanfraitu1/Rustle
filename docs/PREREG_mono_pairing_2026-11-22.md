# Row 1897: pairing single-exon models with the multi-exon transcripts around them (step 0, pre-registered 2026-11-22)

User goal (2026-11-22): "find how to either pair single exon transcripts with multi-exon ones in a sensible way or ... better
rules for the ones only covered by 1 long read or the junctions shared by 1 read". This row takes the first part.

State: a single-exon model that lies under a spliced gene is dropped whole by the shadow rule (`--polish-mono-shadow`), whatever
its position. Rows 1831 / 1832: the rule raises single-exon precision on real libraries and costs 13-48% of the true
single-exon transcripts that are found without it, equally on mouse and apes. It asks one question ("is there a spliced gene
here") where the relation to that gene can be read: a single-exon model inside the LAST EXON of a spliced transcript and ending
where it ends is a 3' fragment of that transcript; one that crosses an exon-intron boundary is unspliced precursor; one that
lies wholly inside an intron, or on the other strand, has no multi-exon transcript it could be a piece of.

This row measures those pairings. Nothing is built. Workspace `r1897_mono_pairing`.

## Runs
Seven datasets: fly, mouse, zebrafish simulated; mouse ES (mouse contigs), mouse heart; OR6737 gorilla testis, WTC11 r2. Per
dataset `assemble` at today's defaults (`base`) and the same with `--polish-mono-shadow` taken out of the driver snapshot
(`noshadow`), `RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_PTC=0`; `x.gtf` of both.

## Measure (`pair.py`)
SHADOWED = single-exon models of `noshadow` that `base` lacks (no `base` single-exon model overlapping half of both lengths).
Each is paired with the multi-exon models of `base` it overlaps, first class that fits:
1. **3' FRAGMENT:** inside the last exon of a same-strand multi-exon model, its 3' end within 50 bp of that model's 3' end.
2. **INSIDE AN EXON, other:** at least 90% of its length inside one exon of a same-strand model.
3. **EXON + INTRON:** overlaps a same-strand exon and has at least 20 bp in an intron of that model.
4. **INTRONIC:** overlaps a same-strand model's span, no exon of any same-strand model.
5. **OPPOSITE STRAND ONLY:** overlaps multi-exon models of the other strand only.
6. OTHER.
Per class: models; models matching a truth single-exon transcript (at least half of both lengths, same strand: row 1831's
measure); truth single-exon transcripts found. Beside it the single-exon models `base` keeps.

## Reading, fixed now
- A class is **WORTH SPARING** on a dataset when it holds at least 20 shadowed models and their share matching the truth is at
  least the share of the single-exon models `base` keeps.
- **LEAD for a pairing rule** = one class (or classes 4 and 5 together) is worth sparing on at least 5 of 7 datasets, at least 3
  of the 5 non-primate ones, AND classes 1-3 together are under half that share on the same datasets (so the pairing, not the
  shadow as a whole, carries the difference). A lead becomes its own row: the shadow rule applied to classes 1-3 only.
- Reported, no bar: how many class-1 models there are per multi-exon transcript (the "paired" 3' fragments), as a possible
  attribute instead of a drop.
- Otherwise the shadow rule stays as it is.

## Limits, stated in advance
- Annotated single-exon genes are few and incomplete on real libraries; simulated reads have no precursor or antisense
  signal, so classes 3-5 will be nearly empty there and the simulated sets can only show classes 1, 2 and 6.
- The pairing is against Rustle's own multi-exon models, not the truth.

## Result (2026-11-22): no LEAD by the bar; the pairing does separate, by a factor of 2.5 to 17

Jobs 4494884-4494890, both arms exit 0 on seven datasets; `r1897_mono_pairing/pair.NAME.txt` (`pair.py`).
Shadowed single-exon models, matching a truth single-exon transcript / all (share):
| dataset | kept by `base` | 1 3' fragment | 2 inside an exon | 3 exon + intron | 4 intronic | 5 opposite strand only | classes 1-3 | classes 4 + 5 | truth found by 4 + 5, new |
|---|---|---|---|---|---|---|---|---|---|
| fly simulated | 94 / 114 (0.825) | 2 / 2 | 0 | 3 / 5 | 8 / 8 | 3 / 3 | 5 / 7 | 11 / 11 | 12 |
| mouse simulated | 30 / 32 (0.938) | 0 | 0 | 2 / 3 | 3 / 3 | 1 / 1 | 2 / 3 | 4 / 4 | 4 |
| zebrafish simulated | 51 / 72 (0.708) | 0 / 4 | 0 | 0 / 1 | 0 | 0 | 0 / 5 | 0 | 0 |
| mouse ES | 266 / 1,937 (0.137) | 6 / 478 (0.013) | 4 / 77 | 29 / 959 (0.030) | 11 / 304 (0.036) | 19 / 163 (0.117) | 39 / 1,514 (0.026) | 30 / 467 (0.064) | 31 |
| mouse heart | 50 / 87 (0.575) | 0 / 94 | 0 / 11 | 6 / 58 (0.103) | 1 / 6 | 7 / 53 (0.132) | 6 / 163 (0.037) | 8 / 59 (0.136) | 8 |
| OR6737 | 76 / 311 (0.244) | 8 / 711 (0.011) | 3 / 93 | 20 / 1,217 (0.016) | 4 / 30 (0.133) | 20 / 95 (0.211) | 31 / 2,021 (0.015) | 24 / 125 (0.192) | 24 |
| WTC11 r2 | 113 / 183 (0.617) | 2 / 278 (0.007) | 1 / 36 | 12 / 272 (0.044) | 5 / 11 | 14 / 33 (0.424) | 15 / 586 (0.026) | 19 / 44 (0.432) | 19 |

- **No class is WORTH SPARING on any dataset** (bar: at least the share of the models `base` keeps): classes 4 + 5 reach 0.47,
  0.24, 0.79, 0.70 of that share on the four real libraries. Registered verdict: the shadow rule stays as it is.
- **The pairing separates all the same.** Single-exon models that are a piece of a spliced transcript (classes 1-3: 3' fragment,
  inside an exon, exon + intron) match the truth 1.5-3.7% of the time on all four real libraries; the ones with no spliced
  transcript they could be a piece of (intronic, or only an opposite-strand gene) match 6.4%, 13.6%, 19.2%, 43.2%: 2.5, 3.7,
  13 and 17 times as often. Opposite-strand-only is the best single class (0.12-0.42).
- Sparing classes 4 + 5 would find 31 / 8 / 24 / 19 more true single-exon transcripts (+11%, +16%, +32%, +17% of what `base`
  finds) for 437 / 51 / 101 / 25 unannotated models. That is a recall trade at a known price, not a rule that pays.
- **The 3' fragments are a countable thing:** 478 / 94 / 711 / 278 shadowed models sit inside the last exon of a spliced model
  and end where it ends, and 0.7-1.3% of them are annotated single-exon genes. They are what the user's "pair them" asks for:
  evidence for the 3' end of their multi-exon partner, not transcripts.
- The simulated sets hold 4-18 shadowed models: no unspliced precursor or antisense reads are simulated, as stated above.

**Candidates, the user's call (none proposed as a default):** (a) spare classes 4 + 5 from the shadow (recall of single-exon
genes, price above); (b) report class-1 fragments as an attribute of their multi-exon partner instead of dropping them
silently. Either needs the pairing inside the Rust polish (the shadow is applied there).
