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
