# PREREG: do near-tied alignments join one tandem copy to the next, and is Rustle fooled (row 1881, 2026-11-21)

User (2026-11-21): "can we try to see if the near tied multimappers also have a lot of mapping artifacts like
readthroughs that join one tandem copy to another and if rustle effectively does not get fooled by that".
Written before any count. Seen: row 1875 (planted reads, primaries only: 0 of 27,000 tandem-pair reads cover both copies).
Workspace `r1881_tandem_joins`.

## Definitions (CHM13 CAT/Liftoff, protein-coding genes not in `chm13_readthrough.ids`)
- TANDEM pair: Ensembl paralogues with last common ancestor Homo sapiens to Catarrhini, same chromosome and strand, not
  overlapping, gap <= 200 kb (all of them, no draw).
- CONTROL pair (ordinary readthrough): a gene and the next gene downstream on its strand, gap <= 200 kb, not
  overlapping, not paralogues at any level of the table.
- A record or a transcript is ON a pair with >= 50 exon bases of either gene on the pair's strand, JOINED with >= 50
  of both.
- Record classes (A119b, AS as in `--seed-pool good`): primary of an untied molecule (second AS < 0.98 x best), primary
  of a near-tied molecule, near-tie secondary (AS >= 0.98 x best). Other secondaries and supplementaries are not counted.
- A joined record has a SINGLE-COPY ALTERNATIVE when the same molecule has another counted record on one copy only of
  that pair: the join is then one of two equally good readings, the artifact the user describes.

## Part A, reads (A119b)
Joined share of records on the pair, by class, TANDEM and CONTROL; joined records with a single-copy alternative; pairs
with >= 2 joined molecules.
- H1 (near ties join more): at TANDEM pairs the joined share of near-tie records (primaries of near-tied molecules +
  near-tie secondaries) is >= 3 x that of primaries of untied molecules. SUPPORTED / NOT.
- H2 (tandem copies join more than ordinary neighbours): joined share of all counted records, TANDEM >= 2 x CONTROL.

## Part B, assemblies on A119b (existing GTFs: Rustle primaries only, Rustle with near-tie secondaries, the mode,
StringTie, FLAIR, IsoQuant, bambu)
Per tool: pairs with a joined transcript and joined transcripts, TANDEM and CONTROL; of the TANDEM pairs with a joined
transcript, how many have an annotated transcript joining both genes. No truth here: a comparison between tools.
- H3a: the mode has a joined transcript at no more TANDEM pairs than StringTie.

## Part C, planted reads (no read is a readthrough: every joined transcript is an error)
Rows 1874 (133 family copies, scenario E) and 1875 (pairs and pseudogenes, scenario E): records of any class joining two
planted loci; then `assemble` with primaries only and with `--seed-pool good`, stage `multicopy` if it runs on such a
BAM (said if not), and `stringtie -L`; transcripts joining two planted loci.
- H3b (not fooled): Rustle writes a joined transcript at no more than 1 planted locus pair in each of the two data
  sets, in every Rustle arm that runs.
"Rustle is not fooled" needs H3b; H3a is the comparison on real data.

## Predictions (mine)
H1 supported. H2 supported. Most joined records at TANDEM pairs have a single-copy alternative. H3b holds (0). H3a: the
mode has more joined pairs than Rustle with primaries only and about as many as StringTie; I do not know which side of
StringTie.

## Result (2026-11-21, job 4494669): H1 NOT SUPPORTED, H2 NOT SUPPORTED, H3b HOLDS, H3a FAILS as written
Pairs: 2,319 TANDEM, 12,447 CONTROL.

Part A, A119b records (joined share of records on the pair):
| class | TANDEM | CONTROL |
|---|---|---|
| primary, untied molecule | 470 of 63,728 = 0.74% | 41,709 of 18,897,967 = 0.22% |
| primary, near-tied molecule | 32 of 19,778 = 0.16% | 978 of 226,355 = 0.43% |
| near-tie secondary | 84 of 129,257 = 0.065% | 920 of 173,996 = 0.53% |
| all counted | 586 of 212,763 = 0.28% | 43,607 of 19,298,318 = 0.23% |
- H1 not supported: at tandem pairs near-tie records join 0.078% (116 of 149,035), untied primaries 0.74%: ten times
  LESS, not three times more. 39 of the 116 near-tie joins have a single-copy alternative.
- H2 not supported: 0.28% against 0.23% (1.2 x, bar 2 x). 83 tandem pairs have a joined molecule, 43 have >= 2.

Part C, planted reads (no readthrough planted):
- Records joining two planted loci: family copies 1 primary, 433 secondaries (276 of 7,980 molecules); pairs and
  pseudogenes 2 primaries, 2,421 secondaries (1,890 of 50,700 molecules). The artifact exists, in secondaries.
- Joined transcripts: Rustle primaries only 0 / 0; Rustle with near-tie secondaries 1 (one pair) / 0; the mode 1 / 0;
  StringTie 0 / 0. H3b holds (at most 1 pair per set). Stage `multicopy` ran on both.

Part B, A119b assemblies, TANDEM pairs with a joined transcript (registered measure) / CONTROL:
Rustle primaries 8 / 398; Rustle with secondaries 198 / 420; mode 198 / 420; StringTie 20 / 654; FLAIR 34 / 1,642;
IsoQuant 6 / 249; bambu 1 / 24. H3a fails: 198 against 20.
Looked at after (post hoc): 189 of the 198 pairs are ONE transcript, `DN_chr5_17456902_1`, a single exon of 113 kb
over the TAF11L array (946 reads). Distinct joined transcripts at tandem pairs: Rustle primaries 13, with secondaries
15, StringTie 49. Without that transcript the mode has 9 pairs against StringTie's 20. So Rustle writes fewer
copy-joining transcripts than StringTie, and near-tie secondaries add two: the 113 kb single-exon model, which is not a
joined read but single-exon coverage merged across an array of single-exon copies, and one at NPIPA7.
Predictions: H1 wrong, H2 wrong, "most joins have a single-copy alternative" wrong (39 of 116), H3b right (I said 0;
it is 1), H3a open in my prediction and it failed.
Limits: `joined` needs 50 exon bases of both genes, so a join into an unannotated copy is not seen; planted reads show
that Rustle does not invent a join, not that it keeps a real one; one sample.
- Added after (the split of the planted joining secondaries by score): family copies 7 near-tie (of 60,457 near-tie
  secondaries), 426 below 0.98 x best (of 78,936); pairs and pseudogenes 0 near-tie (of 180,765), 2,421 below (of
  407,735). The 0.98 filter of `--seed-pool good` removes all but 7 of 2,854 joining secondaries.
