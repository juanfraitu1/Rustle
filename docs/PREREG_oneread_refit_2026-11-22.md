# Row 1899: the one-read ranking refitted at today's defaults, judged with each species held out (pre-registered 2026-11-22)

Row 1898 (seen before this registration): `tools/data/oneread_r1806.json` was fitted on 2026-11-18 under that day's defaults;
a refit on the single library WTC11 r2 at today's defaults beats it on all thirteen primate datasets at StringTie's size
(Sumatran testis from 686 behind StringTie to 1,168 ahead) and loses on mouse ES and mouse cortex. Neither the stale model nor
a one-library model is right. This row refits the same 23 features and the same design on all datasets at today's defaults
and judges the fit only on data of a species it did not see.

## Fit
- Features and design exactly as `oneread_tier.py` (dumps of row 1898: `r1898_replication_score/feat/NAME.npz`, one-read chains
  the default lacks, standardised inside each library). Label = intron chain in the dataset's truth (complete truth on the
  simulated sets; annotation elsewhere). Logistic fit `chain_rules.fit`.
- Balance: at most 20,000 one-read chains per dataset, drawn with a fixed seed (numpy `default_rng(1899)`).
- **Held-out models:** for each of nine species groups (fly; zebrafish; mouse: simulated + ES + heart + adrenal + cortex; human:
  WTC11 r2 + A119b + HG002; gorilla; chimpanzee; bonobo; Bornean orangutan; Sumatran orangutan) one model fitted on the other
  eight groups. A dataset is always scored by the model that never saw its species.
- **Final model** (only if the verdict allows): the fit on all twenty datasets.

## Measures (`r1899_oneread_refit/lodo.py`)
Per dataset, tier 2 ordered by the frozen score and by its held-out model:
- correct chains of the ranked list at StringTie's size (the 17 datasets where StringTie is larger than `final.gtf`);
- annotated chains in the best 10% of the tier (the `list_sensitive` cut), all 20;
- row 1880's AHEAD / BEHIND against every rival under the held-out order.

## Reading, fixed now
- **ADOPTABLE** = the held-out order gives at least the frozen order's correct chains at StringTie's size on at least 15 of 17
  datasets and at least 3 of the 4 non-primate ones among them; no comparison that row 1880 found AHEAD turns to BEHIND; and on
  the three simulated sets the best 10% of the tier holds at least 95% of the frozen order's correct chains on at least 2 of 3.
- ADOPTABLE makes the all-data model a candidate replacement of `oneread_r1806.json` (the user's call). It changes the order of
  `PREFIX.oneread.ranked.gtf` and of tier 2 of `PREFIX.ranked.gtf`, and with partners the promoted set; never `final.gtf`.
- Not ADOPTABLE: the frozen model stays; the table says which group fails.

## Limits, stated in advance
- Every one of the twenty datasets has been used in this project; "held out" means held out of THIS fit, by species.
- The label on real libraries is the annotation, with the non-human ape sets being largely computed (row 1883).
- A species group with one dataset (fly, zebrafish) is judged by a model fitted on vertebrates or on everything else.
