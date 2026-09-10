# PREREG — the HARD-LOCUS bakeoff (advisor's metric, Part 0h; 2026-09-09, before any restricted number is computed)

**The advisor's instruction (ADVISOR_QUESTIONS Part 0h):** the deliverable is a GTF judged only on the
transcripts other assemblers struggle to produce; "score only molecules that are AS-tied, and count transcripts
emitted per copy where the competitor emitted none", and whether they are attributed to the right copy.
§6hb scored ALL family-region molecules (ours placed-in-one-copy 0.606 vs flair 0.346 / StringTie 0.306).

**Hard-locus molecule set** = the molecules the AS-tied gate admits with a placement in a catalog copy = the
rows of `bakeoff/human/ours_final.assignments.tsv` (7,657; every one AS-tied by construction; the gate drops
33,547 unique/clear-best molecules). Strata from the same file: `contested` (origin-pass ∧ ≥ 2 candidates,
1,143 = 230 assigned / 531 tied / 382 ambiguous), `tie_outside_catalog`, origin-rejected. The easy set = the
family-region molecules NOT in the file (the §6hb population minus the hard set).

**Rule** = `bench/tool_bakeoff.py` unchanged (exact intron chain → transcript → copy; `--fuzz 5` for isoseq,
reported symmetrically) with a new `--restrict <names>` that scores only the listed molecules; a tool
"carries" a molecule when its state is `derived_one` or `derived_multi`. `bench/hard_locus_bakeoff.py`
compares the per-tool `.calls.tsv` files on the hard set and per stratum. Tools: ours (`ours.gtf`, the GTF
sees every read), flair (`flair_family.gtf`), StringTie (`stringtie_family.gtf`), isoseq (`isoseq_family.gff`).

| # | prediction | refuted by |
|---|---|---|
| P1 | hard set: **ours carries ≥ 85 %**; **each competitor ≤ 70 %** | ours < 75 %, or any competitor ≥ ours |
| P2 | discordance, every competitor X: molecules carried by ours-not-X **≥ 3 ×** X-not-ours | < 1.5 × for any X |
| P3 | hard-vs-easy: each competitor's `derived_none` on the hard set is **≥ 1.5 ×** its rate on the easy set (the hard loci are where they fail); ours changes by < 1.5 × | any competitor < 1.2 ×, or ours ≥ 1.5 × |
| P4 | per copy: ≥ 1 copy where ours carries ≥ 1 hard molecule and X carries none, for every X; **at most 1** copy the other way | 0 such copies for some X, or ≥ 2 the other way |
| P5 | the contested stratum (1,143) repeats P1–P2 | P1 or P2 fails on it |
| P6 | copy attribution of the 230 O2-assigned: our GTF transcript's `copy_index` = O2's copy ≥ 95 % (report; near-tautological, same pipeline); competitors' derived copy = O2's copy — reported, not predicted (their location is the aligner's primary) | — |
Human only (the gorilla flair/StringTie/isoseq arms are TESTIS, `GGO_mm.bam`; never compared to fibroblast).
Fuzz 0 for flair/StringTie/ours, fuzz 5 for isoseq (and fuzz 5 for all, as a secondary table).
