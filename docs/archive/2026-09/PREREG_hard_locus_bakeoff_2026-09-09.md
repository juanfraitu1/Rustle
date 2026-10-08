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

## Outcome (2026-09-09) — `bakeoff/human/hard/`, `bench/hard_locus_bakeoff.py`; our GTF regenerated on the final pipeline (`ours_final_gtf.gtf`, assignments md5 `0590d544`; derived numbers identical to the 09-08 GTF — the GTF is assignment-independent)
⚠ The hard set inside the bakeoff's molecule universe is **4,115**, not 7,657: 2,710 gate rows have their primary in a copy's LOCUS EXTENT but outside its copy span (the EIF3C-class tie-outside molecules) and 569 have no primary in the `-M -L` subset BAM at all.
| # | fuzz 0 (isoseq at 5) | verdict |
|---|---|---|
| P1 | hard set carried: **ours 0.852**, flair 0.473, StringTie 0.441, **isoseq 0.732** | ✓ ours ≥ 85 % (boundary); flair/StringTie ≤ 70 % ✓; isoseq 73 % — missed, not refuted (< ours) |
| P2 | ours-not-X / X-not-ours: flair **15.7 ×**, StringTie **19.2 ×**, isoseq **2.45 ×** | ✓ ✓; isoseq missed (< 3), not refuted (> 1.5) |
| P3 | `derived_none` hard/easy: ours 0.15/0.48, flair 0.53/0.70, StringTie 0.56/0.74, isoseq 0.27/0.45 — **every tool does BETTER on the hard set** | ⛔ refuted for all: the "easy" set is not a matched control — it holds the singleton and poorly-supported molecules no tool carries (metric trap: an unmatched control) |
| P4 | copies with a hard molecule: ours 21; ours-only vs flair [18] / flair-only [22]; vs StringTie ours-only [0,13,15,17,18] / [22]; **vs isoseq ours-only [] / isoseq-only [14]** | ✓ flair, ✓ StringTie, ⛔ isoseq |
| P5 | contested (1,070): **ours 0.729**, flair 0.363, StringTie 0.321, **isoseq 0.804**; isoseq-not-ours 208 vs ours-not-isoseq 126 | ⛔ ours < 75 % and isoseq > ours |
| P6 | the 230 O2-assigned: ours carries 153, 140/153 = 92 % at O2's copy; flair 75 (99 %), StringTie 103 (97 %), isoseq 173 (91 %) | report |
Fuzz 5 for every tool changes nothing by more than 0.01.

### Why P5 fails: a support policy, not resolution (post-hoc control, labelled as such)
Of the 290 contested molecules our GTF does not carry, **220 have a singleton intron chain** (no other molecule
in the region shares it) and 7 are unspliced; isoseq collapse carries 206 of the 290, **175 of them singletons**
— it emits single-read transcripts, our collapse requires `min_reads` 3 (flair's default is also 3). The
42 copy-22 O2-assigned molecules are ALL singleton chains (derived_none for ours and, 40/42, for isoseq) — the
~490 bp insertion is realized inconsistently by the aligner, so no two reads share an intron chain — which is
also why copy 22 shows in no tool's per-copy coverage. **Restricting every tool to chains carried by ≥ 2
molecules** (`--min-mult 2 --bam`; 12,168 of 15,922 molecules): hard set **ours 0.950** / flair 0.507 /
StringTie 0.488 / isoseq 0.735; contested **ours 0.918 / isoseq 0.806** (ours-not-isoseq 126 vs 31, **4.1 ×**);
per-copy: isoseq-only none. The pre-registered P1/P5 hold on the support-matched comparison and fail on the
unmatched one; both are reported.
