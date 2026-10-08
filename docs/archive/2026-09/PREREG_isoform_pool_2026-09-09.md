# PREREG — the ISOFORM-level certificate: pooling a chain's reads before the pairwise test (2026-09-09, before the prototype runs)

**Why.** O2 decides read by read; the GTF's unit is the isoform (intron chain). 53 human family transcripts sit
on abstaining reads only (§6hl). A chain's reads are independent observations of the same molecule species: if
the isoform comes from copy A, every read carries A's allele at every PSV it spans, so the pairwise test can be
run on the POOLED trials — Σ decisive columns, Σ matches — and can decide where single reads cannot.

**Rule (prototype `bench/isoform_pool.py` on `ours_final2_dump.star_reads.tsv`, the shipped pipeline's per-read
columns; the same arithmetic the binary uses).** Group the AS-tied, origin-pass molecules by (chrom, intron
chain) from their primary records — the collapse key. For a group with ≥ 2 reads and candidate set C (union of
the reads' candidates): for every pair (A, B), n_AB = Σ_r decisive columns (obs ∈ {A allele, B allele}, A ≠ B),
k_AB = Σ_r columns where the read matches A; p_AB = P(Bin(n_AB, e/3) ≥ k_AB) with e = 0.003 (the read-level
`psv_qual` is empty, so every column has the same error); LLR_AB = (2 k_AB − n_AB) · ln((1−e)/(e/3)).
bk = maximin_A min_B LLR_AB (twins tie at 0 → the group stays tied); **isoform assigned** iff every B ≠ bk has
p_bk,B < alpha/(|C|−1) and LLR_bk,B > 0 (alpha 0.001); **tied** if some B has n = 0; else **ambiguous**.
Single-read groups are the read-level verdicts and are excluded from P2–P4.

| # | prediction | refuted by |
|---|---|---|
| P1 | of the 856 human abstaining contested molecules (512 tied + 344 ambiguous), **≥ 40 %** belong to a chain group with ≥ 2 such reads | < 25 % |
| P2 | of those groups, **20–60 % become isoform-assigned** by pooling | < 10 % or > 80 % |
| P3 | pooled verdicts never contradict a member read's own certificate: for every group containing an `assigned` read, the isoform copy equals that read's copy | any group where they differ |
| P4 | **excision control**: remove the isoform's certified copy from the catalog (the §6hg catalogs), re-dump, re-pool — **≥ 90 %** of certified isoforms abstain | < 80 % |
| P5 | where an isoform's UNIQUE mappers sit at copy A and the pooled abstaining reads certify a copy, that copy is A in **≥ 70 %** (report; a certified copy ≠ A can be a second expressing copy) | — |
| P6 | the 53 evidence-less transcripts (§6hl): pooling gives **≥ 10** of them a certified address | < 4 |
No shipped behaviour changes in this step; the prototype's arithmetic is later ported to the binary (step 3).

## Outcome (2026-09-09) — `bench/isoform_pool.py`, `bench/isoform_pool_excision.py`; dump `ours_final2_dump`, excised dumps `excise_all/pool/no{2,6,8,12,19,22}`
| # | verdict |
|---|---|
| P1 | ⚠ missed, not refuted: **320/856 = 37.4 %** of the abstaining contested reads share a chain with ≥ 1 other contested read; 63 % are chain singletons — nothing to pool |
| P2 | ✓ 78 groups → **35 isoform-assigned (45 %)**, 18 tied, 25 ambiguous — but 25 of the 35 already contained a read-level `assigned` read; **only 27 abstaining reads (3.2 %) gain an address**, 10 groups from abstaining reads alone |
| P3 | ✓ 0 contradictions between the pooled copy and a member read's own certificate |
| P4 | ⛔ **REFUTED: 25/35 = 71.4 % abstain when the certified copy is removed.** The 10 that do not: 2 → 8 (×4, margins 7–62), 2 → 1, 12 → 19 (**margin 366**), 12 → 9, 19 → 21 (**297**), 8 → 6 (×2) — the twins again, and pooling makes the wrong call MORE confident: the pooled margin is measured against the rivals that remain, all of which the reads beat at the columns the twin shares with them; the only guard is the origin certificate, blind below ≈ 0.6 % over a read |
| P5 | 6/9 certified copies coincide with the isoform's unique-mapper copy (67 %; report) |
| P6 | ⚠ missed by one: 9 of the 42 evidence-less groups certified (copies 12 ×4, 8 ×2, 6 ×2, 19) — exactly the twin copies whose excision fails |
**Reading.** Pooling adds power only inside the certificate's sibling limit, and at the twins that is where
the wrong copy lives: read-level excision precision 94 % becomes **71 % at the isoform level**. The pooled
certificate is NOT a placer. What survives for step 3: pooled evidence as a reported attribute, placement by
read-level evidence only (unique mappers, read-level certificates), and one transcript with a copy SET where
neither exists.
