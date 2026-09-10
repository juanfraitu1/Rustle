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
