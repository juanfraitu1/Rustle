# PREREG — `--best-by-duel`: the read-star best chosen by PAIRWISE duels, not the column-count score (2026-09-09, before the code change)

**Finding it follows from (§6he):** `psv_score(k)` sums ±1 over EVERY column candidate k carries, so a candidate
that covers columns its rival lacks collects matches there (copies 6/7/8 carry the 57 bp segment reference copy 2
lacks) — while the verdict is PAIRWISE, on the columns both carry. When the two disagree, `bk` is the wrong
candidate and the certificate compares everyone against it: in the shipped human output **65 of 1,143 contested
rows have margin < 0** (bk loses a duel; 38 ambiguous, 27 tied, median 1 decisive column), gorilla MCL1 0/33.

**Rule (behind `--best-by-duel`, default OFF; OFF = byte-identical).** Among the non-partner candidates, bk =
argmax_k min_{k' ≠ k} LLR(k vs k'), where LLR is the SAME pairwise sum the certificate uses (±ln((1−e)/(e/3)) per
column both carry where they differ and the read matches one of them; junction terms where both cover, when
`read_star_junctions`). A candidate that beats every rival has min LLR > 0 and is the unique maximin winner; twins
(LLR 0 both ways) tie and fall through to the existing key `(psv_score, matches, Reverse(k))`. Everything after
bk — the certificate, statuses, columns — is untouched.

| # | prediction | refuted by |
|---|---|---|
| P0 | flag OFF byte-identical: human `0590d544`, MCL1 `8f8b4f42`, MCL58 `1934cbf4`; escapes `91081887` / `ff0b8f16` | any md5 differs |
| P1 | **no `assigned` row changes status or copy** in any family (230 / 4 / 2): an assigned bk already beats every rival, so it is the maximin winner | any change |
| P2 | human: of the 65 contested rows with margin < 0, **between 10 and 65 become `assigned`**; none of the others do | 0, or an assignment from a row with margin ≥ 0 |
| P3 | rows with margin ≥ 0 keep bk and status **≥ 99.5 %** (maximin ties at 0 fall through to the old key) | < 99 % |
| P4 | the new assignments: n_decisive median ≥ 4, margin median ≥ 14; excision of each copy receiving ≥ 3 of them (catalogs from §6hg) — **≥ 90 % abstain** | either evidence floor missed by half, or < 80 % abstain |
| P5 | gorilla MCL1 and MCL58 outputs **byte-identical** to the flag-off runs (no margin < 0 rows to move) | any difference |
| P6 | the flag composes: with `--best-by-alignment` also set, the duel still governs (documented); escapes without it unchanged | — |
Human / gorilla never pooled. Default = the user's decision after the measurement.
