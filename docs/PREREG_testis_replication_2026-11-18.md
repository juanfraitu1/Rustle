# Row 1792: are the "wrong" chains the score keeps in testis real transcripts RefSeq lacks? (pre-registered 2026-11-18, diagnostic)

Rows 1770, 1771 and 1789 each raise precision in cell lines and HG002 and hardly at all in ape testis. Reading offered to the
user: testis precision against RefSeq is capped by real unannotated transcripts. This row tests that reading. It is a
diagnostic: it uses the other individual of the species, so it cannot become a rule. Registered before any count was made.

**Data.** Six species, each with a testis sample T and a cell-line sample C on the same genome (human A119b / HG002; gorilla
OR6737 / KB3781; chimpanzee, bonobo, Bornean and Sumatran orangutan testis / reference individual). Default multi-exon chains
with the register 1769 score (`r1769_multisample_rules/p_r1769.npy`); kept = score >= 0.074.

**Measure.** A chain of one sample is replicated when the other sample of the species assembled the same intron chain (any
score). For the kept chains of a sample, by decile of read support (`l_reads`, deciles of the kept chains):
- replication of kept correct chains (RefSeq match) per decile
- observed replication of kept wrong chains; expected = what correct chains of the same support deciles would give
- ratio R = observed / expected. Also reported: replication of dropped wrong chains.

**Predictions, fixed now.**
- If kept wrong testis chains are mostly real: R(T) >= 0.5 in at least 4 of 6 species.
- If the reading explains why features work in cell lines only (cell-line errors are artifacts, testis "errors" are real):
  R(T) > R(C) in at least 5 of 6 species.
- Both hold: reading supported. First fails: reading not supported by this measure. Only second fails: wrong chains look
  equally real in both tissues, so the tissue difference needs another explanation.

**Limits.** Testis-specific real transcripts are not expected in fibroblasts or lymphoblastoid cells, which lowers R(T) even
if they are real; correct chains of equal support are the reference for that reason, but RefSeq transcripts are likely more
broadly expressed than unannotated ones, so R is biased down. A shared artifact (same aligner, same genome) can replicate too,
which biases R up. One number per species, no test statistics.

## Result (2026-11-18, `r1792_testis_replication/repl.txt`): reading not supported by this measure
| species | R testis | R cell line | kept wrong replicated, testis / cell line | kept correct replicated, testis / cell line |
|---|---|---|---|---|
| human | 0.33 | 0.52 | 0.170 / 0.243 | 0.610 / 0.657 |
| gorilla | 0.61 | 0.67 | 0.348 / 0.369 | 0.724 / 0.734 |
| chimpanzee | 0.39 | 0.49 | 0.206 / 0.122 | 0.604 / 0.431 |
| bonobo | 0.32 | 0.38 | 0.176 / 0.152 | 0.626 / 0.545 |
| Bornean orangutan | 0.29 | 0.69 | 0.162 / 0.234 | 0.632 / 0.547 |
| Sumatran orangutan | 0.23 | 0.61 | 0.123 / 0.116 | 0.560 / 0.372 |

- R(testis) >= 0.5 in 1 of 6 species (prediction: at least 4). R(testis) > R(cell line) in 0 of 6 (prediction: at least 5).
- Kept wrong chains replicate more than dropped wrong chains in every sample (for example A119b 0.170 against 0.101), so the
  score does rank by something that replicates; but relative to correct chains of equal support, testis wrong chains replicate
  less than cell-line wrong chains, the opposite of the reading.
- The stated limit stands: a testis-specific real transcript is not expected in a fibroblast or lymphoblastoid line, so a low
  R(testis) does not show the chains are artifacts. What the row shows is that this measure gives no support for "testis errors
  are mostly real", and that statement should not be made on the present evidence.
