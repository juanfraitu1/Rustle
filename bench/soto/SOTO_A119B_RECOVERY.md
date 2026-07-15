# Soto families recovered in A119b (human Iso-Seq → CHM13v2.0)

Soto 2025 defined these human segmental-duplication gene families on T2T-CHM13. **No liftover was
needed** — the BED accessions (`NC_060925.1`…) are CHM13v2.0 RefSeq, byte-identical in length to the
A119b BAM `@SQ` (chr3=201,105,948, chr15=99,753,195, …); we only renamed the chromosome column
(`NC_0609xx.1`→`chrN`). A119b = human Iso-Seq aligned `minimap2 -ax splice:hq -N 50 -p 0.1` to CHM13v2.0.

## Expression recovery (each member = a Soto paralog locus; reads = primary Iso-Seq reads)

- **80/83 families (96%) have EVERY Soto member expressed** in A119b
- **76/83 families (92%) recovered as multi-copy** (≥2 members expressed)
- **357/362 members (98.6%) expressed**

| family | Soto ID | members (Soto) | expressed in A119b | reads |
|---|---|--:|--:|--:|
| RGPD | ID_395 | 6 | 6 | 17580 |
| GUSBP | ID_163 | 6 | 6 | 16835 |
| NOTCH | ID_400 | 6 | 6 | 9908 |
| FAR | ID_65 | 9 | 9 | 9552 |
| ANKRD | ID_280 | 3 | 3 | 8703 |
| BCRP | ID_283 | 5 | 5 | 8558 |
| NPIPB | ID_154 | 14 | 14 | 6793 |
| LRRC | ID_14 | 7 | 7 | 6726 |
| GOLGA | ID_113 | 16 | 16 | 6537 |
| UBE | ID_78 | 5 | 5 | 5803 |
| HERC | ID_188 | 4 | 4 | 5379 |
| POM | ID_12 | 4 | 4 | 5078 |
| ANKRD | ID_251 | 2 | 2 | 4501 |
| GTF | ID_208 | 5 | 5 | 4192 |
| AC | ID_209 | 3 | 3 | 3579 |
| AL | ID_35 | 10 | 10 | 3471 |
| FAM | ID_182 | 4 | 4 | 3348 |
| PMS | ID_8 | 9 | 9 | 3332 |
| SHLD | ID_443 | 3 | 3 | 3317 |
| GOLGA | ID_88 | 4 | 4 | 3254 |
| SRGAP | ID_462 | 4 | 4 | 3199 |
| SPDYE | ID_207 | 17 | 17 | 3059 |
| CNTNAP | ID_245 | 4 | 4 | 2987 |
| AC | ID_215 | 7 | 7 | 1855 |
| ANAPC | ID_63 | 8 | 8 | 1839 |
| TBC | ID_468 | 9 | 9 | 1747 |
| NSUN | ID_407 | 4 | 4 | 1703 |
| LIMS | ID_386 | 2 | 2 | 1652 |
| AC | ID_147 | 4 | 4 | 1618 |
| GOLGA | ID_116 | 6 | 6 | 1287 |
| AL | ID_104 | 3 | 3 | 1252 |
| SEC | ID_448 | 3 | 3 | 1236 |
| NAIPP | ID_399 | 4 | 4 | 1178 |
| AC | ID_156 | 5 | 5 | 1116 |
| TP | ID_474 | 3 | 3 | 1116 |
| WASH | ID_391 | 5 | 5 | 1061 |
| AC | ID_22 | 6 | 6 | 988 |
| FAM | ID_354 | 4 | 4 | 975 |
| BMS | ID_300 | 2 | 2 | 965 |
| TCAF | ID_127 | 3 | 3 | 901 |
| AC | ID_146 | 2 | 2 | 883 |
| ULK | ID_179 | 3 | 3 | 860 |
| AC | ID_212 | 5 | 5 | 856 |
| GTF | ID_374 | 3 | 3 | 727 |
| AC | ID_211 | 4 | 4 | 707 |
| DDX | ID_49 | 3 | 3 | 680 |
| AMY | ID_131 | 6 | 6 | 594 |
| TEKT | ID_175 | 2 | 2 | 538 |
| FGF | ID_359 | 5 | 5 | 446 |
| BOLA | ID_302 | 1 | 1 | 445 |
| CSPG | ID_324 | 5 | 5 | 431 |
| CDH | ID_313 | 2 | 2 | 413 |
| UBE | ID_481 | 5 | 5 | 393 |
| MST | ID_393 | 2 | 2 | 357 |
| LSP | ID_24 | 3 | 3 | 355 |
| H | ID_226 | 5 | 5 | 343 |
| CTSLP | ID_327 | 4 | 4 | 235 |
| ZNF | ID_490 | 3 | 3 | 203 |
| DDX | ID_332 | 6 | 5 ⚠ | 173 |
| CNTNAP | ID_261 | 5 | 5 | 152 |
| FCGR | ID_357 | 3 | 3 | 146 |
| AL | ID_267 | 1 | 1 | 146 |
| AC | ID_141 | 1 | 1 | 133 |
| CR | ID_142 | 1 | 1 | 108 |
| NF | ID_403 | 3 | 3 | 102 |
| MEP | ID_260 | 3 | 3 | 101 |
| DNM | ID_68 | 4 | 4 | 99 |
| AC | ID_71 | 5 | 5 | 97 |
| AC | ID_148 | 2 | 2 | 80 |
| SYT | ID_240 | 3 | 3 | 77 |
| NCF | ID_402 | 3 | 3 | 76 |
| TRIM | ID_43 | 4 | 4 | 58 |
| DEFB | ID_338 | 3 | 3 | 57 |
| AC | ID_167 | 3 | 3 | 54 |
| AC | ID_213 | 4 | 4 | 36 |
| SPAG | ID_458 | 2 | 2 | 25 |
| OR | ID_411 | 3 | 3 | 20 |
| AC | ID_26 | 2 | 2 | 18 |
| PPIAL | ID_431 | 5 | 2 ⚠ | 7 |
| AC | ID_222 | 1 | 1 | 6 |
| AC | ID_92 | 2 | 2 | 5 |
| DUX | ID_348 | 1 | 1 | 3 |
| DEFB | ID_334 | 1 | 0 ⚠ | 0 |

⚠ = not all members expressed. The only gaps are DEFB (ID_334, single tissue-specific defensin, 0
reads), PPIAL4 (ID_431, 2/5 — very low expression, 7 reads), DDX (ID_332, 5/6). Everything else: full.

**This is expression-level recovery (reads present at the loci).** The copy-resolution step — running
our de-novo family detection (`gw_family_catalog --cross-chrom --homology-primary`) to group these
into families and resolve members as distinct copies — is the next panel.
