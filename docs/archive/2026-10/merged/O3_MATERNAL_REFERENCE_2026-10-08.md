# What happens to the reads of a copy the reference haplotype lacks, and what a reference-free step recovers (2026-10-08)

Study of `docs/PREREG_o3_maternal_reference_2026-10-08.md` (body, Amendment 1, Amendment 2); plan `docs/superpowers/plans/2026-10-08-o3-maternal-reference.md`; code `bench/o3_maternal/` (86 unit tests);
work dir `/mnt/linuxdisk/tmp/o3_mat/`; artifact https://claude.ai/artifact/R1QoJa9GLfxnt51haKG6k8 (private). Substrate: KB3781 (male) fibroblast Iso-Seq against its own maternal (`mat`) and paternal (`pat`) assemblies.
Every number is copied from `W/<ref>/fate/fate.tsv`, `side.json`, `W/<ref>/score/*.json`, `W/compare_*.json` and `W/unm.txt`. A fresh reviewer recomputed them from the BAMs and PAFs; the interpretation below was rewritten after that review.
Per prereg Q1 and S7, nothing is pooled across loci.

## 1. Results

**What the "missing" copies are.** 19 of the 20 large judged copies (7 with the mother's genome as reference, 13 with the father's) have a relative on the reference genome within the allele cutoff: nearest-relative identity 0.9911 to 1.0000 (table below). They are absent by the registered synteny rule and mostly present by sequence. GWFAM175_B0 is the only large copy whose nearest relative on the reference is farther than the cutoff (identity 0.9126), so it is the only test of whether a missing copy hides.

**Q1, where the reads of a missing copy go (per copy; mother as reference, then father).**
- No judged copy has an unmapped fate (0 reads unmapped; 1 and 8 reads partial in the two directions), but this is observed, not guaranteed: every read of these sets had already mapped to the combined primary assembly (section 3).
- GWFAM175_B0 (304 reads): absorbed (285 on another locus, 17 on the nearest relative, 2 tied) at median divergence 0.0662 against 0.0007 on the father's genome, MAPQ 60 on both: a divergence pile, not a hole.
- The other 6 judged copies on the mother's genome: REFUTED at 4 (GWFAM208_B0, 214_B0, 64_B17, 227_B0; the first three have identical or near-identical counterparts, identity 1.0000, 0.9995, 0.9999, de about 0.001 on both genomes), SUPPORTED at GWFAM205_B0 (81 of 82 reads tied; the reads' primary sits on a mat copy at 0.9913 with a second equidistant one at 0.9935), between bars at GWFAM64_B11 (3 of 21 tied).
- On the father's genome: REFUTED at 11 of 13 judged copies, PARTLY GWFAM382_B0, SUPPORTED GWFAM390_B1 (20 of 20 tied, two equidistant relatives at 0.9919 and 0.9942). All 13 have relatives at 0.9927 or closer.
- The registered verdict function uses fractions only (its clause "with de ≈ divergence" was not implemented), so REFUTED here includes copies whose reads sit on an identical copy at divergence near zero.
- Smaller expressed copies beyond the cutoff (reported, not judged): on `mat`, GWFAM175_B1 (7 reads: 1 tied, 6 absorbed elsewhere, de 0.0719) and GWFAM175_B2 (4 reads, absorbed); on `pat`, GWFAM26_B2 (17 reads, identity 0.9672: 9 tied, 8 absorbed, de 0.0102), GWFAM26_B3 (15 reads, 0.9898: 5 tied), GWFAM26_B6 (8 reads, 0.9862: 7 tied) and GWFAM26_B1 (4 reads, 0.9879: none tied). The tied-heavy ones are the closest to the advisor's mechanism; their n is 8 to 17 reads.
- Control, all other reads: 6.4% tied on `mat` (1,910 of 29,807), 5.0% on `pat` (1,461 of 29,148). A further 3,996 (`mat`) and 4,560 (`pat`) of the 35,094 reads, 11% and 13%, are tied or unplaced on the truth haplotype and are excluded from both the loci and the control.

**Unmapped reads exist, in one place.** Of the 959 reads unmapped on the combined primary assembly, 132 map on the mother's genome and 3 on the father's. 125 of the 132 lie in one 71-kb stretch of mat chromosome 12 (95.96 to 96.03 Mb, median MAPQ 60, median read length 3,607 bp; none of the 378 family copies is there): an expressed sequence the mother's genome has and the primary assembly and the father's genome do not, whose reads are unmapped without it. This is the unmapped fate, in the one read set that was not pre-filtered by the primary alignment, and it is a single locus.

**LRPAP1.** By the registered rule the only LRPAP1 locus absent from the mother's genome is the chrY copy (sex control, 1 read). p12 is not absent: it lifts as class `T?` to a mat locus 97.2% identical (Amendment 2). The chr12 cluster has 4 loci on `pat` and 3 on `mat`, and the registered lifts are many-to-one there. All 83 p12 reads are TIED with p14 even on `pat` (4 bp in 2,446); on `mat` 66 are absorbed on the p14 ortholog (identity 0.9959), 16 tied, median divergence 0.0025 against 0.0008 on `pat`. Separately, the 9 flagged IsoCon transcripts of the LRPAP1 family match pat chr12 22.55 Mb (c01, the highest-expressed chr12 copy) at identity 0.998 to 1.000, while their best mat match is 0.987 to 0.990 (mat chr3 and chr16 copies): c01's transcripts have no mat sequence within the allele cutoff, although the registered lift calls c01 present (it sends it to mat 34.29 Mb with 222 mismatches in 20 kb). So the extra chr12 locus may be c01, not p12; this study did not test it.

**Q2, recovery from the mother's genome alone.** IsoCon recovers GWFAM175_B0 (R1 PASS) and flags nothing within the cutoff (R2 PASS, 13 copies). The registered R3 FAILS (6 of 24 families without a missing copy, 25%, bar 20%), but 4 of those 6 families (GWFAM309, GWFAM318, GWFAM70, LRPAP1) are sequences found on the paternal genome at identity ≥ 0.998 whose best maternal match is 0.980 to 0.990 and that no truth locus lists, because the truth table (`bonly`) only lists loci on chromosomes the primary assembly took from the reference haplotype. Only GWFAM169 (a partial consensus at 0.95 on both genomes) and GWFAM382 (0.92 on both) match nothing. Both readings are given; the registered verdict is FAIL. The in-house stage (with the Amendment 15/15b fix) flags 1 candidate and recovers nothing (R1 FAIL, R2 PASS, R3 PASS 1 of 24). With the father's genome as reference (held-out direction) R1 is not testable (no judged copy beyond the cutoff), R2 passes for both arms (26 copies), R3 passes (IsoCon 2 of 21, in-house 0 of 21), and no candidate falls in a truth locus; flagged IsoCon candidates in GWFAM162, GWFAM256 and GWFAM70 are again sequences present on the mother's genome at ≥ 0.998 with the best paternal match below 0.99.

**Side by side.** On the genome that has the copy the reads sit on it at median divergence within the allele cutoff (≤ 0.00958) for 20 of the 21 large copies; the exception is GWFAM382_B0, whose reads are 5 to 6% divergent from both genomes (section 6, item 7). The chain's view of GWFAM175_B0: in the mother's-genome run IsoCon has 1 transcript of the copy, kept as a new copy and recovered by a candidate with 4 transcripts; in the father's-genome run the same transcript matches the reference and is not flagged.

## 2. Fate tables (LARGE = at least 20 reads; per locus)

### Reference = mother's genome (`mat`); truth = `pat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | 1 − identity | de (reference) | de (other haplotype) | verdict (fractions only) | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM208_B0 | catalog | 551 | 0 | 1 | 550 | 0 | 1.0000 | 0.0000 | 0.0008 | 0.0008 | REFUTED | held |
| GWFAM175_B0 | catalog | 304 | 0 | 2 | 17 | 285 | 0.9126 | 0.0874 | 0.0662 | 0.0007 | REFUTED | held |
| GWFAM214_B0 | catalog | 122 | 0 | 0 | 122 | 0 | 0.9995 | 0.0005 | 0.0010 | 0.0010 | REFUTED | held |
| LRPAP1_p12 | lrpap1_desc | 83 | 1 | 16 | 66 | 0 | 0.9959 | 0.0041 | 0.0025 | 0.0008 | descriptive |  |
| GWFAM205_B0 | catalog | 82 | 1 | 81 | 0 | 0 | 0.9935 | 0.0065 | – | 0.0006 | SUPPORTED | FAILED |
| GWFAM64_B17 | catalog | 31 | 0 | 0 | 31 | 0 | 0.9999 | 0.0001 | 0.0015 | 0.0015 | REFUTED | held |
| GWFAM227_B0 | catalog | 21 | 0 | 0 | 21 | 0 | 0.9911 | 0.0089 | 0.0057 | 0.0013 | REFUTED | held |
| GWFAM64_B11 | catalog | 21 | 0 | 3 | 18 | 0 | 0.9980 | 0.0020 | 0.0008 | 0.0024 | between bars | FAILED |

Control, all other reads (29,807): 193 partial, 1,910 tied, 27,704 mapped untied, median de 0.0016.

### Reference = father's genome (`pat`); truth = `mat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | 1 − identity | de (reference) | de (other haplotype) | verdict (fractions only) | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM309_B1 | catalog | 482 | 2 | 0 | 480 | 0 | 0.9991 | 0.0009 | 0.0010 | 0.0011 | REFUTED | held |
| GWFAM70_B2 | catalog | 195 | 1 | 0 | 194 | 0 | 0.9999 | 0.0001 | 0.0007 | 0.0007 | REFUTED | held |
| GWFAM402_B0 | catalog | 100 | 2 | 1 | 97 | 0 | 0.9991 | 0.0009 | 0.0018 | 0.0019 | REFUTED | held |
| GWFAM70_B4 | catalog | 97 | 0 | 0 | 97 | 0 | 0.9927 | 0.0073 | 0.0013 | 0.0090 | REFUTED | held |
| GWFAM70_B3 | catalog | 65 | 0 | 0 | 65 | 0 | 0.9998 | 0.0002 | 0.0012 | 0.0012 | REFUTED | held |
| GWFAM390_B0 | catalog | 59 | 0 | 1 | 58 | 0 | 0.9996 | 0.0004 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM169_B1 | catalog | 46 | 0 | 0 | 46 | 0 | 1.0000 | 0.0000 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM162_B0 | catalog | 41 | 0 | 0 | 41 | 0 | 0.9993 | 0.0007 | 0.0008 | 0.0011 | REFUTED | held |
| GWFAM309_B4 | catalog | 35 | 0 | 0 | 35 | 0 | 1.0000 | 0.0000 | 0.0009 | 0.0009 | REFUTED | held |
| GWFAM70_B1 | catalog | 29 | 0 | 0 | 29 | 0 | 0.9967 | 0.0033 | 0.0043 | 0.0026 | REFUTED | held |
| GWFAM169_B0 | catalog | 26 | 0 | 0 | 1 | 25 | 0.9992 | 0.0008 | 0.0008 | 0.0012 | REFUTED | held |
| GWFAM382_B0 | catalog | 22 | 3 | 2 | 17 | 0 | 0.9948 | 0.0052 | 0.0565 | 0.0570 | PARTLY | FAILED |
| GWFAM390_B1 | catalog | 20 | 0 | 20 | 0 | 0 | 0.9942 | 0.0058 | – | 0.0013 | SUPPORTED | FAILED |

Control, all other reads (29,148): 189 partial, 1,461 tied, 27,498 mapped untied, median de 0.0017.

In the control row there is no nearest-relative split, so "elsewhere" there means "mapped untied on a primary". The verdict uses the fractions only.

## 3. Selection caveats

- Every read of the 34-family and LRPAP1 sets was mapped on the combined primary assembly before this test, so UNMAPPED on `mat` or `pat` in the fate table means "mapped there, lost here" and was 0 in practice; reads unmapped on the primary assembly are the separate set of 959 (see above), and they enter no chain net (prereg S6 asked for unmapped reads in the net; this was not done).
- One animal and one tissue. The judged copies are 7 (mother) and 13 (father); no number here is a rate.
- 19 of the 20 judged copies have a relative within the allele cutoff on the reference, which limits what the verdict REFUTED says.
- The R4 read-level O2 readout (reads of a recovered copy before and after) was not computed.

## 4. Amendment 1 predictions

- **S1** (median de lower on the haplotype that has the copy, within the allele range): holds for 4 of 8 large loci with the mother as reference (GWFAM175_B0, 205_B0, 227_B0, p12) and 2 of 13 with the father as reference; the second clause (≤ 0.00958 on the present haplotype) holds for 20 of 21. Failures on the mother's side are copies with an identical or near-identical counterpart (equal medians) and GWFAM64_B11 (higher on `pat`). On the father's side many failures come from the label, not the biology: GWFAM70_B4's 97 reads sit at de 0.0013 on pat (the reference) and 0.0090 on the mat locus they are labelled to, so they belong to a copy the reference has.
- **S2** (the transcripts exist in both runs; kept as new and recovered only in the absent run): holds for GWFAM175_B0 (1 output on the locus in each run; kept new and recovered only in the `mat` run). For the other loci IsoCon has outputs on the locus in both runs (for example GWFAM208_B0 9 and 5) and none is kept as new, because they link to a reference locus within the cutoff.
- **S3** (the `pat` run applies the registered bars unchanged): reported above; R1 not testable, R2 and R3 pass, no candidate in a truth locus.

## 5. Arms

### Reference = `mat`

| arm | R1 | R2 | R3 (bar ≤ 20%, registered rule) | R3 false-flag families, post hoc split | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|---|
| IsoCon | PASS (1 of 1) | PASS (0 of 13 flagged) | FAIL (6 of 24) | other-genome sequence only: GWFAM309, GWFAM318, GWFAM70, LRPAP1; matches nothing: GWFAM169, GWFAM382 | 21 / 9 | 1 of 16 | {'c_unmatched': 2, 'b_other': 6, 'a_recovered': 1} |
| in-house stage | FAIL (0 of 1) | PASS (0 of 13 flagged) | PASS (1 of 24) | other-genome sequence only: none; matches nothing: GWFAM382 | 20 / 1 | 0 of 16 | {'c_unmatched': 1} |

p12 line: IsoCon outputs with identity × coverage ≥ 0.999 on p12 (paternal) 1, on p14 (paternal) 1, on p14 (maternal) 1, on a p12 counterpart in `mat` 0. The p12 and p14 outputs differ by 4 bp, below the allele cutoff, so the link step treats p12 as an allele of p14 (designed-undetectable). The output counted on p14 (maternal) is a 2,562-bp transcript that is about 6% from both p14 copies (a different isoform) and is assigned to the maternal copy on 2 bp; no claim of three sequence types at one locus is made.

### Reference = `pat` (held-out direction)

| arm | R1 | R2 | R3 (bar ≤ 20%, registered rule) | R3 false-flag families, post hoc split | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|---|
| IsoCon | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (2 of 21) | other-genome sequence only: GWFAM256; matches nothing: GWFAM175 | 22 / 7 | 0 of 30 | {'b_other': 4, 'c_unmatched': 3} |
| in-house stage | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (0 of 21) | other-genome sequence only: none; matches nothing: none | 16 / 1 | 0 of 30 | {'c_unmatched': 1} |

## 6. What did not go as predicted

1. The registered prediction "REFUTED at every LARGE locus" fails at GWFAM205_B0 (SUPPORTED) and GWFAM64_B11 (between bars) on `mat`, and at GWFAM390_B1 (SUPPORTED) and GWFAM382_B0 (PARTLY) on `pat`; and 19 of 20 judged copies are within the cutoff of a reference relative, so most of the REFUTED verdicts test nothing about hidden copies.
2. p12 is not a mother-absent locus (Amendment 2); the registered rule leaves the chr12 LRPAP1 difference unresolved, and the chain's transcripts point to c01.
3. The synteny-based truth rule misstates absence both ways: it calls copies with identical counterparts absent (GWFAM208_B0, 214_B0, 64_B17) and calls copies with no counterpart within the cutoff present (c01, GWFAM309/318/70 candidates), and `bonly` cannot list copies on chromosomes the primary assembly took from the reference haplotype.
4. S1 fails for 4 of 8 (mat) and 11 of 13 (pat) loci.
5. IsoCon R3 FAILS on `mat` by the registered rule (25%); the in-house stage fails R1.
6. The in-house arm is not the registered one: the binary contains the Amendment 15/15b fix (commit a13b817f), contrary to the prereg's "not applied" (Amendment 2); its acceptance re-run was not run.
7. Observation, not a result: GWFAM382_B0's 22 reads are 5 to 6% divergent from both haplotype assemblies, so they come from a copy present in neither genome, the class of the 2026-08-13 reference-absent candidates; this study did not test it.
8. The bar clause "with de ≈ divergence" is not in the verdict function.

## 7. Provenance

minimap2 2.30-r1287 (baseline BAM 2.31); reads aligned with `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M`; nearest-relative identities asm20 `-N 50 -p 0.1`; IsoCon 0.3.3 (`IsoCon pipeline --nr_cores 4`, defaults); `o3_candidates` binary built 2026-10-06 09:09 (sha256 prefix 1e01f514e387cadf), after commit a13b817f; delta 0.00958; the `mat` verdicts were written to `W/mat/VERDICTS.txt` at 2026-10-08 08:50:06 and the `pat` files were created after it, the `pat` run applying the same code and constants. Label cross-check against Amendment 10 (`express`, both haplotypes jointly): same order for GWFAM175_B0 304/281, GWFAM205_B0 82/77, GWFAM227_B0 21/13, but GWFAM208_B0 551/0, GWFAM214_B0 122/0, GWFAM64_B17 31/0, GWFAM64_B11 21/0 (those reads tie across haplotypes in Amendment 10 and are untied on `pat` here). Execution rulings and deferred minors are listed in section 9 and in Amendment 2.

## 8. Not covered

Between-individual differences (one animal); the testis library; copies absent from both haplotypes (GWFAM382_B0 above); the two `_pri` copies with a lift below 50% (GWFAM175:2, GWFAM491:1); which chr12 LRPAP1 locus is the extra one (c01 vs p12); a truth table that covers chromosomes the primary assembly took from the reference haplotype; unmapped reads in any chain net; the R4 read-level O2 readout; the Amendment 15 acceptance re-run of the in-house stage; any threshold tuning (none); independent replication.

## 9. Execution rulings and deferred minors (from the execution ledger, kept because the ledger is deleted)

### Rulings
- Task 2: Ruling: tools/rlock.sh has no exec bit — invoke as 'bash tools/rlock.sh' everywhere (plan says tools/rlock.sh) — do not chmod a tracked file — cost if wrong: none
- Task 4: Ruling: p12 lifts as truth_lift class T? to mat chr12:27.43M (97.2% id, cov 1.0); chr12 cluster 4 pat loci vs 3 mat with many-to-one lifts (c01,c03 -> mat 34.29M), so the registered rule finds only chrY (sex) absent from mat for LRPAP1 — p12 kept as DESCRIPTIVE kind lrpap1_desc (no bar, not in R4 expressed), prereg Amendment 2 appended BEFORE labelling — cost if wrong: p12 would be reported as mat-absent when its mat counterpart may be an unrelated diverged copy; the user's LRPAP1 headline shrinks to 'chr12 has one more locus on pat'
- Task 4: Ruling: all 83 reads whose pat primary is on p12 are TIED on pat (p12 vs p14 pat, 99.84%) so the registered untied-placement label gives p12 0 reads — descriptive loci (lrpap1_desc) are labelled by the read's PRIMARY record on the truth haplotype (ties allowed), restricted to the LRPAP1 net; applies to no catalog/bar locus — cost if wrong: p12 reads may include reads that truly belong to p14 (the two are indistinguishable by score, which is itself reported)
- Task 6: Ruling: the mat FASTA order differs from the index/BAM-header order (name_map order-pairing assertion fired: JAQQLJ020000027.1 17,248 bp vs chr3_mat_hsa4_random_utig4-32 194,496 bp) but the length multisets are identical (225, one duplicated length) — pair by length (unique), duplicated lengths by relative order; refuse if the multisets differ — cost if wrong: two same-length scaffolds could swap names (no catalog locus or LRPAP1 copy lies on a scaffold; checked by the nonexistence of such hits in the loci tables)
- Task 6: Ruling: refabsent/copies.<hap>.paf targets are FASTA accessions (CM..), LRPAP1 PAFs use index names — panel chromosomes must be index names (BAM/FASTA/mmi): rename accessions through the same length-based name_map before building the panel (rename_hits) — cost if wrong: none (panel_to_copies fails loudly on an unknown name)
- Task 8: Ruling: o3_candidates on disk (built 2026-10-06 09:09, sha256 prefix 1e01f514e387cadf) post-dates commit a13b817f (2026-10-05 11:27 'Amendment 15/15b consensus fix — majority test, cs normalisation') so it INCLUDES the Amendment 15/15b fix, contradicting the prereg/plan text ('as shipped at b29afa55, Amendment 15 not applied'); use it as is (the fix repairs the documented A14 false-flag defect; a b29afa55 rebuild is a full cargo build and would test a known-broken stage), disclose in prereg Amendment 2 + artifact + write-up as 'o3_candidates with the Amendment 15/15b fix, acceptance (A13/A14/held-out re-run) NOT run' — cost if wrong: the in-house arm's false-flag behaviour is unvalidated rather than known-bad
- Task 12: Ruling: artifact page checked with ONE desktop-width light render (artifact-design rule: one look, then publish); dark theme and 400-px phone width were NOT rendered (palette validated in both modes with validate_palette.js; CSS built mobile-first with wrapping grids and overflow containers) — cost if wrong: a cosmetic issue on phone/dark that the user reports. Artifact: https://claude.ai/artifact/R1QoJa9GLfxnt51haKG6k8 (private). In-house caption states the Amendment 15/15b fix, not the plan's 'known defect'.
- Final: Ruling: the registered verdict function stays as implemented (fractions only) and all post hoc breakdowns are labelled post hoc — changing the verdict after seeing results would be tuning — cost if wrong: REFUTED labels stay generous for loci with identical counterparts, which the write-up now says
- Final: Ruling: reviewer 'Declined to judge' — (a) Amendment 2 ordering stands: it was appended before labels were computed (file edit precedes the labels run in this session; the same commit carries it with truth.py) ; (b) real divergence vs consensus error of the GWFAM309/318 candidates: unknown, stated as unknown; (c) in-house binary internals: accepted as is and disclosed; (d) artifact scan: stands

### Deferred minors from the final review
- '8.7% away' for 175_B0 — 277/285 reads sit on the second asm20 hit (id .9108), not the registered nearest (.9126; 17 reads)
- 205_B0/390_B1 ties are between two equidistant within-cutoff copies (mechanism not stated)
- 'UNMAPPED = 0 by construction' is 0 observed (the chr12 reads show the converse is possible)
- GWFAM382_B0 de on pat 0.0565 (absorbed only) vs 0.0526 (all primaries) unexplained
- section 4 S1 'holds for exactly four' is post hoc
- label_reads family restriction not in prereg S5 (10 mat / 19 pat reads, none at a LARGE locus)
- score.recovered_loci does not filter kind (chrY/p12 hits count as recovered, never false); test pins it
- Review Focus tests 3 and 5 weak (verdict(None) unreachable; sex exclusion not distinguished from non-catalog)
- section 7 'first pat file 08:50:20' (lrpap1.genes.tsv 08:50:13; gate holds)

(Some minors were corrected in passing when the affected sentences were rewritten after the review: the unmapped-by-construction wording and the post hoc S1 wording.)
